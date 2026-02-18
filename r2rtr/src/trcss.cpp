#include "trcss.h"
#include <algorithm>
#include "edlib.h"

int get_trcss(const std::string &qn_seq, const int &pos, const int &tr_l, std::string &tr_css)
{

    if (tr_l <= 0 || pos < 0)
    {
        return 0;
    }

    size_t current_pos = static_cast<size_t>(pos);
    size_t step = static_cast<size_t>(tr_l);
    size_t total_len = qn_seq.length();
    
    int max_sub_n = 10; // Subsequences number to generate ccs, it can be changed
    if (step <= 10)
        max_sub_n = 5;

    if (current_pos + step > total_len)
    {
        return 0;
    }

    // size_t estimated_count = (total_len - current_pos) / step;

    std::vector<std::string> sub_seq_vec;
    // sub_seq_vec.reserve(estimated_count);
    sub_seq_vec.reserve(max_sub_n);

    int t = 0;
    while (current_pos + step <= total_len)
    {
        if (t >= max_sub_n)
            break;

        sub_seq_vec.emplace_back(qn_seq.substr(current_pos, step));
        current_pos += step;
        ++t;
    }

    if (sub_seq_vec.empty())
    {
        return 0;
    }

    if (sub_seq_vec.size() == 1)
    {
        tr_css = sub_seq_vec[0];
        return 1;
    }
    else
    {
        spoa_css(sub_seq_vec, tr_css);
        int up_num = refine_trcss(sub_seq_vec, tr_css, 0.90); // 0.90 is max divergence between css and seqs, it can be changed
        return up_num;
    }

    // return sub_seq_vec.size();
}

// int get_trcss(const std::string &qn_seq, const int &pos, const int &tr_l, std::string &tr_css)
// {
//     if (tr_l <= 0 || pos < 0)
//         return 0;

//     const size_t total_len = qn_seq.size();
//     size_t current_pos = static_cast<size_t>(pos);
//     size_t step = static_cast<size_t>(tr_l);

//     const int max_sub_n = 20;
//     const int k = 5;
//     const int w = 5;

//     if (current_pos + step + k > total_len)
//         return 0;

//     const std::string kmer = qn_seq.substr(current_pos, k);

//     std::vector<std::string> sub_seq_vec;
//     sub_seq_vec.reserve(max_sub_n);

//     int t = 0;
//     size_t max_offset = step + 2 * w + k;

//     while (current_pos + max_offset <= total_len && t < max_sub_n)
//     {
//         int best_pos = -1;
//         int best_diff = k;

//         int start = std::max(1, static_cast<int>(step - w));
//         int end = static_cast<int>(step + 2 * w);

//         for (int i = start; i <= end; ++i)
//         {
//             int cur_diff = 0;
//             size_t offset = current_pos + i;
//             for (int j = 0; j < k; ++j)
//             {
//                 if (qn_seq[offset + j] != kmer[j])
//                     ++cur_diff;
//             }
//             if (cur_diff < best_diff)
//             {
//                 best_diff = cur_diff;
//                 best_pos = i;
//             }
//         }

//         if (best_pos <= 0)
//             best_pos = step;

//         sub_seq_vec.emplace_back(qn_seq.substr(current_pos, best_pos));
//         current_pos += best_pos;
//         ++t;
//     }

//     if (sub_seq_vec.empty())
//         return 0;

//     if (sub_seq_vec.size() == 1)
//     {
//         tr_css = sub_seq_vec[0];
//         return 1;
//     }

//     spoa_css(sub_seq_vec, tr_css);
//     return refine_trcss(sub_seq_vec, tr_css, 0.85);
// }

void spoa_css(const std::vector<std::string> &seq_vec, std::string &css)
{
    auto alignment_engine = spoa::AlignmentEngine::Create(
        spoa::AlignmentType::kNW, 3, -5, -3); // linear gaps

    spoa::Graph graph{};

    for (const std::string &it : seq_vec)
    {
        auto alignment = alignment_engine->Align(it, graph);
        graph.AddAlignment(alignment, it);
    }

    css = graph.GenerateConsensus();
}

int refine_trcss(std::vector<std::string> &sub_seq_vec, std::string &tr_css, const float &max_div)
{
    if (sub_seq_vec.empty())
        return 0;
    if (tr_css.empty())
        return 0;

    bool converged = false;
    int iteration = 0;
    const int MAX_ITER = 5; // Prevent infinite loops

    while (!converged && iteration < MAX_ITER)
    {
        iteration++;
        std::vector<std::string> kept_vec;
        kept_vec.reserve(sub_seq_vec.size());

        // Configure Edlib: Global alignment (NW), calculate distance only (fast)
        EdlibAlignConfig config = edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_DISTANCE, NULL, 0);

        bool outliers_found = false;

        // terate through all units to calculate divergence
        for (const auto &seq : sub_seq_vec)
        {
            // Calculate edit distance between the unit (query) and current consensus (target)
            EdlibAlignResult result = edlibAlign(seq.c_str(), seq.length(),
                                                 tr_css.c_str(), tr_css.length(),
                                                 config);

            if (result.status == EDLIB_STATUS_OK)
            {
                // Calculate divergence: EditDistance / ConsensusLength
                // Note: Guard against division by zero
                size_t min_l = std::min(tr_css.length(), seq.length());
                float divergence = (tr_css.length() > 0)
                                       ? (static_cast<float>(result.editDistance) / min_l)
                                       : 1.0f;

                if (divergence < max_div)
                {
                    kept_vec.push_back(seq);
                }
                else
                {
                    outliers_found = true;
                    // std::cerr << "[Info] Dropped outlier. Div: " << divergence << "\n";
                }
            }
            else
            {
                // If alignment fails (rare), discard or keep based on policy. Here we discard.
                outliers_found = true;
            }

            edlibFreeAlignResult(result);
        }

        // Check convergence conditions
        // Condition A: No outliers found in this round -> Converged
        if (!outliers_found)
        {
            converged = true;
        }
        // Condition B: All sequences were filtered out (Bad threshold or bad initial consensus)
        else if (kept_vec.empty())
        {
            // Strategy: Stop filtering, keep the last valid set, or return 0.
            // Here we choose to stop and keep the previous state (or return 0 if strict).
            // std::cerr << "[Warning] All sequences filtered out. Reverting...\n";
            return 0;
        }
        // Condition C: Outliers removed, update vector and regenerate consensus
        else
        {
            // Update the vector
            sub_seq_vec = std::move(kept_vec);

            // Regenerate consensus with the cleaner dataset
            // Optimization: If only 1 seq left, it is the consensus
            if (sub_seq_vec.size() == 1)
            {
                tr_css = sub_seq_vec[0];
                converged = true;
            }
            else
            {
                // Call SPOA again
                spoa_css(sub_seq_vec, tr_css);
            }
        }
    }

    return sub_seq_vec.size();
}