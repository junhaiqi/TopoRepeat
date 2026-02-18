#include "aln2tr.h"
#include "trcss.h"
#include <cmath>
#include <iostream>
#include <omp.h>
#include <zlib.h>
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

void get_trs(const char *fa_path, trs_h &r_trs, tr_p_aln_info &tr_aln_lib, score_h &s_lib, const bool &hc)
{
    gzFile fp;
    kseq_t *ks;

    if ((fp = gzopen(fa_path, "r")) == 0)
    {
        fprintf(stderr, "Error: %s not found!\n", fa_path);
        exit(1);
    }

    float max_len_diff = 0.05;
    ks = kseq_init(fp);

    while (kseq_read(ks) >= 0)
    {
        std::string seq = ks->seq.s;
        std::string name = ks->name.s;

        if (r_trs.find(name) == r_trs.end())
            continue;
        if (tr_aln_lib.find(name) == tr_aln_lib.end())
            continue;

        auto &inner_map = tr_aln_lib[name];
        if (inner_map.empty())
            continue;

        auto it_inner = inner_map.begin();
        if (it_inner->second.empty())
            continue;

        int pos = it_inner->second[0].qs;
        int tr_len = r_trs[name];

        if (pos < 0 || pos >= (int)seq.length())
            continue;

        if (!hc)
        {
            std::transform(seq.begin() + pos, seq.end(), seq.begin() + pos, ::toupper);

            std::string tr_css;
            int cp_num_to_css = get_trcss(seq, pos, tr_len, tr_css);

            if (cp_num_to_css < 2 ||
                tr_css.size() > (1 + max_len_diff) * tr_len ||
                tr_css.size() < (1 - max_len_diff) * tr_len)
                continue;

            if (s_lib.find(name) != s_lib.end())
            {
                std::cout << ">" << name << ":" << std::to_string(tr_len) << ":"
                          << tr_css.length() << ":" << std::to_string(s_lib[name]) << "\n";
                std::cout << tr_css << "\n";
            }
        }
        else
        {

            if (pos + tr_len > (int)seq.length())
                continue;

            std::transform(seq.begin() + pos, seq.begin() + pos + tr_len, seq.begin() + pos, ::toupper);
            std::string tr_css = seq.substr(pos, tr_len);

            if (s_lib.find(name) != s_lib.end())
            {
                std::cout << ">" << name << ":" << std::to_string(tr_len) << ":"
                          << tr_css.length() << ":" << std::to_string(s_lib[name]) << "\n";
                std::cout << tr_css << "\n";
            }
        }
    }
    kseq_destroy(ks);
    gzclose(fp);
}

// Define a simple struct to buffer sequence data
struct SeqRecord
{
    std::string name;
    std::string seq;
};

void para_get_trs(const char *fa_path, trs_h &r_trs, tr_p_aln_info &tr_aln_lib, score_h &s_lib, const int &batch_size, const int &td_n, const bool &hc)
{
    gzFile fp;
    kseq_t *ks;

    if ((fp = gzopen(fa_path, "r")) == 0)
    {
        fprintf(stderr, "Error: %s not found!\n", fa_path);
        exit(1);
    }

    ks = kseq_init(fp);
    float max_len_diff = 0.05;

    std::vector<SeqRecord> batch;
    batch.reserve(batch_size);

    auto process_batch = [&](const std::vector<SeqRecord> &current_batch)
    {
#pragma omp parallel for schedule(dynamic) num_threads(td_n)
        for (size_t i = 0; i < current_batch.size(); ++i)
        {

            const std::string &name = current_batch[i].name;
            std::string seq = current_batch[i].seq;

            if (r_trs.find(name) == r_trs.end())
                continue;

            auto it_aln = tr_aln_lib.find(name);
            if (it_aln == tr_aln_lib.end() || it_aln->second.empty())
                continue;

            auto it_inner = it_aln->second.begin();
            if (it_inner->second.empty())
                continue;

            auto it_score = s_lib.find(name);
            if (it_score == s_lib.end())
                continue;

            int pos = it_inner->second[0].qs;
            int tr_len = r_trs.at(name);

            if (pos < 0 || pos >= static_cast<int>(seq.length()))
                continue;

            if (pos + tr_len > static_cast<int>(seq.length()))
                continue;

            std::string output;

            if (!hc)
            {

                int end_transform = std::min(static_cast<int>(seq.length()), static_cast<int>(seq.length()));
                std::transform(seq.begin() + pos, seq.end(), seq.begin() + pos, ::toupper);

                std::string tr_css;
                int cp_num_to_css = get_trcss(seq, pos, tr_len, tr_css);

                if (cp_num_to_css < 2 ||
                    tr_css.size() > (1 + max_len_diff) * tr_len ||
                    tr_css.size() < (1 - max_len_diff) * tr_len)
                    continue;

                output = ">" + name + ":" + std::to_string(tr_len) + ":" +
                         std::to_string(tr_css.length()) + ":" +
                         std::to_string(it_score->second) + "\n" + tr_css + "\n";
            }
            else
            {
                std::transform(seq.begin() + pos, seq.begin() + pos + tr_len, seq.begin() + pos, ::toupper);

                std::string tr_css = seq.substr(pos, tr_len);

                output = ">" + name + ":" + std::to_string(tr_len) + ":" +
                         std::to_string(tr_css.length()) + ":" +
                         std::to_string(it_score->second) + "\n" + tr_css + "\n";
            }

#pragma omp critical
            {
                std::cout << output;
            }
        }
    };

    while (kseq_read(ks) >= 0)
    {
        batch.push_back({ks->name.s, ks->seq.s});
        if (batch.size() >= static_cast<size_t>(batch_size))
        {
            process_batch(batch);
            batch.clear();
        }
    }

    if (!batch.empty())
    {
        process_batch(batch);
    }

    kseq_destroy(ks);
    gzclose(fp);
}

void from_alns_to_trs(const tr_p_aln_info &tr_aln_lib, trs_h &r_trs, const int &min_sp_num, const int &min_tr_l, score_h &s_lib)
{
    for (const auto &qn_kv : tr_aln_lib)
    {
        std::vector<int> dist_vec;
        for (const auto &rn_kv : qn_kv.second)
        {
            for (int i = 0; i < rn_kv.second.size() - 1; ++i)
            {
                for (int j = i + 1; j < rn_kv.second.size(); ++j)
                {
                    int dist = -1;
                    /*
                    Case 1:
                                        qs     qe
                    ------------------------------ qn
                                        rs     re
                                    ----------------------- tn
                    */
                    int tv1 = rn_kv.second[i].qs - rn_kv.second[i].ts;
                    int tv2 = rn_kv.second[j].qs - rn_kv.second[j].ts;
                    int tv3 = (rn_kv.second[i].ql - rn_kv.second[i].qe) - (rn_kv.second[i].tl - rn_kv.second[i].te);
                    int tv4 = (rn_kv.second[j].ql - rn_kv.second[j].qe) - (rn_kv.second[j].tl - rn_kv.second[j].te);
                    if (tv1 >= 0 && tv2 >= 0)
                    {
                        dist = std::abs(tv1 - tv2);
                        // std::cerr << "Case1:" << qn_kv.first << "\t" << dist << "\n";
                        dist_vec.emplace_back(dist);
                    }
                    // /*
                    // Case 2:
                    //                     qs     qe
                    //                 ------------------------------ qn
                    //                     rs     re
                    // -------------------------------------- tn
                    // */
                    else if (tv3 >= 0 && tv4 >= 0)
                    {
                        dist = std::abs(tv3 - tv4);
                        // std::cerr << "Case2:" << qn_kv.first << "\t" << dist << "\n";
                        dist_vec.emplace_back(dist);
                    }
                }
            }
        }
        std::sort(dist_vec.begin(), dist_vec.end());
        // std::cerr << "**" << dist_vec.size() << "\n";
        int score = 0;
        int tr_l = dist_vec_to_tr_l(dist_vec, min_sp_num, score);

        // if (tr_l != -1)
        //     r_trs[ qn_kv.first ] = tr_l;

        // std::cerr << tr_l << "****\n";

        if (tr_l >= min_tr_l)
        {
            r_trs[qn_kv.first] = tr_l;
            s_lib[qn_kv.first] = score;
        }

        // std::cerr << "success\n";
    }
}

int dist_vec_to_tr_l(const std::vector<int> &d_vec, const int &min_sp_num, int &score)
{
    if (d_vec.empty())
    {
        return -1;
    }

    // Note: If the vector is too large, the O(N^2) complexity will be slow. So, we use a cutoff.
    int cutoff = 500;
    // Note: If the distance is less than $min_dist bp, two alignments are considered the same.
    int min_dist = 5;
    // Note: Length should be supported by $min_l_supp_n distance values
    // int min_l_supp_n = 50;

    if (d_vec.size() == 1)
        return -1;

    std::vector<int> tmp_vec;
    int t = 0;
    for (const int &it : d_vec)
    {
        if (it < min_dist)
            continue;

        else
        {
            if (t >= cutoff)
                break;

            tmp_vec.emplace_back(it);
            ++t;
        }
    }

    std::vector<int> score_vec(cutoff, 0);

    int count = 0;
    for (size_t i = 0; i < tmp_vec.size(); ++i)
    {
        // if (count >= cutoff)
        // {
        //     break;
        // }

        // if (d_vec[i] < min_dist)
        // {
        //     // std::cerr << "**%" << d_vec[i] << "\n";
        //     continue;
        // }
        // else
        // {
        //     count++;
        // }

        for (size_t j = i + 1; j < tmp_vec.size(); ++j)
        {
            int divisor = tmp_vec[i];
            int dividend = tmp_vec[j];
            int remainder = dividend % divisor;

            // std::cerr << d_vec[i] << "\t" << d_vec[j] << "\t" << remainder << "\n";

            // Check if d_vec[j] is roughly an integer multiple of d_vec[i].
            // We use integer arithmetic instead of floating-point division for better performance.
            // Condition: remainder / divisor <= 0.1  OR  remainder / divisor >= 0.9

            bool is_match = false;

            // Case A: The remainder is very small (close to exact division)
            // Equivalent to: remainder < 0.01 * divisor
            if (static_cast<long long>(remainder) * 100 < divisor)
            {
                is_match = true;
            }
            // Case B: The remainder is very large (close to the next multiple)
            // Equivalent to: remainder > 0.99 * divisor
            else if (static_cast<long long>(remainder) * 100 > static_cast<long long>(divisor) * 99)
            {
                is_match = true;
            }

            if (is_match)
            {
                score_vec[i]++;
            }
        }
    }

    // Find the index with the maximum score
    auto max_it = std::max_element(score_vec.begin(), score_vec.end());
    int max_idx = std::distance(score_vec.begin(), max_it);
    score = *max_it;
    // for (const int &it : score_vec)
    // {
    //     std::cerr << "*****" << it << "****" << "\n";
    // }

    if (score < min_sp_num) // Note: Need add a parameter
        return -1;

    return tmp_vec[max_idx];
}