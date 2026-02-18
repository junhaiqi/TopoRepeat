#include "mahit.h"
#include <omp.h>
void f_mul_aln_ps(const char *fn, r_pair_vec &mul_aln_ps)
{
    paf_file_t *pf = paf_open(fn);
    if (pf == 0)
    {
        fprintf(stderr, "Error: Cannot open file %s\n", fn);
        return;
    }

    paf_rec_t r; // Note: In future, developling HOR length prediction method direct use this data structrue.

    while (paf_read(pf, &r) >= 0)
    {
        mul_aln_ps[r.qn][r.tn]++;
    }
    paf_close(pf);
}

void f_aln_ps_tr(const r_pair_vec &mul_aln_ps, r_pair &tr_ps)
{
    for (const auto &it : mul_aln_ps)
    {
        auto max_it = std::max_element(it.second.begin(), it.second.end(),
                                       [](const auto &a, const auto &b)
                                       {
                                           return a.second < b.second;
                                       }); // Select a refence read with max value
        if (tr_ps.find(max_it->first) != tr_ps.end() && tr_ps[max_it->first] == it.first)
        {
            continue;
        }
        tr_ps[it.first] = max_it->first;
    }
}

void f_tr_aln_info(const char *fn, tr_p_aln_info &tr_aln_lib, const int &l_cutoff)
{
    paf_file_t *pf = paf_open(fn);
    if (pf == 0)
    {
        fprintf(stderr, "Error: Cannot open file %s\n", fn);
        return;
    }

    paf_rec_t r;
    tn_alns cur_tn_alns; // All target alignments for the current Query
    std::string last_qn = "";
    bool is_first = true; // Flag to handle the first record

    // Lambda: Encapsulate saving logic to avoid code duplication
    auto save_best_target = [&](const std::string &qn)
    {
        if (cur_tn_alns.empty())
            return; // Prevent crash on empty map

        // Find the vector with the maximum size (best alignment set)
        auto max_it = std::max_element(cur_tn_alns.begin(), cur_tn_alns.end(),
                                       [](const auto &a, const auto &b)
                                       {
                                           return a.second.size() < b.second.size();
                                       });

        // Save to the main library
        tr_aln_lib[qn][max_it->first] = max_it->second;
    };

    while (paf_read(pf, &r) >= 0)
    {
        // Unified filtering logic: Process positive strand only
        // Comment out this line if strand/length filtering is not required
        if (r.rev != 0 || r.ql < l_cutoff || r.tl < l_cutoff)
            continue;

        // Initialize for the first record
        if (is_first)
        {
            last_qn = r.qn;
            is_first = false;
        }

        // New Query Name detected, save the previous batch first
        if (r.qn != last_qn)
        {
            save_best_target(last_qn); // Save the previous group

            // Reset state
            tn_alns().swap(cur_tn_alns); // Use this logic to completely release memory if needed
            // cur_tn_alns.clear(); // Standard clear is sufficient here
            last_qn = r.qn;
        }

        // Add current record (add here if it passed the filter, regardless of whether it's a new group)
        // Note: Using initializer list
        cur_tn_alns[r.tn].push_back({r.ql, r.qs, r.qe, r.tl, r.ts, r.te});
    }

    // Loop finished, ensure the last batch of data is saved
    if (!last_qn.empty())
    {
        save_best_target(last_qn);
    }

    paf_close(pf);
}

void f_tr_aln_info_omp(const char *fn,
                       tr_p_aln_info &tr_aln_lib,
                       int &num_threads)
{
    if (num_threads > 0)
        omp_set_num_threads(num_threads);

    paf_file_t *pf = paf_open(fn);
    if (!pf)
    {
        fprintf(stderr, "Error: Cannot open file %s\n", fn);
        return;
    }

    paf_rec_t r;
    tn_alns cur_tn_alns;
    std::string last_qn;
    bool is_first = true;

#pragma omp parallel
    {
#pragma omp single
        {
            while (paf_read(pf, &r) >= 0)
            {
                if (r.rev != 0)
                    continue;

                if (is_first)
                {
                    last_qn = r.qn;
                    is_first = false;
                }

                if (r.qn != last_qn)
                {
                    std::string qn = last_qn;
                    tn_alns block;
                    block.swap(cur_tn_alns);

#pragma omp task firstprivate(qn, block)
                    {
                        if (!block.empty())
                        {
                            auto max_it = std::max_element(
                                block.begin(), block.end(),
                                [](const auto &a, const auto &b)
                                {
                                    return a.second.size() < b.second.size();
                                });

#pragma omp critical
                            {
                                tr_aln_lib[qn][max_it->first] = max_it->second;
                            }
                        }
                    }

                    last_qn = r.qn;
                }

                cur_tn_alns[r.tn].push_back(
                    {r.ql, r.qs, r.qe, r.tl, r.ts, r.te});
            }

            if (!last_qn.empty())
            {
                std::string qn = last_qn;
                tn_alns block;
                block.swap(cur_tn_alns);

#pragma omp task firstprivate(qn, block)
                {
                    if (!block.empty())
                    {
                        auto max_it = std::max_element(
                            block.begin(), block.end(),
                            [](const auto &a, const auto &b)
                            {
                                return a.second.size() < b.second.size();
                            });

#pragma omp critical
                        {
                            tr_aln_lib[qn][max_it->first] = max_it->second;
                        }
                    }
                }
            }

#pragma omp taskwait
        }
    }

    paf_close(pf);
}

void f_tr_aln_info_hifiasm(const char *fn, tr_p_aln_info &tr_aln_lib) // Note: Incomplete development!!!
{
    paf_file_t *pf = paf_open(fn);
    if (pf == 0)
    {
        fprintf(stderr, "Error: Cannot open file %s\n", fn);
        return;
    }

    paf_rec_t r;
    tn_alns cur_qn_alns; // All query alignments for the current target
    std::string last_tn = "";
    bool is_first = true; // Flag to handle the first record

    // Lambda: Encapsulate saving logic to avoid code duplication
    auto save_best_target = [&](const std::string &tn)
    {
        if (cur_qn_alns.empty())
            return; // Prevent crash on empty map

        // Find the vector with the maximum size (best alignment set)
        auto max_it = std::max_element(cur_qn_alns.begin(), cur_qn_alns.end(),
                                       [](const auto &a, const auto &b)
                                       {
                                           return a.second.size() < b.second.size();
                                       });

        // Save to the main library
        tr_aln_lib[tn][max_it->first] = max_it->second;
    };

    while (paf_read(pf, &r) >= 0)
    {
        // Unified filtering logic: Process positive strand only
        // Comment out this line if strand filtering is not required
        if (r.rev != 0)
            continue;

        // Initialize for the first record
        if (is_first)
        {
            last_tn = r.tn;
            is_first = false;
        }

        // New Query Name detected, save the previous batch first
        if (r.tn != last_tn)
        {
            save_best_target(last_tn); // Save the previous group

            // Reset state
            tn_alns().swap(cur_qn_alns); // Use this logic to completely release memory if needed
            // cur_tn_alns.clear(); // Standard clear is sufficient here
            last_tn = r.tn;
        }

        // Add current record (add here if it passed the filter, regardless of whether it's a new group)
        // Note: Using initializer list
        // cur_qn_alns[r.qn].push_back({r.ql, r.qs, r.qe, r.tl, r.ts, r.te});
        cur_qn_alns[r.tn].push_back({r.ql, r.qs, r.qe, r.tl, r.ts, r.te});
    }

    // Loop finished, ensure the last batch of data is saved
    if (!last_tn.empty())
    {
        save_best_target(last_tn);
    }

    paf_close(pf);
}