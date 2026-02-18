#ifndef ALN2TR_H
#define ALN2TR_H

#include "mahit.h"

// typedef robin_hood::unordered_map<std::string, uint32_t[2]> trs_h;    // Read name -> start positon and end position.
typedef robin_hood::unordered_map<std::string, int> trs_h;                                                                                               // Read name -> tandem repeat length.
typedef robin_hood::unordered_map<std::string, int> score_h;                                                                                             // Read name -> tandem repeat length score.
void from_alns_to_trs(const tr_p_aln_info &tr_aln_lib, trs_h &r_trs, const int &min_sp_num, const int &min_tr_l, score_h &s_lib);                        // From muti-alignment infomation to get tandem repeats in a read by calculating distance between two alignments.
int dist_vec_to_tr_l(const std::vector<int> &d_vec, const int &min_sp_num, int &score);                                                                  // From a sorted vector to infer true tandem repeat length.
void get_trs(const char *fa_path, trs_h &r_trs, tr_p_aln_info &tr_aln_lib, score_h &s_lib, const bool &hc);                                              // Get tandem repeat sequences based on inferred length.
void para_get_trs(const char *fa_path, trs_h &r_trs, tr_p_aln_info &tr_aln_lib, score_h &s_lib, const int &batch_size, const int &td_n, const bool &hc); // Muti-thread version of `get_trs`

#endif // !ALN2TR_H./