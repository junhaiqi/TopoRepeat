#ifndef TRCSS_H
#define TRCSS_H

#include <iostream>
#include "spoa.h"

int get_trcss(const std::string &qn_seq, const int &pos, const int &tr_l, std::string &tr_css);
void spoa_css(const std::vector<std::string> &seq_vec, std::string &css);
int refine_trcss(std::vector<std::string> &sub_seq_vec, std::string &tr_css, const float &max_div);
#endif // !TRCSS_H