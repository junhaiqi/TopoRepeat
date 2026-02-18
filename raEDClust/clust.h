#ifndef CLUST_H
#define CLUST_H

#include "edlib.h"
#include <iostream>
#include <omp.h>
#include <atomic>
#include <vector>

char comp_base(char c);
std::string revcomp(const std::string &s);
int rotation_aware_edit_distance(
    const std::string &s,
    const std::string &t,
    int max_dist = -1);

int approximate_rotation_aware_edit_distance(
    const std::string &s,
    const std::string &t,
    int max_dist = -1);

int rotation_aware_edit_distance_cutoff(
    const std::string &s,
    const std::string &t,
    int &cutoff,
    int max_dist = -1);

struct Read
{
    std::string name;
    std::string seq;
};
std::vector<std::vector<int>>
rotation_aware_clustering(
    const std::vector<Read> &seqs,
    float length_ratio,
    float sim_threshold);

inline bool length_ratio_ok(
    const std::string &a,
    const std::string &b,
    double ratio)
{
    double l1 = a.size(), l2 = b.size();
    return std::min(l1, l2) / std::max(l1, l2) >= ratio;
}

inline double rotation_similarity(
    const std::string &a,
    const std::string &b)
{
    int d = rotation_aware_edit_distance(a, b);
    return 1.0 - double(d) / std::max(a.size(), b.size());
}

inline double rotation_similarity_cutoff(
    const std::string &a,
    const std::string &b)
{
    // int d = rotation_aware_edit_distance_cutoff(a, b, cutoff);
    int d = approximate_rotation_aware_edit_distance(a, b);
    return 1.0 - double(d) / std::max(a.size(), b.size());
}

int select_rep_read_r2rtr(
    const std::vector<int> &cluster,
    const std::vector<Read> &reads);

int select_rep_read_by_length(
    const std::vector<int> &cluster,
    const std::vector<Read> &reads);

void write_clusters_fasta(
    const std::vector<std::vector<int>> &clusters,
    const std::vector<Read> &reads,
    const char *out_fa,
    bool &r2rtr);

void write_clusters_txt(
    const std::vector<std::vector<int>> &clusters,
    const std::vector<Read> &reads,
    const char *clus_file);
#endif // !CLUST_H
