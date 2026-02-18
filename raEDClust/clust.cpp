

#include "clust.h"
#include <mutex>

char comp_base(char c)
{
    switch (c)
    {
    case 'A':
        return 'T';
    case 'T':
        return 'A';
    case 'C':
        return 'G';
    case 'G':
        return 'C';
    default:
        return 'N';
    }
}

std::string revcomp(const std::string &s)
{
    std::string rc(s.size(), 'N');
    for (size_t i = 0; i < s.size(); ++i)
        rc[s.size() - 1 - i] = comp_base(s[i]);
    return rc;
}

int approximate_rotation_aware_edit_distance(
    const std::string &s,
    const std::string &t,
    int max_dist)
{
    int best = -1;

    auto check_one_orientation = [&](const std::string &base_t)
    {
        std::string doubled = base_t + base_t;

        EdlibAlignConfig config = edlibNewAlignConfig(
            max_dist,
            EDLIB_MODE_HW,
            EDLIB_TASK_DISTANCE,
            NULL, 0);

        EdlibAlignResult r = edlibAlign(
            s.c_str(), (int)s.size(),
            doubled.c_str(), (int)doubled.size(),
            config);

        if (r.status == EDLIB_STATUS_OK && r.editDistance >= 0)
        {
            if (best == -1 || r.editDistance < best)
            {
                best = r.editDistance;
            }
        }
        edlibFreeAlignResult(r);
    };

    check_one_orientation(t);
    if (best == 0)
        return 0;
    std::string rc = revcomp(t);
    check_one_orientation(rc);

    return best;
}

int rotation_aware_edit_distance(
    const std::string &s,
    const std::string &t,
    int max_dist)
{
    int best = (max_dist >= 0)
                   ? max_dist
                   : std::max(s.size(), t.size());

    const int L = t.size();
    std::string t2 = t + t;

    auto try_all_rotations = [&](const std::string &base)
    {
        std::string doubled = base + base;
        for (int i = 0; i < L; ++i)
        {
            EdlibAlignConfig config = edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_DISTANCE, NULL, 0);

            EdlibAlignResult r = edlibAlign(
                s.c_str(), (int)s.size(),
                doubled.c_str() + i, (int)L,
                config);

            if (r.editDistance >= 0 && r.editDistance < best)
            {
                best = r.editDistance;
                if (best == 0)
                {
                    edlibFreeAlignResult(r);
                    return;
                }
            }
            edlibFreeAlignResult(r);
        }
    };

    // forward
    try_all_rotations(t);

    // reverse complement
    std::string rc = revcomp(t);
    try_all_rotations(rc);

    return best;
}

int rotation_aware_edit_distance_cutoff(
    const std::string &s,
    const std::string &t,
    int &cutoff,
    int max_dist)
{
    int best = (max_dist >= 0)
                   ? max_dist
                   : std::max(s.size(), t.size());

    const int L = t.size();
    std::string t2 = t + t;

    auto try_all_rotations = [&](const std::string &base)
    {
        std::string doubled = base + base;
        for (int i = 0; i < L; ++i)
        {
            EdlibAlignConfig config = edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_DISTANCE, NULL, 0);

            EdlibAlignResult r = edlibAlign(
                s.c_str(), (int)s.size(),
                doubled.c_str() + i, (int)L,
                config);

            if (r.editDistance >= 0 && r.editDistance < best)
            {
                best = r.editDistance;
                if (best == 0 || best <= cutoff)
                {
                    edlibFreeAlignResult(r);
                    return;
                }
            }
            edlibFreeAlignResult(r);
        }
    };

    // forward
    try_all_rotations(t);

    // reverse complement
    std::string rc = revcomp(t);
    try_all_rotations(rc);

    return best;
}

std::vector<std::vector<int>>
rotation_aware_clustering(
    const std::vector<Read> &reads,
    float length_ratio,
    float sim_threshold)
{
    int N = reads.size();
    std::vector<char> assigned(N, 0);
    std::vector<std::vector<int>> clusters;

    int assigned_cnt = 0;

    for (int i = 0; i < N; ++i)
    {
        if (assigned[i])
            continue;

        std::vector<int> cluster;
        cluster.push_back(i);
        assigned[i] = 1;
        assigned_cnt++;

        const std::string &s1 = reads[i].seq;

#pragma omp parallel for schedule(dynamic)
        for (int j = i + 1; j < N; ++j)
        {
            if (assigned[j])
                continue;

            const std::string &s2 = reads[j].seq;

            if (!length_ratio_ok(s1, s2, length_ratio))
                continue;

            int cutoff =
                std::max(s1.size(), s2.size()) * (1.0 - sim_threshold);

            double sim = rotation_similarity_cutoff(s1, s2);

            if (sim >= sim_threshold)
            {
#pragma omp critical
                {
                    if (!assigned[j])
                    {
                        assigned[j] = 1;
                        cluster.push_back(j);
                        assigned_cnt++;
                    }
                }
            }
        }

        clusters.push_back(std::move(cluster));
        if (clusters.size() % 10 == 0 || cluster.size() > 1)
        {
            std::cerr
                << "[Cluster " << clusters.size() << "] "
                << "seed=" << i
                << " size=" << cluster.size()
                << " assigned=" << assigned_cnt
                << " remaining=" << (N - assigned_cnt)
                << "\n";
        }
    }

    return clusters;
}

int select_rep_read_r2rtr(
    const std::vector<int> &cluster,
    const std::vector<Read> &reads)
{
    int best = cluster[0];
    int best_score = 0;

    auto parse_score = [](const std::string &name)
    {
        auto pos = name.find_last_of(':');
        if (pos == std::string::npos)
            return 0;
        return atoi(name.c_str() + pos + 1);
    };

    best_score = parse_score(reads[best].name);

    for (int idx : cluster)
    {
        int s = parse_score(reads[idx].name);
        if (s > best_score)
        {
            best_score = s;
            best = idx;
        }
    }
    return best;
}

int select_rep_read_by_length(
    const std::vector<int> &cluster,
    const std::vector<Read> &reads)
{
    int best = cluster[0];
    size_t best_len = reads[best].seq.size();

    for (int idx : cluster)
    {
        size_t len = reads[idx].seq.size();
        if (len > best_len)
        {
            best_len = len;
            best = idx;
        }
    }
    return best;
}

void write_clusters_fasta(
    const std::vector<std::vector<int>> &clusters,
    const std::vector<Read> &reads,
    const char *out_fa,
    bool &r2rtr)

{
    FILE *fp = fopen(out_fa, "w");
    if (!fp)
    {
        fprintf(stderr, "Error: cannot write %s\n", out_fa);
        exit(1);
    }

    for (size_t cid = 0; cid < clusters.size(); ++cid)
    {
        const auto &cluster = clusters[cid];

        // Select representative by longest sequence
        int rep;
        if (r2rtr)
            rep = select_rep_read_r2rtr(cluster, reads);
        else
            rep = select_rep_read_by_length(cluster, reads);

        fprintf(fp,
                ">cluster%zu|%s\n%s\n",
                cid + 1,
                reads[rep].name.c_str(),
                reads[rep].seq.c_str());
    }

    fclose(fp);
}

void write_clusters_txt(
    const std::vector<std::vector<int>> &clusters,
    const std::vector<Read> &reads,
    const char *clus_file)
{
    FILE *fp = fopen(clus_file, "w");
    if (!fp)
    {
        fprintf(stderr, "Error: cannot write %s\n", clus_file);
        exit(1);
    }

    for (size_t cid = 0; cid < clusters.size(); ++cid)
    {
        const auto &cluster = clusters[cid];
        fprintf(fp, ">cluster%zu\n", cid + 1);
        for (size_t j = 0; j < cluster.size(); ++j)
        {
            fprintf(fp, "%s\n", reads[cluster[j]].name.c_str());
        }
    }

    fclose(fp);
}