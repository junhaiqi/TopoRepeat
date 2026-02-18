#include <iostream>
#include "clust.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <unordered_map>
#include <vector>

#include "kseq.h"

#include "ketopt.h"
#include "log.h"

#include <zlib.h>
KSEQ_INIT(gzFile, gzread)
#ifdef _OPENMP
#include <omp.h>
#endif

// typedef std::unordered_map<std::string, std::string> SeqDict;

std::vector<Read> read_fasta_vec(const char *fn)
{
    std::vector<Read> reads;

    gzFile fp = gzopen(fn, "r");
    if (!fp)
    {
        fprintf(stderr, "Error: cannot open %s\n", fn);
        exit(1);
    }

    kseq_t *ks = kseq_init(fp);
    while (kseq_read(ks) >= 0)
    {
        reads.push_back({ks->name.s, ks->seq.s});
    }

    kseq_destroy(ks);
    gzclose(fp);

    return reads;
}

void print_usage(FILE *fp)
{
    fprintf(fp,
            "Usage: rotCluster <in.fa> <out.fa> [options]\n\n"
            "Options:\n"
            "  -l FLOAT   Minimum length ratio (default: 0.85)\n"
            "  -s FLOAT   Similarity threshold (default: 0.90)\n"
            "  -c STR     Output clustering results to this file (default: None)\n"
            "  -x BOOL    in.fa is in r2rtr format (default: 0)\n"
            "  -t INT     Number of threads (default: max available)\n"
            "  -h         Show this help message\n");
}

int main(int argc, char **argv)
{
    // -----------------------------
    // Defaults
    // -----------------------------
    float length_ratio = 0.85f;
    float similarity = 0.90f;
    int n_threads = -1;
    bool r2rtr_out = 0;
    ketopt_t opt = KETOPT_INIT;
    char *clus_file = "";
    // -----------------------------
    // Parse options
    // -----------------------------
    int c;
    while ((c = ketopt(&opt, argc, argv, 1, "l:s:c:t:x:h", 0)) >= 0)
    {
        if (c == 'l')
            length_ratio = atof(opt.arg);
        else if (c == 's')
            similarity = atof(opt.arg);
        else if (c == 'c')
            clus_file = opt.arg;
        else if (c == 't')
            n_threads = atoi(opt.arg);
        else if (c == 'x')
            r2rtr_out = atoi(opt.arg);
        else if (c == 'h')
        {
            print_usage(stdout);
            return 0;
        }
    }

    // -----------------------------
    // Positional arguments (I/O placeholders)
    // -----------------------------
    if (argc - opt.ind != 2)
    {
        print_usage(stderr);
        return 1;
    }

    const char *in_fa = argv[opt.ind];
    const char *out_fa = argv[opt.ind + 1];

    // -----------------------------
    // Validate parameters
    // -----------------------------
    if (length_ratio <= 0.0 || length_ratio > 1.0)
    {
        fprintf(stderr, "Error: -l must be in (0, 1]\n");
        return 1;
    }

    if (similarity <= 0.0 || similarity > 1.0)
    {
        fprintf(stderr, "Error: -s must be in (0, 1]\n");
        return 1;
    }

#ifdef _OPENMP
    if (n_threads > 0)
    {
        omp_set_num_threads(n_threads);
    }
#endif

    // -----------------------------
    // Load sequences
    // -----------------------------
    start_main_timer();
    std::stringstream ss;
    fprintf(stderr, "[INFO] Reading FASTA...\n");
    std::vector<Read> seqs = read_fasta_vec(in_fa);
    fprintf(stderr, "[INFO] Loaded %zu sequences\n", seqs.size());

    // -----------------------------
    // Run clustering
    // -----------------------------
    fprintf(stderr, "[INFO] Rotation-aware clustering started\n");

    auto clusters = rotation_aware_clustering(
        seqs,
        length_ratio,
        similarity);

    fprintf(stderr, "[INFO] Clustering finished\n");
    fprintf(stderr, "[INFO] Total clusters: %zu\n", clusters.size());

    // -----------------------------
    // Write output
    // -----------------------------
    write_clusters_fasta(clusters, seqs, out_fa, r2rtr_out);

    if (clus_file != "")
    {
        write_clusters_txt(clusters, seqs, clus_file);
    }

    // -----------------------------
    // Summary
    // -----------------------------
    fprintf(stderr,
            "[SUMMARY]\n"
            "  Input FASTA     : %s\n"
            "  Output FASTA    : %s\n"
            "  Sequences       : %zu\n"
            "  Clusters        : %zu\n"
            "  length_ratio    : %.3f\n"
            "  similarity      : %.3f\n",
            in_fa, out_fa,
            seqs.size(),
            clusters.size(),
            length_ratio,
            similarity);
    print_resource_usage();
    return 0;
}
