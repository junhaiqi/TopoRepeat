#!/usr/bin/env python3
import sys
import random
import argparse
import pyfastx
import os


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Filter FASTA/FASTQ reads by minimum length and randomly subsample reads.\n\n"
            "Required arguments:\n"
            "  -i / -o / -p"
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument("-i", "--input", required=True,
                        help="Input FASTA/FASTQ file (plain or .gz)")

    parser.add_argument("-o", "--output", required=True,
                        help="Output FASTA/FASTQ file (plain or .gz)")

    parser.add_argument("-p", "--percent", type=float, required=True,
                        help="Sampling percentage (0–100)")

    parser.add_argument("-l", "--length", type=int, default=10000,
                        help="Minimum read length")

    parser.add_argument("-s", "--seed", type=int, default=0,
                        help="Random seed (0 means random)")

    return parser.parse_args()


def detect_format(filename):
    """
    Detect whether input file is FASTA or FASTQ.
    """
    opener = open
    if filename.endswith(".gz"):
        import gzip
        opener = gzip.open

    with opener(filename, "rt") as f:
        first_char = f.read(1)

    if first_char == ">":
        return "fasta"
    elif first_char == "@":
        return "fastq"
    else:
        raise ValueError("Cannot detect file format (not FASTA/FASTQ)")


def process_file(args):

    if args.seed != 0:
        random.seed(args.seed)

    probability = args.percent / 100.0
    min_len = args.length

    total_reads = 0
    passed_length_filter = 0
    retained_reads = 0

    try:
        filetype = detect_format(args.input)

        if filetype == "fastq":
            reader = pyfastx.Fastq(args.input, build_index=False)
        else:
            reader = pyfastx.Fasta(args.input, build_index=False)

        with open(args.output, "w", buffering=1024*1024) as fout:

            if filetype == "fastq":

                for name, seq, qual in reader:
                    total_reads += 1

                    if len(seq) >= min_len:
                        passed_length_filter += 1

                        if random.random() < probability:
                            fout.write(f"@{name}\n{seq}\n+\n{qual}\n")
                            retained_reads += 1

            else:  # FASTA

                for name, seq in reader:
                    total_reads += 1

                    if len(seq) >= min_len:
                        passed_length_filter += 1

                        if random.random() < probability:
                            fout.write(f">{name}\n{seq}\n")
                            retained_reads += 1

    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

    # Statistics
    print("-" * 40)
    print("Processing Complete")
    print(f"Input type:                    {filetype.upper()}")
    print(f"Total Raw Reads:               {total_reads}")
    print(f"Reads >= {min_len}bp:                {passed_length_filter}")
    print(f"Retained Reads ({args.percent}% target): {retained_reads}")

    if total_reads > 0:
        actual_ratio = (retained_reads / total_reads) * 100
        print(f"Actual Retention Ratio:        {actual_ratio:.2f}%")
    print("-" * 40)


if __name__ == "__main__":
    args = parse_args()
    process_file(args)
