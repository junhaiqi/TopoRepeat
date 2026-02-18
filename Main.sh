#!/bin/bash
set -euo pipefail

# ==============================================================================
# Integrated Pipeline: Subsampling -> Minimap2 -> r2rtr -> SRF Analysis
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Configuration & Default Parameters
# ------------------------------------------------------------------------------

# --- Tool Paths (Please adjust these absolute paths or ensure relative paths exist) ---
SCRIPT_DIR=$(dirname "$(readlink -f "$0")")

# Python script for sampling
PY_SCRIPT="${SCRIPT_DIR}/sample_reads.py"

# Minimap2 (Assumes in PATH, or specify absolute path)
MINIMAP2_BIN="minimap2"

# r2rtr
R2RTR_BIN="${SCRIPT_DIR}/r2rtr/r2rtr"

# SRF Tools
RAEDCLUST_BIN="${SCRIPT_DIR}/raEDClust/raEDClust"
SRFUTILS_BIN="${SCRIPT_DIR}/srfutils.js"

# --- Default Values ---
THREADS=16
OUTPUT_DIR=""
INPUT_FILE=""

# Step 1: Sampling defaults
PERCENT=10
MM2_PRESET="ava-pb" # ava-pb or ava-ont
SEED=0
MIN_READ_LEN=10000

# Step 2: r2rtr defaults
R2RTR_MIN_SUPPORT=10  # -n
R2RTR_MIN_LEN=10      # -l
R2RTR_ACCURATE=0      # -c (0 for noisy, 1 for HiFi/ONT-R10)

# Step 3: SRF/Clustering defaults
CLUST_LEN_RATIO=0.85
CLUST_SIM=0.90
R2RTR_MODE=1          # Default 1 for r2rtr output

# Colors
GREEN='\033[0;32m'
BLUE='\033[0;34m'
RED='\033[0;31m'
NC='\033[0m'

# ------------------------------------------------------------------------------
# 2. Help / Usage
# ------------------------------------------------------------------------------
usage() {
    echo -e "${BLUE}Usage: $0 -i <input_fastq> -o <output_dir> [options]${NC}"
    echo ""
    echo "Description: Automated workflow for TR analysis using Subsampling -> r2rtr -> SRF."
    echo ""
    echo "Required:"
    echo "  -i <path>   Input FASTQ file (raw reads)"
    echo "  -o <path>   Output directory"
    echo ""
    echo "Options [Global]:"
    echo "  -t <int>    Threads (default: $THREADS)"
    echo ""
    echo "Options [Step 1: Sampling & Alignment]:"
    echo "  -p <float>  Sampling percentage (default: $PERCENT)"
    echo "  -x <str>    Minimap2 preset for self-alignment (ava-pb/ava-ont) (default: $MM2_PRESET)"
    echo "  -L <int>    Min read length for sampling (default: $MIN_READ_LEN)"
    echo ""
    echo "Options [Step 2: r2rtr Inference]:"
    echo "  -n <int>    r2rtr min support supports (default: $R2RTR_MIN_SUPPORT)"
    echo "  -k <int>    r2rtr min TR length (default: $R2RTR_MIN_LEN)"
    echo "  -c <0|1>    r2rtr accurate reads mode (0=Noisy, 1=HiFi/R10) (default: $R2RTR_ACCURATE, the recommendation is 0.)"
    echo ""
    echo "Options [Step 3: SRF Analysis]:"
    echo "  -S <float>  Clustering similarity threshold (default: $CLUST_SIM)"
    echo "  -R <float>  Clustering length ratio (default: $CLUST_LEN_RATIO)"
    echo ""
    echo "  -h          Show this help message"
    exit 1
}

# ------------------------------------------------------------------------------
# 3. Argument Parsing
# ------------------------------------------------------------------------------
while getopts "i:o:t:p:x:L:n:k:c:S:R:h" opt; do
    case $opt in
        i) INPUT_FILE="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        p) PERCENT="$OPTARG" ;;
        x) MM2_PRESET="$OPTARG" ;;
        L) MIN_READ_LEN="$OPTARG" ;;
        n) R2RTR_MIN_SUPPORT="$OPTARG" ;;
        k) R2RTR_MIN_LEN="$OPTARG" ;;
        c) R2RTR_ACCURATE="$OPTARG" ;;
        S) CLUST_SIM="$OPTARG" ;;
        R) CLUST_LEN_RATIO="$OPTARG" ;;
        h) usage ;;
        *) usage ;;
    esac
done

# Check requirements
if [[ -z "$INPUT_FILE" || -z "$OUTPUT_DIR" ]]; then
    echo -e "${RED}Error: Input file (-i) and Output directory (-o) are required.${NC}"
    usage
fi

# ------------------------------------------------------------------------------
# 4. Preparation & Dependency Check
# ------------------------------------------------------------------------------
echo -e "${GREEN}=== Initializing Pipeline ===${NC}"

# Check input file
if [[ ! -f "$INPUT_FILE" ]]; then
    echo -e "${RED}Error: Input file does not exist: $INPUT_FILE${NC}"
    exit 1
fi
INPUT_ABS=$(readlink -f "$INPUT_FILE")

# Check tools
for tool in "$PY_SCRIPT" "$R2RTR_BIN" "$RAEDCLUST_BIN" "$SRFUTILS_BIN"; do
    if [[ ! -f "$tool" ]]; then
        echo -e "${RED}Error: Tool not found at: $tool${NC}"
        echo "Please edit the 'Tool Paths' section in the script."
        exit 1
    fi
done

if ! command -v "$MINIMAP2_BIN" &> /dev/null; then
    echo -e "${RED}Error: minimap2 not found in PATH.${NC}"
    exit 1
fi

# Create Directory Structure
mkdir -p "$OUTPUT_DIR"
SUBDIR_1="${OUTPUT_DIR}/01_Sampling"
SUBDIR_2="${OUTPUT_DIR}/02_r2rtr"
SUBDIR_3="${OUTPUT_DIR}/03_raEDClust_SRF"
mkdir -p "$SUBDIR_1" "$SUBDIR_2" "$SUBDIR_3"

# Log parameters
echo "Input:    $INPUT_ABS"
echo "Output:   $OUTPUT_DIR"
echo "Threads:  $THREADS"
echo "----------------------------"

# ------------------------------------------------------------------------------
# 5. Step 1: Sampling & Self-Alignment
# ------------------------------------------------------------------------------
echo -e "${GREEN}[Step 1] Sampling Reads & Minimap2 Self-Alignment${NC}"

# Define Step 1 Filenames
BASENAME=$(basename "$INPUT_FILE" | sed 's/.gz$//' | sed 's/.fastq$//' | sed 's/.fq$//')
SAMPLED_FQ="${SUBDIR_1}/${BASENAME}_${PERCENT}pct.fq.gz"
SELF_PAF="${SUBDIR_1}/${BASENAME}_${PERCENT}pct.paf"
LOG_MM2="${SUBDIR_1}/minimap2.log"

# 1.1 Run Python Sampling
if [[ -f "$SAMPLED_FQ" ]]; then
    echo "  > Sampled file exists, skipping sampling..."
else
    echo "  > Running sample_reads.py..."
    python3 "$PY_SCRIPT" \
        -i "$INPUT_ABS" \
        -o "$SAMPLED_FQ" \
        -p "$PERCENT" \
        -l "$MIN_READ_LEN" -s "$SEED"
fi

# 1.2 Run Minimap2 (Self-alignment)
if [[ -f "$SELF_PAF" ]]; then
    echo "  > PAF file exists, skipping alignment..."
else
    echo "  > Running Minimap2 ($MM2_PRESET)..."
    "$MINIMAP2_BIN" -x "$MM2_PRESET" \
        -t "$THREADS" \
        "$SAMPLED_FQ" "$SAMPLED_FQ" \
        > "$SELF_PAF" 2> "$LOG_MM2"
fi

echo "  > Step 1 Output: $SAMPLED_FQ"
echo "  > Step 1 Output: $SELF_PAF"

# ------------------------------------------------------------------------------
# 6. Step 2: r2rtr Unit Inference
# ------------------------------------------------------------------------------
echo -e "${GREEN}[Step 2] Inferring Units with r2rtr${NC}"

# Define Step 2 Filenames
UNITS_FASTA="${SUBDIR_2}/${BASENAME}_r2rtr_units.fasta"
LOG_R2RTR="${SUBDIR_2}/r2rtr.log"

# r2rtr command construction
# Usage: ./r2rtr [Options] <in.paf> -f <read.fq>
# Note: r2rtr output typically goes to stdout, redirecting to file.

echo "  > Running r2rtr..."
echo "    Parameters: -n $R2RTR_MIN_SUPPORT -l $R2RTR_MIN_LEN -c $R2RTR_ACCURATE"

"$R2RTR_BIN" \
    "$SELF_PAF" \
    -f "$SAMPLED_FQ" \
    -n "$R2RTR_MIN_SUPPORT" \
    -l "$R2RTR_MIN_LEN" \
    -c "$R2RTR_ACCURATE" \
    -t "$THREADS" \
    > "$UNITS_FASTA" 2> "$LOG_R2RTR"

if [[ ! -s "$UNITS_FASTA" ]]; then
    echo -e "${RED}Error: r2rtr produced empty output! Check log: $LOG_R2RTR${NC}"
    exit 1
fi

echo "  > Step 2 Output: $UNITS_FASTA"

# ------------------------------------------------------------------------------
# 7. Step 3: SRF Analysis (Clustering -> Enlong -> Abundance)
# ------------------------------------------------------------------------------
echo -e "${GREEN}[Step 3] SRF Analysis (Clustering & Abundance)${NC}"

# Define Step 3 Filenames
PREFIX="${SUBDIR_3}/${BASENAME}"
CLUSTER_FA="${PREFIX}.clustered.fa"
CLUSTER_TXT="${PREFIX}.clusters.txt"
ENLONG_FA="${PREFIX}.enlong.fa"
FINAL_PAF="${PREFIX}.mapping.paf"
BED_FILE="${PREFIX}.bed"
ABUN_FILE="${PREFIX}_abundance.txt"

# 3.1 Clustering (raEDClust)
echo "  > Clustering units (raEDClust)..."
"$RAEDCLUST_BIN" \
    "$UNITS_FASTA" \
    "$CLUSTER_FA" \
    -l "$CLUST_LEN_RATIO" \
    -s "$CLUST_SIM" \
    -c "$CLUSTER_TXT" \
    -x "$R2RTR_MODE" \
    -t "$THREADS"

# 3.2 Enlong (srfutils)
echo "  > Extending clusters (enlong)..."
"$SRFUTILS_BIN" enlong "$CLUSTER_FA" > "$ENLONG_FA"

# 3.3 Mapping Reads back to Units (Minimap2)
# NOTE: Mapping the SAMPLED reads to the units. 
# If you want to map original reads, change $SAMPLED_FQ to $INPUT_ABS below.
echo "  > Mapping sampled reads to units (Minimap2)..."
"$MINIMAP2_BIN" -c \
    -N1000000 \
    -f1000 \
    -r100,100 \
    -t "$THREADS" \
    "$ENLONG_FA" \
    "$SAMPLED_FQ" \
    > "$FINAL_PAF"

# 3.4 PAF -> BED
echo "  > Converting PAF to BED..."
"$SRFUTILS_BIN" paf2bed "$FINAL_PAF" > "$BED_FILE"

# 3.5 BED -> Abundance
echo "  > Calculating Abundance..."
"$SRFUTILS_BIN" bed2abun "$BED_FILE" > "$ABUN_FILE"

# ------------------------------------------------------------------------------
# 8. Completion
# ------------------------------------------------------------------------------
echo -e "${GREEN}=== Pipeline Completed Successfully ===${NC}"
echo "Final Abundance File: $ABUN_FILE"
echo "Inferred Units:       $UNITS_FASTA"
echo "Clustered Units:      $CLUSTER_FA"
echo ""