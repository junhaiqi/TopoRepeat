#!/bin/bash

# ==============================================================================
# Automated Pipeline Script: FASTQ Subsampling + Minimap2 Alignment
# ==============================================================================

# Default Parameters
PERCENT=10
THREADS=16
SCRIPT_DIR=$(dirname "$(readlink -f "$0")")
PY_SCRIPT="${SCRIPT_DIR}/sample_reads.py"
MM2P=ava-pb

# Color Definitions
GREEN='\033[0;32m'
RED='\033[0;31m'
NC='\033[0m' # No Color

# ------------------------------------------------------------------------------
# Print Usage/Help Function
# ------------------------------------------------------------------------------
usage() {
    echo -e "Description: Sample reads from a FastQ file and generate a PAF file based on Minimap2"
    echo -e "Usage: $0 -i <input_fastq> -d <output_dir> [-p <percent>] [-x <ava-pb/ava-ont>] [-t <threads>]"
    echo ""
    echo "Required:"
    echo "  -i <path>   Input FASTQ file path (supports .fastq, .fq, .gz)"
    echo "  -d <path>   Output directory path (will be created if not exists)"
    echo ""
    echo "Optional:"
    echo "  -p <float>  Sampling percentage (Default: 10)"
    echo "  -x <string> minimap2 -x parameter (Default: ava-pb)"
    echo "  -t <int>    Minimap2 threads (Default: 16)"
    echo "  -h          Show this help message"
    echo ""
    echo "Example:"
    echo "  bash $0 -i /data/input.fq -d /data/output -p 10 -t 24"
    exit 1
}

# ------------------------------------------------------------------------------
# Parse Command Line Arguments
# ------------------------------------------------------------------------------
while getopts "i:d:p:t:h:x:" opt; do
    case $opt in
        i) INPUT_FILE="$OPTARG" ;;
        d) OUTPUT_DIR="$OPTARG" ;;
        p) PERCENT="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        x) MM2P="$OPTARG" ;;
        h) usage ;;
        *) usage ;;
    esac
done

# Check required arguments
if [[ -z "$INPUT_FILE" || -z "$OUTPUT_DIR" ]]; then
    echo -e "${RED}Error: Input file (-i) and Output directory (-d) are required.${NC}"
    usage
fi

# Check if Python script exists
if [[ ! -f "$PY_SCRIPT" ]]; then
    echo -e "${RED}Error: Sampling script not found: $PY_SCRIPT${NC}"
    echo "Please ensure 'sample_reads.py' is in the same directory as this script."
    exit 1
fi

# ------------------------------------------------------------------------------
# Prepare Paths and Filenames
# ------------------------------------------------------------------------------

# 1. Create output directory
if [[ ! -d "$OUTPUT_DIR" ]]; then
    echo "Creating output directory: $OUTPUT_DIR"
    mkdir -p "$OUTPUT_DIR"
fi

# 2. Extract basename (remove path and extension)
FILENAME=$(basename "$INPUT_FILE")
# Remove common extensions
BASENAME="${FILENAME%.gz}"        # Remove .gz
BASENAME="${BASENAME%.fastq}"     # Remove .fastq
BASENAME="${BASENAME%.fq}"        # Remove .fq

# 3. Construct output file paths
# Format: Name_10%.fq.gz
OUT_FQ_NAME="${BASENAME}_${PERCENT}%.fq.gz"
OUT_FQ_PATH="${OUTPUT_DIR}/${OUT_FQ_NAME}"

# Format: paf_Name_10%.paf
OUT_PAF_NAME="paf_${BASENAME}_${PERCENT}%.paf"
OUT_PAF_PATH="${OUTPUT_DIR}/${OUT_PAF_NAME}"

# Log file
LOG_FILE="${OUTPUT_DIR}/${BASENAME}_minimap2.log"

# ------------------------------------------------------------------------------
# Execute Pipeline
# ------------------------------------------------------------------------------

echo -e "${GREEN}=== Pipeline Started ===${NC}"
echo "Input File:      $INPUT_FILE"
echo "Output Dir:      $OUTPUT_DIR"
echo "Sampling Rate:   $PERCENT%"
echo "MM2 Parameter:   $MM2P"
echo "Threads:         $THREADS"
echo "----------------------------------------"

# Step 1: Run Python Sampling Script
echo -e "${GREEN}[Step 1] Running sample_reads.py ...${NC}"

python3 "$PY_SCRIPT" \
  -i "$INPUT_FILE" \
  -o "$OUT_FQ_PATH" \
  -p "$PERCENT" \
  -l 10000 -s 0  # Minimum read length cutoff and random seed, they can be adjusted

# Check if Step 1 was successful
if [[ $? -ne 0 ]]; then
    echo -e "${RED}Error: Sampling script execution failed!${NC}"
    exit 1
fi

echo -e "Sampling complete. File saved to: $OUT_FQ_PATH"
echo "----------------------------------------"

# Step 2: Run Minimap2
echo -e "${GREEN}[Step 2] Running Minimap2 (ava-pb) ...${NC}"
echo "Log file: $LOG_FILE"

# Note: Running sequentially in foreground to ensure dependency
minimap2 -x "$MM2P" \
    "$OUT_FQ_PATH" \
    "$OUT_FQ_PATH" \
    -t "$THREADS" \
    > "$OUT_PAF_PATH" 2> "$LOG_FILE"

# Check if Step 2 was successful
if [[ $? -ne 0 ]]; then
    echo -e "${RED}Error: Minimap2 execution failed! Check log: $LOG_FILE${NC}"
    exit 1
fi

echo -e "${GREEN}=== Pipeline Completed Successfully ===${NC}"
echo "Result PAF File: $OUT_PAF_PATH"