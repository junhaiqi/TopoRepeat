#!/bin/bash
set -euo pipefail

# ==============================================================================
# SRF Analysis Pipeline:
# Clustering (raEDClust) -> Enlong -> Minimap2 -> BED -> Abundance
# ==============================================================================

# ---------------------- Configuration Section ----------------------

# Tool Paths
SRF_DIR="$(pwd)/srf"
MINIMAP2_DIR="$(pwd)/minimap2"
RAEDCLUST_BIN="$(pwd)/raEDClust/raEDClust"

SRFUTILS_BIN="$SRF_DIR/srfutils.js"
MINIMAP2_BIN="$MINIMAP2_DIR/minimap2"

# Default Parameters
THREADS=16

# raEDClust default parameters
CLUST_LEN_RATIO=0.85
CLUST_SIM=0.90
R2RTR_MODE=1   # 0: normal fasta, 1: r2rtr format

# Color Definitions
GREEN='\033[0;32m'
RED='\033[0;31m'
NC='\033[0m'

# ------------------------------------------------------------------
# Help / Usage
# ------------------------------------------------------------------
usage() {
    echo -e "Usage: $0 -i <input_fasta> -r <ref_fasta> -o <output_dir> [options]"
    echo ""
    echo "Required:"
    echo "  -i <path>   Input query FASTA (e.g., raw sequencing reads / reference genome)"
    echo "  -r <path>   Reference FASTA file containing inferred TR units (e.g., r2rtr/srf inferred units)"
    echo "  -o <path>   Output directory"
    echo ""
    echo "Optional:"
    echo "  -t <int>    Threads (default: ${THREADS})"
    echo "  -l <float>  raEDClust min length ratio (default: ${CLUST_LEN_RATIO})"
    echo "  -s <float>  raEDClust similarity threshold (default: ${CLUST_SIM})"
    echo "  -x <0|1>    Input query FASTA is r2rtr format (default: ${R2RTR_MODE})"
    echo "  -h          Show this help message"
    echo ""
    echo "Example:"
    echo "  bash $0 -i reads.fa -r units.fa -o results -t 24 -l 0.9 -s 0.95"
    exit 1
}

# ------------------------------------------------------------------
# Parse Arguments
# ------------------------------------------------------------------
while getopts "i:r:o:t:l:s:x:h" opt; do
    case $opt in
        i) INPUT_FILE="$OPTARG" ;;
        r) REF_FILE="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        l) CLUST_LEN_RATIO="$OPTARG" ;;
        s) CLUST_SIM="$OPTARG" ;;
        x) R2RTR_MODE="$OPTARG" ;;
        h) usage ;;
        *) usage ;;
    esac
done

if [[ -z "${INPUT_FILE:-}" || -z "${REF_FILE:-}" || -z "${OUTPUT_DIR:-}" ]]; then
    echo -e "${RED}Error: -i, -r and -o are required.${NC}"
    usage
fi

# ------------------------------------------------------------------
# Environment Checks
# ------------------------------------------------------------------
for f in "$INPUT_FILE" "$REF_FILE"; do
    [[ -f "$f" ]] || { echo -e "${RED}File not found: $f${NC}"; exit 1; }
done

for tool in "$SRFUTILS_BIN" "$MINIMAP2_BIN" "$RAEDCLUST_BIN"; do
    [[ -x "$tool" ]] || { echo -e "${RED}Tool not executable: $tool${NC}"; exit 1; }
done

INPUT_ABS=$(readlink -f "$INPUT_FILE")
REF_ABS=$(readlink -f "$REF_FILE")

mkdir -p "$OUTPUT_DIR"
cd "$OUTPUT_DIR"

INPUT_BASE=$(basename "$INPUT_FILE" .fa)
REF_BASE=$(basename "$REF_FILE" .fa)
PREFIX="${INPUT_BASE}_${REF_BASE}"

# Output files
CLUSTER_FA="${PREFIX}.clustered.fa"
CLUSTER_TXT="${PREFIX}.clusters.txt"
ENLONG_FA="${PREFIX}.enlong.fa"
PAF_FILE="${PREFIX}.paf"
BED_FILE="${PREFIX}.bed"
ABUN_FILE="${PREFIX}_abundance.txt"

# ------------------------------------------------------------------
# Execute Pipeline
# ------------------------------------------------------------------

echo -e "${GREEN}=== SRF Pipeline Started ===${NC}"
echo "Threads: ${THREADS}"
echo "Clustering: length_ratio=${CLUST_LEN_RATIO}, sim=${CLUST_SIM}, r2rtr=${R2RTR_MODE}"
echo "------------------------------------------------------------"

# --- Step 0: raEDClust clustering ---
echo -e "${GREEN}[Step 0] Clustering reference sequences (raEDClust)...${NC}"

"$RAEDCLUST_BIN" \
    "${REF_ABS}" \
    "${CLUSTER_FA}" \
    -l "${CLUST_LEN_RATIO}" \
    -s "${CLUST_SIM}" \
    -c "${CLUSTER_TXT}" \
    -x "${R2RTR_MODE}" \
    -t "${THREADS}" 
    
if [[ $? -ne 0 ]]; then
    echo -e "${RED}Error: raEDClust failed!${NC}"
    exit 1
fi

# --- Step 1: Enlong ---
echo -e "${GREEN}[Step 1] Extending clustered reference (srfutils enlong)...${NC}"
"$SRFUTILS_BIN" enlong "${CLUSTER_FA}" > "${ENLONG_FA}"

# --- Step 2: Minimap2 ---
echo -e "${GREEN}[Step 2] Running Minimap2...${NC}"
"$MINIMAP2_BIN" -c \
    -N1000000 \
    -f1000 \
    -r100,100 \
    -t"${THREADS}" \
    "${ENLONG_FA}" \
    "${INPUT_ABS}" \
    > "${PAF_FILE}"

# --- Step 3: PAF to BED ---
echo -e "${GREEN}[Step 3] PAF → BED...${NC}"
"$SRFUTILS_BIN" paf2bed "${PAF_FILE}" > "${BED_FILE}"

# --- Step 4: BED to Abundance ---
echo -e "${GREEN}[Step 4] BED → Abundance...${NC}"
"$SRFUTILS_BIN" bed2abun "${BED_FILE}" > "${ABUN_FILE}"

# ------------------------------------------------------------------
# Finalize
# ------------------------------------------------------------------
echo -e "${GREEN}=== Pipeline Completed Successfully ===${NC}"
echo "Clustered ref: ${CLUSTER_FA}"
echo "Cluster table: ${CLUSTER_TXT}"
echo "Abundance:     ${ABUN_FILE}"
