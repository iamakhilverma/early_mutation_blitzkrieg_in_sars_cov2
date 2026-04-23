#!/bin/bash
#===============================================================================
# Run Nextclade on pre-filtered GISAID sequences via Apptainer/Singularity
#===============================================================================
#
# Usage:
#   bash scripts/01_run_nextclade.sh [options]
#
# Options:
#   --jobs N         Number of parallel workers (default: all CPUs)
#   --fasta PATH     Path to FASTA (default: filtered_data/filtered_sequences.fasta)
#   --outdir PATH    Nextclade output directory (default: nextclade_output)
#
# Minerva (96 cores):
#   bsub -Is -n 96 -R "rusage[mem=2000]" -W 24:00 -q premium bash
#   bash scripts/01_run_nextclade.sh --jobs 96
#
# Parallelism: Nextclade uses --jobs for embarrassingly parallel processing.
# Each worker needs ~1-2GB RAM. 96 cores ≈ 45-90 min for ~16M sequences.
#===============================================================================

set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
FASTA_INPUT="${PROJECT_DIR}/filtered_data/filtered_sequences.fasta"
NEXTCLADE_DIR="${PROJECT_DIR}/nextclade_output"
DATASET_DIR="${NEXTCLADE_DIR}/dataset_sars-cov-2"
CONTAINER_DIR="${PROJECT_DIR}/containers"
CONTAINER_SIF="${CONTAINER_DIR}/nextclade_latest.sif"
JOBS=0

while [[ $# -gt 0 ]]; do
    case $1 in
        --jobs)    JOBS="$2"; shift 2 ;;
        --fasta)   FASTA_INPUT="$2"; shift 2 ;;
        --outdir)  NEXTCLADE_DIR="$2"; shift 2 ;;
        *)         echo "Unknown: $1"; exit 1 ;;
    esac
done

if [ "$JOBS" -eq 0 ] 2>/dev/null; then
    JOBS=$(nproc 2>/dev/null || echo 16)
fi

echo "================================================================"
echo "Nextclade Analysis (Step 1 of 3)"
echo "================================================================"
echo "  FASTA:   ${FASTA_INPUT}"
echo "  Workers: ${JOBS}"
echo ""

# --- Container runtime ---
# Load module if not already available (Minerva HPC)
if ! command -v apptainer &>/dev/null && ! command -v singularity &>/dev/null; then
    echo "  Loading apptainer module..."
    module load apptainer/1.3.6 2>/dev/null || module load apptainer 2>/dev/null || true
fi

if command -v apptainer &>/dev/null; then
    CMD="apptainer"
elif command -v singularity &>/dev/null; then
    CMD="singularity"
else
    echo "ERROR: apptainer/singularity not found."
    echo "  Try: module load apptainer/1.3.6"
    exit 1
fi

mkdir -p "${NEXTCLADE_DIR}" "${CONTAINER_DIR}"

if [ ! -f "${CONTAINER_SIF}" ]; then
    echo "[1/3] Pulling Nextclade container (requires internet)..."
    ${CMD} pull "${CONTAINER_SIF}" docker://nextstrain/nextclade:latest
    if [ ! -f "${CONTAINER_SIF}" ]; then
        echo "ERROR: Container pull failed. Compute nodes may not have internet."
        echo "  Run on login node first:"
        echo "    ml apptainer/1.3.6"
        echo "    mkdir -p containers"
        echo "    apptainer pull containers/nextclade_latest.sif docker://nextstrain/nextclade:latest"
        exit 1
    fi
else
    echo "[1/3] Container exists: ${CONTAINER_SIF}"
fi

# --- Dataset ---
if [ -d "${DATASET_DIR}" ] && [ -f "${DATASET_DIR}/reference.fasta" ]; then
    echo "[2/3] Dataset already exists: ${DATASET_DIR}"
else
    echo "[2/3] Getting SARS-CoV-2 reference dataset..."
    echo ""
    echo "  NOTE: This requires internet. If on a compute node without internet,"
    echo "  run this on the login node first:"
    echo "    ml apptainer/1.3.6"
    echo "    apptainer exec containers/nextclade_latest.sif nextclade dataset get \\"
    echo "        --name sars-cov-2 --output-dir nextclade_output/dataset_sars-cov-2"
    echo ""
    ${CMD} exec "${CONTAINER_SIF}" nextclade dataset get \
        --name sars-cov-2 --output-dir "${DATASET_DIR}"
fi

# --- Run ---
if [ ! -f "${FASTA_INPUT}" ]; then
    echo "ERROR: FASTA not found: ${FASTA_INPUT}"
    echo "Run 00_prefilter_gisaid.py first."
    exit 1
fi

echo "[3/3] Running Nextclade with ${JOBS} workers..."
echo "  Started: $(date)"

OUTPUT_TSV="${NEXTCLADE_DIR}/nextclade.tsv"

# Note: On Minerva, /sc/arion/, $HOME, and $PWD are auto-mounted
# into Apptainer containers, so no --bind needed.
${CMD} exec \
    "${CONTAINER_SIF}" \
    nextclade run \
        --input-dataset "${DATASET_DIR}" \
        --output-tsv "${OUTPUT_TSV}" \
        --jobs "${JOBS}" \
        "${FASTA_INPUT}"

echo "  Finished: $(date)"

if [ -f "${OUTPUT_TSV}" ]; then
    NROWS=$(( $(wc -l < "${OUTPUT_TSV}") - 1 ))
    echo "  ✓ Output: ${OUTPUT_TSV} (${NROWS} sequences)"
else
    echo "ERROR: Output not found."
    exit 1
fi
