#!/bin/bash
# =============================================================================
# End-to-end pipeline: GISAID raw → figures
# =============================================================================
#
# Chains the four stages of the analysis. Each stage can be skipped if its
# output is already present. Typical use: if you downloaded our pre-processed
# tables in processed_data/, you only need stage 4 (submit/submit_figures.sh).
#
# Usage:
#   bash scripts/run_pipeline.sh [--skip-prefilter] [--skip-nextclade] \
#                                [--skip-process]   [--only-figures] \
#                                [--jobs N]
#
# Stages:
#   1. Prefilter GISAID          scripts/00_prefilter_gisaid.py      (~1 h)
#   2. Nextclade alignment       scripts/01_run_nextclade.sh         (~3.5 h on 96 cores)
#   3. Process Nextclade output  scripts/02_process_nextclade.R      (~30 min)
#   4. Generate figures          scripts/03_generate_figures.R       (~1 min)
# =============================================================================

set -euo pipefail

SKIP_PREFILTER=0
SKIP_NEXTCLADE=0
SKIP_PROCESS=0
ONLY_FIGURES=0
N_JOBS=96

while [[ $# -gt 0 ]]; do
  case "$1" in
    --skip-prefilter) SKIP_PREFILTER=1 ;;
    --skip-nextclade) SKIP_NEXTCLADE=1 ;;
    --skip-process)   SKIP_PROCESS=1 ;;
    --only-figures)   ONLY_FIGURES=1; SKIP_PREFILTER=1; SKIP_NEXTCLADE=1; SKIP_PROCESS=1 ;;
    --jobs)           shift; N_JOBS="$1" ;;
    *) echo "unknown arg: $1" >&2; exit 1 ;;
  esac
  shift
done

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
cd "${PROJECT_DIR}"

say() { echo ""; echo "=========================================="; echo "$1"; echo "=========================================="; }

# --- Stage 1: Prefilter GISAID ---
if [[ $SKIP_PREFILTER -eq 0 ]]; then
  say "Stage 1/4: Prefilter GISAID"
  python3 "${SCRIPT_DIR}/00_prefilter_gisaid.py" \
    --metadata   raw_data/metadata.tsv \
    --sequences  raw_data/sequences.fasta \
    --exclude    reference_data/nextstrain_exclude.txt \
    --out-meta   filtered_data/filtered_metadata.tsv \
    --out-fasta  filtered_data/filtered_sequences.fasta
fi

# --- Stage 2: Nextclade ---
if [[ $SKIP_NEXTCLADE -eq 0 ]]; then
  say "Stage 2/4: Nextclade alignment"
  bash "${SCRIPT_DIR}/01_run_nextclade.sh" --jobs "${N_JOBS}"
fi

# --- Stage 3: Process Nextclade → processed_data/ ---
if [[ $SKIP_PROCESS -eq 0 ]]; then
  say "Stage 3/4: Process Nextclade output"
  Rscript "${SCRIPT_DIR}/02_process_nextclade.R" "${PROJECT_DIR}"
fi

# --- Stage 4: Figures + source data ---
say "Stage 4/4: Generate figures + source data"
Rscript "${SCRIPT_DIR}/03_generate_figures.R" "${PROJECT_DIR}"

say "Pipeline complete."
echo "Figures:    ${PROJECT_DIR}/figures/"
echo "Source data: ${PROJECT_DIR}/source_data/"
