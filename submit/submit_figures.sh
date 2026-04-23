#!/bin/bash
#BSUB -J figures
#BSUB -P acc_YOUR_PROJECT
#BSUB -q premium
#BSUB -n 4
#BSUB -W 01:00
#BSUB -R rusage[mem=32000]
#BSUB -R span[hosts=1]
#BSUB -o logs/figures_%J.log
#BSUB -e logs/figures_%J.log
#BSUB -L /bin/bash

# Regenerate all four manuscript figures + per-panel source data.
# Edit `#BSUB -P acc_YOUR_PROJECT` above to match your LSF account, then:
#   bsub < submit/submit_figures.sh
#
# Inputs:  processed_data/  (pre-populated — see REPRODUCIBILITY.md step 3)
# Outputs: figures/, source_data/
# Runtime: ~45-60 s on 4 cores / 32 GB.

PROJECT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
cd "${PROJECT_DIR}"
mkdir -p logs

# R environment — adapt module / env name to your cluster
ml anaconda3/2025.06
source activate av3env

echo "Figure generation started: $(date)"
echo "Host: $(hostname)"
echo "Project: ${PROJECT_DIR}"

Rscript scripts/03_generate_figures.R "${PROJECT_DIR}"

echo ""
echo "Figure generation finished: $(date)"
