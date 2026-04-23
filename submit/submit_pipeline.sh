#!/bin/bash
#BSUB -J pipeline
#BSUB -P acc_YOUR_PROJECT
#BSUB -q premium
#BSUB -n 96
#BSUB -W 24:00
#BSUB -R rusage[mem=192000]
#BSUB -R span[hosts=1]
#BSUB -o logs/pipeline_%J.log
#BSUB -e logs/pipeline_%J.log
#BSUB -L /bin/bash

# Full pipeline from GISAID raw data → figures.
# Only needed if you're re-running from scratch. If you just want the figures
# from our pre-processed tables, use submit/submit_figures.sh instead.
#
# Required inputs (see REPRODUCIBILITY.md):
#   raw_data/metadata.tsv            (14 GB, from GISAID download)
#   raw_data/sequences.fasta         (487 GB, from GISAID download)
#   reference_data/pipeline_config.txt
#   reference_data/nextstrain_exclude.txt
#   reference_data/variant_emergence_dates.csv
#   reference_data/pango_alias_key.json
#   reference_data/nejm_early_cases.csv
#   reference_data/owid-covid-data.csv   (via scripts/fetch_owid_data.py)
#
# Wall time: ~4-8 h depending on node count.

PROJECT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
cd "${PROJECT_DIR}"
mkdir -p logs filtered_data nextclade_output processed_data

ml anaconda3/2025.06
ml apptainer/1.3.6
source activate mutation-tracking

echo "Pipeline started: $(date)"
echo "Host: $(hostname)"
echo "Project: ${PROJECT_DIR}"

bash scripts/run_pipeline.sh --jobs 96 2>&1

echo ""
echo "Pipeline finished: $(date)"
