# REPRODUCIBILITY.md

Two paths depending on what you have:

- **Path A** (fast, ~1 min): Regenerate the four final figures and their per-panel
  source CSVs from the pre-processed tables already in `processed_data/`.
- **Path B** (slow, ~4-8 h): Rebuild `processed_data/` from the original GISAID
  bulk download plus auxiliary sources.

If you just want to see or verify the figures, use Path A.

## 1. Software prerequisites

R + Python in one shot — [`environment.yml`](environment.yml):

```bash
conda env create -f environment.yml
conda activate mutation-tracking
```

| Tool | Version | Used for | How installed |
|---|---|---|---|
| R | 4.3.3 | Processing + figures | `environment.yml` |
| R packages | `data.table`, `tidyverse`, `patchwork`, `ggpubr`, `ggrepel` | — | `environment.yml` |
| Python | 3.10 | OWID fetch (scripts 00 and 08 use only stdlib) | `environment.yml` |
| Python packages | `pandas`, `owid-catalog` | — | `environment.yml` |
| Nextclade | v3 | Mutation calling | `apptainer pull docker://nextstrain/nextclade` |
| Apptainer | ≥ 1.3 | Container runtime | Cluster module or system package |
| LSF (optional) | — | HPC submit templates | Cluster-provided |

## 2. Dataset

- **Primary data:** GISAID EpiCoV bulk download (metadata + FASTA), accessed
  ~2026-02-21. GISAID's terms prohibit redistribution, so the raw files are not
  in this repository. Users with GISAID credentials can reproduce the download.
- **Secondary data:**
  - Our World in Data `owid-covid-data.csv` — fetched programmatically
    (`scripts/fetch_owid_data.py`, writes to `reference_data/`).
  - `reference_data/nejm_early_cases.csv` — daily symptom-onset counts for
    Dec 1, 2019 – Jan 21, 2020 from Li et al. (NEJM 2020) Figure 1. Included.
  - `reference_data/variant_emergence_dates.csv` — curated WHO/GISAID/Nextstrain
    emergence dates used for the mis-dating audit. Included.
  - `reference_data/pango_alias_key.json` — Pango alias mapping. Included.
  - `reference_data/nextstrain_exclude.txt` — Nextstrain problematic-sequences
    list. Included (update from https://github.com/nextstrain/ncov if older
    than a few weeks).

## 3. Path A — regenerate figures from processed tables

### 3a. On a workstation / login node

```bash
# If .csv.gz compression was applied, either decompress or let fread() handle it directly:
#   gunzip processed_data/substitutions_only/*.gz        (optional)

Rscript scripts/03_generate_figures.R .
```

Runtime: ~45-60 s on 4 cores / 32 GB. Outputs:

- `figures/figure1.{pdf,jpeg}`, `figure2.{pdf,jpeg}`, `figure3.{pdf,jpeg}`,
  `supplementary_figure1.{pdf,jpeg}`
- `source_data/fig{1..3,s1}*.csv` (18 files, one per panel / element)

### 3b. On LSF

```bash
cp submit/submit_figures.sh submit/submit_figures.local.sh
sed -i 's/acc_YOUR_PROJECT/acc_<your_account>/' submit/submit_figures.local.sh
bsub < submit/submit_figures.local.sh
```

## 4. Path B — rebuild from GISAID raw

### 4a. Populate inputs

```
raw_data/
├── metadata.tsv              # 14 GB, from GISAID bulk download
└── sequences.fasta           # 487 GB, from GISAID bulk download

reference_data/
├── owid-covid-data.csv       # python3 scripts/fetch_owid_data.py  (~450 MB)
├── nejm_early_cases.csv      # included
├── variant_emergence_dates.csv
├── pango_alias_key.json
├── nextstrain_exclude.txt
└── pipeline_config.txt
```

### 4b. Pull the Nextclade container (login node, needs internet)

```bash
mkdir -p containers
apptainer pull containers/nextclade_latest.sif docker://nextstrain/nextclade:latest
```

### 4c. Run the pipeline

```bash
# Full pipeline
bash scripts/run_pipeline.sh --jobs 96

# Or on LSF (edit -P YOUR_PROJECT in submit_pipeline.sh first)
bsub < submit/submit_pipeline.sh
```

Stages and wall-time expectations (on 96 cores / 192 GB):

| Stage | Script | Time |
|---|---|---|
| 1. Prefilter GISAID | `00_prefilter_gisaid.py` | ~1 h |
| 2. Nextclade alignment | `01_run_nextclade.sh` | ~3.5 h |
| 3. Process Nextclade → `processed_data/` | `02_process_nextclade.R` | ~30 min |
| 4. Figures + source data | `03_generate_figures.R` | ~1 min |

Each stage is idempotent; `--skip-prefilter`, `--skip-nextclade`,
`--skip-process`, or `--only-figures` resume from later stages.

### 4d. (Optional) audit mis-dated VOC/VOI sequences

```bash
python3 scripts/audit_variant_dates.py \
  --metadata         filtered_data/filtered_metadata.tsv \
  --emergence-dates  reference_data/variant_emergence_dates.csv \
  --alias-key        reference_data/pango_alias_key.json \
  --output           processed_data/audit_misdated_variants.csv
```

This produces the 11-row sub-lineage audit used to generate the 1,883-sequence
manual-exclusion list. Re-running is only necessary if the GISAID metadata has
been refreshed or the emergence-dates table has been updated.

## 5. Filter and cutoff settings (final)

| Figure | First-seen cutoff | Second-occurrence heuristic | Frequency bins | Notes |
|---|---|---|---|---|
| Fig 1 A-F | ≥ 0.01% pandemic-wide | — | [0.01%, 1%) / [1%, 10%) / ≥10% | 4-color scatter + pie inset |
| Fig 2 A, B | First seen by **2020-01-15** | — | Same as Fig 1 | NT (A) / AA (B) scatter + Spearman; y-axis `limits = c(NA, 100)` |
| Fig 2 C | First seen by 2020-01-15 | — | 3 bins: [0.01%, 0.1%) / [0.1%, 10%) / ≥10% | NT + AA median ± SD, connecting lines, Wilcoxon brackets |
| Fig 2 D | First seen by 2020-01-15 | **≤ 28 days** vs > 28 days | ≥ 1% only | Dodged bar (NT vs AA × ≤28d vs >28d) |
| Fig 3 A, B | Full pandemic | — | 4 bins [0.01,0.1,1,10,100]% | Viridis green-yellow |
| Fig 3 C | Full pandemic | — | 3 tiers: ≥0.01%, ≥0.1%, ≥1% | No reference line |
| Fig 3 D | Full pandemic | — | Spike vs rest (AA ≥ 0.1%) | No reference line |
| Supp Fig 1 A-D | — | — | — | Overview only (4 panels) |

## 6. Key numbers (sanity-check after re-running)

| Check | Expected value |
|---|---|
| Final retained sequences (`processed_data/nseqs.txt`) | 5,628,868 |
| Unique NT substitutions ≥ 0.01% | 10,147 |
| Unique AA substitutions ≥ 0.01% | 5,665 |
| Manual exclusions | 1,883 (1,846 wastewater/Vero + 26 misdated VOC + 11 misdated VOI) |
| Fig 2A (NT scatter) points | 101 |
| Fig 2B (AA scatter) points | 51 |
| Fig 2 Spearman ρ (NT) | −0.47, p = 6.2 × 10⁻⁷ |
| Fig 2 Spearman ρ (AA) | −0.62, p = 1.3 × 10⁻⁶ |

If all of these match after a fresh run, the pipeline is reproducing correctly.

## 7. Contact

For questions: akhil.kumar@mssm.edu (or open an issue on the repository).
