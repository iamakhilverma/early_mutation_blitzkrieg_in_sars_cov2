# Early mutation dynamics of SARS-CoV-2

Code, processed data, and figures for the manuscript on early-pandemic SARS-CoV-2
mutation dynamics (5,633,141 high-coverage GISAID sequences, December 2019 –
December 2025).

**Central finding.** Of nucleotide and amino acid substitutions first detected by
January 15, 2020, the speed of second observation is a strong negative predictor
of eventual pandemic-wide frequency (Spearman ρ = −0.47 for NT, −0.62 for AA;
p < 10⁻⁶). Among those early mutations that ultimately reached ≥1% frequency,
84.6% (NT) and 78.9% (AA) recurred within 28 days of first detection.

## Repository layout

```
.
├── README.md                      This file
├── LICENSE                        MIT
├── REPRODUCIBILITY.md             How to reproduce the figures from scratch or from processed tables
├── .gitignore
│
├── scripts/                       End-to-end analysis pipeline
│   ├── 00_prefilter_gisaid.py     (Py) GISAID raw → filtered FASTA + metadata
│   ├── 01_run_nextclade.sh        (bash) Nextclade alignment via Apptainer
│   ├── 02_process_nextclade.R     (R)  Nextclade output → processed_data/
│   ├── 03_generate_figures.R      (R)  processed_data/ → figures/ + source_data/
│   ├── audit_variant_dates.py     (Py) flags mis-dated VOC/VOI sequences for exclusion
│   ├── fetch_owid_data.py         (Py) OWID case & vaccination CSV download
│   └── run_pipeline.sh            (bash) orchestrator with --skip-* flags
│
├── submit/                        LSF submit templates (edit account, adapt for your cluster)
│   ├── submit_figures.sh          (stage 4 only — ~1 min; uses pre-populated processed_data/)
│   └── submit_pipeline.sh         (full stages 1-4 — ~4-8 h; needs GISAID raw)
│
├── reference_data/                Small curated inputs (< 1 MB)
│   ├── variant_emergence_dates.csv   27 WHO/GISAID/Nextstrain lineage emergence dates
│   ├── nejm_early_cases.csv          Daily case counts for Dec 2019 - Jan 21, 2020 (Li et al.)
│   ├── pango_alias_key.json          Pango lineage alias mapping
│   ├── nextstrain_exclude.txt        Nextstrain "problematic sequences" list
│   └── pipeline_config.txt           START_DATE/END_DATE for this run
│
├── processed_data/                Tables produced by 02_process_nextclade.R (~120 MB)
│   ├── README.md                  Per-file description
│   ├── substitutions_only/        First-occurrence + by-date (heaviest: ~50 MB gzipped each)
│   ├── second_occurrence/         Full + plot-subset second-occurrence tables
│   ├── world_infections.csv       Daily global case counts (NEJM + OWID)
│   ├── vaccination_data.csv       OWID "fully vaccinated" %
│   ├── gisaid_daily_sequence_counts.csv
│   ├── sequence_counts_by_date.csv
│   ├── exclude_additional.txt     1,883 manually-excluded accessions
│   ├── excluded_sequences_detail.csv
│   ├── audit_misdated_variants.csv
│   ├── nseqs.txt                  5633141
│   └── processing_summary.txt
│
├── source_data/                   Exact data plotted in each figure panel (≈1 MB)
│   ├── README.md
│   ├── fig1ab_daily_first_occurrences.csv
│   ├── fig1_pie_inset_data.csv
│   ├── fig1cd_monthly_rates.csv
│   ├── fig1e_nt_scatter.csv
│   ├── fig1f_aa_scatter.csv
│   ├── fig2a_nt_scatter.csv       101 points
│   ├── fig2b_aa_scatter.csv        51 points
│   ├── fig2c_median_sd.csv          6 rows (3 bins × NT/AA)
│   ├── fig2c_wilcoxon_brackets.csv  6 rows (3 pairs × NT/AA)
│   ├── fig2d_bar.csv                4 rows
│   ├── fig3a_nt_freq_by_period.csv
│   ├── fig3b_aa_freq_by_period.csv
│   ├── fig3c_nonsyn_by_tier.csv
│   ├── fig3d_spike_vs_rest.csv
│   ├── suppfig1a_first4m_vs_rest.csv
│   ├── suppfig1b_infections_by_period.csv
│   ├── suppfig1c_daily_seqs_vs_cases.csv
│   └── suppfig1d_vaccination.csv
│
├── figures/                       Final figures — one PDF and one JPEG per figure
│   ├── figure1.{pdf,jpeg}
│   ├── figure2.{pdf,jpeg}
│   ├── figure3.{pdf,jpeg}
│   └── supplementary_figure1.{pdf,jpeg}
│
└── docs/                          Manuscript-supporting documents
    ├── methods.md                 Camera-ready Methods section
    ├── figure_legends.md          Camera-ready figure legends
    ├── summary_statistics.md      All numbers that appear in text and tables
    └── table_guide.md             Per-CSV descriptions, key columns, quick reference
```

## Quickstart

### Option A — regenerate figures from our processed tables (fast, ~1 min)

```bash
git clone <this-repo-url>
cd <repo>

# Decompress the two heavy per-mutation-by-date tables (optional; fread reads .gz directly)
gunzip processed_data/substitutions_only/*.gz

# Run
Rscript scripts/03_generate_figures.R .
# … or on LSF:  bsub < submit/submit_figures.sh
```

Outputs land in `figures/` (PDFs + JPEGs) and `source_data/` (CSVs).

### Option B — rebuild everything from GISAID raw

1. Get GISAID credentials and download the bulk `metadata.tsv` and `sequences.fasta` into `raw_data/`.
2. Fetch OWID data: `python3 scripts/fetch_owid_data.py` (writes `reference_data/owid-covid-data.csv`).
3. Pull the Nextclade Apptainer container (see `scripts/01_run_nextclade.sh` for the `apptainer pull` command).
4. Run the full pipeline: `bash scripts/run_pipeline.sh --jobs 96` (or `bsub < submit/submit_pipeline.sh`).

Full recipe with HPC details, software versions, and expected outputs is in [REPRODUCIBILITY.md](REPRODUCIBILITY.md).

## Software

- **R ≥ 4.3.3** with packages: `data.table`, `stringr`, `tidyverse`, `ggplot2`, `patchwork`, `ggpubr`, `ggrepel`
- **Python ≥ 3.10** with `pandas`, `biopython`, `pyarrow` (for prefilter), `owid-catalog` (for OWID fetch)
- **Nextclade v3** (via Apptainer image; `docker://nextstrain/nextclade`)
- **LSF** submit templates are provided; adapt `-P acc_YOUR_PROJECT` and module names for your cluster.

## Data

- GISAID (EpiCoV) bulk download, accessed February 2026. GISAID data is not
  redistributed here — users must request access through gisaid.org. The
  `processed_data/` tables in this repo are derived summary counts and are
  shareable.
- Daily case counts: Our World in Data (OWID) from 2020-01-22 onward, plus
  Li et al. (NEJM 2020) for the preceding four weeks.
- Variant emergence dates curated from WHO, GISAID, and Nextstrain.

## Citation

If you use this code or processed data, please cite the manuscript (to be added
once accepted) and acknowledge the originating labs of the GISAID sequences as
described at https://gisaid.org.

## License

MIT — see [LICENSE](LICENSE).
