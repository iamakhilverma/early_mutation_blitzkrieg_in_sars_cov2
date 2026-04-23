# processed_data/ — reproduce the analysis without re-running Nextclade

Tables produced by `scripts/02_process_nextclade.R`. If you clone the repo and
run `scripts/03_generate_figures.R` directly, the figures are built from these
files (no GISAID access needed).

Strict filter is canonical: **5,628,868 sequences** (high-coverage flag, human
host, complete date, ≥29 kb, ≤5% N, Nextstrain exclude, Nextclade QC
good/mediocre, plus 1,883 manual exclusions — wastewater, cell-passaged, misdated).

## Mutations

| File | Rows | Description |
|---|---:|---|
| `substitutions_only/nt_first_occurrence.csv` | 10,147 | Every unique NT substitution ≥ 0.01%. Columns: `mutation_name`, `pos`, `ref`, `alt`, `collection_date` (first), `count`, `fraction`, `per_abundance`, `geq_1pct`, `geq_pt1pct`, `geq_pt01pct`, `mutation_type` |
| `substitutions_only/aa_first_occurrence.csv` | 49,995 | Every unique AA substitution (all frequencies). Use `geq_pt01pct == TRUE` for the 5,665 rows ≥ 0.01%. Adds `aa_label`, `region`, `is_spike` |
| `substitutions_only/nt_mutations_by_date.csv.gz` | ~7.5 M | One row per (mutation, collection_date) with daily counts — input for second-occurrence computation. Gzipped to fit the GitHub 100 MB file limit. `fread()` reads `.gz` transparently |
| `substitutions_only/aa_mutations_by_date.csv.gz` | ~5 M | Same for AA |
| `second_occurrence/nt_second_occurrence.csv` | 71,050 | Every unique NT mutation with `first_date`, `second_date`, `n_dates`, `gap_days`, `count`, `per_abundance`. Singletons (n_dates = 1) have blank `second_date`/`gap_days` |
| `second_occurrence/aa_second_occurrence.csv` | 49,995 | Same for AA (8,663 singletons) |

## Contextual data

| File | Source | Description |
|---|---|---|
| `world_infections.csv` | NEJM (Dec 19 – Jan 21, 2020) + OWID (Jan 22 onward) | Daily global new-case counts |
| `vaccination_data.csv` | OWID | Daily `pct_fully_vaccinated` (World) |
| `gisaid_daily_sequence_counts.csv` | GISAID | Daily `n_seqs`, `cum_seqs` — full GISAID |
| `gisaid_monthly_sequence_counts.csv` | GISAID | Same aggregated monthly |
| `sequence_counts_by_date.csv` | Our cohort | Daily retained-sequence counts in our 5,628,868-sequence cohort |

## Exclusions / audit trail

| File | Description |
|---|---|
| `exclude_additional.txt` | 1,883 GISAID accession IDs manually excluded |
| `excluded_sequences_detail.csv` | Per-sequence reason (wastewater / cell-passaged / misdated) |
| `audit_misdated_variants.csv` | 11 newly flagged misdated VOI / sub-lineage sequences |
| `audit_new_exclusions.txt` | Summary of the latest exclusion audit |

## Pipeline metadata

| File | Description |
|---|---|
| `nseqs.txt` | Final retained-sequence count: `5628868` |
| `processing_summary.txt` | Phase-1 processing log (substitutions counts per threshold) |

## Reproducing these files

See [REPRODUCIBILITY.md](../REPRODUCIBILITY.md). In short:

- `nt_first_occurrence.csv`, `aa_first_occurrence.csv`, `nt_mutations_by_date.csv`, `aa_mutations_by_date.csv`, `{nt,aa}_second_occurrence.csv`, `sequence_counts_by_date.csv`, `nseqs.txt`, `processing_summary.txt` are all produced by `scripts/02_process_nextclade.R`.
- `world_infections.csv`, `vaccination_data.csv`, `gisaid_daily/monthly_sequence_counts.csv` are produced by `scripts/fetch_owid_data.py` (OWID) plus the Nextclade prefilter summary (GISAID counts).
- `exclude_additional.txt` and `excluded_sequences_detail.csv` come from the wastewater/Vero keyword search plus `scripts/audit_variant_dates.py`.
