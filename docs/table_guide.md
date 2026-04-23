# Table Guide

Every data file shipped with the manuscript, what figure/table it backs, and where it lives.

## Directory layout

| Directory | Purpose | Audience |
|---|---|---|
| `communication/source_data_final/` | Exact points/lines/bars plotted in each figure (one CSV per panel) | Journal submission |
| `communication/supplementary_data_final/` | Upstream processed tables sufficient to reproduce the full analysis | Journal submission |
| `communication/tables_for_journal/` | Curated tables on data-quality decisions (exclusions, audits) | Journal submission |
| `communication/tables_for_review/` | Digested summaries used during revision for teammate review | Record — not in submission |

No file duplication across directories. `source_data_final/` is small (≈1 MB) and matches plotted numbers exactly; `supplementary_data_final/` is ~420 MB including the per-mutation-by-date tables used upstream.

---

## `source_data_final/` — exactly what the plots show

### Figure 1

| File | Panel | Rows | Key columns |
|---|---|---:|---|
| `fig1ab_daily_first_occurrences.csv` | 1A, 1B | 2,223 | `date`, `nt_new_geq_pt01pct`, `aa_new_geq_pt01pct` |
| `fig1_pie_inset_data.csv` | 1A/1B insets | — | `TimelineLabel`, `name`, `bins`, `countMutations` |
| `fig1cd_monthly_rates.csv` | 1C, 1D | — | `month_date`, `nt_count`, `aa_count`, `newCases`, `nt_per_1000`, `aa_per_1000`, `flag` |
| `fig1e_nt_scatter.csv` | 1E | 10,147 | `mutation_name`, `pos`, `ref`, `alt`, `collection_date`, `count`, `per_abundance` |
| `fig1f_aa_scatter.csv` | 1F | 5,665 | `aa_label`, `region`, `mutation_name`, `pos`, `ref`, `alt`, `collection_date`, `count`, `per_abundance` |

### Figure 2 — Jan 15, 2020 cutoff for every panel; D also filtered to ≥ 1% and split at 28 d

| File | Panel | Rows | Key columns |
|---|---|---:|---|
| `fig2a_nt_scatter.csv` | 2A | 101 | `mutation_name`, `first_date`, `second_date`, `days_to_second`, `per_abundance`, `count`, `abundance_bin` |
| `fig2b_aa_scatter.csv` | 2B | 51 | `aa_label`, `first_date`, `second_date`, `days_to_second`, `per_abundance`, `count`, `abundance_bin` |
| `fig2c_median_sd.csv` | 2C | 6 | `type`, `freq_group`, `n`, `median_days`, `sd_days` |
| `fig2c_wilcoxon_brackets.csv` | 2C | 6 | `type`, `grp1`, `grp2`, `p`, `stars` |
| `fig2d_bar.csv` | 2D | 4 | `type`, `days_group`, `n`, `pct` |

### Figure 3

| File | Panel | Description |
|---|---|---|
| `fig3a_nt_freq_by_period.csv` | 3A | NT new-mutations per period × frequency bin (with `pct`) |
| `fig3b_aa_freq_by_period.csv` | 3B | AA new-mutations per period × frequency bin (with `pct`) |
| `fig3c_nonsyn_by_tier.csv` | 3C | Per period × tier: `nt_new`, `aa_new`, `pct_nonsyn`, `tier` |
| `fig3d_spike_vs_rest.csv` | 3D | Per period × `gene_group`: `n`, `pct` (AA ≥ 0.1% only) |

### Supplementary Figure 1

| File | Panel | Description |
|---|---|---|
| `suppfig1a_first4m_vs_rest.csv` | S1A | Type (NT/AA) × period (Dec-19–Mar-20 vs rest), `n`, `pct` |
| `suppfig1b_infections_by_period.csv` | S1B | Same two periods, `total_cases` (log-scale bar) |
| `suppfig1c_daily_seqs_vs_cases.csv` | S1C | Daily: `collection_date`, `cum_seqs`, `cum_cases`, `cum_seqs_M`, `cum_cases_M` |
| `suppfig1d_vaccination.csv` | S1D | OWID daily `pct_fully_vaccinated` from 2021-01-01 onward |

---

## `supplementary_data_final/` — reproduce the analysis

| File | Rows | Purpose |
|---|---:|---|
| `nt_first_occurrence.csv` | 10,147 | Every unique NT substitution ≥ 0.01% with first date, count, frequency |
| `aa_first_occurrence.csv` | 49,995 | Every unique AA substitution (all frequencies). Filter `geq_pt01pct == TRUE` for 5,665 rows matching the ≥0.01% threshold |
| `nt_mutations_by_date.csv` | ~7.5 M | Per-mutation-per-date observation counts (NT) — input for second-occurrence computation |
| `aa_mutations_by_date.csv` | ~5 M | Same for AA |
| `nt_second_occurrence_full.csv` | 71,050 | Every unique NT mutation with first/second date (includes 7,739 singletons with blank `second_date` / `gap_days`) |
| `aa_second_occurrence_full.csv` | 49,995 | Every unique AA mutation with first/second date (8,663 singletons) |
| `nt_second_occurrence_jan15_plotted.csv` | 101 | The NT subset first-seen ≤ Jan 15 with a second occurrence — matches `fig2a_nt_scatter.csv` |
| `aa_second_occurrence_jan15_plotted.csv` | 51 | Analogous for AA |
| `world_infections.csv` | 2,223 | Daily global new-case counts (NEJM Dec '19 – Jan 21, 2020 + OWID from Jan 22 onward) |
| `vaccination_data.csv` | — | OWID daily `pct_fully_vaccinated` (World) |
| `gisaid_daily_sequence_counts.csv` | ~2,170 | Daily GISAID sequence deposits (complete set before our filter) |
| `gisaid_monthly_sequence_counts.csv` | — | Monthly aggregate of the above |
| `sequence_counts_by_date.csv` | — | Daily retained-sequence counts in our 5,628,868-sequence cohort |
| `exclude_additional.txt` | 1,883 | Accession IDs of manual exclusions |
| `excluded_sequences_detail.csv` | 1,883 | Per-sequence reason (wastewater / cell-passaged / misdated) |
| `audit_misdated_variants.csv` | 11 | Per-sequence flags for mis-dated VOI / sub-lineage |
| `audit_new_exclusions.txt` | — | Summary of the latest exclusion audit |
| `variant_emergence_dates.csv` | 27 | Curated WHO/GISAID/Nextstrain variant emergence dates |
| `nejm_early_cases.csv` | — | Daily symptom-onset counts (Li et al., NEJM 2020 Fig 1) |
| `nseqs.txt` | 1 | Final retained-sequence count: 5628868 |
| `processing_summary.txt` | — | Brief Phase 1 processing log |

---

## `tables_for_journal/`

Curated tables on data-quality decisions (not directly plotted).

| File | Rows | Description |
|---|---:|---|
| `strict_excluded_sequences.csv` | 1,883 | Every excluded sequence: accession, reason (wastewater/cell-passaged, Alpha/Beta/Gamma/Delta/Omicron before emergence), collection date, lineage |
| `strict_misdated_variants_audit.csv` | 11 | Sequences with collection dates preceding their Pango lineage's emergence — comprehensive audit of 27 lineage prefixes against the curated emergence table |

---

## `tables_for_review/` (not in submission — record of revision)

Digested summaries from the revision cycle. Kept for our records and for reviewer questions; not part of the journal package.

| File | Description |
|---|---|
| `strict_ppv_summary.csv` | Predictive-evaluation PPV at 4 date cutoffs × 2 gap thresholds × 3 freq thresholds × 2 types (search grid) |
| `lenient_ppv_summary.csv` | Same on the lenient branch (abandoned for final) |
| `sequence_count_summary.csv` | Sequence counts by date cutoff (pre- and post-exclusion) |

The search grid itself (methods in `predictive_evaluation_methods.md`, results in `docs/predictive_evaluation_memory.md`) was dropped from the final submission. The summary CSVs are retained here for reviewer reference only.

---

## Quick reference — which file answers what?

| Question | File |
|---|---|
| What points are on the Fig 2A NT scatter? | `source_data_final/fig2a_nt_scatter.csv` (101 rows) |
| What are the medians and SDs behind the Fig 2C line? | `source_data_final/fig2c_median_sd.csv` |
| What's the p-value for the `***` bracket in Fig 2C? | `source_data_final/fig2c_wilcoxon_brackets.csv` |
| How many NT/AA mutations ≥ 1% reached <=28d after 2nd obs? | `source_data_final/fig2d_bar.csv` |
| Full list of NT mutations ≥ 0.01% with first date | `supplementary_data_final/nt_first_occurrence.csv` |
| Per-mutation-per-date observation counts | `supplementary_data_final/{nt,aa}_mutations_by_date.csv` |
| Daily infection counts (Fig 1C/D denominator, Supp Fig 1B) | `supplementary_data_final/world_infections.csv` |
| Vaccination timeline for Supp Fig 1D | `source_data_final/suppfig1d_vaccination.csv` |
| What sequences were excluded and why? | `tables_for_journal/strict_excluded_sequences.csv` |
| How many sequences after each filter step? | `summary_statistics.md` (Dataset section) |
