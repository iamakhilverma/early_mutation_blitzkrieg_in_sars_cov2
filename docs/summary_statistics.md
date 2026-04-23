# Summary Statistics — Final Manuscript

Strict filter is canonical for the manuscript. The lenient branch and the predictive-evaluation search grid were explored during revision but are not part of the submission; see `docs/PROJECT_HISTORY.md` and `docs/predictive_evaluation_memory.md` for record.

## Dataset

| Step | N sequences |
|---|---:|
| GISAID bulk download (Feb 2026) | 17,564,626 |
| After prefilter (high-coverage, human, complete date, ≥29 kb, ≤5% N, Nextstrain exclude) | 5,732,294 |
| After Nextclade QC (good / mediocre) | 5,635,025 |
| After manual exclusions (wastewater, Vero-passaged, misdated VOC/VOI) | **5,633,141** |

### Exclusions breakdown

| Reason | Count |
|---|---:|
| Wastewater or Vero cell-passaged | 1,846 |
| Misdated VOC (collection before variant emergence) | 26 |
| Misdated VOI / sub-lineage | 11 |
| **Total manual exclusions** | **1,883** |
| GISAID "under investigation" in our analysis cohort | **0** (of 42,419 flagged) |

## Mutation counts

| Threshold | NT | AA |
|---|---:|---:|
| Total unique (≥ 0.01%) | 10,135 | 5,653 |
| ≥ 0.1% | 1,911 | 1,054 |
| ≥ 1% | 246 | 166 |
| ≥ 10% | 61 | 48 |
| ≥ 50% | 32 | 27 |

## Figure 2 — Second-occurrence analysis (first-date cutoff: January 15, 2020)

### Scatter panels (A = NT, B = AA)

| Metric | NT (panel A) | AA (panel B) |
|---|---:|---:|
| N mutations plotted | 101 | 51 |
| Spearman ρ (days to 2nd vs per_abundance) | **−0.47** | **−0.62** |
| p-value | 6.2 × 10⁻⁷ | 1.3 × 10⁻⁶ |

### Panel C — Median days to second occurrence ± 1 SD, per frequency bin

| Frequency bin | NT median (SD), n | AA median (SD), n |
|---|---|---|
| [0.01%, 0.1%) | 68 d (49.6), n = 37 | 69 d (32.0), n = 14 |
| [0.1%, 10%) | 30 d (24.9), n = 56 | 26.5 d (24.0), n = 32 |
| ≥ 10% | 9 d (25.1), n = 8 | 21 d (29.3), n = 5 |

### Panel C — Pairwise Wilcoxon rank-sum tests

| Type | Comparison | p | Significance |
|---|---|---:|---|
| NT | [0.01%, 0.1%) vs [0.1%, 10%) | 2.8 × 10⁻⁴ | *** |
| NT | [0.1%, 10%) vs ≥ 10% | 5.2 × 10⁻² | ns |
| NT | [0.01%, 0.1%) vs ≥ 10% | 3.1 × 10⁻³ | ** |
| AA | [0.01%, 0.1%) vs [0.1%, 10%) | 5.9 × 10⁻⁴ | *** |
| AA | [0.1%, 10%) vs ≥ 10% | 4.1 × 10⁻¹ | ns |
| AA | [0.01%, 0.1%) vs ≥ 10% | 2.9 × 10⁻² | * |

### Panel D — Mutations reaching ≥ 1% pandemic-wide frequency (Jan 15 cutoff), by days-to-second-occurrence bucket

| Type | ≤ 28 days | > 28 days |
|---|---|---|
| NT | 22 (84.6%) | 4 (15.4%) |
| AA | 15 (78.9%) | 4 (21.1%) |

## Figure inventory (final)

| File | Description |
|---|---|
| `figures/figure1.{pdf,jpeg}` | Timeline of first occurrence; 6 panels (A–F) |
| `figures/figure2.{pdf,jpeg}` | Second-occurrence analysis; 4 panels (A–D); Jan 15 cutoff for all panels; panel D ≥1%, ≤28 d heuristic |
| `figures/figure3.{pdf,jpeg}` | Frequency distribution and non-synonymous enrichment; 4 panels (A–D), viridis green-yellow palette, no reference lines |
| `figures/supplementary_figure1.{pdf,jpeg}` | Overview: mutations/infections/sequencing/vaccination; 4 panels (A–D) |

## Source-data and supplementary-data inventory (final)

| Path | Files | Purpose |
|---|---:|---|
| `communication/source_data_final/` | 18 CSVs + README | Exact points/lines plotted (one per panel or bracket set) |
| `communication/supplementary_data_final/` | 21 files + README | Processed tables sufficient to reproduce the analysis without rerunning Nextclade |

## Tables for the journal / review

| Path | File | Rows | Description |
|---|---|---:|---|
| `communication/tables_for_journal/` | `strict_excluded_sequences.csv` | 1,883 | Per-sequence exclusion reasons |
| `communication/tables_for_journal/` | `strict_misdated_variants_audit.csv` | 11 | Sub-lineage audit flags |
| `communication/tables_for_review/` | `strict_ppv_summary.csv` | — | Predictive-evaluation search-grid summary (for reviewer reference; not in submission) |
| `communication/tables_for_review/` | `lenient_ppv_summary.csv` | — | Same on lenient branch (for reviewer reference) |
| `communication/tables_for_review/` | `sequence_count_summary.csv` | — | Sequence counts by date cutoff |
