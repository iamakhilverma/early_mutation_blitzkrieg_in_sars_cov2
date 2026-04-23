## Methods

### Data Acquisition

We obtained 17,564,626 SARS-CoV-2 consensus genome sequences and associated metadata from the GISAID EpiCoV database (accessed February 2026). Global COVID-19 case counts and vaccination data were obtained from Our World in Data (OWID; accessed via the owid-catalog Python package). For the period preceding OWID coverage (December 1, 2019 to January 21, 2020), daily case counts were derived from symptom onset data reported in Figure 1 of Li et al. (NEJM, 2020), which documented the first 425 confirmed cases in Wuhan, China.

### Quality Control and Filtering

Sequences were subjected to stringent quality filters: (1) GISAID "Is high coverage?" flag set to True; (2) "Is low coverage?" flag set to False; (3) human host only; (4) complete collection date in YYYY-MM-DD format within the range December 1, 2019 to December 31, 2025; (5) sequence length >= 29,000 bp; (6) ambiguous base (N) content <= 5%; (7) not present on the Nextstrain exclusion list (10,837 known problematic sequences). These filters retained 5,732,294 sequences (32.6% of the total). The primary exclusion criterion was the high-coverage requirement, which removed 11,219,881 sequences (63.9%).

### Mutation Calling

Retained sequences were aligned to the Wuhan-Hu-1 reference genome (NC_045512.2) using Nextclade v3 (https://clades.nextstrain.org), which performs independent pairwise alignment of each sequence. The Nextclade dataset for SARS-CoV-2 was used to provide the reference sequence, gene map, and quality control parameters. Sequences receiving a Nextclade overall QC status of "bad" were excluded, retaining 5,635,025 sequences with "good" or "mediocre" QC status. We further removed 1,883 sequences identified on manual audit as wastewater samples, cell-passaged (Vero) isolates, or sequences with collection dates preceding their Pango lineage's earliest documented emergence (see "Audit of problematic sequences" below). The final analysis cohort comprised **5,633,141 sequences**.

From the Nextclade output, we extracted nucleotide substitutions, amino acid substitutions, and insertions/deletions. For each mutation, we computed the pandemic-wide frequency as the number of sequences carrying that mutation divided by the total number of sequences (5,633,141), expressed as a percentage. Mutations were retained for analysis if their pandemic-wide frequency exceeded 0.01% (approximately 564 or more sequences). This threshold yielded 10,135 unique nucleotide substitutions and 5,653 unique amino acid substitutions.

Nextclade reports amino acid mutations using ORF-level coordinates (ORF1a, ORF1b) rather than individual protein-level coordinates (NSP1-NSP16) used in GISAID annotations. For example, the well-known NSP12:P323L mutation (nucleotide C14408T) is reported as ORF1b:P314L by Nextclade. All non-ORF1ab proteins (Spike, N, M, E, and accessory proteins) use identical naming conventions in both systems.

### First and Second Occurrence Analysis

For each mutation, we identified the earliest collection date at which it was observed (first occurrence) and the second-earliest collection date (second occurrence). The interval between these dates, termed "days to second occurrence," was computed as the difference in calendar days. This metric serves as a proxy for the mutation's early detection and recurrence rate.

To assess whether early recurrence predicts eventual pandemic-wide success, we computed the Spearman rank correlation between days to second occurrence and pandemic-wide frequency for all mutations first observed by January 15, 2020 — the earliest informative window, when only ~100 nucleotide and ~50 amino acid substitutions at ≥0.01% frequency had been seen globally. We additionally performed pairwise Wilcoxon rank-sum tests across three frequency groups: [0.01%, 0.1%), [0.1%, 10%), and ≥10%. To quantify the enrichment of rapidly-recurring mutations among those that reached appreciable pandemic-wide frequency, we tabulated mutations first observed by January 15, 2020 that ultimately reached ≥1% pandemic-wide frequency by whether their second observation occurred within 28 days of the first (Figure 2D).

### Frequency Distribution and Selection Analysis

New mutations were assigned to time periods (monthly for December 2019 through June 2020, then semi-annual periods) and classified by pandemic-wide frequency bins. The proportion of mutations in each frequency bin was computed per time period to assess temporal trends in mutation discovery at different frequency thresholds.

The non-synonymous fraction was computed as the ratio of new amino acid substitutions to new nucleotide substitutions first observed in each time period, at three frequency thresholds (>= 0.01%, >= 0.1%, >= 1%). Under neutral evolution, approximately 75% of random single-nucleotide substitutions in coding regions are expected to be non-synonymous due to codon degeneracy. Deviations below this expectation indicate purifying selection; enrichment at higher frequency thresholds indicates positive selection.

The proportion of new amino acid mutations occurring in the Spike protein versus the remainder of the genome was computed per time period for mutations at >= 0.1% frequency. The null expectation was derived from the ratio of Spike amino acid positions (1,273) to total coding positions across all ORFs (9,744), yielding 13.1%. Enrichment above this baseline indicates preferential accumulation of Spike mutations, consistent with immune-driven positive selection.

### Infection Normalization

To account for the relationship between the number of infections and the number of mutations discovered, we computed the rate of new mutation discoveries per thousand documented infections for each time period. Daily COVID-19 case counts were obtained from OWID for the period January 22, 2020 onward. For the preceding period (December 1, 2019 to January 21, 2020), we used symptom onset counts from Li et al. (NEJM, 2020). Cumulative case counts were computed as a running sum, with the transition between data sources occurring at January 22, 2020 without discontinuity.

### Statistical Tests

Spearman rank correlations were used to assess monotonic associations between days to second occurrence and pandemic-wide frequency. Pairwise Wilcoxon rank-sum tests compared days to second occurrence across three frequency groups ([0.01%, 0.1%), [0.1%, 10%), ≥10%); brackets in Figure 2C mark significance at `***` p < 0.001, `**` p < 0.01, `*` p < 0.05, `ns` not significant. All statistical analyses were performed in R v4.3.3 using the data.table, ggpubr, and tidyverse packages.

### Audit of problematic sequences

To remove sequences with suspect lineage calls or non-natural origins, we performed an additional manual audit on top of Nextclade QC. A total of 1,883 sequences were excluded, comprising (i) 1,846 wastewater or Vero-cell-passaged isolates identified by free-text searches of GISAID metadata (keywords: "wastewater", "Vero", "cell"); (ii) 26 sequences whose Pango lineage corresponded to a WHO-designated Variant of Concern but whose collection date preceded that VOC's earliest documented community emergence (e.g., Alpha-labelled sequences predating October 2020, Omicron BA.1-labelled sequences predating November 2021); and (iii) 11 sequences flagged similarly for Variants of Interest or lineage sub-variants, using a curated table of 27 lineage-emergence dates derived from WHO, GISAID, and Nextstrain reports. The exclusion list (`exclude_additional.txt`) and per-sequence reasons (`excluded_sequences_detail.csv`) are provided as supplementary data.

### Verification Against Flagged Sequences

To further ensure data integrity, we cross-referenced our filtered dataset against the 42,419 GISAID sequences flagged as "under investigation" (accessed March 22, 2026; human host, collection date December 1, 2019 to January 1, 2026). Of these, 41,690 were present in our raw metadata download (17,564,626 sequences). None passed our quality control filters: zero flagged sequences were present in our filtered dataset (5,732,294 sequences) or final analysis cohort (5,633,141 sequences). All flagged sequences were excluded at the prefilter stage, primarily by the high-coverage requirement.

### Software and Reproducibility

All analyses were performed on the Minerva high-performance computing cluster at the Icahn School of Medicine at Mount Sinai. Nextclade was executed via Apptainer (v1.3.6) containers. R (v4.3.3) and Python (v3.10) scripts were run in a single Conda environment whose exact package set is pinned in the `environment.yml` file accompanying the code. Code and processed data are available at https://github.com/iamakhilverma/early_mutation_blitzkrieg_in_sars_cov2.

---
