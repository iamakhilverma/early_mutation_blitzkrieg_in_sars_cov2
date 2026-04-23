#!/usr/bin/env Rscript
#' =============================================================================
#' Phase 1: Process Nextclade Output into Analysis-Ready Tables
#' =============================================================================
#'
#' Author: Akhil Kumar
#'
#' Produces only the tables needed for Figures 1, 2, 3, and Supplementary Figure
#' 1 in the final manuscript. Deletion/insertion calls and per-lineage
#' aggregations (which the previous pipeline wrote for Phase 2 analyses that are
#' not part of the final paper) are intentionally omitted.
#'
#' Outputs (in processed_data/):
#'   substitutions_only/
#'     nt_first_occurrence.csv      one row per unique NT substitution (≥0.01%)
#'     aa_first_occurrence.csv      one row per unique AA substitution (all freqs; flag `geq_pt01pct`)
#'     nt_mutations_by_date.csv     one row per (mutation, collection_date) — NT
#'     aa_mutations_by_date.csv     one row per (mutation, collection_date) — AA
#'   second_occurrence/
#'     nt_second_occurrence.csv     full NT second-occurrence table (incl. singletons with blank second_date)
#'     aa_second_occurrence.csv     full AA second-occurrence table
#'   sequence_counts_by_date.csv
#'   world_infections.csv           (if OWID CSV available)
#'   vaccination_data.csv           (if OWID CSV available)
#'   nseqs.txt
#'   processing_summary.txt
#'
#' Usage:  Rscript scripts/02_process_nextclade.R [project_dir]
#' Memory: ~30-50 GB (with ~5.7 M sequences)
#' =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
})

msg <- function(...) message(sprintf(...))

# =============================================================================
# Configuration
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
PROJECT_DIR <- if (length(args) >= 1) args[1] else "."

NEXTCLADE_TSV <- file.path(PROJECT_DIR, "nextclade_output", "nextclade.tsv")
METADATA_TSV  <- file.path(PROJECT_DIR, "filtered_data", "filtered_metadata.tsv")
CONFIG_FILE   <- file.path(PROJECT_DIR, "filtered_data", "pipeline_config.txt")
OWID_FILE     <- file.path(PROJECT_DIR, "reference_data", "owid-covid-data.csv")
if (!file.exists(OWID_FILE)) OWID_FILE <- file.path(PROJECT_DIR, "raw_data", "owid-covid-data.csv")
OUT_DIR       <- file.path(PROJECT_DIR, "processed_data")

SUBS_DIR    <- file.path(OUT_DIR, "substitutions_only")
SECOND_DIR  <- file.path(OUT_DIR, "second_occurrence")
for (d in c(OUT_DIR, SUBS_DIR, SECOND_DIR)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

START_DATE <- as.Date("2019-12-01")
END_DATE   <- as.Date("2025-12-31")
if (file.exists(CONFIG_FILE)) {
  config <- readLines(CONFIG_FILE)
  for (line in config) {
    parts <- strsplit(line, "=")[[1]]
    if (length(parts) == 2) {
      if (parts[1] == "START_DATE") START_DATE <- as.Date(parts[2])
      if (parts[1] == "END_DATE")   END_DATE   <- as.Date(parts[2])
    }
  }
}
MIN_FREQ <- 0.0001

gene_map <- c(
  "S" = "Spike", "N" = "N", "M" = "M", "E" = "E",
  "ORF1a" = "ORF1a", "ORF1b" = "ORF1b",
  "ORF3a" = "NS3", "ORF6" = "NS6", "ORF7a" = "NS7a",
  "ORF7b" = "NS7b", "ORF8" = "NS8", "ORF9b" = "ORF9b"
)

msg(paste(rep("=", 70), collapse = ""))
msg("Phase 1: Process Nextclade (substitutions + second-occurrence only)")
msg(paste(rep("=", 70), collapse = ""))
msg("Project:  %s", PROJECT_DIR)
msg("Window:   %s to %s", START_DATE, END_DATE)

# =============================================================================
# 1. Read Nextclade TSV (selected columns only)
# =============================================================================

msg("\n[1/7] Reading Nextclade TSV: %s", NEXTCLADE_TSV)
nc <- fread(NEXTCLADE_TSV, sep = "\t", header = TRUE, quote = "",
            select = c("seqName", "substitutions", "aaSubstitutions",
                       "qc.overallStatus", "Nextclade_pango"))
nc_total_rows <- nrow(nc)
msg("  Rows: %s, Cols: %d", format(nc_total_rows, big.mark = ","), ncol(nc))
nc[, AccessionID := str_extract(seqName, "EPI_ISL_\\d+")]

# =============================================================================
# 2. Read pre-filtered metadata
# =============================================================================

msg("\n[2/7] Reading pre-filtered metadata: %s", METADATA_TSV)
meta <- fread(METADATA_TSV)
msg("  Rows: %s", format(nrow(meta), big.mark = ","))
meta[, collection_date := as.Date(collection_date)]

# =============================================================================
# 3. Join + QC filter + manual exclusions
# =============================================================================

msg("\n[3/7] Joining with metadata + QC filter")
nc2 <- merge(nc, meta[, .(accession_id, collection_date, lineage, location)],
             by.x = "AccessionID", by.y = "accession_id")
msg("  After join: %s", format(nrow(nc2), big.mark = ","))
rm(nc, meta); gc()

qc_counts <- table(nc2$qc.overallStatus, useNA = "always")
msg("  QC distribution:")
for (nm in names(qc_counts)) {
  msg("    %s: %s", ifelse(is.na(nm), "NA", nm), format(qc_counts[nm], big.mark = ","))
}
nc2 <- nc2[qc.overallStatus %in% c("good", "mediocre")]
nc2 <- nc2[!is.na(AccessionID) & !is.na(collection_date)]
nc2 <- nc2[collection_date >= START_DATE & collection_date <= END_DATE]
nc2[is.na(lineage) | lineage == "", lineage := "unassigned"]

# Manual exclusion list (wastewater, Vero-cell-passaged, misdated VOC/VOI)
exclude_file <- file.path(OUT_DIR, "exclude_additional.txt")
if (file.exists(exclude_file)) {
  exclude_ids <- readLines(exclude_file)
  n_before <- nrow(nc2)
  nc2 <- nc2[!AccessionID %in% exclude_ids]
  n_removed <- n_before - nrow(nc2)
  msg("  Excluded additional problematic sequences: %d (wastewater, cell-passaged, misdated)", n_removed)
}

# De-duplicate by AccessionID (Nextclade occasionally emits >1 row per sequence;
# metadata is already unique). Keep first. This ensures every downstream count
# — NSEQS, mutation frequencies, sequence_list — treats each sequence exactly once.
n_before_dedup <- nrow(nc2)
nc2 <- unique(nc2, by = "AccessionID")
n_dedup_removed <- n_before_dedup - nrow(nc2)
if (n_dedup_removed > 0) {
  msg("  Removed %s duplicate-AccessionID rows (Nextclade-side dups)",
      format(n_dedup_removed, big.mark = ","))
}

NSEQS <- nrow(nc2)
msg("  After QC + exclusions + dedup: %s", format(NSEQS, big.mark = ","))
msg("\n  >>> TOTAL RETAINED: %s <<<", format(NSEQS, big.mark = ","))

# =============================================================================
# 4. Save sequence counts
# =============================================================================

msg("\n[4/7] Saving sequence counts")
seqs_by_date <- unique(nc2[, .(AccessionID, collection_date)], by = "AccessionID")[
  , .(n_sequences = .N), by = collection_date]
setorder(seqs_by_date, collection_date)
seqs_by_date[, cumulative := cumsum(n_sequences)]
fwrite(seqs_by_date, file.path(OUT_DIR, "sequence_counts_by_date.csv"))
msg("  Dates: %d", nrow(seqs_by_date))
rm(seqs_by_date); gc()

# =============================================================================
# 5. NT substitutions
# =============================================================================

msg("\n[5/7] Processing NT substitutions")
nt_data <- nc2[substitutions != "" & !is.na(substitutions),
               .(AccessionID, collection_date, lineage, substitutions)]
nt_long <- nt_data[, .(mutation_name = unlist(strsplit(substitutions, ","))),
                   by = .(AccessionID, collection_date)]
rm(nt_data); gc()
nt_long[, `:=`(ref = str_sub(mutation_name, 1, 1),
               pos = as.integer(str_extract(mutation_name, "\\d+")),
               alt = str_sub(mutation_name, -1, -1),
               mutation_type = "substitution")]
nt_npairs <- nrow(nt_long)
msg("  Pairs: %s", format(nt_npairs, big.mark = ","))

# By date (input for second occurrence + Figure 2 timing analysis)
nt_by_date <- nt_long[, .(count = .N), by = .(collection_date, mutation_name, pos, ref, alt)]
fwrite(nt_by_date, file.path(SUBS_DIR, "nt_mutations_by_date.csv"))

# Global + first occurrence
nt_global <- nt_long[, .(count = .N), by = .(mutation_name, pos, ref, alt)]
nt_global[, `:=`(fraction = count / NSEQS, per_abundance = (count / NSEQS) * 100)]
nt_first <- nt_long[, .(collection_date = min(collection_date)), by = .(mutation_name, pos, ref, alt)]
nt_first <- merge(nt_first, nt_global[, .(mutation_name, count, fraction, per_abundance)], by = "mutation_name")
setorder(nt_first, collection_date)
nt_subs_out <- nt_first[per_abundance >= MIN_FREQ * 100]
nt_subs_out[, `:=`(geq_1pct = per_abundance >= 1, geq_pt1pct = per_abundance >= 0.1,
                    geq_pt01pct = per_abundance >= 0.01, mutation_type = "substitution")]
fwrite(nt_subs_out, file.path(SUBS_DIR, "nt_first_occurrence.csv"))

nt_counts <- list(pairs = nt_npairs, total = nrow(nt_global),
                  geq01 = nrow(nt_subs_out), geq1 = sum(nt_subs_out$geq_pt1pct),
                  geq10 = sum(nt_subs_out$geq_1pct))
msg("  >=0.01%%: %d | >=0.1%%: %d | >=1%%: %d", nt_counts$geq01, nt_counts$geq1, nt_counts$geq10)

rm(nt_long, nt_by_date, nt_global, nt_first, nt_subs_out); gc()
msg("  Memory freed")

# =============================================================================
# 6. AA substitutions
# =============================================================================

msg("\n[6/7] Processing AA substitutions")
aa_data <- nc2[aaSubstitutions != "" & !is.na(aaSubstitutions),
               .(AccessionID, collection_date, lineage, aaSubstitutions)]
aa_long <- aa_data[, .(aa_raw = unlist(strsplit(aaSubstitutions, ","))),
                   by = .(AccessionID, collection_date)]
rm(aa_data); gc()
aa_long[, `:=`(gene_short = str_split_fixed(aa_raw, ":", 2)[, 1],
               mutation_name = str_split_fixed(aa_raw, ":", 2)[, 2])]
aa_long[, region := fifelse(gene_short %in% names(gene_map), gene_map[gene_short], gene_short)]
aa_long[, aa_label := paste0(region, "_", mutation_name)]
aa_long[, `:=`(ref = sub("^([A-Za-z*]+)(\\d+)([A-Za-z*]+)$", "\\1", mutation_name),
               pos = as.integer(str_extract(mutation_name, "\\d+")),
               alt = sub("^([A-Za-z*]+)(\\d+)([A-Za-z*]+)$", "\\3", mutation_name))]
aa_long <- aa_long[!is.na(pos) & ref != "" & alt != ""]
aa_long <- aa_long[ref != "ins" & ref != "del" & alt != "ins" & alt != "del"]
aa_long[, mutation_type := "substitution"]
aa_npairs <- nrow(aa_long)
msg("  Pairs: %s", format(aa_npairs, big.mark = ","))

# By date
aa_by_date <- aa_long[, .(count = .N), by = .(collection_date, aa_label, region, mutation_name, pos, ref, alt)]
fwrite(aa_by_date, file.path(SUBS_DIR, "aa_mutations_by_date.csv"))

# Global + first occurrence
aa_global <- aa_long[, .(count = .N), by = .(aa_label, region, mutation_name, pos, ref, alt)]
aa_global[, `:=`(fraction = count / NSEQS, per_abundance = (count / NSEQS) * 100)]
aa_first <- aa_long[, .(collection_date = min(collection_date)),
                    by = .(aa_label, region, mutation_name, pos, ref, alt)]
aa_first <- merge(aa_first, aa_global[, .(aa_label, count, fraction, per_abundance)], by = "aa_label")
setorder(aa_first, collection_date)
aa_first[, `:=`(geq_1pct = per_abundance >= 1, geq_pt1pct = per_abundance >= 0.1,
                geq_pt01pct = per_abundance >= 0.01, is_spike = (region == "Spike"),
                mutation_type = "substitution")]
fwrite(aa_first, file.path(SUBS_DIR, "aa_first_occurrence.csv"))

aa_counts <- list(pairs = aa_npairs, total = nrow(aa_global),
                  geq01 = sum(aa_first$geq_pt01pct), geq1 = sum(aa_first$geq_pt1pct),
                  geq10 = sum(aa_first$geq_1pct), spike = sum(aa_first$is_spike))
msg("  >=0.01%%: %d | >=0.1%%: %d | >=1%%: %d | Spike: %d",
    aa_counts$geq01, aa_counts$geq1, aa_counts$geq10, aa_counts$spike)

rm(aa_long, aa_by_date, aa_global, aa_first); gc()
msg("  Memory freed")

rm(nc2); gc()

# =============================================================================
# 7. Second occurrence + OWID + summary
# =============================================================================

msg("\n[7/7] Second occurrence + contextual data")

# NT second occurrence
nt_by_date <- fread(file.path(SUBS_DIR, "nt_mutations_by_date.csv"))
nt_first   <- fread(file.path(SUBS_DIR, "nt_first_occurrence.csv"))
nt_second <- nt_by_date[, .(dates = list(sort(unique(collection_date)))), by = mutation_name]
nt_second[, first_date  := sapply(dates, function(x) x[1])]
nt_second[, second_date := sapply(dates, function(x) if (length(x) >= 2) x[2] else NA)]
nt_second[, n_dates     := sapply(dates, length)]
nt_second[, dates := NULL]
nt_second[, first_date  := as.Date(first_date,  origin = "1970-01-01")]
nt_second[, second_date := as.Date(second_date, origin = "1970-01-01")]
nt_second[, gap_days    := as.integer(second_date - first_date)]
nt_second <- merge(nt_second, nt_first[, .(mutation_name, count, per_abundance)],
                   by = "mutation_name", all.x = TRUE)
fwrite(nt_second, file.path(SECOND_DIR, "nt_second_occurrence.csv"))
msg("  NT second-occurrence rows: %d", nrow(nt_second))
rm(nt_by_date, nt_first, nt_second); gc()

# AA second occurrence
aa_by_date <- fread(file.path(SUBS_DIR, "aa_mutations_by_date.csv"))
aa_first   <- fread(file.path(SUBS_DIR, "aa_first_occurrence.csv"))
aa_second <- aa_by_date[, .(dates = list(sort(unique(collection_date)))), by = aa_label]
aa_second[, first_date  := sapply(dates, function(x) x[1])]
aa_second[, second_date := sapply(dates, function(x) if (length(x) >= 2) x[2] else NA)]
aa_second[, n_dates     := sapply(dates, length)]
aa_second[, dates := NULL]
aa_second[, first_date  := as.Date(first_date,  origin = "1970-01-01")]
aa_second[, second_date := as.Date(second_date, origin = "1970-01-01")]
aa_second[, gap_days    := as.integer(second_date - first_date)]
aa_second <- merge(aa_second, aa_first[, .(aa_label, count, per_abundance)],
                   by = "aa_label", all.x = TRUE)
fwrite(aa_second, file.path(SECOND_DIR, "aa_second_occurrence.csv"))
msg("  AA second-occurrence rows: %d", nrow(aa_second))
rm(aa_by_date, aa_first, aa_second); gc()

# Build world_infections.csv as NEJM (Dec 1 2019 – Jan 21 2020) + OWID (Jan 22+).
# Prefer reference_data/world_infections_owid_raw.csv (a daily-granularity OWID
# snapshot kept in the repo) — current owid-covid-data.csv has weekly-aggregated
# values for early 2020 which would distort Fig 1C/D.
NEJM_FILE     <- file.path(PROJECT_DIR, "reference_data", "nejm_early_cases.csv")
OWID_RAW_FILE <- file.path(PROJECT_DIR, "reference_data", "world_infections_owid_raw.csv")

if (file.exists(NEJM_FILE) && file.exists(OWID_RAW_FILE)) {
  nejm <- fread(NEJM_FILE)[, .(collection_date = as.Date(collection_date),
                                new_cases = as.numeric(new_cases))]
  owid_raw <- fread(OWID_RAW_FILE)[, .(collection_date = as.Date(collection_date),
                                        new_cases = as.numeric(new_cases))]
  owid_post <- owid_raw[collection_date >= as.Date("2020-01-22")]
  world <- rbindlist(list(nejm, owid_post), use.names = TRUE)
  setorder(world, collection_date)
  world[, total_cases := cumsum(new_cases)]
  world[, new_cases_smoothed := round(frollmean(new_cases, 7, na.rm = TRUE,
                                                align = "right"), 2)]
  world[is.na(new_cases_smoothed), new_cases_smoothed := 0]
  world <- world[collection_date >= START_DATE & collection_date <= END_DATE]
  fwrite(world, file.path(OUT_DIR, "world_infections.csv"))
  msg("  Infections (NEJM + OWID-snapshot merge): %d rows", nrow(world))
} else if (file.exists(OWID_FILE)) {
  owid_full <- fread(OWID_FILE)
  world <- owid_full[location == "World", .(
    collection_date = as.Date(date),
    new_cases = fifelse(is.na(new_cases), 0, new_cases),
    new_cases_smoothed = fifelse(is.na(new_cases_smoothed), 0, new_cases_smoothed),
    total_cases = total_cases)]
  world <- world[collection_date >= START_DATE & collection_date <= END_DATE]
  fwrite(world, file.path(OUT_DIR, "world_infections.csv"))
  msg("  Infections (OWID only — fallback, NEJM/OWID-snapshot missing): %d rows", nrow(world))
  rm(owid_full)
}

# Vaccination data still comes from full OWID file
if (file.exists(OWID_FILE)) {
  owid_full <- fread(OWID_FILE)
  vacc <- owid_full[location == "World" & !is.na(people_fully_vaccinated_per_hundred),
                    .(collection_date = as.Date(date),
                      pct_fully_vaccinated = people_fully_vaccinated_per_hundred)]
  vacc <- vacc[collection_date >= START_DATE & collection_date <= END_DATE]
  fwrite(vacc, file.path(OUT_DIR, "vaccination_data.csv"))
  msg("  Vaccination: %d rows", nrow(vacc))
  rm(owid_full, vacc)
} else {
  msg("  OWID file not found at %s — skipping vaccination_data.csv", OWID_FILE)
}
if (exists("world")) rm(world)
if (exists("owid"))  rm(owid)

# Summary
summary_text <- sprintf(
"Phase 1 Processing Summary
===================================================
Date: %s | Window: %s to %s
Sequences: %s retained of %s input

SUBSTITUTIONS:
  NT pairs: %s | Unique >=0.01%%: %d | >=0.1%%: %d | >=1%%: %d
  AA pairs: %s | Unique >=0.01%%: %d | >=0.1%%: %d | >=1%%: %d | Spike: %d
",
  Sys.time(), START_DATE, END_DATE,
  format(NSEQS, big.mark = ","), format(nc_total_rows, big.mark = ","),
  format(nt_counts$pairs, big.mark = ","), nt_counts$geq01, nt_counts$geq1, nt_counts$geq10,
  format(aa_counts$pairs, big.mark = ","), aa_counts$geq01, aa_counts$geq1, aa_counts$geq10, aa_counts$spike)

writeLines(summary_text, file.path(OUT_DIR, "processing_summary.txt"))
writeLines(as.character(NSEQS), file.path(OUT_DIR, "nseqs.txt"))
cat(summary_text)

msg("\n%s", paste(rep("=", 70), collapse = ""))
msg("Phase 1 complete! Outputs: %s", OUT_DIR)
