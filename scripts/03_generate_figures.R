#!/usr/bin/env Rscript
#' =============================================================================
#' Generate all manuscript figures + per-panel source data
#' =============================================================================
#'
#' Author: Akhil Kumar
#'
#' Reads the tables produced by scripts/02_process_nextclade.R and writes the
#' four final figures plus one CSV per figure panel (exact points/lines/bars
#' plotted). No intermediate or exploratory output.
#'
#' Inputs (from <project_dir>/processed_data/):
#'   substitutions_only/nt_first_occurrence.csv
#'   substitutions_only/aa_first_occurrence.csv
#'   substitutions_only/nt_mutations_by_date.csv{,.gz}
#'   substitutions_only/aa_mutations_by_date.csv{,.gz}
#'   world_infections.csv
#'   vaccination_data.csv
#'   gisaid_daily_sequence_counts.csv
#'   nseqs.txt
#'
#' Outputs:
#'   <project_dir>/figures/figure1.{pdf,jpeg}            Timeline (6 panels A-F)
#'   <project_dir>/figures/figure2.{pdf,jpeg}            Second-occurrence analysis
#'                                                       (A=NT scatter, B=AA scatter,
#'                                                        C=joint median±SD,
#'                                                        D=>=1% bar @ 28 d)
#'   <project_dir>/figures/figure3.{pdf,jpeg}            Frequency + non-synonymous (4 panels)
#'   <project_dir>/figures/supplementary_figure1.{pdf,jpeg}  Overview (panels A-D)
#'   <project_dir>/source_data/fig{1..3,s1*}.csv         One CSV per panel
#'
#' Usage:
#'   Rscript scripts/03_generate_figures.R [project_dir]
#'
#' Memory: ~16-32 GB
#' Time:   ~1 minute on 4 cores / 32 GB
#' =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(ggpubr)
  library(ggrepel)
  library(data.table)
})

set.seed(1)

# =============================================================================
# Configuration
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
PROJECT_DIR <- if (length(args) >= 1) args[1] else "."

DATA_DIR   <- file.path(PROJECT_DIR, "processed_data")
FIG_DIR    <- file.path(PROJECT_DIR, "figures")
# Exactly the data points / lines plotted in the final figures, one CSV per panel.
SOURCE_DATA_DIR <- file.path(PROJECT_DIR, "source_data")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(SOURCE_DATA_DIR, showWarnings = FALSE, recursive = TRUE)

# Read a CSV that may be gzipped — uses system `gunzip -c` on .gz so we don't
# need the R.utils package. Returns a data.table.
read_csv_maybe_gz <- function(path) {
  gzpath <- paste0(path, ".gz")
  if (file.exists(path)) {
    fread(path)
  } else if (file.exists(gzpath)) {
    fread(cmd = paste("gunzip -c", shQuote(gzpath)))
  } else {
    stop("Neither ", path, " nor ", gzpath, " exists", call. = FALSE)
  }
}

# Read dates from config if available
config_file <- file.path(PROJECT_DIR, "filtered_data", "pipeline_config.txt")
if (file.exists(config_file)) {
  config <- readLines(config_file)
  for (line in config) {
    parts <- strsplit(line, "=")[[1]]
    if (length(parts) == 2) {
      if (parts[1] == "START_DATE") START_DATE <- as.Date(parts[2])
      if (parts[1] == "END_DATE")   END_DATE   <- as.Date(parts[2])
    }
  }
} else {
  START_DATE <- as.Date("2019-12-01")
  END_DATE   <- as.Date("2025-12-31")
}

# Colors
COLOR_NT <- "#7C9082"
COLOR_AA <- "#846267"
COLOR_VALS  <- c("[0.01%, 1%)" = "#EC746C", "[1%, 10%)" = "#59C0C5", ">=10%" = "#C77BFE")
SHAPE_VALS  <- c("[0.01%, 1%)" = 16, "[1%, 10%)" = 17, ">=10%" = 15)

msg <- function(...) message(sprintf(...))

# Month factor levels for plotting (dynamic based on date range)
MONTH_LEVELS <- format(seq.Date(
  as.Date(format(START_DATE, "%Y-%m-01")),
  as.Date(format(END_DATE, "%Y-%m-01")),
  by = "month"
), "%b-%y")

msg("=" |> rep(70) |> paste(collapse = ""))
msg("Generating Manuscript Figures")
msg("=" |> rep(70) |> paste(collapse = ""))

# =============================================================================
# Load Data
# =============================================================================

msg("\n[1/6] Loading data...")

NSEQS <- as.integer(readLines(file.path(DATA_DIR, "nseqs.txt")))
msg("  Total sequences (NSEQS): %s", format(NSEQS, big.mark = ","))

SUBS_DIR <- file.path(DATA_DIR, "substitutions_only")

# NT data
nt_first <- fread(file.path(SUBS_DIR, "nt_first_occurrence.csv"))
nt_first[, collection_date := as.Date(collection_date)]
msg("  NT first occurrence: %s mutations", format(nrow(nt_first), big.mark = ","))

# AA data
aa_first <- fread(file.path(SUBS_DIR, "aa_first_occurrence.csv"))
aa_first[, collection_date := as.Date(collection_date)]
msg("  AA first occurrence: %s mutations", format(nrow(aa_first), big.mark = ","))

# NT by date (for second occurrence analysis)
nt_by_date <- read_csv_maybe_gz(file.path(SUBS_DIR, "nt_mutations_by_date.csv"))
nt_by_date[, collection_date := as.Date(collection_date)]

# AA by date (for second occurrence analysis)
aa_by_date <- read_csv_maybe_gz(file.path(SUBS_DIR, "aa_mutations_by_date.csv"))
aa_by_date[, collection_date := as.Date(collection_date)]

# Infections
infections <- fread(file.path(DATA_DIR, "world_infections.csv"))
infections[, collection_date := as.Date(collection_date)]

# Dec-2019 and early Jan 2020 infections now included in world_infections.csv
# (merged from NEJM Li et al. 2020 Figure 1 data + OWID from Jan 22 onward)
msg("  Infections data: %d rows, %s to %s",
    nrow(infections), min(infections$collection_date), max(infections$collection_date))

infections[, TimelineLabel := format(collection_date, "%b-%y")]
infections_monthly <- infections[, .(newCases = sum(new_cases, na.rm = TRUE)),
                                 by = .(TimelineLabel)]
infections_monthly[, TimelineLabel := factor(TimelineLabel, levels = MONTH_LEVELS)]

# Vaccination (load if exists)
vacc_file <- file.path(DATA_DIR, "vaccination_data.csv")
has_vacc <- file.exists(vacc_file)
if (has_vacc) {
  vacc <- fread(vacc_file)
  vacc[, collection_date := as.Date(collection_date)]
}

# =============================================================================
# Helper: create threshold subsets as tibbles for ggplot
# =============================================================================

# NT subsets
nt.geq1   <- as_tibble(nt_first[geq_1pct == TRUE])
nt.geqpt1 <- as_tibble(nt_first[geq_pt1pct == TRUE])
nt.geqpt01<- as_tibble(nt_first[geq_pt01pct == TRUE])

# AA subsets
aa.geq1   <- as_tibble(aa_first[geq_1pct == TRUE])
aa.geqpt1 <- as_tibble(aa_first[geq_pt1pct == TRUE])
aa.geqpt01<- as_tibble(aa_first[geq_pt01pct == TRUE])

# Top 100
nt.top100 <- as_tibble(nt_first[order(-per_abundance)][1:min(100, .N)])
aa.top100 <- as_tibble(aa_first[order(-per_abundance)][1:min(100, .N)])

msg("  NT >= 0.01%%: %d, >= 0.1%%: %d, >= 1%%: %d",
    nrow(nt.geqpt01), nrow(nt.geqpt1), nrow(nt.geq1))
msg("  AA >= 0.01%%: %d, >= 0.1%%: %d, >= 1%%: %d",
    nrow(aa.geqpt01), nrow(aa.geqpt1), nrow(aa.geq1))

# =============================================================================
# Build Timeline Dataframe (daily counts of first occurrences)
# =============================================================================

msg("\n[2/6] Building timeline dataframe...")

all_days <- tibble(Timeline = seq(START_DATE, END_DATE, by = "1 day"))

count_first <- function(df, col_name) {
  df |>
    filter(collection_date >= START_DATE, collection_date <= END_DATE) |>
    count(collection_date, name = col_name)
}

DF <- all_days
for (info in list(
  list(nt.geq1, "nt.geq1"), list(nt.geqpt1, "nt.geqpt1"), list(nt.geqpt01, "nt.geqpt01"),
  list(aa.geq1, "aa.geq1"), list(aa.geqpt1, "aa.geqpt1"), list(aa.geqpt01, "aa.geqpt01")
)) {
  counts <- count_first(info[[1]], info[[2]])
  DF <- left_join(DF, counts, by = c("Timeline" = "collection_date")) |>
    mutate(!!info[[2]] := replace_na(!!sym(info[[2]]), 0))
}
DF$TimelineLabel <- format(DF$Timeline, "%b-%y")

msg("  Timeline rows: %d", nrow(DF))

# =============================================================================
# FIGURE 1: Timeline of first occurrence of all unique mutations
# =============================================================================

msg("\n[3/6] Generating Figure 1...")

# --- Pie chart data ---
DF_pie <- DF |>
  pivot_longer(cols = matches("^(nt|aa)\\.(geq|top)")) |>
  mutate(
    TimelineLabel = factor(TimelineLabel, levels = MONTH_LEVELS),
    bins = ifelse(Timeline <= as.Date("2020-03-31"), "Dec-19 to Mar-20", "Rest")
  ) |>
  group_by(TimelineLabel, name, bins) |>
  summarize(countMutations = sum(value), .groups = "drop")

# Helper to make timeline + pie inset
make_timeline_panel <- function(data, col, color, title_text, n_total) {
  p_line <- data |>
    pivot_longer(cols = all_of(col)) |>
    mutate(flag = ifelse(Timeline <= as.Date("2020-03-31"), "Dec-19 to Mar-20", "Rest")) |>
    ggplot(aes(x = Timeline, y = value, group = interaction(name, flag))) +
    geom_path(aes(color = flag)) +
    scale_x_date(breaks = seq.Date(START_DATE, END_DATE, by = "3 months"), date_labels = "%b-%y") +
    theme_classic(base_size = 14) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
          legend.position = "none",
          plot.title = element_text(hjust = 0.5)) +
    scale_color_manual(values = c("Dec-19 to Mar-20" = color, "Rest" = "darkgrey")) +
    labs(x = "Timeline",
         y = paste0("Number of new mutations\nin a day (n = ", format(n_total, big.mark = ","), ")"),
         title = title_text)

  p_pie <- DF_pie |>
    group_by(name) |>
    mutate(norm = countMutations / sum(countMutations)) |>
    filter(name == col) |>
    ggplot(aes(x = "", y = norm, fill = bins)) +
    geom_bar(stat = "identity") + coord_polar(theta = "y") +
    scale_fill_manual(values = c("Dec-19 to Mar-20" = color, "Rest" = "darkgrey")) +
    theme_void(base_size = 14) +
    theme(legend.position = "left", legend.title = element_blank())

  p_line + inset_element(p_pie, left = 0.45, bottom = 0.3, right = 0.98, top = 1, align_to = "plot")
}

# Panel A: NT >= 0.01% timeline
p1A <- make_timeline_panel(DF, "nt.geqpt01", COLOR_NT, "Nucleotide", sum(DF$nt.geqpt01))

# Panel B: AA >= 0.01% timeline
p1B <- make_timeline_panel(DF, "aa.geqpt01", COLOR_AA, "Amino acid", sum(DF$aa.geqpt01))

# --- Panels C & D: Mutations per thousand infections (line plot with date axis) ---
# Aggregate by month using actual dates for proper x-axis
DF_inf <- DF |>
  mutate(month_date = as.Date(format(Timeline, "%Y-%m-01"))) |>
  group_by(month_date) |>
  summarize(
    nt_count = sum(nt.geqpt01),
    aa_count = sum(aa.geqpt01),
    .groups = "drop"
  )

# Join with monthly infections using date-based grouping
inf_monthly_date <- infections |>
  mutate(month_date = as.Date(format(collection_date, "%Y-%m-01"))) |>
  group_by(month_date) |>
  summarize(newCases = sum(new_cases, na.rm = TRUE), .groups = "drop")

DF_inf <- left_join(DF_inf, inf_monthly_date, by = "month_date") |>
  mutate(
    newCases = replace_na(newCases, 0),
    nt_per_1000 = ifelse(newCases > 0, (nt_count / newCases) * 1000, 0),
    aa_per_1000 = ifelse(newCases > 0, (aa_count / newCases) * 1000, 0),
    flag = ifelse(month_date <= as.Date("2020-03-31"), "Dec-19 to Mar-20", "Rest")
  ) |>
  filter(newCases > 0 | nt_count > 0 | aa_count > 0)  # drop months with no data at all

make_per_1000_panel <- function(data, y_col, color, title_text, y_label) {
  early_data <- data |> filter(flag == "Dec-19 to Mar-20")
  ggplot(data, aes(x = month_date, y = .data[[y_col]])) +
    geom_line(aes(group = 1), color = "darkgrey", linewidth = 0.7) +
    geom_line(data = early_data, aes(group = 1), color = color, linewidth = 0.7) +
    geom_point(aes(color = flag), size = 1.5) +
    scale_x_date(breaks = seq.Date(START_DATE, END_DATE, by = "3 months"),
                 date_labels = "%b-%y", limits = c(START_DATE, END_DATE)) +
    scale_color_manual(values = c("Dec-19 to Mar-20" = color, "Rest" = "darkgrey")) +
    theme_classic(base_size = 14) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
          legend.position = "none",
          plot.title = element_text(hjust = 0.5)) +
    labs(x = "Timeline", y = y_label, title = title_text)
}

p1C <- make_per_1000_panel(DF_inf, "nt_per_1000", COLOR_NT, "Nucleotide",
                            "Number of new nucleotide mutations\nper thousand infections")
p1D <- make_per_1000_panel(DF_inf, "aa_per_1000", COLOR_AA, "Amino acid",
                            "Number of new amino acid mutations\nper thousand infections")

# --- Panels E & F: Scatter (frequency vs first occurrence date) ---
make_scatter <- function(df, top100, label_col, title_text, y_label, n_labels = 15) {
  plotDF <- df |>
    mutate(bins = cut(per_abundance, breaks = c(0.01, 0.1, 1, 10, 100),
                      labels = c("(0.01%,0.1%]", "(0.1%,1%]", "(1%,10%]", "(10%,100%]"),
                      include.lowest = TRUE))
  countDF <- plotDF |> group_by(bins) |> summarise(Count = n(), .groups = "drop")

  ggplot(plotDF, aes(x = collection_date, y = per_abundance, color = bins, alpha = 0.05)) +
    geom_point(size = 1.5) +
    theme_classic(base_size = 14) +
    ggtitle(title_text) +
    scale_x_date(breaks = seq.Date(START_DATE, END_DATE, by = "3 months"), date_labels = "%b-%y",
                 limits = c(START_DATE, END_DATE)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
    geom_text_repel(
      aes(label = !!sym(label_col)),
      data = filter(plotDF, !!sym(label_col) %in% head(top100[[label_col]], n_labels)),
      size = 2.8, min.segment.length = 0, seed = 42, box.padding = 0.5,
      max.overlaps = Inf, force = 3, force_pull = 0.3,
      arrow = arrow(length = unit(0.008, "npc")),
      color = "grey10", segment.color = "grey50"
    ) +
    guides(colour = guide_legend(override.aes = list(alpha = 1)), alpha = "none") +
    scale_color_manual(
      values = c("(0.01%,0.1%]" = "#EC746C", "(0.1%,1%]" = "#59C0C5",
                 "(1%,10%]" = "#7CAE00", "(10%,100%]" = "#C77BFE"),
      labels = paste0("n(", levels(countDF$bins), ") = ", countDF$Count)
    ) +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "inside",
          legend.position.inside = c(0.98, 0.98),
          legend.justification = c("right", "top"),
          legend.title = element_blank(), legend.margin = margin(6, 6, 6, 6),
          legend.background = element_rect(fill = "white", color = NA)) +
    labs(x = "Timeline", y = y_label) +
    scale_y_continuous(labels = function(x) paste0(x, "%"),
                       limits = c(0, 100), expand = expansion(mult = c(0, 0.02)))
}

p1E <- make_scatter(nt.geqpt01, nt.top100, "mutation_name", "Nucleotide",
                    "Frequency of nucleotide mutations\nin the population", n_labels = 15)
p1F <- make_scatter(aa.geqpt01, aa.top100, "aa_label", "Amino acid",
                    "Frequency of amino acid mutations\nin the population", n_labels = 20)

# Assemble Figure 1
fig1 <- (p1A | p1B) / (p1C | p1D) / (p1E | p1F) +
  plot_annotation(tag_levels = list(c("A", "", "B", "", "C", "D", "E", "F")))

ggsave(file.path(FIG_DIR, "figure1.pdf"), plot = fig1, width = 16, height = 16)
ggsave(file.path(FIG_DIR, "figure1.jpeg"), plot = fig1, dpi = 300, width = 16.5, height = 16.5)
msg("  Figure 1 saved!")

# --- Source data for Fig 1: exactly the points/lines plotted ---
# DF and DF_inf are tibbles; use dplyr-style selection before fwrite.
fwrite(DF |>
         dplyr::transmute(date = Timeline,
                          nt_new_geq_pt01pct = nt.geqpt01,
                          aa_new_geq_pt01pct = aa.geqpt01),
       file.path(SOURCE_DATA_DIR, "fig1ab_daily_first_occurrences.csv"))
fwrite(DF_inf,
       file.path(SOURCE_DATA_DIR, "fig1cd_monthly_rates.csv"))
fwrite(nt.geqpt01[, c("mutation_name", "pos", "ref", "alt", "collection_date",
                      "count", "per_abundance")],
       file.path(SOURCE_DATA_DIR, "fig1e_nt_scatter.csv"))
fwrite(aa.geqpt01[, c("aa_label", "region", "mutation_name", "pos", "ref", "alt",
                      "collection_date", "count", "per_abundance")],
       file.path(SOURCE_DATA_DIR, "fig1f_aa_scatter.csv"))
fwrite(DF_pie, file.path(SOURCE_DATA_DIR, "fig1_pie_inset_data.csv"))

# =============================================================================
# FIGURE 2: Mutation abundance & second occurrence
#   Panel A: NT scatter (Jan 15 cutoff)
#   Panel B: Joint NT+AA median ± SD by frequency bin (3 bins), connected
#   Panel C: AA scatter (Jan 15 cutoff)
#   Panel D: >= 10% mutations with <= 15-day second-occurrence heuristic
# =============================================================================

msg("\n[4/6] Generating Figure 2...")

# Cutoffs shared across Figure 2 panels
FIG2_FIRST_DATE <- as.Date("2020-01-15")   # "first seen by Jan 15"
FIG2_DAYS_CUTOFF <- 15                      # <= 15 days second-occurrence heuristic

# --- Calculate second occurrence for NT ---
calc_second_occurrence <- function(by_date_dt) {
  # For each mutation, find first and second date it was observed
  by_date_dt[, .(collection_date = sort(unique(collection_date))),
             by = .(mutation_name = if ("mutation_name" %in% names(by_date_dt)) mutation_name
                    else if ("aa_label" %in% names(by_date_dt)) aa_label
                    else stop("No mutation name column found"))
  ][, .(first_date = collection_date[1],
        second_date = if (.N > 1) collection_date[2] else as.Date(NA)),
    by = mutation_name
  ][, days_to_second := as.integer(second_date - first_date)]
}

# NT second occurrence
nt_intervals <- nt_by_date[, .(collection_date = sort(unique(collection_date))),
                           by = .(mutation_name)
][, .(first_date = collection_date[1],
      second_date = if (.N > 1) collection_date[2] else as.Date(NA)),
  by = mutation_name
][, days_to_second := as.integer(second_date - first_date)]

# Merge with global counts, filter to subs >= 0.01%
nt_second <- merge(nt_intervals, as.data.table(nt.geqpt01)[, .(mutation_name, per_abundance, count)],
                   by = "mutation_name")

# AA second occurrence
aa_intervals <- aa_by_date[, .(collection_date = sort(unique(collection_date))),
                           by = .(aa_label)
][, .(first_date = collection_date[1],
      second_date = if (.N > 1) collection_date[2] else as.Date(NA)),
  by = aa_label
][, days_to_second := as.integer(second_date - first_date)]

aa_second <- merge(aa_intervals, as.data.table(aa.geqpt01)[, .(aa_label, per_abundance, count)],
                   by = "aa_label")

# Restrict to mutations first seen by FIG2_FIRST_DATE (Jan 15, 2020)
nt_early <- nt_second[first_date <= FIG2_FIRST_DATE & !is.na(days_to_second)]
nt_early[, abundance := cut(per_abundance, breaks = c(0.01, 1, 10, 100),
                            labels = c("[0.01%, 1%)", "[1%, 10%)", ">=10%"),
                            include.lowest = TRUE, right = FALSE)]
nt_early[, freq_group := cut(per_abundance,
                              breaks = c(0.01, 0.1, 10, 100),
                              labels = c("[0.01%, 0.1%)", "[0.1%, 10%)", ">=10%"),
                              include.lowest = TRUE, right = FALSE)]
msg("  NT mutations with second occurrence (first seen by Jan 15): %d", nrow(nt_early))

aa_early <- aa_second[first_date <= FIG2_FIRST_DATE & !is.na(days_to_second)]
aa_early[, abundance := cut(per_abundance, breaks = c(0.01, 1, 10, 100),
                            labels = c("[0.01%, 1%)", "[1%, 10%)", ">=10%"),
                            include.lowest = TRUE, right = FALSE)]
aa_early[, freq_group := cut(per_abundance,
                              breaks = c(0.01, 0.1, 10, 100),
                              labels = c("[0.01%, 0.1%)", "[0.1%, 10%)", ">=10%"),
                              include.lowest = TRUE, right = FALSE)]
msg("  AA mutations with second occurrence (first seen by Jan 15): %d", nrow(aa_early))

# Save second occurrence data for reproducibility
fwrite(nt_early, file.path(DATA_DIR, "second_occurrence", "second_occurrence_nt.csv"))
fwrite(aa_early, file.path(DATA_DIR, "second_occurrence", "second_occurrence_aa.csv"))

# --- Figure 2 configuration ---
# Jan 15, 2020 first-seen cutoff for every panel; panel D uses a 28-day second-
# occurrence heuristic applied to mutations reaching >= 1% pandemic-wide frequency.
FIG2_CUTOFF      <- list(date = as.Date("2020-01-15"), label = "Jan 15, 2020", short = "jan15")
FIG2_HEURISTIC   <- 28
FIG2_ABUNDANCE   <- 1
FIG2B_GRP_LEVELS <- c("[0.01%, 0.1%)", "[0.1%, 10%)", ">=10%")
FIG2_TYPE_COLORS <- c("Nucleotide" = "#1F78B4", "Amino acid" = "#E31A1C")

# Consistent typography across all Figure 2 panels (A, B, C, D).
# All text black; plot.title left-aligned; legend on top, horizontal.
fig2_common_theme <- function() {
  theme(
    plot.title       = element_text(size = 16, hjust = 0, face = "plain", color = "black"),
    axis.title.x     = element_text(size = 13, color = "black"),
    axis.title.y     = element_text(size = 13, color = "black"),
    axis.text.x      = element_text(size = 12, color = "black"),
    axis.text.y      = element_text(size = 12, color = "black"),
    legend.title     = element_text(size = 13, color = "black"),
    legend.text      = element_text(size = 12, color = "black"),
    legend.position  = "top",
    legend.direction = "horizontal"
  )
}

build_median_sd <- function(dt, type_label) {
  counts <- dt[, .N, by = freq_group]
  all_counts <- data.table(freq_group = factor(FIG2B_GRP_LEVELS, levels = FIG2B_GRP_LEVELS))
  all_counts <- merge(all_counts, counts, by = "freq_group", all.x = TRUE)
  all_counts[is.na(N), N := 0L]
  stats <- dt[!is.na(freq_group), .(median_days = as.numeric(median(days_to_second, na.rm = TRUE)),
                                     sd_days     = as.numeric(stats::sd(days_to_second, na.rm = TRUE)),
                                     n           = .N),
              by = freq_group]
  stats <- merge(all_counts[, .(freq_group)], stats, by = "freq_group", all.x = TRUE)
  stats[, type := factor(type_label, levels = c("Nucleotide", "Amino acid"))]
  stats[, freq_group := factor(freq_group, levels = FIG2B_GRP_LEVELS)]
  stats
}

pairwise_wilcox <- function(dt, bins) {
  pairs <- list(c(bins[1], bins[2]), c(bins[2], bins[3]), c(bins[1], bins[3]))
  rbindlist(lapply(pairs, function(pair) {
    x <- dt[freq_group == pair[1], days_to_second]
    y <- dt[freq_group == pair[2], days_to_second]
    if (length(x) < 2 || length(y) < 2) {
      return(data.table(grp1 = pair[1], grp2 = pair[2], p = NA_real_, stars = "ns"))
    }
    p <- suppressWarnings(wilcox.test(x, y)$p.value)
    stars <- ifelse(is.na(p), "ns",
             ifelse(p < 0.001, "***",
             ifelse(p < 0.01, "**",
             ifelse(p < 0.05, "*", "ns"))))
    data.table(grp1 = pair[1], grp2 = pair[2], p = p, stars = stars)
  }))
}

build_fig2 <- function(cutoff_date, cutoff_label, heuristic_days, min_abundance = 10) {
  nt_f <- nt_second[first_date <= cutoff_date & !is.na(days_to_second)]
  nt_f[, abundance := cut(per_abundance, breaks = c(0.01, 1, 10, 100),
                          labels = c("[0.01%, 1%)", "[1%, 10%)", ">=10%"),
                          include.lowest = TRUE, right = FALSE)]
  nt_f[, freq_group := cut(per_abundance, breaks = c(0.01, 0.1, 10, 100),
                           labels = FIG2B_GRP_LEVELS,
                           include.lowest = TRUE, right = FALSE)]

  aa_f <- aa_second[first_date <= cutoff_date & !is.na(days_to_second)]
  aa_f[, abundance := cut(per_abundance, breaks = c(0.01, 1, 10, 100),
                          labels = c("[0.01%, 1%)", "[1%, 10%)", ">=10%"),
                          include.lowest = TRUE, right = FALSE)]
  aa_f[, freq_group := cut(per_abundance, breaks = c(0.01, 0.1, 10, 100),
                           labels = FIG2B_GRP_LEVELS,
                           include.lowest = TRUE, right = FALSE)]

  pA <- ggscatter(as.data.frame(nt_f), x = "days_to_second", y = "per_abundance",
                  color = "abundance", shape = "abundance", size = 3,
                  add = "reg.line", add.params = list(color = "grey40", fill = "lightgray"),
                  conf.int = TRUE, cor.coef = TRUE,
                  cor.coeff.args = list(method = "spearman", label.sep = "\n",
                                        size = 5, color = "black")) +
    scale_color_manual(values = COLOR_VALS, drop = FALSE) +
    scale_shape_manual(values = SHAPE_VALS, drop = FALSE) +
    scale_y_continuous(labels = function(x) paste0(x, "%"),
                       limits = c(NA, 100), expand = expansion(mult = c(0.02, 0.02))) +
    labs(title = sprintf("Nucleotide (first seen by %s, N=%d)", cutoff_label, nrow(nt_f)),
         x = "Days to Second Occurrence", y = "Abundance at 5 years (%)",
         color = "Abundance", shape = "Abundance") +
    fig2_common_theme()

  pC <- ggscatter(as.data.frame(aa_f), x = "days_to_second", y = "per_abundance",
                  color = "abundance", shape = "abundance", size = 3,
                  add = "reg.line", add.params = list(color = "grey40", fill = "lightgray"),
                  conf.int = TRUE, cor.coef = TRUE,
                  cor.coeff.args = list(method = "spearman", label.sep = "\n",
                                        size = 5, color = "black")) +
    scale_color_manual(values = COLOR_VALS, drop = FALSE) +
    scale_shape_manual(values = SHAPE_VALS, drop = FALSE) +
    scale_y_continuous(labels = function(x) paste0(x, "%"),
                       limits = c(NA, 100), expand = expansion(mult = c(0.02, 0.02))) +
    labs(title = sprintf("Amino acid (first seen by %s, N=%d)", cutoff_label, nrow(aa_f)),
         x = "Days to Second Occurrence", y = "Abundance at 5 years (%)",
         color = "Abundance", shape = "Abundance") +
    fig2_common_theme()

  js <- rbind(build_median_sd(nt_f, "Nucleotide"),
              build_median_sd(aa_f, "Amino acid"))

  nt_n_l <- nt_f[!is.na(freq_group), .(n_nt = .N), by = freq_group]
  aa_n_l <- aa_f[!is.na(freq_group), .(n_aa = .N), by = freq_group]
  bc <- data.table(freq_group = factor(FIG2B_GRP_LEVELS, levels = FIG2B_GRP_LEVELS))
  bc <- merge(bc, nt_n_l, by = "freq_group", all.x = TRUE)
  bc <- merge(bc, aa_n_l, by = "freq_group", all.x = TRUE)
  bc[is.na(n_nt), n_nt := 0L]
  bc[is.na(n_aa), n_aa := 0L]
  bc[, xlab := sprintf("%s\nNT=%d\nAA=%d", freq_group, n_nt, n_aa)]

  nt_sig_l <- pairwise_wilcox(nt_f, FIG2B_GRP_LEVELS); nt_sig_l[, type := "Nucleotide"]
  aa_sig_l <- pairwise_wilcox(aa_f, FIG2B_GRP_LEVELS); aa_sig_l[, type := "Amino acid"]
  sdf <- rbind(nt_sig_l, aa_sig_l)
  sdf[, type := factor(type, levels = c("Nucleotide", "Amino acid"))]
  sdf[, grp1_idx := match(grp1, FIG2B_GRP_LEVELS)]
  sdf[, grp2_idx := match(grp2, FIG2B_GRP_LEVELS)]
  sdf[, pair_kind := ifelse(grp2_idx - grp1_idx == 1, "adj", "across")]

  js[, x_num := as.numeric(freq_group)]
  js[type == "Nucleotide", x_offset := x_num - 0.12]
  js[type == "Amino acid",  x_offset := x_num + 0.12]

  ymax_l <- max(js$median_days + js$sd_days, na.rm = TRUE)
  # Band y-positions: NT and AA are separated slightly more than before for readability.
  Y_NT_ADJ_L    <- ymax_l * 1.25
  Y_NT_ACROSS_L <- ymax_l * 1.55
  Y_AA_ADJ_L    <- ymax_l * 2.05
  Y_AA_ACROSS_L <- ymax_l * 2.55
  Y_TOP_L       <- ymax_l * 3.15

  sdf[, y_pos := dplyr::case_when(
    type == "Nucleotide" & pair_kind == "adj"    ~ Y_NT_ADJ_L,
    type == "Nucleotide" & pair_kind == "across" ~ Y_NT_ACROSS_L,
    type == "Amino acid" & pair_kind == "adj"    ~ Y_AA_ADJ_L,
    type == "Amino acid" & pair_kind == "across" ~ Y_AA_ACROSS_L,
    TRUE ~ NA_real_)]

  js$x_plot <- js$x_offset
  sig_local <- copy(sdf)
  sig_local[type == "Nucleotide", `:=`(x_start = grp1_idx - 0.12, x_end = grp2_idx - 0.12)]
  sig_local[type == "Amino acid",  `:=`(x_start = grp1_idx + 0.12, x_end = grp2_idx + 0.12)]
  sig_local <- sig_local[!is.na(y_pos)]
  TIP <- 0.96

  pB <- ggplot(js[!is.na(median_days)],
               aes(x = x_plot, y = median_days, color = type, group = type)) +
    geom_line(linewidth = 0.9, na.rm = TRUE) +
    geom_errorbar(aes(ymin = pmax(median_days - sd_days, 1),
                      ymax = median_days + sd_days),
                  width = 0.1, linewidth = 0.6, na.rm = TRUE) +
    geom_point(size = 3.5, na.rm = TRUE) +
    geom_segment(data = sig_local,
                 aes(x = x_start, xend = x_end, y = y_pos, yend = y_pos),
                 color = "black", linewidth = 0.45, inherit.aes = FALSE) +
    geom_segment(data = sig_local,
                 aes(x = x_start, xend = x_start, y = y_pos, yend = y_pos * TIP),
                 color = "black", linewidth = 0.45, inherit.aes = FALSE) +
    geom_segment(data = sig_local,
                 aes(x = x_end, xend = x_end, y = y_pos, yend = y_pos * TIP),
                 color = "black", linewidth = 0.45, inherit.aes = FALSE) +
    # Asterisks render visually higher than "ns" at the same vjust because they
    # extend to the cap height, so we pull them down slightly while keeping "ns"
    # at its original offset.
    geom_text(data = sig_local[stars == "ns"],
              aes(x = (x_start + x_end) / 2, y = y_pos, label = stars),
              color = "black", size = 3.8, vjust = -0.35,
              inherit.aes = FALSE) +
    geom_text(data = sig_local[stars != "ns"],
              aes(x = (x_start + x_end) / 2, y = y_pos, label = stars),
              color = "black", size = 3.8, vjust = -0.05,
              inherit.aes = FALSE) +
    scale_x_continuous(breaks = seq_along(FIG2B_GRP_LEVELS),
                       labels = bc$xlab,
                       limits = c(0.5, length(FIG2B_GRP_LEVELS) + 0.5)) +
    scale_y_log10(breaks = c(1, 3, 10, 30, 100, 300),
                  labels = c("1", "3", "10", "30", "100", "300"),
                  expand = expansion(mult = c(0.02, 0.02))) +
    coord_cartesian(ylim = c(1, Y_TOP_L), clip = "off") +
    scale_color_manual(values = FIG2_TYPE_COLORS) +
    theme_classic(base_size = 12) +
    fig2_common_theme() +
    theme(legend.title = element_blank()) +
    labs(x = "Pandemic-wide frequency",
         y = "Days to second occurrence (median ± SD)",
         title = sprintf("NT & AA median days to second occurrence\n(first seen by %s)",
                         cutoff_label))

  dlbl_lo <- sprintf("<=%d days", heuristic_days)
  dlbl_hi <- sprintf(">%d days",  heuristic_days)
  build_hi <- function(dt, type_label) {
    dt |>
      as.data.frame() |>
      dplyr::filter(per_abundance >= min_abundance) |>
      dplyr::mutate(type = type_label,
                    days_group = factor(ifelse(days_to_second <= heuristic_days,
                                               dlbl_lo, dlbl_hi),
                                        levels = c(dlbl_lo, dlbl_hi)))
  }
  dD <- dplyr::bind_rows(build_hi(nt_f, "NT mutations"),
                         build_hi(aa_f, "AA mutations")) |>
    dplyr::count(type, days_group) |>
    dplyr::group_by(type) |>
    dplyr::mutate(pct = n / sum(n)) |>
    dplyr::ungroup() |>
    dplyr::mutate(type = factor(type, levels = c("NT mutations", "AA mutations")))

  pD <- ggplot(dD, aes(x = type, y = n, fill = days_group)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.5),
             color = "black", linewidth = 0.3, width = 0.4) +
    geom_text(aes(label = sprintf("%.1f%%", pct * 100)),
              position = position_dodge(width = 0.5),
              vjust = -0.3, size = 4.5, color = "black") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
    scale_x_discrete(expand = expansion(add = 0.5)) +
    scale_fill_manual(values = setNames(c("coral", "cornflowerblue"),
                                        c(dlbl_lo, dlbl_hi))) +
    theme_classic(base_size = 12) +
    fig2_common_theme() +
    labs(x = "", y = "Number of Mutations",
         fill = "Days to Second Occurrence:",
         title = sprintf(">=%g%% mutations, <=%dd heuristic\n(first seen by %s)",
                         min_abundance, heuristic_days, cutoff_label))

  # Summary tables used by panels B (median±SD) and its Wilcoxon brackets
  median_sd_df <- copy(js)[!is.na(median_days),
                            .(type, freq_group, n, median_days, sd_days)]
  wilcoxon_df  <- copy(sdf)[, .(type, grp1, grp2, p, stars)]

  list(pA = pA, pB = pB, pC = pC, pD = pD,
       # Exact points/bars plotted, for source_data_final/
       nt_scatter_df = as.data.table(nt_f)[, .(mutation_name, first_date, second_date,
                                                days_to_second, per_abundance, count,
                                                abundance_bin = abundance)],
       aa_scatter_df = as.data.table(aa_f)[, .(aa_label, first_date, second_date,
                                                days_to_second, per_abundance, count,
                                                abundance_bin = abundance)],
       median_sd_df = median_sd_df,
       wilcoxon_df  = wilcoxon_df,
       panelD_df    = as.data.table(dD),
       n_nt = nrow(nt_f), n_aa = nrow(aa_f))
}

# --- Build single Figure 2 (grid approach disabled; see commented block above) ---
# Panel layout: (A=NT scatter, B=AA scatter, C=NT+AA median±SD, D=>=1% bar, <=28d).
panels <- build_fig2(FIG2_CUTOFF$date, FIG2_CUTOFF$label, FIG2_HEURISTIC,
                     min_abundance = FIG2_ABUNDANCE)
fig2 <- (panels$pA | panels$pC) / (panels$pB | panels$pD) +
  plot_annotation(tag_levels = "A")
ggsave(file.path(FIG_DIR, "figure2.pdf"),  plot = fig2, width = 14, height = 12)
ggsave(file.path(FIG_DIR, "figure2.jpeg"), plot = fig2, dpi = 300, width = 14, height = 12)
fwrite(panels$panelD_df,
       file.path(DATA_DIR, "second_occurrence", "figure2d_panel_data.csv"))

# --- Source data for Fig 2: exactly the points/bars plotted ---
fwrite(panels$nt_scatter_df, file.path(SOURCE_DATA_DIR, "fig2a_nt_scatter.csv"))
fwrite(panels$aa_scatter_df, file.path(SOURCE_DATA_DIR, "fig2b_aa_scatter.csv"))
fwrite(panels$median_sd_df,  file.path(SOURCE_DATA_DIR, "fig2c_median_sd.csv"))
fwrite(panels$wilcoxon_df,   file.path(SOURCE_DATA_DIR, "fig2c_wilcoxon_brackets.csv"))
fwrite(panels$panelD_df,     file.path(SOURCE_DATA_DIR, "fig2d_bar.csv"))

msg("  Figure 2 saved (cutoff N: NT=%d, AA=%d; heuristic=%dd; panel D >=%g%%)",
    panels$n_nt, panels$n_aa, FIG2_HEURISTIC, FIG2_ABUNDANCE)

# =============================================================================
# FIGURE 3: Frequency of new mutations & non-synonymous mutation percentages
# =============================================================================

msg("\n[5/6] Generating Figure 3...")

# --- Panels A & B: Frequency distribution of new mutations by time period ---
# Group mutations by frequency bin and time period

# (Dead code removed — freq_by_period had incorrect period labels. Use make_freq_barplot instead.)

# Shared period definitions for Figure 3 (monthly for first 7 months, then 6-month bins)
FIG3_BREAKS <- c(
  as.Date("2019-12-01"), as.Date("2020-01-01"), as.Date("2020-02-01"),
  as.Date("2020-03-01"), as.Date("2020-04-01"), as.Date("2020-05-01"),
  as.Date("2020-06-01"), as.Date("2020-07-01"),
  as.Date("2021-01-01"), as.Date("2021-07-01"), as.Date("2022-01-01"),
  as.Date("2022-07-01"), as.Date("2023-01-01"), as.Date("2023-07-01"),
  as.Date("2024-01-01"), as.Date("2024-07-01"), as.Date("2025-01-01"),
  as.Date("2026-01-01")
)
FIG3_LABELS <- c(
  "Dec-19", "Jan-20", "Feb-20", "Mar-20", "Apr-20", "May-20", "Jun-20",
  "Jul-20 to Dec-20", "Jan-21 to Jun-21", "Jul-21 to Dec-21",
  "Jan-22 to Jun-22", "Jul-22 to Dec-22", "Jan-23 to Jun-23",
  "Jul-23 to Dec-23", "Jan-24 to Jun-24", "Jul-24 to Dec-24",
  "Jan-25 to Dec-25"
)

make_freq_barplot <- function(df, label_col, title_text, source_csv = NULL) {
  plotDF <- df |>
    mutate(
      freq_bin = cut(per_abundance, breaks = c(0.01, 0.1, 1, 10, 100),
                     labels = c("(0.01%,0.1%]", "(0.1%,1%]", "(1%,10%]", "(10%,100%]"),
                     include.lowest = TRUE),
      period = cut(collection_date, breaks = FIG3_BREAKS, labels = FIG3_LABELS,
                   include.lowest = TRUE, right = FALSE)
    ) |>
    filter(!is.na(period), !is.na(freq_bin)) |>
    count(period, freq_bin) |>
    group_by(period) |>
    mutate(pct = n / sum(n)) |>
    ungroup()

  if (!is.null(source_csv)) fwrite(plotDF, source_csv)

  ggplot(plotDF, aes(x = period, y = pct, fill = freq_bin)) +
    geom_bar(stat = "identity", width = 0.7, color = "black", linewidth = 0.3) +
    coord_flip() +
    scale_x_discrete(drop = FALSE, limits = rev(FIG3_LABELS)) +
    scale_y_continuous(labels = scales::percent,
                       expand = expansion(mult = c(0, 0.02))) +
    scale_fill_manual(values = c("(0.01%,0.1%]" = "#FDE725",
                                 "(0.1%,1%]" = "#5DC863",
                                 "(1%,10%]" = "#21908C",
                                 "(10%,100%]" = "#440154")) +
    theme_classic(base_size = 11) +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "bottom",
          legend.title = element_blank()) +
    labs(x = "", y = "", title = paste("Frequency of new", title_text, "mutations"))
}

p3A <- make_freq_barplot(nt.geqpt01, "mutation_name", "nucleotide",
                          source_csv = file.path(SOURCE_DATA_DIR, "fig3a_nt_freq_by_period.csv"))
p3B <- make_freq_barplot(aa.geqpt01, "aa_label", "amino acid",
                          source_csv = file.path(SOURCE_DATA_DIR, "fig3b_aa_freq_by_period.csv"))

# --- Panel C: Non-synonymous fraction by frequency tier ---
# Shows that high-frequency mutations are enriched for non-synonymous changes
# (evidence of positive selection), while at low frequency the ratio reflects
# the ~55% baseline (neutral + purifying selection).
#
# NOTE: Both NT and AA must use the same >= 0.01% threshold for valid comparison.
# aa_first has 50,008 total mutations but only 5,656 at >= 0.01%.

nt_coding <- nt.geqpt01  # Already filtered to >= 0.01%

# Compute non-syn fraction at different frequency tiers for each time period
nonsyn_by_tier <- function(nt_df, aa_df, breaks, labels, tier_name, tier_filter_nt, tier_filter_aa) {
  nt_period <- nt_df |>
    filter(tier_filter_nt(per_abundance)) |>
    mutate(period = cut(collection_date, breaks = breaks, labels = labels,
                        include.lowest = TRUE, right = FALSE)) |>
    filter(!is.na(period)) |>
    count(period, name = "nt_new")

  aa_period <- aa_df |>
    filter(tier_filter_aa(per_abundance)) |>
    mutate(period = cut(collection_date, breaks = breaks, labels = labels,
                        include.lowest = TRUE, right = FALSE)) |>
    filter(!is.na(period)) |>
    count(period, name = "aa_new")

  full_join(nt_period, aa_period, by = "period") |>
    mutate(
      period = factor(period, levels = labels),
      nt_new = replace_na(nt_new, 0),
      aa_new = replace_na(aa_new, 0),
      pct_nonsyn = ifelse(nt_new >= 20, pmin(aa_new / nt_new, 1), NA_real_),
      tier = tier_name
    )
}

nonsyn_all <- nonsyn_by_tier(nt_coding, aa.geqpt01, FIG3_BREAKS, FIG3_LABELS,
                             ">= 0.01%", function(x) x >= 0.01, function(x) x >= 0.01)
nonsyn_pt1 <- nonsyn_by_tier(nt_coding, aa.geqpt01, FIG3_BREAKS, FIG3_LABELS,
                             ">= 0.1%", function(x) x >= 0.1, function(x) x >= 0.1)
nonsyn_1   <- nonsyn_by_tier(nt_coding, aa.geqpt01, FIG3_BREAKS, FIG3_LABELS,
                             ">= 1%", function(x) x >= 1, function(x) x >= 1)

nonsyn_combined <- bind_rows(nonsyn_all, nonsyn_pt1, nonsyn_1) |>
  mutate(tier = factor(tier, levels = c(">= 0.01%", ">= 0.1%", ">= 1%")))

p3C <- ggplot(nonsyn_combined |> filter(!is.na(pct_nonsyn)),
              aes(x = period, y = pct_nonsyn, fill = tier)) +
  geom_bar(stat = "identity", position = position_dodge2(width = 0.8, preserve = "single"),
           width = 0.7, color = "black", linewidth = 0.3) +
  coord_flip() +
  scale_x_discrete(drop = FALSE, limits = rev(FIG3_LABELS)) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1),
                     expand = expansion(mult = c(0, 0.02))) +
  scale_fill_manual(values = c(">= 0.01%" = "darkgoldenrod3",
                                ">= 0.1%" = "darkorange",
                                ">= 1%" = "firebrick")) +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "bottom") +
  labs(x = "", y = "", title = "Non-synonymous fraction\nby frequency tier",
       fill = "Threshold")

# --- Panel D: Spike vs rest of genome (AA mutations >= 0.1%) ---
# At >= 0.1%, Spike enrichment during vaccine era is clearly visible.
# Null expectation: Spike = 1,273 / 9,744 = 13.1% of total AA positions.

SPIKE_NULL <- 1273 / 9744  # Spike = 1,273 AA positions / 9,744 total AA positions = 13.1%

aa_spike_vs_rest <- aa.geqpt01 |>
  filter(per_abundance >= 0.1) |>
  mutate(
    gene_group = ifelse(region == "Spike", "Spike", "Rest"),
    period = cut(collection_date, breaks = FIG3_BREAKS,
                 labels = FIG3_LABELS, include.lowest = TRUE, right = FALSE)
  ) |>
  filter(!is.na(period)) |>
  count(period, gene_group) |>
  mutate(period = factor(period, levels = FIG3_LABELS)) |>
  group_by(period) |>
  mutate(pct = n / sum(n)) |>
  ungroup()

p3D <- ggplot(aa_spike_vs_rest, aes(x = period, y = pct, fill = gene_group)) +
  geom_bar(stat = "identity", width = 0.7, color = "black", linewidth = 0.3) +
  coord_flip() +
  scale_x_discrete(drop = FALSE, limits = rev(FIG3_LABELS)) +
  scale_y_continuous(labels = scales::percent,
                     expand = expansion(mult = c(0, 0.02))) +
  scale_fill_manual(values = c("Spike" = "steelblue", "Rest" = "grey70")) +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "bottom") +
  labs(x = "", y = "",
       title = "Spike vs rest\n(AA mutations >= 0.1%)",
       fill = "")

# Assemble Figure 3: single row, aligned timelines on y-axis
fig3 <- (p3A | p3B | p3C | p3D) +
  plot_annotation(tag_levels = "A")

ggsave(file.path(FIG_DIR, "figure3.pdf"), plot = fig3, width = 20, height = 8)
ggsave(file.path(FIG_DIR, "figure3.jpeg"), plot = fig3, dpi = 300, width = 20, height = 8)
msg("  Figure 3 saved!")

# --- Source data for Fig 3 C and D (A and B are written inside make_freq_barplot above) ---
fwrite(nonsyn_combined,  file.path(SOURCE_DATA_DIR, "fig3c_nonsyn_by_tier.csv"))
fwrite(aa_spike_vs_rest, file.path(SOURCE_DATA_DIR, "fig3d_spike_vs_rest.csv"))

# (figure3_alt removed — 4-bin version is the final figure 3)

# =============================================================================
# SUPPLEMENTARY FIGURE 1
# =============================================================================

msg("\n[6/6] Generating Supplementary Figure 1...")

# --- Panel A: Percentage of mutations in first 4 months vs rest ---
nt_3m <- nt.geqpt01 |>
  mutate(period = ifelse(collection_date <= as.Date("2020-03-31"),
                         "Dec-19 to Mar-20", "Apr-20 to Dec-25")) |>
  count(period, name = "n") |>
  mutate(type = factor("Nucleotide", levels = c("Nucleotide", "Amino acid")),
         pct = n / sum(n))

aa_3m <- aa.geqpt01 |>
  mutate(period = ifelse(collection_date <= as.Date("2020-03-31"),
                         "Dec-19 to Mar-20", "Apr-20 to Dec-25")) |>
  count(period, name = "n") |>
  mutate(type = factor("Amino acid", levels = c("Nucleotide", "Amino acid")),
         pct = n / sum(n))

combined_3m <- bind_rows(nt_3m, aa_3m)
combined_3m$period <- factor(combined_3m$period, levels = c("Dec-19 to Mar-20", "Apr-20 to Dec-25"))

pS1A <- ggplot(combined_3m, aes(x = period, y = pct, fill = type)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9),
           width = 0.8, color = "black", linewidth = 0.3) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1),
                     expand = expansion(mult = c(0, 0.02))) +
  scale_x_discrete(expand = expansion(add = 0.5),
                   labels = function(x) gsub(" to ", "\nto ", x)) +
  scale_fill_manual(values = c("Nucleotide" = "#7C9082", "Amino acid" = "#C0C0C0")) +
  theme_classic(base_size = 14) +
  theme(legend.position = "inside",
        legend.position.inside = c(0.98, 0.98),
        legend.justification = c(1, 1),
        axis.text.x = element_text(size = 10)) +
  labs(x = "", y = "Percentage of mutations", fill = "",
       title = "Mutations: first 4 months vs rest")

# --- Panel B: Comparison of documented infections ---
inf_3m <- infections |>
  mutate(period = ifelse(collection_date <= as.Date("2020-03-31"),
                         "Dec-19 to Mar-20", "Apr-20 to Dec-25")) |>
  group_by(period) |>
  summarise(total_cases = sum(new_cases, na.rm = TRUE), .groups = "drop") |>
  mutate(period = factor(period, levels = c("Dec-19 to Mar-20", "Apr-20 to Dec-25")))

pS1B <- ggplot(inf_3m, aes(x = period, y = total_cases, fill = period)) +
  geom_bar(stat = "identity", width = 0.4, color = "black", linewidth = 0.3) +
  scale_x_discrete(expand = expansion(add = 0.5),
                   labels = function(x) gsub(" to ", "\nto ", x)) +
  scale_y_log10(labels = scales::label_log(), expand = expansion(mult = c(0, 0.05))) +
  scale_fill_manual(values = c("Dec-19 to Mar-20" = "#7C9082", "Apr-20 to Dec-25" = "#C0C0C0")) +
  theme_classic(base_size = 14) +
  theme(legend.position = "none") +
  labs(x = "", y = "Number of documented infections (log scale)",
       title = "Infections by time period")

# --- Panel C: Correlation between cumulative cases and sequences ---
# Use FULL GISAID daily counts (17.1M sequences, 2,195 daily data points)
gisaid_daily_file <- file.path(DATA_DIR, "gisaid_daily_sequence_counts.csv")
gisaid_monthly_file <- file.path(DATA_DIR, "gisaid_monthly_sequence_counts.csv")

if (file.exists(gisaid_daily_file)) {
  seq_daily <- fread(gisaid_daily_file)
  seq_daily[, collection_date := as.Date(collection_date)]
  msg("  Using full GISAID daily counts (%s total sequences, %d days)",
      format(max(seq_daily$cum_seqs), big.mark = ","), nrow(seq_daily))
} else if (file.exists(gisaid_monthly_file)) {
  # Fallback to monthly
  seq_monthly <- fread(gisaid_monthly_file)
  msg("  WARNING: Daily GISAID counts not found, using monthly (%s total)",
      format(max(seq_monthly$cum_seqs), big.mark = ","))
}

# Daily cumulative infections
inf_daily <- copy(infections)[order(collection_date)]
inf_daily[, cum_cases := cumsum(new_cases)]

# Merge daily sequences with daily infections by date
daily_corr <- merge(seq_daily, inf_daily[, .(collection_date, cum_cases)],
                    by = "collection_date", all.x = TRUE)
# Forward-fill missing infection days (weekends/gaps in OWID)
daily_corr[, cum_cases := nafill(cum_cases, type = "locf")]
daily_corr <- daily_corr[!is.na(cum_cases)]
daily_corr[, cum_seqs_M := cum_seqs / 1e6]
daily_corr[, cum_cases_M := cum_cases / 1e6]

# Fit linear model on daily data
lm_fit <- lm(cum_cases_M ~ cum_seqs_M, data = daily_corr)
lm_coefs <- coef(lm_fit)
r_sq <- summary(lm_fit)$r.squared
eq_label <- sprintf("Linear Trendline: y = %g + %g x",
                     round(lm_coefs[1], 1), round(lm_coefs[2], 1))
rsq_label <- sprintf("R\u00b2 = %.2f", r_sq)

# Also keep monthly for backward compat (used by other sections)
if (file.exists(gisaid_monthly_file)) {
  seq_monthly <- fread(gisaid_monthly_file)
}
inf_monthly_cum <- infections[, .(cases = sum(new_cases, na.rm = TRUE)),
                              by = .(month = format(collection_date, "%Y-%m"))
][order(month)][, cum_cases := cumsum(cases)]
corr_data <- merge(seq_monthly, inf_monthly_cum, by = "month")
corr_data[, cum_seqs_M := cum_seqs / 1e6]
corr_data[, cum_cases_M := cum_cases / 1e6]

pS1C <- ggplot(daily_corr, aes(x = cum_seqs_M, y = cum_cases_M)) +
  geom_point(color = "steelblue", size = 0.8, alpha = 0.4) +
  geom_smooth(method = "lm", color = "red", se = TRUE, fill = "pink", alpha = 0.2) +
  annotate("text", x = max(daily_corr$cum_seqs_M) * 0.05,
           y = max(daily_corr$cum_cases_M) * 0.95,
           label = rsq_label, hjust = 0, size = 5, fontface = "bold") +
  annotate("text", x = max(daily_corr$cum_seqs_M) * 0.05,
           y = max(daily_corr$cum_cases_M) * 0.85,
           label = eq_label, hjust = 0, size = 4, color = "grey30") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.02)),
                     limits = c(0, NA)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02)),
                     limits = c(0, NA)) +
  theme_classic(base_size = 14) +
  labs(x = "Cumulative sequences deposited in GISAID\n(in millions)",
       y = "Cumulative infections worldwide (in millions)",
       title = "Cases vs sequences")

# --- Panel D: Vaccination timeline ---
# Note: OWID World aggregate starts at ~10% due to Chinese emergency-use backfill.
# Filter to start from 2021 when global rollout meaningfully began.
if (has_vacc && nrow(vacc) > 0) {
  vacc_plot <- vacc[collection_date >= as.Date("2021-01-01")]
  pS1D <- ggplot(vacc_plot, aes(x = collection_date, y = pct_fully_vaccinated)) +
    geom_line(color = "darkgreen", linewidth = 1) +
    scale_x_date(date_labels = "%b-%y",
                 breaks = seq.Date(as.Date("2021-01-01"), as.Date("2024-07-01"), by = "6 months")) +
    scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0, 0.02))) +
    theme_classic(base_size = 14) +
    labs(x = "Timeline", y = "% fully vaccinated (World)",
         title = "Percentage of world population\nfully vaccinated")
} else {
  pS1D <- ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Vaccination data\nnot available") +
    theme_void()
}

# --- Assemble Supplementary Figure 1 (panels A-D only) ---
# Former panels E and F moved to Figure 2 (panel D) / dropped; Fisher variants removed.
# S1C and S1D both get width 2 (equal); S1A and S1B slightly wider than 1 so their
# x-axis labels don't crowd; add per-panel plot.margin for visible whitespace.
pS1A <- pS1A + theme(plot.margin = margin(t = 5, r = 12, b = 5, l = 12))
pS1B <- pS1B + theme(plot.margin = margin(t = 5, r = 12, b = 5, l = 12))
pS1C <- pS1C + theme(plot.margin = margin(t = 5, r = 12, b = 5, l = 12))
pS1D <- pS1D + theme(plot.margin = margin(t = 5, r = 12, b = 5, l = 12))

supp_fig1 <- (pS1A + pS1B + pS1C + pS1D + plot_layout(widths = c(1.3, 1.3, 2, 2))) +
  plot_annotation(tag_levels = "A")

ggsave(file.path(FIG_DIR, "supplementary_figure1.pdf"), plot = supp_fig1, width = 22, height = 7)
ggsave(file.path(FIG_DIR, "supplementary_figure1.jpeg"), plot = supp_fig1, dpi = 300, width = 22, height = 7)
msg("  Supplementary Figure 1 (A-D only) saved!")

# --- Source data for Supp Fig 1 ---
fwrite(as.data.table(combined_3m),
       file.path(SOURCE_DATA_DIR, "suppfig1a_first4m_vs_rest.csv"))
fwrite(as.data.table(inf_3m),
       file.path(SOURCE_DATA_DIR, "suppfig1b_infections_by_period.csv"))
fwrite(daily_corr[, .(collection_date, cum_seqs, cum_cases, cum_seqs_M, cum_cases_M)],
       file.path(SOURCE_DATA_DIR, "suppfig1c_daily_seqs_vs_cases.csv"))
if (has_vacc && exists("vacc_plot") && nrow(vacc_plot) > 0) {
  fwrite(vacc_plot, file.path(SOURCE_DATA_DIR, "suppfig1d_vaccination.csv"))
}


# =============================================================================
# Print Fisher's Exact Test Results (for manuscript text)
# =============================================================================

msg("\n" |> rep(1) |> paste(collapse = ""))
msg("=" |> rep(70) |> paste(collapse = ""))
msg("Statistical Tests for Manuscript")
msg("=" |> rep(70) |> paste(collapse = ""))

# Fisher's test at the Figure 2 heuristic (>= 10% vs <10%, <=FIG2_DAYS_CUTOFF vs >FIG2_DAYS_CUTOFF)
tryCatch({
  nt_fisher_df <- nt_early |>
    as.data.frame() |>
    mutate(
      freq_cat = ifelse(per_abundance >= 10, ">=10%", "<10%"),
      days_cat = ifelse(days_to_second <= FIG2_DAYS_CUTOFF,
                        sprintf("<=%d days", FIG2_DAYS_CUTOFF),
                        sprintf(">%d days", FIG2_DAYS_CUTOFF))
    )
  fisher_nt <- fisher.test(table(nt_fisher_df$freq_cat, nt_fisher_df$days_cat))
  msg("  NT Fisher's exact test (>=10%% vs <10%%, %dd heuristic): p = %g",
      FIG2_DAYS_CUTOFF, fisher_nt$p.value)
}, error = function(e) msg("  NT Fisher's test error: %s", e$message))

tryCatch({
  aa_fisher_df <- aa_early |>
    as.data.frame() |>
    mutate(
      freq_cat = ifelse(per_abundance >= 10, ">=10%", "<10%"),
      days_cat = ifelse(days_to_second <= FIG2_DAYS_CUTOFF,
                        sprintf("<=%d days", FIG2_DAYS_CUTOFF),
                        sprintf(">%d days", FIG2_DAYS_CUTOFF))
    )
  fisher_aa <- fisher.test(table(aa_fisher_df$freq_cat, aa_fisher_df$days_cat))
  msg("  AA Fisher's exact test (>=10%% vs <10%%, %dd heuristic): p = %g",
      FIG2_DAYS_CUTOFF, fisher_aa$p.value)
}, error = function(e) msg("  AA Fisher's test error: %s", e$message))

# Spearman correlations
tryCatch({
  sp_nt <- cor.test(nt_early$days_to_second, nt_early$per_abundance, method = "spearman")
  msg("  NT Spearman rho = %.3f, p = %g", sp_nt$estimate, sp_nt$p.value)
}, error = function(e) msg("  NT Spearman error: %s", e$message))

tryCatch({
  sp_aa <- cor.test(aa_early$days_to_second, aa_early$per_abundance, method = "spearman")
  msg("  AA Spearman rho = %.3f, p = %g", sp_aa$estimate, sp_aa$p.value)
}, error = function(e) msg("  AA Spearman error: %s", e$message))

# =============================================================================
# Summary
# =============================================================================

msg("\n" |> rep(1) |> paste(collapse = ""))
msg("=" |> rep(70) |> paste(collapse = ""))
msg("All figures generated!")
msg("=" |> rep(70) |> paste(collapse = ""))
msg("Outputs saved to: %s", FIG_DIR)
msg("\nFigures:")
msg("  - figure1.pdf / .jpeg")
msg("  - figure2.pdf / .jpeg  (Jan 15 cutoff; panel D >=1%%, <=28d heuristic)")
msg("  - figure3.pdf / .jpeg")
msg("  - supplementary_figure1.pdf / .jpeg")

msg("\nSession Info:")
sessionInfo()
