#!/usr/bin/env Rscript
# Rebuild processed_data/world_infections.csv as the NEJM (Dec 2019 – Jan 21 2020)
# + OWID daily-snapshot (Jan 22 2020+) merge. This is the same merge that
# scripts/02_process_nextclade.R Step 9 performs; this standalone helper recreates
# the file without re-running Phase 1, e.g. if you already have processed_data/
# from an older pipeline run that used a different infections source.
suppressPackageStartupMessages(library(data.table))

PROJECT_DIR <- if (length(commandArgs(trailingOnly = TRUE)) >= 1)
  commandArgs(trailingOnly = TRUE)[1] else "."

nejm <- fread(file.path(PROJECT_DIR, "reference_data", "nejm_early_cases.csv"))[
  , .(collection_date = as.Date(collection_date), new_cases = as.numeric(new_cases))]

# Use the historical daily-granularity OWID snapshot (kept in reference_data/).
# The current owid-covid-data.csv has weekly-aggregated values for early 2020
# that would distort Fig 1C/D.
owid_post <- fread(file.path(PROJECT_DIR, "reference_data", "world_infections_owid_raw.csv"))[
  , .(collection_date = as.Date(collection_date), new_cases = as.numeric(new_cases))]
owid_post <- owid_post[collection_date >= as.Date("2020-01-22")]

world <- rbindlist(list(nejm, owid_post), use.names = TRUE)
setorder(world, collection_date)
world[, total_cases := cumsum(new_cases)]
world[, new_cases_smoothed := round(frollmean(new_cases, 7, na.rm = TRUE,
                                              align = "right"), 2)]
world[is.na(new_cases_smoothed), new_cases_smoothed := 0]
world <- world[collection_date >= as.Date("2019-12-01") &
                collection_date <= as.Date("2025-12-31")]

out <- file.path(PROJECT_DIR, "processed_data", "world_infections.csv")
fwrite(world, out)
cat(sprintf("Wrote %s (%d rows)\n", out, nrow(world)))
cat("Sanity check (Dec 2019 / early Jan 2020):\n")
print(world[collection_date <= as.Date("2020-01-25")])
