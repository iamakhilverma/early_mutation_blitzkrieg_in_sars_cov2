#!/usr/bin/env python3
"""
Fetch COVID-19 case data and vaccination data from the OWID catalog.
Saves raw extracts as CSVs for the mutation tracking pipeline.

Usage:
    python3 scripts/fetch_owid_data.py [out_dir]

Defaults to writing into ./processed_data relative to the current working
directory (i.e., run this from the project root).
"""

import os
import sys

import pandas as pd
from owid.catalog import search

OUT_DIR = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.getcwd(), "processed_data")
os.makedirs(OUT_DIR, exist_ok=True)

# ── 1. Search for COVID tables ────────────────────────────────────────────────
print("=" * 70)
print("Searching OWID catalog for COVID tables...")
print("=" * 70)
results = search(namespace='covid', kind='table')
print(f"\nFound {len(results)} tables. Listing all:\n")
for i, r in enumerate(results):
    print(f"  [{i}] {r.path}")

# ── 2. Fetch cases_deaths table ──────────────────────────────────────────────
print("\n" + "=" * 70)
print("Fetching cases_deaths table...")
print("=" * 70)

# Find the cases_deaths table by path
cases_idx = None
vacc_idx = None
for i, r in enumerate(results):
    if "cases_deaths/cases_deaths" in r.path:
        cases_idx = i
    if "vaccinations_global/vaccinations_global" in r.path:
        vacc_idx = i

if cases_idx is None:
    raise RuntimeError("Could not find cases_deaths table in OWID catalog")
if vacc_idx is None:
    raise RuntimeError("Could not find vaccinations_global table in OWID catalog")

print(f"  cases_deaths index: {cases_idx} -> {results[cases_idx].path}")
print(f"  vaccinations_global index: {vacc_idx} -> {results[vacc_idx].path}")

tb_cases = results[cases_idx].fetch()
print(f"\nRaw table shape: {tb_cases.shape}")
print(f"Index names: {tb_cases.index.names}")
print(f"Columns: {list(tb_cases.columns)}")

# Filter for World
if "country" in tb_cases.index.names:
    tb_world = tb_cases.xs("World", level="country")
else:
    tb_world = tb_cases[tb_cases["country"] == "World"]

tb_world = tb_world.reset_index()
print(f"\nAfter filtering for 'World': {tb_world.shape}")
print(f"Columns available: {list(tb_world.columns)}")

# Extract needed columns
keep_cols_cases = ["date", "new_cases", "total_cases"]
# Check which columns actually exist
for c in keep_cols_cases:
    if c not in tb_world.columns:
        print(f"  WARNING: column '{c}' not found. Available: {list(tb_world.columns)}")

df_cases = tb_world[keep_cols_cases].copy()

# Convert date to datetime if needed, then filter range
df_cases["date"] = pd.to_datetime(df_cases["date"])
df_cases = df_cases[(df_cases["date"] >= "2020-01-22") & (df_cases["date"] <= "2025-12-31")]
df_cases = df_cases.dropna(subset=["new_cases", "total_cases"], how="all")
df_cases = df_cases.sort_values("date").reset_index(drop=True)

# Rename columns
df_cases.rename(columns={"date": "collection_date"}, inplace=True)

# Save
out_cases = f"{OUT_DIR}/world_infections_owid_raw.csv"
df_cases.to_csv(out_cases, index=False)
print(f"\nSaved to: {out_cases}")
print(f"Total rows: {len(df_cases)}")
print(f"Date range: {df_cases['collection_date'].min()} to {df_cases['collection_date'].max()}")

print("\nFirst 30 rows:")
print(df_cases.head(30).to_string(index=False))
print("\nLast 10 rows:")
print(df_cases.tail(10).to_string(index=False))

# ── 3. Fetch vaccinations_global table ────────────────────────────────────────
print("\n" + "=" * 70)
print("Fetching vaccinations_global table...")
print("=" * 70)

tb_vacc = results[vacc_idx].fetch()
print(f"\nRaw table shape: {tb_vacc.shape}")
print(f"Index names: {tb_vacc.index.names}")
print(f"Columns: {list(tb_vacc.columns)}")

# Filter for World
if "country" in tb_vacc.index.names:
    tb_vacc_world = tb_vacc.xs("World", level="country")
else:
    tb_vacc_world = tb_vacc[tb_vacc["country"] == "World"]

tb_vacc_world = tb_vacc_world.reset_index()
print(f"\nAfter filtering for 'World': {tb_vacc_world.shape}")
print(f"Columns available: {list(tb_vacc_world.columns)}")

# Extract needed columns
if "people_fully_vaccinated_per_hundred" in tb_vacc_world.columns:
    vacc_col = "people_fully_vaccinated_per_hundred"
else:
    # Try to find the right column
    candidates = [c for c in tb_vacc_world.columns if "fully" in c.lower() and "hundred" in c.lower()]
    if candidates:
        vacc_col = candidates[0]
        print(f"  Using column: {vacc_col}")
    else:
        print(f"  Available columns: {list(tb_vacc_world.columns)}")
        raise RuntimeError("Could not find people_fully_vaccinated_per_hundred column")

df_vacc = tb_vacc_world[["date", vacc_col]].copy()
df_vacc["date"] = pd.to_datetime(df_vacc["date"])
df_vacc = df_vacc.dropna(subset=[vacc_col])
df_vacc = df_vacc.sort_values("date").reset_index(drop=True)

# Rename columns
df_vacc.rename(columns={"date": "collection_date", vacc_col: "pct_fully_vaccinated"}, inplace=True)

# Save
out_vacc = f"{OUT_DIR}/vaccination_data_owid_raw.csv"
df_vacc.to_csv(out_vacc, index=False)
print(f"\nSaved to: {out_vacc}")
print(f"Total rows: {len(df_vacc)}")
print(f"Date range: {df_vacc['collection_date'].min()} to {df_vacc['collection_date'].max()}")

print("\nFirst 30 rows:")
print(df_vacc.head(30).to_string(index=False))
print("\nLast 10 rows:")
print(df_vacc.tail(10).to_string(index=False))

print("\n" + "=" * 70)
print("DONE. Both files saved.")
print("=" * 70)
