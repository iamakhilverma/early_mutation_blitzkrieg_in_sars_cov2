#!/usr/bin/env python3
"""
Audit: Find sequences with collection dates before their variant's emergence.

Uses variant_emergence_dates.csv (curated from WHO/GISAID/Nextstrain sources)
to flag sequences that claim to be a variant before that variant could have
existed — strong evidence of misdated collection or misassignment.

Handles Pango aliases:
  Q.*   = B.1.1.7.* (Alpha sub-lineages)
  AY.*  = B.1.617.2.* (Delta sub-lineages)
  BA.*  = B.1.1.529.* (Omicron sub-lineages, handled via sub-lineage entries)

Output: list of accessions + reasons for exclusion
"""

import csv
import sys
from collections import defaultdict
from datetime import datetime

PROJECT_DIR = "/sc/arion/projects/Tsankov_Normal_Lung/users/akhil/projects/av/projects/mutation_tracking"

# ── Load variant emergence dates ──
emergence = {}  # lineage_prefix → (impossible_before_date, label)

with open(f"{PROJECT_DIR}/raw_data/variant_emergence_dates.csv") as f:
    reader = csv.DictReader(f)
    for row in reader:
        label = row["who_label"]
        impossible = row["impossible_before_date"]
        pango = row["pango_lineage"]

        # Handle semicolon-separated lineages (e.g., Epsilon: B.1.427;B.1.429)
        for lin in pango.split(";"):
            lin = lin.strip()
            emergence[lin] = (impossible, label)

# Add Pango aliases manually
# Q.* = B.1.1.7.* (Alpha sub-lineages)
emergence["Q"] = (emergence["B.1.1.7"][0], "Alpha (Q alias)")
# AY.* = B.1.617.2.* (Delta sub-lineages)
emergence["AY"] = (emergence["B.1.617.2"][0], "Delta (AY alias)")

print(f"Loaded {len(emergence)} lineage prefixes to check")
for lin, (date, label) in sorted(emergence.items()):
    print(f"  {lin:<15s} impossible before {date}  ({label})")

# ── Load existing exclusion list ──
existing_excluded = set()
try:
    with open(f"{PROJECT_DIR}/processed_data/exclude_additional.txt") as f:
        for line in f:
            line = line.strip()
            if line:
                existing_excluded.add(line)
    print(f"\nExisting exclusion list: {len(existing_excluded)} accessions")
except FileNotFoundError:
    print("\nNo existing exclusion list found")

# ── Scan sequence_list.csv ──
print(f"\nScanning sequence_list.csv...")

flagged = []
n_scanned = 0
lineage_counts = defaultdict(int)

with open(f"{PROJECT_DIR}/processed_data/sequence_list.csv") as f:
    reader = csv.DictReader(f)
    for row in reader:
        n_scanned += 1
        if n_scanned % 1_000_000 == 0:
            print(f"  ... {n_scanned:,} sequences scanned")

        acc = row["AccessionID"]
        date_str = row["collection_date"]
        lineage = row["lineage"]

        if not lineage or lineage == "unassigned" or not date_str:
            continue

        # Check each emergence rule
        for prefix, (impossible_before, label) in emergence.items():
            # Match: lineage == prefix exactly, or lineage starts with prefix.
            if lineage == prefix or lineage.startswith(prefix + "."):
                coll_date = datetime.strptime(date_str, "%Y-%m-%d").date()
                cutoff = datetime.strptime(impossible_before, "%Y-%m-%d").date()

                if coll_date < cutoff:
                    already = acc in existing_excluded
                    flagged.append({
                        "accession": acc,
                        "collection_date": date_str,
                        "lineage": lineage,
                        "variant_label": label,
                        "impossible_before": impossible_before,
                        "already_excluded": already,
                    })
                    lineage_counts[label] += 1
                break  # Only need to match one rule per sequence

print(f"\nScanned {n_scanned:,} sequences total")
print(f"Flagged {len(flagged)} sequences with impossible dates")

# ── Summary by variant ──
print(f"\n{'Variant':<30s} {'Count':>6s}  {'New':>4s}  Impossible before")
print("-" * 80)
new_total = 0
for label in sorted(lineage_counts.keys()):
    count = lineage_counts[label]
    new_count = sum(1 for f in flagged if f["variant_label"] == label and not f["already_excluded"])
    new_total += new_count
    # Find the impossible_before for this label
    imp_date = [f["impossible_before"] for f in flagged if f["variant_label"] == label][0]
    print(f"  {label:<28s} {count:>6d}  {new_count:>4d}  {imp_date}")

print(f"\n  {'TOTAL':<28s} {len(flagged):>6d}  {new_total:>4d}")

# ── Write new flagged accessions ──
new_flagged = [f for f in flagged if not f["already_excluded"]]

if new_flagged:
    outpath = f"{PROJECT_DIR}/processed_data/audit_misdated_variants.csv"
    with open(outpath, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["accession", "collection_date", "lineage",
                                                "variant_label", "impossible_before", "already_excluded"])
        writer.writeheader()
        for row in flagged:
            writer.writerow(row)
    print(f"\nFull details written to: {outpath}")

    new_ids_path = f"{PROJECT_DIR}/processed_data/audit_new_exclusions.txt"
    with open(new_ids_path, "w") as f:
        for row in new_flagged:
            f.write(row["accession"] + "\n")
    print(f"New accession IDs to exclude: {new_ids_path} ({len(new_flagged)} IDs)")
else:
    print("\nNo NEW misdated sequences found beyond existing exclusion list.")

# ── Show some examples ──
if flagged:
    print(f"\nExample flagged sequences:")
    for f in flagged[:20]:
        status = " [ALREADY EXCLUDED]" if f["already_excluded"] else " [NEW]"
        print(f"  {f['accession']}  {f['collection_date']}  {f['lineage']:<15s}  "
              f"{f['variant_label']:<20s}  before {f['impossible_before']}{status}")
    if len(flagged) > 20:
        print(f"  ... and {len(flagged) - 20} more")
