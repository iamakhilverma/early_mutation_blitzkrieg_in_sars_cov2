#!/usr/bin/env python3
"""
Audit: flag sequences whose collection date precedes their Pango lineage's
earliest documented emergence — strong evidence of a misdated or misassigned
sample.

This script produced the `audit_misdated_variants.csv` and the 11 sub-lineage
accessions appended to `exclude_additional.txt` that ship with this repository.
You do NOT need to re-run it unless you refresh the GISAID metadata — the
bundled exclusion list is already applied by scripts/02_process_nextclade.R.

If you do want to re-run it, you need a `sequence_list.csv` that has columns
(AccessionID, collection_date, lineage). Pass its path via --sequence-list.
That file is not included in this repo (it's per-sequence metadata derived
from GISAID, whose terms prohibit redistribution of per-sequence data).

Usage
-----
    python3 scripts/audit_variant_dates.py \\
        --sequence-list    path/to/sequence_list.csv \\
        --emergence-dates  reference_data/variant_emergence_dates.csv \\
        --existing-exclude processed_data/exclude_additional.txt \\
        --out-csv          processed_data/audit_misdated_variants.csv \\
        --out-txt          processed_data/audit_new_exclusions.txt

Pango alias handling
--------------------
Q.*   → B.1.1.7.*    (Alpha sub-lineages)
AY.*  → B.1.617.2.*  (Delta sub-lineages)
BA.*  → B.1.1.529.*  (Omicron sub-lineages — covered via the sub-lineage rows
                      in variant_emergence_dates.csv)
"""

import argparse
import csv
import sys
from collections import defaultdict
from datetime import datetime


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--sequence-list",    required=True,
                   help="CSV with AccessionID, collection_date, lineage columns.")
    p.add_argument("--emergence-dates",  default="reference_data/variant_emergence_dates.csv")
    p.add_argument("--existing-exclude", default="processed_data/exclude_additional.txt",
                   help="Optional — accessions already in this list are flagged as 'already_excluded'.")
    p.add_argument("--out-csv", default="processed_data/audit_misdated_variants.csv")
    p.add_argument("--out-txt", default="processed_data/audit_new_exclusions.txt")
    return p.parse_args()


def load_emergence(path):
    emergence = {}  # lineage_prefix -> (impossible_before_date, label)
    with open(path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            label = row["who_label"]
            impossible = row["impossible_before_date"]
            pango = row["pango_lineage"]
            for lin in pango.split(";"):
                lin = lin.strip()
                if lin:
                    emergence[lin] = (impossible, label)
    # Manual aliases
    if "B.1.1.7" in emergence:
        emergence["Q"] = (emergence["B.1.1.7"][0], "Alpha (Q alias)")
    if "B.1.617.2" in emergence:
        emergence["AY"] = (emergence["B.1.617.2"][0], "Delta (AY alias)")
    return emergence


def load_existing(path):
    try:
        with open(path) as f:
            return {line.strip() for line in f if line.strip()}
    except FileNotFoundError:
        return set()


def main():
    args = parse_args()

    emergence = load_emergence(args.emergence_dates)
    print(f"Loaded {len(emergence)} lineage prefixes to check")
    for lin, (date, label) in sorted(emergence.items()):
        print(f"  {lin:<15s} impossible before {date}  ({label})")

    existing_excluded = load_existing(args.existing_exclude)
    print(f"\nExisting exclusion list: {len(existing_excluded)} accessions")

    print(f"\nScanning {args.sequence_list} ...")
    flagged = []
    n_scanned = 0
    lineage_counts = defaultdict(int)

    with open(args.sequence_list) as f:
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

            for prefix, (impossible_before, label) in emergence.items():
                if lineage == prefix or lineage.startswith(prefix + "."):
                    coll_date = datetime.strptime(date_str, "%Y-%m-%d").date()
                    cutoff = datetime.strptime(impossible_before, "%Y-%m-%d").date()
                    if coll_date < cutoff:
                        flagged.append({
                            "accession": acc,
                            "collection_date": date_str,
                            "lineage": lineage,
                            "variant_label": label,
                            "impossible_before": impossible_before,
                            "already_excluded": acc in existing_excluded,
                        })
                        lineage_counts[label] += 1
                    break

    print(f"\nScanned {n_scanned:,} sequences total")
    print(f"Flagged {len(flagged)} sequences with impossible dates")

    print(f"\n{'Variant':<30s} {'Count':>6s}  {'New':>4s}  Impossible before")
    print("-" * 80)
    new_total = 0
    for label in sorted(lineage_counts.keys()):
        count = lineage_counts[label]
        new_count = sum(1 for f in flagged if f["variant_label"] == label and not f["already_excluded"])
        new_total += new_count
        imp_date = next(f["impossible_before"] for f in flagged if f["variant_label"] == label)
        print(f"  {label:<28s} {count:>6d}  {new_count:>4d}  {imp_date}")
    print(f"\n  {'TOTAL':<28s} {len(flagged):>6d}  {new_total:>4d}")

    new_flagged = [f for f in flagged if not f["already_excluded"]]

    if not new_flagged:
        print("\nNo NEW misdated sequences found beyond existing exclusion list.")
        return

    with open(args.out_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["accession", "collection_date", "lineage",
                                               "variant_label", "impossible_before", "already_excluded"])
        writer.writeheader()
        for row in flagged:
            writer.writerow(row)
    print(f"\nFull details written to: {args.out_csv}")

    with open(args.out_txt, "w") as f:
        for row in new_flagged:
            f.write(row["accession"] + "\n")
    print(f"New accession IDs to exclude: {args.out_txt} ({len(new_flagged)} IDs)")

    if flagged:
        print("\nExample flagged sequences:")
        for row in flagged[:20]:
            status = " [ALREADY EXCLUDED]" if row["already_excluded"] else " [NEW]"
            print(f"  {row['accession']}  {row['collection_date']}  {row['lineage']:<15s}  "
                  f"{row['variant_label']:<20s}  before {row['impossible_before']}{status}")
        if len(flagged) > 20:
            print(f"  ... and {len(flagged) - 20} more")


if __name__ == "__main__":
    sys.exit(main())
