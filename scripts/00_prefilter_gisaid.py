#!/usr/bin/env python3
"""
=============================================================================
Pre-filter GISAID FASTA Using Metadata TSV
=============================================================================

Reads the GISAID bulk metadata TSV to identify clean sequences, then
filters the FASTA to keep only those. Also downloads the Nextstrain
exclusion list to remove contaminated/problematic sequences (including
the "ML" multiple-lineage sequences that GISAID flags in its web UI
but does NOT include in the bulk metadata download).

GISAID metadata columns (verified Feb 2026 bulk download):
  [0]  Virus name
  [1]  Last vaccinated
  [2]  Passage details/history
  [3]  Type
  [4]  Accession ID                 ← used
  [5]  Collection date              ← used (filter)
  [6]  Location
  [7]  Additional location information
  [8]  Sequence length              ← used (filter)
  [9]  Host                         ← used (filter)
  [10] Patient age
  [11] Gender
  [12] Clade
  [13] Pango lineage
  [14] Pango version
  [15] Variant
  [16] AA Substitutions             (not used — we get AA from Nextclade)
  [17] Submission date
  [18] Is reference?
  [19] Is complete?                 ← used (filter)
  [20] Is high coverage?
  [21] Is low coverage?             ← used (filter)
  [22] N-Content                    ← used (filter)
  [23] GC-Content

Filters applied:
  1. Host = Human
  2. Is complete? = True
  3. Is low coverage? ≠ True
  4. Collection date is YYYY-MM-DD (complete, not partial)
  5. Collection date within specified range
  6. N-Content ≤ 5%
  7. Sequence length ≥ 29000 bp
  8. Nextstrain exclusion list (removes known contaminants, duplicates,
     and sequences flagged for multiple lineages / under investigation)

Inputs:
  - metadata_tsv_*.tar.xz  OR  extracted metadata.tsv
  - sequences_fasta_*.tar.xz  OR  extracted .fasta
  - Nextstrain exclude list (auto-downloaded)

Outputs (in --outdir):
  - filtered_sequences.fasta    (clean sequences only)
  - filtered_metadata.tsv       (metadata for kept sequences)
  - prefilter_report.txt        (filtering summary)
  - pipeline_config.txt         (date range for downstream scripts)

Usage:
  python3 scripts/00_prefilter_gisaid.py \\
      --metadata raw_data/metadata.tsv \\
      --fasta raw_data/sequences.fasta \\
      --start-date 2019-12-01 \\
      --end-date 2025-12-31 \\
      --outdir filtered_data

Author: Akhil Kumar
=============================================================================
"""

import argparse
import gzip
import io
import lzma
import os
import re
import subprocess
import sys
import tarfile
import time
from collections import Counter
from datetime import datetime


def parse_args():
    p = argparse.ArgumentParser(
        description="Pre-filter GISAID FASTA using metadata TSV"
    )
    p.add_argument("--metadata", required=True,
                   help="Path to metadata .tar.xz or .tsv")
    p.add_argument("--fasta", required=True,
                   help="Path to FASTA .tar.xz or .fasta")
    p.add_argument("--start-date", default="2019-12-01",
                   help="Collection date start (default: 2019-12-01)")
    p.add_argument("--end-date", default="2025-12-31",
                   help="Collection date end (default: 2025-12-31)")
    p.add_argument("--max-n-content", type=float, default=0.05,
                   help="Max N-content fraction (default: 0.05 = 5%%)")
    p.add_argument("--min-length", type=int, default=29000,
                   help="Min sequence length in bp (default: 29000)")
    p.add_argument("--outdir", default="filtered_data",
                   help="Output directory (default: filtered_data)")
    p.add_argument("--nextstrain-exclude", default=None,
                   help="Path to Nextstrain exclude.txt (auto-downloaded if not provided)")
    p.add_argument("--skip-nextstrain-exclude", action="store_true",
                   help="Skip Nextstrain exclusion list entirely")
    p.add_argument("--inspect-only", action="store_true",
                   help="Print metadata columns + sample rows and exit")
    return p.parse_args()


def open_from_tarxz_or_plain(filepath):
    """Open a .tar.xz, .gz, .xz, or plain text file. Returns text file handle."""
    if filepath.endswith(".tar.xz") or filepath.endswith(".tar.gz"):
        mode = "r:xz" if filepath.endswith(".tar.xz") else "r:gz"
        tar = tarfile.open(filepath, mode)
        for member in tar.getmembers():
            name = member.name.lower()
            if name.endswith(".tsv") or name.endswith(".fasta") or name.endswith(".fa"):
                print(f"  Extracted from archive: {member.name} "
                      f"({member.size / 1e9:.2f} GB)")
                f = tar.extractfile(member)
                return io.TextIOWrapper(f, encoding="utf-8", errors="replace")
        raise FileNotFoundError(f"No TSV/FASTA in archive: {filepath}")
    elif filepath.endswith(".gz"):
        return io.TextIOWrapper(
            gzip.open(filepath, "rb"), encoding="utf-8", errors="replace")
    elif filepath.endswith(".xz"):
        return io.TextIOWrapper(
            lzma.open(filepath, "rb"), encoding="utf-8", errors="replace")
    else:
        return open(filepath, "r", encoding="utf-8", errors="replace")


def download_nextstrain_exclude(outdir):
    """Download Nextstrain exclusion list. Returns path or None."""
    url = "https://raw.githubusercontent.com/nextstrain/ncov/master/defaults/exclude.txt"
    outpath = os.path.join(outdir, "nextstrain_exclude.txt")

    for cmd in [["curl", "-sL", "-o", outpath, url],
                ["wget", "-q", "-O", outpath, url]]:
        try:
            subprocess.run(cmd, check=True, timeout=60)
            lines = open(outpath).readlines()
            entries = [l.strip() for l in lines
                       if l.strip() and not l.startswith("#")]
            print(f"  ✓ Nextstrain exclude list: {len(entries)} entries")
            return outpath
        except (subprocess.CalledProcessError, FileNotFoundError, OSError):
            continue
    print("  ⚠ Could not download Nextstrain exclude list (non-critical)")
    return None


def parse_nextstrain_exclude(filepath):
    """Parse exclude.txt → set of virus names (strain IDs)."""
    names = set()
    if filepath and os.path.exists(filepath):
        for line in open(filepath):
            line = line.strip()
            if line and not line.startswith("#"):
                # Some lines have comments after the name
                name = line.split("#")[0].strip()
                if name:
                    names.add(name)
    return names


def is_valid_complete_date(date_str, start_date, end_date):
    """Check date is YYYY-MM-DD and within range. Returns (ok, reason)."""
    if not date_str or date_str in ("", "?", "unknown", "NA"):
        return False, "missing_date"
    if not re.match(r"^\d{4}-\d{2}-\d{2}$", date_str):
        return False, "incomplete_date"
    try:
        d = datetime.strptime(date_str, "%Y-%m-%d").date()
        if d < start_date or d > end_date:
            return False, "out_of_date_range"
        return True, None
    except ValueError:
        return False, "invalid_date"


def main():
    args = parse_args()
    start_date = datetime.strptime(args.start_date, "%Y-%m-%d").date()
    end_date = datetime.strptime(args.end_date, "%Y-%m-%d").date()
    os.makedirs(args.outdir, exist_ok=True)

    # =========================================================================
    # INSPECT MODE
    # =========================================================================
    if args.inspect_only:
        print("=" * 70)
        print("INSPECT MODE — showing metadata structure")
        print("=" * 70)
        fh = open_from_tarxz_or_plain(args.metadata)
        header = fh.readline().rstrip("\n\r")
        cols = header.split("\t")
        print(f"\nColumns ({len(cols)}):")
        for i, c in enumerate(cols):
            print(f"  [{i:2d}] {c}")
        print(f"\nFirst 3 data rows:")
        for row_num in range(3):
            line = fh.readline()
            if not line:
                break
            fields = line.rstrip("\n\r").split("\t")
            print(f"\n  --- Row {row_num + 1} ---")
            for i, c in enumerate(cols):
                val = fields[i] if i < len(fields) else "(missing)"
                # Truncate long values
                if len(val) > 80:
                    val = val[:77] + "..."
                print(f"  [{i:2d}] {c}: {val}")
        fh.close()
        return

    # =========================================================================
    # STEP 0: Download Nextstrain exclusion list
    # =========================================================================
    print("=" * 70)
    print("STEP 0: Nextstrain exclusion list")
    print("=" * 70)
    if args.skip_nextstrain_exclude:
        print("  Skipped (--skip-nextstrain-exclude)")
        exclude_names = set()
    elif args.nextstrain_exclude:
        print(f"  Using provided file: {args.nextstrain_exclude}")
        exclude_names = parse_nextstrain_exclude(args.nextstrain_exclude)
        print(f"  Entries: {len(exclude_names)}")
    else:
        # Check for pre-downloaded file first (compute nodes have no internet)
        local_exclude = os.path.join(args.outdir, "nextstrain_exclude.txt")
        raw_exclude = os.path.join(os.path.dirname(args.metadata), "nextstrain_exclude.txt")
        if os.path.exists(local_exclude):
            print(f"  Using existing file: {local_exclude}")
            exclude_names = parse_nextstrain_exclude(local_exclude)
            print(f"  Entries: {len(exclude_names)}")
        elif os.path.exists(raw_exclude):
            print(f"  Using existing file: {raw_exclude}")
            exclude_names = parse_nextstrain_exclude(raw_exclude)
            print(f"  Entries: {len(exclude_names)}")
        else:
            print("  Attempting download from GitHub...")
            exclude_path = download_nextstrain_exclude(args.outdir)
            exclude_names = parse_nextstrain_exclude(exclude_path)
            if not exclude_names:
                print("  ⚠ No exclude list available. ML/contaminated sequences")
                print("    will NOT be filtered. To fix, run on login node:")
                print(f"    curl -sL https://raw.githubusercontent.com/nextstrain/ncov/master/defaults/exclude.txt -o {local_exclude}")

    # =========================================================================
    # STEP 1: Read and filter metadata
    # =========================================================================
    print()
    print("=" * 70)
    print("STEP 1: Filtering metadata")
    print("=" * 70)
    print(f"  File: {args.metadata}")
    print(f"  Date range: {args.start_date} to {args.end_date}")
    print(f"  Min length: {args.min_length} bp")
    print(f"  Max N-content: {args.max_n_content}")
    print()

    t0 = time.time()
    fh = open_from_tarxz_or_plain(args.metadata)

    # Read and validate header
    header = fh.readline().rstrip("\n\r")
    cols = header.split("\t")
    print(f"  Columns: {len(cols)}")

    # Build column index by name (exact match)
    col_map = {c.strip(): i for i, c in enumerate(cols)}

    # Define required columns with their expected names
    # (based on verified GISAID bulk download Feb 2026)
    REQUIRED = {
        "accession_id":    "Accession ID",
        "collection_date": "Collection date",
        "host":            "Host",
        "virus_name":      "Virus name",
    }
    OPTIONAL = {
        "is_complete":      "Is complete?",
        "is_high_coverage": "Is high coverage?",
        "is_low_coverage":  "Is low coverage?",
        "n_content":        "N-Content",
        "seq_length":       "Sequence length",
        "location":         "Location",
        "lineage":          "Pango lineage",
        "submission_date":  "Submission date",
        "variant":          "Variant",
        "clade":            "Clade",
    }

    idx = {}
    for key, col_name in REQUIRED.items():
        if col_name not in col_map:
            print(f"  ✗ REQUIRED column not found: '{col_name}'")
            print(f"    Available: {list(col_map.keys())}")
            sys.exit(1)
        idx[key] = col_map[col_name]
        print(f"  ✓ {key}: [{idx[key]}] '{col_name}'")

    for key, col_name in OPTIONAL.items():
        if col_name in col_map:
            idx[key] = col_map[col_name]
            print(f"  ✓ {key}: [{idx[key]}] '{col_name}'")
        else:
            print(f"  - {key}: not found (skipping this filter)")

    print()

    # ── Filter loop ──
    keep_ids = {}        # accession_id → {virus_name, date, location, lineage, ...}
    reasons = Counter()  # exclusion reason → count
    total = 0
    report_every = 2_000_000

    for line in fh:
        total += 1
        if total % report_every == 0:
            elapsed = time.time() - t0
            rate = total / elapsed
            print(f"    {total:>12,} rows | {len(keep_ids):>10,} kept | "
                  f"{elapsed:.0f}s | {rate:,.0f} rows/s")

        f = line.rstrip("\n\r").split("\t")
        # Pad short rows
        while len(f) <= max(idx.values()):
            f.append("")

        # --- Accession ID ---
        acc = f[idx["accession_id"]].strip()
        if not acc.startswith("EPI_ISL_"):
            reasons["bad_accession_id"] += 1
            continue

        # --- Host ---
        host = f[idx["host"]].strip()
        if host.lower() != "human":
            reasons["non_human"] += 1
            continue

        # --- Collection date ---
        date_str = f[idx["collection_date"]].strip()
        ok, reason = is_valid_complete_date(date_str, start_date, end_date)
        if not ok:
            reasons[reason] += 1
            continue

        # --- Is complete? ---
        if "is_complete" in idx:
            val = f[idx["is_complete"]].strip()
            if val and val != "True":
                reasons["not_complete"] += 1
                continue

        # --- Is high coverage? (REQUIRED) ---
        if "is_high_coverage" in idx:
            val = f[idx["is_high_coverage"]].strip()
            if val != "True":
                reasons["not_high_coverage"] += 1
                continue

        # --- Is low coverage? (additional safety) ---
        if "is_low_coverage" in idx:
            val = f[idx["is_low_coverage"]].strip()
            if val == "True":
                reasons["low_coverage"] += 1
                continue

        # --- N-Content ---
        if "n_content" in idx:
            val = f[idx["n_content"]].strip()
            if val:
                try:
                    if float(val) > args.max_n_content:
                        reasons["high_n_content"] += 1
                        continue
                except ValueError:
                    pass

        # --- Sequence length ---
        if "seq_length" in idx:
            val = f[idx["seq_length"]].strip()
            if val:
                try:
                    if int(val) < args.min_length:
                        reasons["short_sequence"] += 1
                        continue
                except ValueError:
                    pass

        # --- Nextstrain exclude list ---
        virus_name = f[idx["virus_name"]].strip()
        # Strip hCoV-19/ prefix for matching
        name_stripped = virus_name
        if name_stripped.startswith("hCoV-19/"):
            name_stripped = name_stripped[8:]
        if name_stripped in exclude_names or virus_name in exclude_names:
            reasons["nextstrain_excluded"] += 1
            continue

        # --- Passed all filters ---
        info = {"virus_name": virus_name, "date": date_str}
        if "location" in idx:
            info["location"] = f[idx["location"]].strip()
        if "lineage" in idx:
            info["lineage"] = f[idx["lineage"]].strip()
        if "submission_date" in idx:
            info["submission_date"] = f[idx["submission_date"]].strip()
        if "seq_length" in idx:
            info["length"] = f[idx["seq_length"]].strip()

        keep_ids[acc] = info

    fh.close()

    # Build virus_name → accession_id lookup for FASTA matching
    # GISAID bulk FASTA headers use virus names, NOT accession IDs
    keep_names = {}  # virus_name → accession_id
    for acc, info in keep_ids.items():
        vname = info["virus_name"]
        keep_names[vname] = acc
        # Also store without hCoV-19/ prefix
        if vname.startswith("hCoV-19/"):
            keep_names[vname[8:]] = acc
    elapsed = time.time() - t0
    excluded = total - len(keep_ids)

    print()
    print(f"  Metadata filtering done ({elapsed:.0f}s)")
    print(f"  ─────────────────────────────────")
    print(f"  Total rows:       {total:>14,}")
    print(f"  Passed filters:   {len(keep_ids):>14,}")
    print(f"  Excluded:         {excluded:>14,}")
    print(f"  ─────────────────────────────────")
    print(f"  Exclusion breakdown:")
    for reason, count in sorted(reasons.items(), key=lambda x: -x[1]):
        pct = 100.0 * count / total
        print(f"    {reason:25s} {count:>12,}  ({pct:.1f}%)")

    # =========================================================================
    # STEP 2: Filter FASTA
    # =========================================================================
    print()
    print("=" * 70)
    print("STEP 2: Filtering FASTA")
    print("=" * 70)
    print(f"  File: {args.fasta}")
    print(f"  Keeping: {len(keep_ids):,} sequences")
    print(f"  Virus name lookup entries: {len(keep_names):,}")
    print()

    t1 = time.time()
    fasta_fh = open_from_tarxz_or_plain(args.fasta)

    out_fasta_path = os.path.join(args.outdir, "filtered_sequences.fasta")
    out_fasta = open(out_fasta_path, "w")

    kept = 0
    skipped = 0
    total_seqs = 0
    writing = False
    current_acc = None

    for line in fasta_fh:
        if line.startswith(">"):
            total_seqs += 1
            if total_seqs % 2_000_000 == 0:
                elapsed = time.time() - t1
                print(f"    {total_seqs:>12,} scanned | {kept:>10,} kept | "
                      f"{elapsed:.0f}s")

            header = line.rstrip("\n\r")
            # GISAID bulk FASTA header format: >hCoV-19/Country/ID/Year
            # (no accession ID, no date — just virus name)
            virus_name = header[1:].strip()  # strip ">"

            # Also try without pipe-delimited extras (some formats have |EPI_ISL|date)
            virus_name_base = virus_name.split("|")[0].strip()

            # Try matching: full name, base name, without prefix
            matched_acc = None
            if virus_name_base in keep_names:
                matched_acc = keep_names[virus_name_base]
            elif virus_name in keep_names:
                matched_acc = keep_names[virus_name]
            else:
                # Try stripping hCoV-19/ prefix
                stripped = virus_name_base
                if stripped.startswith("hCoV-19/"):
                    stripped = stripped[8:]
                if stripped in keep_names:
                    matched_acc = keep_names[stripped]

            if matched_acc:
                writing = True
                current_acc = matched_acc
                kept += 1
                # Write header with accession ID appended for Nextclade
                out_fasta.write(f">{virus_name_base}|{matched_acc}\n")
            else:
                writing = False
                current_acc = None
                skipped += 1
        else:
            if writing:
                out_fasta.write(line)

    out_fasta.close()
    fasta_fh.close()

    elapsed = time.time() - t1
    print()
    print(f"  FASTA filtering done ({elapsed:.0f}s)")
    print(f"  Total in FASTA:   {total_seqs:>14,}")
    print(f"  Kept:             {kept:>14,}")
    print(f"  Skipped:          {skipped:>14,}")

    # =========================================================================
    # STEP 3: Write outputs
    # =========================================================================
    print()
    print("=" * 70)
    print("STEP 3: Writing outputs")
    print("=" * 70)

    # Filtered metadata
    meta_path = os.path.join(args.outdir, "filtered_metadata.tsv")
    with open(meta_path, "w") as f:
        f.write("accession_id\tvirus_name\tcollection_date\t"
                "location\tlineage\tsubmission_date\tlength\n")
        for acc, info in sorted(keep_ids.items()):
            f.write(f"{acc}\t{info.get('virus_name','')}\t"
                    f"{info.get('date','')}\t{info.get('location','')}\t"
                    f"{info.get('lineage','')}\t"
                    f"{info.get('submission_date','')}\t"
                    f"{info.get('length','')}\n")
    print(f"  ✓ {meta_path} ({len(keep_ids):,} rows)")

    # Report
    report_path = os.path.join(args.outdir, "prefilter_report.txt")
    with open(report_path, "w") as f:
        f.write("GISAID Pre-Filter Report\n")
        f.write(f"Generated: {datetime.now():%Y-%m-%d %H:%M:%S}\n")
        f.write(f"{'=' * 50}\n\n")
        f.write(f"Metadata:      {args.metadata}\n")
        f.write(f"FASTA:         {args.fasta}\n")
        f.write(f"Date range:    {args.start_date} to {args.end_date}\n")
        f.write(f"Min length:    {args.min_length} bp\n")
        f.write(f"Max N-content: {args.max_n_content}\n")
        f.write(f"Nextstrain exclude: {len(exclude_names)} entries\n\n")
        f.write(f"Metadata rows:        {total:>14,}\n")
        f.write(f"Passed filters:       {len(keep_ids):>14,}\n")
        f.write(f"FASTA total:          {total_seqs:>14,}\n")
        f.write(f"FASTA kept:           {kept:>14,}\n\n")
        f.write(f"Exclusion breakdown:\n")
        for reason, count in sorted(reasons.items(), key=lambda x: -x[1]):
            pct = 100.0 * count / total
            f.write(f"  {reason:25s} {count:>12,}  ({pct:.1f}%)\n")
    print(f"  ✓ {report_path}")

    # Config for downstream R scripts
    config_path = os.path.join(args.outdir, "pipeline_config.txt")
    with open(config_path, "w") as f:
        f.write(f"START_DATE={args.start_date}\n")
        f.write(f"END_DATE={args.end_date}\n")
        f.write(f"METADATA_ROWS={total}\n")
        f.write(f"PASSED_FILTERS={len(keep_ids)}\n")
        f.write(f"FASTA_KEPT={kept}\n")
    print(f"  ✓ {config_path}")

    print()
    print("=" * 70)
    print("PRE-FILTERING COMPLETE")
    print("=" * 70)
    print(f"  Output: {out_fasta_path}")
    print(f"  → Next: bash scripts/01_run_nextclade.sh --jobs 96")
    print()


if __name__ == "__main__":
    main()
