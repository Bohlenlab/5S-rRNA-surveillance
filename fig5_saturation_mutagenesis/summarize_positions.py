#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# summarize_positions.py — Sum per-base A/C/G/T coverage across per-position
# mismatch CSVs into per-position and overall totals.
# Author: Jonathan Bohlen
# (c) 2026 Jonathan Bohlen. Code accompanying Sengl et al. (2026),
# "Surveillance and selection of 5S ribosomal RNA genes in the human genome."
# Released under the MIT License; see LICENSE at the repository root.
# -----------------------------------------------------------------------------
import os, csv, glob, sys
from collections import defaultdict
import argparse

def parse_args():
    p = argparse.ArgumentParser(description="Sum A,C,G,T,total_coverage across *_mismatches.csv files.")
    p.add_argument("--pattern", default=os.path.join(os.environ.get("FIVES_DATA", "data"), "*_mismatches.csv"),
                   help="Glob for input files (default: *_mismatches.csv)")
    p.add_argument("--out-by-pos", default=os.path.join(os.environ.get("FIVES_OUT", "output"), "summary_by_position.csv"),
                   help="Output CSV with per-position sums")
    p.add_argument("--out-total", default=os.path.join(os.environ.get("FIVES_OUT", "output"), "summary_totals.csv"),
                   help="Output CSV with one overall totals row")
    return p.parse_args()

def main():
    a = parse_args()
    files = sorted(glob.glob(a.pattern))
    if not files:
        print(f"No files matched: {a.pattern}", file=sys.stderr)
        return 1

    # Aggregate per position
    per_pos = defaultdict(lambda: {"ref_base": None, "A":0, "C":0, "G":0, "T":0, "total_coverage":0})
    totals  = {"A":0, "C":0, "G":0, "T":0, "total_coverage":0}

    required = {"position","ref_base","A","C","G","T","total_coverage"}

    for path in files:
        with open(path, newline="") as fh:
            r = csv.DictReader(fh)
            if not required.issubset(r.fieldnames or []):
                missing = required - set(r.fieldnames or [])
                raise RuntimeError(f"{path} missing columns: {missing}")
            for row in r:
                try:
                    pos = int(row["position"])
                except (ValueError, TypeError):
                    continue  # skip malformed lines
                slot = per_pos[pos]
                # keep the first ref_base we see; warn on inconsistency
                rb = (row.get("ref_base") or "").upper()
                if slot["ref_base"] is None:
                    slot["ref_base"] = rb
                elif rb and rb != slot["ref_base"]:
                    # Positions should agree; keep the first
                    pass
                # add counts
                for key in ("A","C","G","T","total_coverage"):
                    val = int(row.get(key) or 0)
                    slot[key] += val
                    totals[key] += val

    # Write per-position sums
    with open(a.out_by_pos, "w", newline="") as out:
        w = csv.writer(out)
        w.writerow(["position","ref_base","total_coverage","A","C","G","T"])
        for pos in sorted(per_pos):
            d = per_pos[pos]
            w.writerow([pos, d["ref_base"] or "", d["total_coverage"], d["A"], d["C"], d["G"], d["T"]])

    # Write one-row overall totals
    with open(a.out_total, "w", newline="") as out2:
        w = csv.writer(out2)
        w.writerow(["A","C","G","T","total_coverage"])
        w.writerow([totals["A"], totals["C"], totals["G"], totals["T"], totals["total_coverage"]])

    print(f"Processed {len(files)} files → {a.out_by_pos} (positions={len(per_pos)}) and {a.out_total}")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
