#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# 03_aggregate_45S.py — Aggregates per-call ONT 45S methylation into per-read,
# per-bin signed distance from array edges.
# Author: Jonathan Bohlen
# (c) 2026 Jonathan Bohlen. Code accompanying Sengl et al. (2026),
# "Surveillance and selection of 5S ribosomal RNA genes in the human genome."
# Released under the MIT License; see LICENSE at the repository root.
# -----------------------------------------------------------------------------
"""
03_aggregate_45S.py

Aggregate per-call ONT methylation (methylation/calls_*.tsv) into per-read, per-bin
methylation expressed as signed distance from the array edge, per NOR and per edge side.

  dkb < 0  = in the unique flank (outside the array)
  dkb > 0  = into the rDNA array
Each call is assigned to the NEARER array edge of its NOR:
  DJ (distal, low-coord edge = array_lo):  dkb = (pos - array_lo)/1000
  PJ (proximal, high-coord edge = array_hi): dkb = (array_hi - pos)/1000

Output: <FIVES_DATA>/methylation/molecule_bin_meth_45S_ont.tsv  (nor, edge, sample, read, dkb, n, meth)
Reads with < MINCALL confident calls on a (NOR,edge) are dropped.

Paths are read from environment variables (FIVES_DATA).
"""
import glob, os
from pathlib import Path
import numpy as np, pandas as pd

H45 = Path(os.environ.get("FIVES_DATA", "data"))
MD  = H45 / "methylation"
NOR = {"chr13": (5_770_548, 9_348_041), "chr14": (2_099_537, 2_817_811),
       "chr15": (2_506_442, 4_707_485), "chr21": (3_108_298, 5_612_715),
       "chr22": (4_793_794, 5_720_650)}
BW, MINCALL = 0.5, 3            # 500 bp bins; >=3 confident calls per read-edge

def main():
    files = sorted(glob.glob(str(MD / "calls_*.tsv")))
    if not files:
        print("no calls_*.tsv yet"); return
    rows = []
    summ = []
    for f in files:
        s = Path(f).name[len("calls_"):-len(".tsv")]
        df = pd.read_csv(f, sep="\t")
        if df.empty:
            summ.append((s, 0, 0)); continue
        df = df[df.chrom.isin(NOR)].copy()
        lo = df.chrom.map(lambda c: NOR[c][0]); hi = df.chrom.map(lambda c: NOR[c][1])
        d_dj = df.ref_position - lo            # signed dist from distal edge (into array = +)
        d_pj = hi - df.ref_position            # signed dist from proximal edge (into array = +)
        near_dj = (df.ref_position - lo).abs() <= (hi - df.ref_position).abs()
        df["edge"] = np.where(near_dj, "DJ", "PJ")
        df["dkb"]  = np.where(near_dj, d_dj, d_pj) / 1000.0
        summ.append((s, df.read_id.nunique(), len(df)))
        for (nor, edge, rid), g in df.groupby(["chrom", "edge", "read_id"]):
            if len(g) < MINCALL: continue
            gb = g.assign(b=(np.floor(g.dkb / BW) * BW)).groupby("b")["meth"].agg(["size", "sum"])
            for b, r in gb.iterrows():
                rows.append((nor, edge, s, str(rid)[:16], round(b + BW/2, 3),
                             int(r["size"]), int(r["sum"])))
    out = MD / "molecule_bin_meth_45S_ont.tsv"
    pd.DataFrame(rows, columns=["nor", "edge", "sample", "read", "dkb", "n", "meth"]
                 ).to_csv(out, sep="\t", index=False)
    print(f"wrote {out}  ({len(rows)} read-bins from {len(files)} sample(s))")
    print(f"{'sample':10s}{'flank_reads':>12}{'calls':>10}")
    for s, nr, nc in summ:
        print(f"{s:10s}{nr:>12}{nc:>10}")

if __name__ == "__main__":
    main()
