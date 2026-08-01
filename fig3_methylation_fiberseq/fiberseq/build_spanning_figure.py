#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# build_spanning_figure.py — Spanning-molecule CpG methylation and m6A accessibility across the 5S array ends (Fiber-seq).
# Author: Jonathan Bohlen
# (c) 2026 Jonathan Bohlen. Code accompanying Sengl et al. (2026),
# "Surveillance and selection of 5S ribosomal RNA genes in the human genome."
# Released under the MIT License; see LICENSE at the repository root.
# -----------------------------------------------------------------------------
"""
Spanning-molecule CpG methylation and m6A accessibility across the 5S array ends,
both flanks resolved separately (telomere-proximal high-coord edge and
centromere-proximal low-coord edge), per Fiber-seq sample. m6A accessibility is
taken from the same single molecules.

Distance-into-array dkb = sign*(refpos-edge)/1000 ; <0 = unique flank, >0 = into array.
Spanning molecules = reads that touch the flank and reach into the array (span>=MINSPAN).
Inputs: data/<sample>/<sample>.mods.tsv.gz (modkit extract: read_id,refpos,code,qual), refs/region_wide.fa
Output: figures/spanning_5S_fiberseq.pdf  +  data/spanning_edge_summary.tsv

Paths are read from environment variables (see repository README):
    FIVES_DATA  input derived-data directory
    FIVES_REFS  reference fasta directory
    FIVES_OUT   output directory
"""
import os, gzip, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt

DATA, REFS, FIG = os.environ.get("FIVES_DATA", "data"), os.environ.get("FIVES_REFS", "refs"), os.environ.get("FIVES_OUT", "output")
SAMPLES = ["HG002", "CHM1", "CHM13", "GM12878", "K562"]
ALO, AHI = 227745360, 228024984          # CHM13 5S array edges (consensus asm20 -N500)
BASE0 = 227695000 - 1                     # 0-based genomic of region_wide.fa index 0
P, GENE_S, GENE_E, ULEN = 2.241, 630, 748, 2168
GW = (GENE_E - GENE_S) / 1000.0
MINSPAN = 8.0                             # Fiber-seq reads ~16 kb -> require >=8 kb span
BINKB, MINCALLS = 0.25, 8                 # fixed dkb bins
FLANK, ARRAY = 12.0, 20.0                 # plot window
ENDS = {"telomere": (-1, AHI, "#b2182b"), "centromere": (+1, ALO, "#2166AC")}

# --- reference CpG cytosine positions (0-based genomic) ---
seq = []
for ln in open(f"{REFS}/region_wide.fa"):
    if not ln.startswith(">"): seq.append(ln.strip())
seq = "".join(seq).upper()
cpg = set()
for i in range(len(seq) - 1):
    if seq[i] == "C" and seq[i + 1] == "G":
        cpg.add(BASE0 + i); cpg.add(BASE0 + i + 1)
print(f"{len(cpg)//2} CpG sites in region")

def load(s):
    p = f"{DATA}/{s}/{s}.mods.tsv.gz"
    df = pd.read_csv(p, sep="\t")
    df["meth"] = (df["qual"] >= 0.8).astype(int)   # confident modified
    return df

def binned(df, sign, edge, code, cpg_only):
    d = df[df["code"] == code].copy()
    if cpg_only: d = d[d["refpos"].isin(cpg)]
    if d.empty: return None, 0
    d["dkb"] = sign * (d["refpos"] - edge) / 1000.0
    # spanning-read selection: read must touch flank (<0) and reach array (>0), span>=MINSPAN
    g = d.groupby("read_id")["dkb"]
    keep = g.agg(lambda x: (x.min() < -0.5) and (x.max() > 2.0) and (x.max() - x.min() >= MINSPAN))
    d = d[d["read_id"].isin(keep[keep].index)]
    nmol = d["read_id"].nunique()
    if d.empty: return None, 0
    d["bin"] = (np.floor(d["dkb"] / BINKB) * BINKB + BINKB / 2)
    gb = d.groupby("bin")["meth"].agg(["size", "sum"]).reset_index()
    gb = gb[gb["size"] >= MINCALLS]
    gb["pct"] = gb["sum"] / gb["size"] * 100
    return gb, nmol

fig, axes = plt.subplots(len(SAMPLES), 2, figsize=(15, 16), sharex="col")
summary = []
for row, s in enumerate(SAMPLES):
    df = load(s)
    for col, (end, (sign, edge, ccol)) in enumerate(ENDS.items()):
        ax = axes[row, col]; ax2 = ax.twinx()
        cpgb, nmol = binned(df, sign, edge, "m", True)
        m6b, _ = binned(df, sign, edge, "a", False)
        ax.axvspan(0, ARRAY, color="gold", alpha=.10, lw=0)
        ax.axvspan(-FLANK, 0, color="#dddddd", alpha=.30, lw=0)
        ax.axvline(0, color="k", lw=.7, ls="--", alpha=.6)
        # 5S gene period marks into the array
        k = 0
        while k * P < ARRAY:
            ax.axvspan(k * P, k * P + GW, color="#1a9850", alpha=.18, lw=0); k += 1
        if cpgb is not None:
            sm = cpgb.set_index("bin")["pct"].rolling(3, center=True, min_periods=1).mean()
            ax.plot(sm.index, sm.values, color=ccol, lw=2.0, zorder=4)
        if m6b is not None:
            sm6 = m6b.set_index("bin")["pct"].rolling(3, center=True, min_periods=1).mean()
            ax2.plot(sm6.index, sm6.values, color="#e08214", lw=1.3, alpha=.8, zorder=3)
        ax.set_xlim(-FLANK, ARRAY); ax.set_ylim(0, 100); ax2.set_ylim(0, 25)
        ax.set_ylabel(f"{s}\n5mCpG %", color=ccol, fontsize=8)
        ax2.set_ylabel("m6A %", color="#b5650a", fontsize=7)
        ax.axhline(65, color="grey", ls=":", lw=.7)
        if row == 0:
            ax.set_title(f"{end.capitalize()}-proximal end  (<0 = flank, gold = array, green = 5S gene)", fontsize=10)
        if row == len(SAMPLES) - 1:
            ax.set_xlabel(f"Distance into array from {end} edge (kb)")
        ax.text(.02, .04, f"n_mol={nmol}", transform=ax.transAxes, fontsize=7, color="#555")
        # edge (0-4kb) vs interior (8-20kb) summary
        if cpgb is not None:
            e = cpgb.loc[(cpgb.bin >= 0) & (cpgb.bin < 4), "pct"].mean()
            it = cpgb.loc[(cpgb.bin >= 8) & (cpgb.bin <= 20), "pct"].mean()
            fl = cpgb.loc[cpgb.bin < 0, "pct"].mean()
            summary.append((s, end, nmol, round(fl, 1), round(e, 1), round(it, 1), round(it - e, 1)))
fig.suptitle("Spanning-molecule 5mCpG (left) + m6A accessibility (right) across the 5S array ends — Fiber-seq, both flanks resolved", y=.997, fontsize=12)
fig.tight_layout()
out = f"{FIG}/spanning_5S_fiberseq.pdf"
fig.savefig(out, dpi=150); print("wrote", out)
sm = pd.DataFrame(summary, columns=["sample", "end", "n_mol", "flank%", "edge0_4kb%", "interior8_20kb%", "int_minus_edge"])
sm.to_csv(f"{DATA}/spanning_edge_summary.tsv", sep="\t", index=False)
print(sm.to_string(index=False))
