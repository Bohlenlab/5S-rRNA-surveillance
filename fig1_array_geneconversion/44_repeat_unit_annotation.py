#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# 44_repeat_unit_annotation.py — Builds a feature annotation of the 2168 bp 5S
# rDNA repeat-unit consensus (feature map, GC% and CpG tracks) and writes the
# feature table.
# Author: Jonathan Bohlen
# (c) 2026 Jonathan Bohlen. Code accompanying Sengl et al. (2026),
# "Surveillance and selection of 5S ribosomal RNA genes in the human genome."
# Released under the MIT License; see LICENSE at the repository root.
# -----------------------------------------------------------------------------
"""
44_repeat_unit_annotation.py

Feature annotation of the 2168 bp human 5S rDNA repeat unit (population consensus),
replacing the coarse NTS-pre / gene / NTS-post split with biologically defined
features:

  5S rRNA gene (Pol III type-1), 630-749 (+1 = 630, 100% id to canonical 5S rRNA)
    internal control region (ICR): Box A, Intermediate Element (IE), Box C
  Pol III terminator: oligo-dT immediately 3' of the gene
  Alu SINE (AluY-related, antisense, ~93% id): 787-1066
  5' / 3' external spacers, with GC-rich low-complexity / CpG microsatellite

Produces a figure (feature map + GC% + CpG tracks) and a feature TSV.

Outputs:
  figures/44_repeat_unit_annotation.pdf
  repeat_unit_features.tsv
"""

import json
import os
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrow, Rectangle, FancyBboxPatch
import numpy as np

HPRC = Path(os.environ.get("FIVES_OUT", "output"))
CONS = Path(os.environ.get("FIVES_DATA", "data")) / "consensus_reference.json"
FIG_DIR = HPRC / "figures" / "07_sequence_annotation"

seq = json.load(open(CONS))["consensus"].upper()
L = len(seq)
GENE_S, GENE_E = 630, 749     # +1 = 630

# ── feature table (1-based, inclusive) ────────────────────────────────────────
# ICR offsets relative to +1 (canonical 5S internal control region)
def g(off):  # gene-relative + -> consensus position
    return GENE_S - 1 + off

FEATURES = [
    ("5'_external_spacer", 1, 629, "spacer", "5' non-transcribed spacer (G-rich)"),
    ("CpG_microsatellite", 9, 21, "lowcomplex", "(CG)n CpG-rich microsatellite"),
    ("GT_rich_lowcomplex", 34, 93, "lowcomplex", "(GT)n / G-rich low-complexity (CpG-island-like 5' junction)"),
    ("5S_rRNA_gene", GENE_S, GENE_E, "gene", "5S rRNA gene, Pol III type-1 (100% id)"),
    ("ICR_BoxA", g(50), g(64), "promoter", "Internal control region: Box A"),
    ("ICR_IE", g(67), g(72), "promoter", "Internal control region: Intermediate Element"),
    ("ICR_BoxC", g(80), g(97), "promoter", "Internal control region: Box C"),
    ("PolIII_terminator", 750, 765, "terminator", "Pol III terminator (oligo-dT)"),
    ("Alu_SINE", 787, 1066, "alu", "Alu SINE, AluY-related, antisense, ~93% id"),
    ("polymorphic_indel", 1181, 1190, "indel", "GC-rich polymorphic indel (84% del)"),
    ("CT_microsatellite", 1916, 1968, "lowcomplex", "CT-rich pyrimidine microsatellite"),
    ("CTGT_microsatellite", 2009, 2155, "lowcomplex", "(CTGT/CTCT)n pyrimidine microsatellite (TRF period-8, 8 copies)"),
    ("3'_external_spacer", 766, 2168, "spacer", "3' non-transcribed spacer (C-rich)"),
]

# write TSV
with open(HPRC / "repeat_unit_features.tsv", "w") as fh:
    fh.write("feature\tstart\tend\tlength\tclass\tdescription\n")
    for name, s, e, cls, desc in FEATURES:
        fh.write(f"{name}\t{s}\t{e}\t{e-s+1}\t{cls}\t{desc}\n")
print("Feature table -> repeat_unit_features.tsv")
for name, s, e, cls, desc in FEATURES:
    print(f"  {name:20s} {s:5d}-{e:<5d} ({e-s+1:4d} bp)  {desc}")

# ── GC% and CpG tracks ────────────────────────────────────────────────────────
WIN = 40
xs = np.arange(WIN//2, L-WIN//2)
gc = np.array([ (seq[i-WIN//2:i+WIN//2].count("G")+seq[i-WIN//2:i+WIN//2].count("C"))/WIN*100 for i in xs])
cpg = np.array([ sum(1 for j in range(i-WIN//2, i+WIN//2-1) if seq[j:j+2]=="CG") for i in xs])

# ── figure ────────────────────────────────────────────────────────────────────
fig, (axm, axg, axc) = plt.subplots(3, 1, figsize=(15, 8),
        gridspec_kw=dict(height_ratios=[2.4, 1, 1], hspace=0.32))

COL = {"spacer":"#e8e8e8","gene":"#f4b942","promoter":"#b8860b","terminator":"#d1495b",
       "alu":"#7b68ee","lowcomplex":"#7fb069","indel":"#ff9505"}

# baseline
axm.plot([0, L], [0, 0], color="black", lw=1.2, zorder=1)
def box(ax, s, e, y0, h, color, ec="black", lw=0.8, z=3, alpha=1):
    ax.add_patch(Rectangle((s, y0), e-s, h, facecolor=color, edgecolor=ec,
                           linewidth=lw, zorder=z, alpha=alpha))

# spacers (thin)
box(axm, 1, 629, -0.18, 0.36, COL["spacer"], z=2)
box(axm, 766, 2168, -0.18, 0.36, COL["spacer"], z=2)
# gene (tall)
box(axm, GENE_S, GENE_E, -0.45, 0.9, COL["gene"], z=3)
# ICR sub-boxes inside gene
for nm in ["ICR_BoxA","ICR_IE","ICR_BoxC"]:
    f=[x for x in FEATURES if x[0]==nm][0]
    box(axm, f[1], f[2], -0.30, 0.6, COL["promoter"], z=4)
# terminator
box(axm, 750, 765, -0.45, 0.9, COL["terminator"], z=4)
# Alu (tall) with antisense arrow
box(axm, 787, 1066, -0.40, 0.8, COL["alu"], z=3)
axm.annotate("", xy=(800, 0), xytext=(1055, 0),
             arrowprops=dict(arrowstyle="-|>", color="white", lw=2), zorder=5)
# low-complexity / microsat ticks (CpG, GT-rich, CT/CTGT microsatellites)
for s_, e_ in [(9, 21), (34, 93), (1916, 1968), (2009, 2155)]:
    box(axm, s_, e_, -0.12, 0.24, COL["lowcomplex"], z=5)
# polymorphic indel marker
axm.plot(1185, 0.55, marker="*", ms=16, color=COL["indel"], markeredgecolor="black",
         markeredgewidth=0.5, zorder=6)

# transcription start arrow (+1, rightward)
axm.annotate("", xy=(GENE_S+70, 0.75), xytext=(GENE_S, 0.75),
             arrowprops=dict(arrowstyle="-|>", color="black", lw=1.6))
axm.text(GENE_S, 0.92, "+1", fontsize=9, ha="left", fontweight="bold")

# labels
labels = [
    (315, 0.78, "5′ external spacer", "center"),
    (689, -0.72, "5S rRNA gene", "center"),
    (689, -0.92, "(120 bp, Pol III type-1)", "center"),
    (757, 1.08, "terminator", "center"),
    (926, -0.70, "Alu SINE", "center"),
    (926, -0.90, "(antisense, ~93% AluY)", "center"),
    (1185, 1.02, "polymorphic\nindel", "center"),
    (1700, 0.78, "3′ external spacer", "center"),
    (60, -0.46, "(CG)n/(GT)n", "center"),
    (2080, -0.52, "(CT/CTGT)n\nmicrosat", "center"),
]
for x,y,t,ha in labels:
    axm.text(x, y, t, fontsize=8.5, ha=ha, va="center")
# ICR label (point to gene, text well below)
axm.annotate("internal control region (ICR):\nBox A · IE · Box C",
             xy=(g(73), -0.30), xytext=(689, -1.30),
             fontsize=7.5, ha="center", color="#7a5c00",
             arrowprops=dict(arrowstyle="-", color="#b8860b", lw=0.8))

axm.set_xlim(-20, L+20); axm.set_ylim(-1.5, 1.3)
axm.axis("off")
axm.set_title("Human 5S rDNA repeat unit (2168 bp consensus) — feature annotation",
              fontsize=13, fontweight="bold", loc="left")

# GC track
axg.fill_between(xs, gc, 50, where=gc>=50, color="#4575b4", alpha=0.5, lw=0)
axg.fill_between(xs, gc, 50, where=gc<50, color="#d73027", alpha=0.5, lw=0)
axg.plot(xs, gc, color="black", lw=0.6)
axg.axhline(50, color="grey", lw=0.6, ls="--")
axg.axvspan(GENE_S, GENE_E, color="#f4b942", alpha=0.25)
axg.axvspan(787, 1066, color="#7b68ee", alpha=0.15)
axg.set_xlim(-20, L+20); axg.set_ylabel("GC %", fontsize=9)
axg.set_ylim(20, 100); axg.set_xticklabels([])
axg.text(0.01, 0.85, f"mean GC {100*(seq.count('G')+seq.count('C'))/L:.0f}%",
         transform=axg.transAxes, fontsize=8)

# CpG track
axc.fill_between(xs, cpg, color="#2c7fb8", alpha=0.7, lw=0)
axc.axvspan(GENE_S, GENE_E, color="#f4b942", alpha=0.25)
axc.axvspan(787, 1066, color="#7b68ee", alpha=0.15)
axc.set_xlim(-20, L+20); axc.set_ylabel(f"CpG / {WIN}bp", fontsize=9)
axc.set_xlabel("Position in repeat unit (bp)", fontsize=10)

plt.savefig(FIG_DIR / "44_repeat_unit_annotation.pdf", bbox_inches="tight", dpi=200)
plt.close()
print(f"\nFigure -> figures/44_repeat_unit_annotation.pdf")
PY_DONE = True
