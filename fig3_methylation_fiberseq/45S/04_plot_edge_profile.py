#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# 04_plot_edge_profile.py — Per-NOR, per-edge 45S array-edge methylation profiles.
# Author: Jonathan Bohlen
# (c) 2026 Jonathan Bohlen. Code accompanying Sengl et al. (2026),
# "Surveillance and selection of 5S ribosomal RNA genes in the human genome."
# Released under the MIT License; see LICENSE at the repository root.
# -----------------------------------------------------------------------------
"""
04_plot_edge_profile.py

Per-NOR, per-EDGE 45S array-edge methylation profile from molecule_bin_meth_45S_ont.tsv.
7 rows (ALL, cross-NOR, 5 NORs) x 2 columns (DJ | PJ).
  DJ = distal junction (telomere-side) — plotted normally (flank left -> into array right).
  PJ = proximal junction (centromere-side) — x-axis INVERTED (array left -> flank right) so the two
       array interiors meet at the centre and the figure reads as the contiguous chr q-arm locus.
Into-array distance is PERIOD-NORMALISED per NOR to a common 44 kb so gene bodies stay phase-aligned
at all depths (flank <0 left unscaled). Style: 8 cm x 4 cm panels, 8-pt font, thin lines; DJ = blue,
PJ = red. Two band options are produced:
  <label>.pdf         band = mean +/- 95% CI (call-weighted)
  <label>_spread.pdf  band = median + IQR (25-75th) + decile (10-90th) ACROSS DONORS

Input is read from <FIVES_DATA>/methylation; figures are written to <FIVES_OUT>.
Paths are read from environment variables (FIVES_DATA, FIVES_OUT).
"""
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"pdf.fonttype": 42, "ps.fonttype": 42, "font.size": 8,
    "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "xtick.major.size": 2, "ytick.major.size": 2})
import glob, os
CM = 1 / 2.54
BASE = os.environ.get("FIVES_DATA", "data")
NORS = ["ALL", "xNOR", "chr13", "chr14", "chr15", "chr21", "chr22"]
# column order = genomic order read left->right: telomere (far left) -> ... -> centromere (far right).
# LEFT col  = DJ (distal/telomere edge), NOT inverted: telomere flank far-left, array to the right.
# RIGHT col = PJ (proximal/centromere edge), INVERTED: array on the left (meets DJ array), centromere flank far-right.
EDGES = [("DJ", "Distal junction (DJ) · telomere side"),
         ("PJ", "Proximal junction (PJ) · centromere side")]
COL = {"chr13": "#1b9e77", "chr14": "#d95f02", "chr15": "#7570b3", "chr21": "#e7298a",
       "chr22": "#66a61e", "xNOR": "#555555", "ALL": "#000000"}   # y-label (row) colour = NOR identity
EDGECOL = {"DJ": "#2166AC", "PJ": "#D6604D"}                       # trace/band colour
NAME = {**{c: c for c in NORS}, "xNOR": "cross-NOR\n(ambiguous)", "ALL": "ALL\n(pooled)"}
# ---- DISPLAY options ----
XLO, BW = -50, 3.0
MINCALLS = 30
XHI_CAP = 500.0
PERIOD = 44.0
GWKB = 13.35
PW_CM, PH_CM = 8.0, 4.0      # panel width / height (cm)
MINCALLS_DONOR = 10          # min calls for a single donor to contribute to a bin's spread
MINDONORS_SPREAD = 8         # min donors with data in a bin to draw a spread band

# pool per-donor edge caches (named NOR = uniquely-mapped MAPQ>=20; xNOR = ambiguous MAPQ<20).
ef = sorted(glob.glob(f"{BASE}/methylation/molbinedge_*.tsv"))
if ef:
    d = pd.concat([pd.read_csv(f, sep="\t") for f in ef], ignore_index=True)
else:
    d = pd.read_csv(f"{BASE}/methylation/molecule_bin_meth_45S_ont.tsv", sep="\t")
rmin = d.groupby(["nor", "edge", "sample", "read"]).dkb.transform("min")
d = d[rmin < 0].copy()

ndonor = d["sample"].nunique()
label = d["sample"].iloc[0] if ndonor == 1 else f"cohort_{ndonor}donors"

NORC = {"chr13": (5_770_548, 9_348_041), "chr14": (2_099_537, 2_817_811), "chr15": (2_506_442, 4_707_485),
        "chr21": (3_108_298, 5_612_715), "chr22": (4_793_794, 5_720_650)}
GENE_IV = {}
try:
    gi = pd.read_csv(f"{BASE}/methylation/gene_intervals.tsv", sep="\t")
    GENE_IV = {c: sorted(zip(g.start, g.end)) for c, g in gi.groupby("chrom")}
except Exception as e:
    print("gene shading skipped:", e)

def nor_period(sn):
    lo, hi = NORC[sn]
    ds = sorted((a - lo) / 1000. for a, b in GENE_IV.get(sn, []))
    return float(np.median(np.diff(ds))) if len(ds) >= 3 else PERIOD
PERIODN = {sn: nor_period(sn) for sn in NORC}
MEANPER = float(np.mean(list(PERIODN.values()))) if PERIODN else PERIOD
def scale_for(nor):
    return PERIOD / PERIODN.get(nor, MEANPER) if nor in PERIODN else PERIOD / MEANPER
d["dkb_n"] = np.where(d.dkb > 0, d.dkb * d["nor"].map(scale_for), d.dkb)

def phi0(edge):
    firsts = []
    for sn, (lo, hi) in NORC.items():
        ds = sorted((a - lo) / 1000. if edge == "DJ" else (hi - b) / 1000. for a, b in GENE_IV.get(sn, []))
        ds = [x for x in ds if x > 0]
        if ds: firsts.append(ds[0] * scale_for(sn))
    return float(np.median(firsts)) if firsts else (4.0 if edge == "DJ" else 38.0)
PHI0 = {e: phi0(e) for e, _ in EDGES}

_named = d[(d.nor != "xNOR") & (d.dkb_n > 0)]
_cov = _named.assign(_b=(np.floor(_named.dkb_n / BW) * BW)).groupby("_b").n.sum()
_cov = _cov[_cov >= MINCALLS]
XHI = float(min(_cov.index.max() + 2 * BW, XHI_CAP)) if len(_cov) else 120.0
key = ["nor", "edge", "sample", "read"]
span = d.groupby(key).agg(mn=("dkb", "min"), mx=("dkb", "max"))
junc = set(span[(span.mn < 0) & (span.mx > 0)].index)

def binprofile(sub):         # call-weighted mean + 95% CI
    sub = sub[(sub.dkb_n >= XLO) & (sub.dkb_n < XHI)].copy()
    if sub.empty: return sub
    sub["b"] = (np.floor(sub.dkb_n / BW) * BW + BW / 2)
    g = sub.groupby("b").agg(m=("meth", "sum"), n=("n", "sum"))
    g = g[g.n >= MINCALLS]
    g["pct"] = 100 * g.m / g.n
    g["ci"] = 100 * 1.96 * np.sqrt((g.pct / 100) * (1 - g.pct / 100) / g.n)
    return g

def binprofile_spread(sub):  # per-donor distribution per bin -> median + IQR (25-75) + decile (10-90)
    sub = sub[(sub.dkb_n >= XLO) & (sub.dkb_n < XHI)].copy()
    if sub.empty: return sub
    sub["b"] = (np.floor(sub.dkb_n / BW) * BW + BW / 2)
    pdn = sub.groupby(["b", "sample"]).agg(m=("meth", "sum"), n=("n", "sum")).reset_index()
    pdn = pdn[pdn.n >= MINCALLS_DONOR]
    pdn["pct"] = 100 * pdn.m / pdn.n
    rows = []
    for b, gg in pdn.groupby("b"):
        v = gg.pct.values
        if len(v) < MINDONORS_SPREAD: continue
        q10, q25, q50, q75, q90 = np.percentile(v, [10, 25, 50, 75, 90])
        rows.append((b, q10, q25, q50, q75, q90))
    return pd.DataFrame(rows, columns=["b", "q10", "q25", "q50", "q75", "q90"]).set_index("b")

def gene_bands(ax, edge):
    k = 0
    while PHI0[edge] + k * PERIOD < XHI:
        x0 = PHI0[edge] + k * PERIOD
        if x0 + GWKB > 0:
            ax.axvspan(max(x0, 0), min(x0 + GWKB, XHI), color="#fdbe85", alpha=0.5, lw=0, zorder=0)
        k += 1

def build(mode):
    figw = 2 * PW_CM * CM / 0.80 + 0.4
    figh = 7 * PH_CM * CM / 0.86 + 0.5
    fig, axes = plt.subplots(7, 2, figsize=(figw, figh), sharex="col", sharey=True,
                             constrained_layout=True)   # sharex per-column so only PJ inverts
    for r, nor in enumerate(NORS):
        for c, (edge, etitle) in enumerate(EDGES):
            ax = axes[r, c]
            sub = d[d.edge == edge] if nor == "ALL" else d[(d.nor == nor) & (d.edge == edge)]
            nrd = sub.groupby(["sample", "read"]).ngroups
            col = EDGECOL[edge]
            ax.axvspan(XLO, 0, color="0.93", zorder=0)               # flanking satellite
            gene_bands(ax, edge)                                     # consensus phase-aligned gene bands
            ax.axvline(0, color="k", lw=0.6, ls="--", zorder=2)
            if mode == "spread":
                g = binprofile_spread(sub)
                if len(g):
                    ax.fill_between(g.index, g.q10, g.q90, color=col, alpha=0.15, lw=0, zorder=3)
                    ax.fill_between(g.index, g.q25, g.q75, color=col, alpha=0.33, lw=0, zorder=3)
                    ax.plot(g.index, g.q50, "-", color=col, lw=1.0, zorder=4)
            else:
                g = binprofile(sub)
                if len(g):
                    ax.fill_between(g.index, g.pct - g.ci, g.pct + g.ci, color=col, alpha=0.20, lw=0)
                    ax.plot(g.index, g.pct, "-", color=col, lw=1.0, zorder=3)
            ax.set_ylim(0, 100)
            ax.set_xlim((XHI, XLO) if c == 1 else (XLO, XHI))        # RIGHT col (DJ) inverted -> arrays meet centre
            ax.tick_params(labelsize=6.5, length=2)
            ax.spines[["top", "right"]].set_visible(False)
            satx = 0.985 if c == 1 else 0.015                        # flank side (satellite label)
            ax.text(satx, 0.06, "satellite", transform=ax.transAxes, fontsize=5, color="0.5",
                    va="bottom", ha=("right" if c == 1 else "left"))
            ax.text(0.5, 0.06, f"{nrd} reads", transform=ax.transAxes, fontsize=5, color="0.4",
                    va="bottom", ha="center")
            if c == 0:
                ax.set_ylabel(f"{NAME[nor]}\n% mCpG", fontsize=7, color=COL[nor], fontweight="bold")
            if r == 0:
                ax.set_title(etitle, fontsize=8, color=col)
    axes[-1, 0].set_xlabel("telomere flank ·  into array from distal (telomere) edge (kb)", fontsize=6.5)
    axes[-1, 1].set_xlabel("into array from proximal (centromere) edge (kb)  · centromere flank  [inverted]", fontsize=6.5)
    bandtxt = ("median + IQR (25-75th) + decile (10-90th) across donors" if mode == "spread"
               else "mean +/- 95% CI (call-weighted)")
    fig.suptitle(f"45S rDNA array-edge methylation, per NOR and per edge  ·  {label} (HPRC ONT)  ·  "
                 f"period-normalised; PJ mirrored; band = {bandtxt}", fontsize=7.5)
    suf = "_spread" if mode == "spread" else ""
    fig.savefig(f'{os.environ.get("FIVES_OUT","output")}/04_edge_methylation_perNOR_perEdge_{label}{suf}.pdf', dpi=300, bbox_inches="tight")
    fig.savefig(f'{os.environ.get("FIVES_OUT","output")}/04_edge_methylation_perNOR_perEdge_{label}{suf}.png', dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote 04_edge_methylation_perNOR_perEdge_{label}{suf}.pdf  (band={mode})")

build("ci")
build("spread")
print("per-NOR periods (kb):", {k: round(v, 1) for k, v in PERIODN.items()}, " mean", round(MEANPER, 1))
print("consensus first-gene phase (kb):", {k: round(v, 1) for k, v in PHI0.items()})
