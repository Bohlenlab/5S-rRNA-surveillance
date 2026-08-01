#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# compute_rnu2_ownasm.py — Per-unit promoter CpG methylation of the RNU2 (U2 snRNA) array on HG002's own assembly.
# Author: Jonathan Bohlen
# (c) 2026 Jonathan Bohlen. Code accompanying Sengl et al. (2026),
# "Surveillance and selection of 5S ribosomal RNA genes in the human genome."
# Released under the MIT License; see LICENSE at the repository root.
# -----------------------------------------------------------------------------
"""RNU2 methylation on HG002's own assembly (maternal 62-unit and paternal
16-unit arrays). Per-CpG 5mC from pysam MM/ML tags. Bulk (all reads) gives the
array-average promoter methylation; MAPQ>=20 restricts to flank-anchored reads.
Reports per-unit promoter 5mC, the edge->interior Spearman trend, and the
maternal vs paternal comparison.

Paths are read from environment variables (see repository README):
    FIVES_DATA  input derived-data directory (BAM)
    FIVES_REFS  reference fasta / unit-coordinate directory
    FIVES_OUT   output directory (per-CpG tables, figures)
"""
import os, numpy as np, pandas as pd, pysam, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy import stats
plt.rcParams.update({"pdf.fonttype":42,"font.family":"Arial","axes.linewidth":0.8})
ML_THR=128; MINCOV=4
fa=pysam.FastaFile(f'{os.environ.get("FIVES_REFS","refs")}/hg002_rnu2_ref.fa')
units=pd.read_csv(f'{os.environ.get("FIVES_REFS","refs")}/hg002_units.txt',sep=r"\s+",header=None,names=["ctg","start"])
seqs={c:fa.fetch(c).upper() for c in fa.references}
cpg={c:set(i for i in range(len(s)-1) if s[i]=="C" and s[i+1]=="G") for c,s in seqs.items()}

def per_cpg(ctg,mapq):
    bam=pysam.AlignmentFile(f'{os.environ.get("FIVES_DATA","data")}/HG002_ownasm.bam',"rb"); s=seqs[ctg]; cg=cpg[ctg]; nm={}; nt={}
    for r in bam.fetch(ctg):
        if r.is_secondary or r.is_supplementary or r.is_unmapped or r.mapping_quality<mapq: continue
        mods=r.modified_bases
        if not mods: continue
        mc=mods.get(('C',0,'m')) or mods.get(('C',1,'m'))
        if not mc: continue
        q2r=dict(r.get_aligned_pairs(matches_only=True))
        for qpos,ml in mc:
            rp=q2r.get(qpos)
            if rp is None: continue
            if s[rp]=="C" and rp in cg: site=rp
            elif rp>0 and s[rp]=="G" and (rp-1) in cg: site=rp-1
            else: continue
            nt[site]=nt.get(site,0)+1
            if ml>=ML_THR: nm[site]=nm.get(site,0)+1
    return pd.DataFrame([(p,nm.get(p,0),nt[p],100*nm.get(p,0)/nt[p]) for p in nt if nt[p]>=MINCOV],
                        columns=["pos","n_m","n_tot","frac"]).sort_values("pos")

def sliding(bulk, W=600, step=150):
    # coverage-weighted sliding-window mean 5mC (+/- 95% CI) along genomic position
    p = bulk.pos.values.astype(float); f = bulk.frac.values; n = bulk.n_tot.values.astype(float)
    if len(p) == 0: return (np.array([]),) * 4
    xs = np.arange(p.min(), p.max() + 1, step); X, Y, LO, HI = [], [], [], []
    for x in xs:
        m = (p >= x - W / 2) & (p <= x + W / 2); k = int(m.sum())
        if k < 3: continue
        w = n[m]; ff = f[m]; mean = np.average(ff, weights=w)
        sem = np.sqrt(np.average((ff - mean) ** 2, weights=w) / k)
        X.append(x); Y.append(mean); LO.append(mean - 1.96 * sem); HI.append(mean + 1.96 * sem)
    return (np.array(X), np.array(Y), np.clip(LO, 0, 100), np.clip(HI, 0, 100))

# ---- publication style ----
CM = 1 / 2.54
plt.rcParams.update({"pdf.fonttype": 42, "ps.fonttype": 42, "font.size": 8, "axes.linewidth": 0.6,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6, "xtick.major.size": 2, "ytick.major.size": 2})
ACT, SIL = "#2166AC", "#D6604D"            # active (low 5mC) blue ; silent (high) red
from matplotlib.lines import Line2D
PANELS = [("HG002_MAT_rnu2", "Maternal (62 units)"), ("HG002_PAT_rnu2", "Paternal (16 units)")]
fig, axes = plt.subplots(1, 2, figsize=(2 * 8 * CM + 1.4, 4.5 * CM + 1.2), sharey=True,
                         constrained_layout=True)
for ax, (ctg, lab) in zip(axes, PANELS):
    u = units[units.ctg == ctg].start.values
    arr_lo, arr_hi = u.min(), u.max() + 188
    bulk = per_cpg(ctg, 0)
    pm = []
    for us in u:
        w = bulk[(bulk.pos >= us - 100) & (bulk.pos <= us + 288)]
        pm.append(np.average(w.frac, weights=w.n_tot) if len(w) else np.nan)
    pm = np.array(pm)
    idx = np.arange(len(u)); edged = np.minimum(idx, len(u) - 1 - idx); ok = ~np.isnan(pm)
    rho = stats.spearmanr(edged[ok], pm[ok]) if ok.sum() > 5 else (np.nan, np.nan)
    print(f"[{lab}] units={len(u)}  median unit-promoter 5mC={np.nanmedian(pm):.1f}%  "
          f"range {np.nanmin(pm):.0f}-{np.nanmax(pm):.0f}%  active(<20%)={int((pm<20).sum())}/{len(u)}  "
          f"edge->interior Spearman={rho[0]:+.3f} p={rho[1]:.1e}")
    for us in u:                                                                             # transcribed U2 gene only (188 bp/unit)
        ax.axvspan(us / 1000, (us + 188) / 1000, color="#fdbe85", alpha=0.7, lw=0, zorder=0)
    sx, sm, slo, shi = sliding(bulk)                                                          # smoothed profile
    ax.fill_between(sx / 1000, slo, shi, color="#6baed6", alpha=0.35, lw=0, zorder=2)         # 95% CI
    ax.plot(sx / 1000, sm, "-", color="#08519c", lw=0.9, zorder=3)                            # sliding-window mean
    ax.axhline(20, color="0.7", ls=":", lw=0.6, zorder=1)                                    # active threshold
    ax.set_ylim(-3, 103); ax.set_xlim((arr_lo - 45000) / 1000, (arr_hi + 45000) / 1000)
    ax.set_title(f"{lab} — median unit-promoter 5mC {np.nanmedian(pm):.1f}%", fontsize=7.5)
    ax.set_xlabel("HG002 chr17 position (kb)", fontsize=7.5)
    ax.tick_params(labelsize=7); ax.spines[["top", "right"]].set_visible(False)
    bulk.to_csv(f'{os.environ.get("FIVES_OUT","output")}/rnu2_ownasm_{ctg}.cpg.tsv', sep="\t", index=False)
axes[0].set_ylabel("5mCpG (%)", fontsize=8)
from matplotlib.patches import Patch
axes[1].legend(handles=[
    Line2D([0], [0], color="#08519c", lw=1.2, label="sliding-window mean 5mC (0.6 kb)"),
    Patch(facecolor="#6baed6", alpha=0.35, edgecolor="none", label="95% CI"),
], fontsize=6, frameon=False, loc='center right', handletextpad=0.5)
fig.suptitle("RNU2 (U2 snRNA) array — per-unit promoter methylation on HG002's own assembly", fontsize=8)
fig.savefig(f'{os.environ.get("FIVES_OUT","output")}/rnu2_ownasm_profile.pdf', dpi=300, bbox_inches="tight")
fig.savefig(f'{os.environ.get("FIVES_OUT","output")}/rnu2_ownasm_profile.png', dpi=200, bbox_inches="tight")
print("wrote rnu2_ownasm_profile.pdf")
