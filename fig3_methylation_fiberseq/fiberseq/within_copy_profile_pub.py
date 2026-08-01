#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# within_copy_profile_pub.py — Within-copy m6A accessibility profiles by copy class (Fiber-seq).
# Author: Jonathan Bohlen
# (c) 2026 Jonathan Bohlen. Code accompanying Sengl et al. (2026),
# "Surveillance and selection of 5S ribosomal RNA genes in the human genome."
# Released under the MIT License; see LICENSE at the repository root.
# -----------------------------------------------------------------------------
"""Within-copy m6A accessibility by copy class, both classes (active vs silent)
overlaid in a single panel. Copies are split at the 25th percentile of per-copy
CpG methylation. Optional argument "own" restricts to own-assembly samples.
Output: figures/within_copy_accessibility_by_class<TAG>_pub.pdf

Paths are read from environment variables (see repository README):
    FIVES_DATA  input derived-data directory
    FIVES_REFS  reference fasta directory
    FIVES_OUT   output directory
"""
import os, sys, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.size": 8, "axes.titlesize": 8, "axes.labelsize": 8,
                     "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 6,
                     "axes.linewidth": 0.6, "pdf.fonttype": 42, "ps.fonttype": 42})
ALLSAMP = ["HG002", "CHM1", "CHM13", "GM12878", "K562"]; OWNASM = ["CHM13", "HG002"]
SAMPLES = OWNASM if (len(sys.argv) > 1 and sys.argv[1] == "own") else ALLSAMP
TAG = "_ownasm" if SAMPLES == OWNASM else ""
ULEN, BINW = 2168, 15
GENE_S, GENE_E, ALU_S, ALU_E = 630, 748, 787, 1066
MIN_A, MIN_CPG = 10, 5
pctf = plt.FuncFormatter(lambda v, _: f"{v:.0f}%")

df = []
for S in SAMPLES:
    d = pd.read_csv(f'{os.environ.get("FIVES_DATA","data")}/{S}/{S}.tandemcalls.tsv.gz', sep="\t"); d["sample"] = S; df.append(d)
df = pd.concat(df, ignore_index=True)
df["meth"] = (df.qual >= 0.8).astype(int)
df["within"] = df.refpos % ULEN
df["copy"] = df["sample"] + "_" + df["read"].astype(str) + "_" + (df.refpos // ULEN).astype(str)

seq = "".join(l.strip() for l in open(f'{os.environ.get("FIVES_REFS","refs")}/5S_t2t_consensus.fa') if not l.startswith(">")).upper()
CPG = set(i for i in range(len(seq) - 1) if seq[i] == "C" and seq[i + 1] == "G")
CPG |= set(i + 1 for i in range(len(seq) - 1) if seq[i] == "C" and seq[i + 1] == "G")

a = df[df.code == "a"].copy()
m = df[(df.code == "m") & (df.within.isin(CPG))].copy()
acc = a.groupby("copy").agg(na=("meth", "size"))
met = m.groupby("copy").agg(m5=("meth", "mean"), nm=("meth", "size"))
g = acc.join(met, how="inner"); g = g[(g.na >= MIN_A) & (g.nm >= MIN_CPG)].copy()
THR = g.m5.quantile(0.25); g["cls"] = np.where(g.m5 < THR, "active", "silent")

a = a[a["copy"].isin(g.index)].copy(); a["bin"] = (a.within // BINW) * BINW
cb = a.groupby(["copy", "bin"]).meth.mean().reset_index().merge(g[["cls"]], left_on="copy", right_index=True)

def prof(cls):
    # mean fraction of molecules accessible per bin (+ bootstrap-free 95% CI of the mean);
    # per-copy percentiles are NOT used because single-molecule m6A is zero-inflated (~75% of
    # copies read 0 at any position) so the median/quartiles collapse to 0 and carry no signal.
    rows = []
    for b, d in cb[cb.cls == cls].groupby("bin"):
        v = d.meth.values * 100
        if len(v) < 3: continue
        se = v.std() / np.sqrt(len(v))
        rows.append((b, v.mean(), v.mean() - 1.96 * se, v.mean() + 1.96 * se))
    p = pd.DataFrame(rows, columns=["bin", "mean", "lo", "hi"]).sort_values("bin")
    for c in ("mean", "lo", "hi"):
        p[c] = p[c].rolling(3, center=True, min_periods=1).mean()
    return p

# colour scheme: methylated/silent = blue, hypomethylated/active = red
CLS = [("silent", "#2166AC", "Silent / high-CpG (top 75%)"),
       ("active", "#D6604D", "Active / low-CpG (bottom 25%)")]
YMAX = 35
fig, axes = plt.subplots(1, 2, figsize=(8, 4), dpi=300, sharey=True, gridspec_kw=dict(wspace=0.08))
for ax, (cls, col, lab) in zip(axes, CLS):
    p = prof(cls); n = int((g.cls == cls).sum())
    ax.axvspan(GENE_S, GENE_E, color="#aec6cf", alpha=0.30, lw=0)
    ax.axvspan(ALU_S, ALU_E, color="#c8a0e8", alpha=0.30, lw=0)
    ax.text((GENE_S+GENE_E)/2, YMAX*1.005, "5S gene", ha="center", va="bottom", fontsize=5.5, color="#555")
    ax.text((ALU_S+ALU_E)/2, YMAX*1.005, "ALU", ha="center", va="bottom", fontsize=5.5, color="#555")
    ax.fill_between(p.bin, p.lo, p.hi, color=col, alpha=0.30, lw=0, label="95% CI of mean")
    ax.plot(p.bin, p["mean"], color=col, lw=1.0, label="Mean accessibility")
    ax.set_xlim(0, ULEN); ax.set_ylim(0, YMAX); ax.yaxis.set_major_formatter(pctf)
    ax.set_xlabel("Position within repeat unit (bp)")
    ax.set_title(f"{lab}  ({n:,} copies)", fontsize=7)
    ax.tick_params(length=2, width=0.6); ax.grid(lw=0.3, alpha=0.35, axis="y")
    ax.legend(loc="upper right", frameon=False, handlelength=1.4, fontsize=5.5)
axes[0].set_ylabel("m6A accessibility (%)")
fig.tight_layout()
out = f'{os.environ.get("FIVES_OUT","output")}/within_copy_accessibility_by_class{TAG}_pub.pdf'
fig.savefig(out, bbox_inches="tight"); plt.close(fig)
print(f"wrote {out}  (4x4in, 300dpi, font8, 1pt; split 25th-pctile CpG = {THR*100:.0f}%; "
      f"active={int((g.cls=='active').sum())} silent={int((g.cls=='silent').sum())})")
