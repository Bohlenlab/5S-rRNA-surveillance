#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# 98_binary_classifications_two_ways.py — per-cohort carrier frequency under two binary variant classifications (incorporation-deficient and expression-deficient, each vs all others).
# Author: Jonathan Bohlen
# (c) 2026 Jonathan Bohlen. Code accompanying Sengl et al. (2026),
# "Surveillance and selection of 5S ribosomal RNA genes in the human genome."
# Released under the MIT License; see LICENSE at the repository root.
# -----------------------------------------------------------------------------
"""
98_binary_classifications_two_ways.py

TWO independent binary classifications, each "deficient vs ALL others", evaluated separately in
each cohort.

  Classification 1 (incorporation): incorp < 0.5  vs  incorp >= 0.5
                                    (over the 284 variants with incorporation measured)
  Classification 2 (expression):    expr  < 0.5   vs  expr  >= 0.5
                                    (over all 342 assayed variants)

Each classification keeps ALL other variants in the "others" group (no third bin,
no expression restriction). Per-individual carrier frequency (VAF >= VAF_CONF).

Outputs (<FIVES_OUT>/06_population_genetics/):
  98_binary_two_ways.pdf       – row1 incorporation contrast, row2 expression contrast
  tables/98_master_per_variant.tsv
  tables/98_incorporation_contrast.tsv
  tables/98_expression_contrast.tsv
"""

import os
import sqlite3
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.transforms as mtransforms
from matplotlib.patches import Rectangle
from scipy.stats import mannwhitneyu, spearmanr, gaussian_kde

plt.rcParams.update({"pdf.fonttype": 42, "font.size": 9})

DB     = Path(os.environ.get("FIVES_DB", "5S_rDNA.db"))
OUTDIR = Path(os.environ.get("FIVES_OUT", "output")) / "06_population_genetics"
TABLES = OUTDIR / "tables"; TABLES.mkdir(parents=True, exist_ok=True)

VAF_CONF   = 0.003
DEF_CUT    = 0.50
N_UKBB     = 490_075
GENE_LO, GENE_HI = 633, 746
COH = ["UKBB", "HPRC", "GTEx"]
COH_COL = {"UKBB": "#33608c", "HPRC": "#2ca25f", "GTEx": "#b8860b"}
COL_INC, COL_EXP, COL_OTH = "#d62728", "#9467bd", "#9aa0a6"

# ── data prep (identical metrics to script 97) ────────────────────────────────
con = sqlite3.connect(DB)
fa = pd.read_sql_query("""SELECT consensus_pos AS pos, ref_base AS ref, alt_base AS alt,
                                 rna_expr_mean AS expr, incorp_60s_mean AS incorp
                          FROM functional_annotation""", con)

# UKBB per-individual carriers
uk = con.execute("""SELECT t2t_pos,t2t_ref,t2t_alt,vaf_array FROM ukbb_population_variants
                    WHERE t2t_pos BETWEEN ? AND ?""", (GENE_LO, GENE_HI)).fetchall()
uk_map = {(p, r, a): int((np.frombuffer(b, np.float32) >= VAF_CONF).sum()) for p, r, a, b in uk}
fa["UKBB"] = [uk_map.get((p, r, a), 0) / N_UKBB for p, r, a in zip(fa.pos, fa.ref, fa.alt)]

# HPRC per-individual (within-individual VAF >= VAF_CONF)
asm_tot = pd.read_sql_query("""SELECT a.assembly_id AS asm, COUNT(*) AS tot
    FROM copy c JOIN haplotype h ON h.haplotype_id=c.haplotype_id
    JOIN assembly a ON a.assembly_id=h.assembly_id
    WHERE a.cohort LIKE 'HPRC%' AND c.array_member=1 GROUP BY a.assembly_id""", con)
N_HPRC = len(asm_tot); tot = dict(zip(asm_tot.asm, asm_tot.tot))
hp = pd.read_sql_query(f"""SELECT v.consensus_pos AS pos, v.alt AS alt, a.assembly_id AS asm,
       COUNT(DISTINCT v.copy_id) AS carr
    FROM variant v JOIN copy c ON c.copy_id=v.copy_id
    JOIN haplotype h ON h.haplotype_id=c.haplotype_id
    JOIN assembly a ON a.assembly_id=h.assembly_id
    WHERE v.alignment_source='consensus_t2t' AND a.cohort LIKE 'HPRC%' AND c.array_member=1
      AND v.consensus_pos BETWEEN {GENE_LO} AND {GENE_HI}
    GROUP BY v.consensus_pos, v.alt, a.assembly_id""", con)
hp["vaf"] = [c / tot[a] for c, a in zip(hp.carr, hp.asm)]
hp_carr = hp[hp.vaf >= VAF_CONF].groupby(["pos", "alt"]).size()
fa["HPRC"] = [int(hp_carr.get((p, a), 0)) / N_HPRC for p, a in zip(fa.pos, fa.alt)]

# GTEx per-individual carriers
N_GTEX = con.execute("SELECT COUNT(*) FROM assembly WHERE cohort='GTEx_v9_WGS'").fetchone()[0]
gt = pd.read_sql_query(f"""SELECT rv.consensus_pos AS pos, rv.alt AS alt,
       COUNT(DISTINCT rv.assembly_id) AS carr
    FROM read_variant rv JOIN assembly a ON a.assembly_id=rv.assembly_id
    WHERE a.cohort='GTEx_v9_WGS' AND rv.modality='illumina' AND rv.vaf >= {VAF_CONF}
      AND rv.consensus_pos BETWEEN {GENE_LO} AND {GENE_HI}
    GROUP BY rv.consensus_pos, rv.alt""", con)
gt_map = {(p, a): n for p, a, n in zip(gt.pos, gt.alt, gt.carr)}
fa["GTEx"] = [gt_map.get((p, a), 0) / N_GTEX for p, a in zip(fa.pos, fa.alt)]
con.close()

# binary class flags
fa["incorp_def"] = fa.incorp < DEF_CUT     # NaN -> False
fa["expr_def"]   = fa.expr   < DEF_CUT
fa.to_csv(TABLES / "98_master_per_variant.tsv", sep="\t", index=False)
print(f"denominators per individual: UKBB={N_UKBB:,} HPRC={N_HPRC} GTEx={N_GTEX} (VAF≥{VAF_CONF*100:.1f}%)")

FLOOR = {c: max(fa[c][fa[c] > 0].min() / 2, 1e-7) for c in COH}
def logf(x, c): return np.log10(np.where(x > 0, x, FLOOR[c]))

def tviolin(ax, x0, y, color, width=0.36):
    """Violin trimmed to the 10th-90th percentile (deciles); thin line = P10-P90,
    thick bar = IQR (Q25-Q75), white dot = median. KDE estimated on full data."""
    y = np.asarray(y, float)
    p10, q25, med, q75, p90 = np.percentile(y, [10, 25, 50, 75, 90])
    if (p90 - p10) > 1e-6 and len(np.unique(y)) >= 3 and y.std() > 1e-6:
        kde = gaussian_kde(y)
        ys = np.linspace(p10, p90, 160)              # body trimmed to deciles
        d = kde(ys); d = d / d.max() * width
        ax.fill_betweenx(ys, x0 - d, x0 + d, color=color, alpha=0.55,
                         lw=0.6, edgecolor=color, zorder=2)
    else:                                  # degenerate (e.g. all at floor)
        ax.add_patch(Rectangle((x0 - width*0.35, p10), width*0.7,
                               max(p90 - p10, 0.03), fc=color, alpha=0.55, lw=0, zorder=2))
    ax.plot([x0, x0], [p10, p90], color="k", lw=0.8, zorder=3)      # 10th-90th decile range
    ax.plot([x0, x0], [q25, q75], color="k", lw=4, solid_capstyle="round", zorder=3)  # IQR
    ax.scatter([x0], [med], s=16, color="w", edgecolor="k", lw=0.8, zorder=4)         # median
def stars(p): return "***" if p < .001 else "**" if p < .01 else "*" if p < .05 else "n.s."
def pstr(p):  return "P<1e-4" if p < 1e-4 else f"P={p:.3g}"

# ── two classifications, each deficient vs ALL others ─────────────────────────
CLASS = [
    ("incorporation", "incorp_def", "incorp", COL_INC,
     fa[fa.incorp.notna()].copy()),    # only incorp-measured variants
    ("expression",    "expr_def",   "expr",   COL_EXP,
     fa.copy()),                        # all variants
]

res = {}   # (contrast, cohort) -> dict
for cname, flag, score, col, dset in CLASS:
    for c in COH:
        defv = dset[dset[flag]][c].values
        oth  = dset[~dset[flag]][c].values
        _, p = mannwhitneyu(defv, oth, alternative="less")
        md, mo = np.median(defv), np.median(oth)
        depl = 100 * (mo - md) / mo if mo > 0 else np.nan
        rho, prho = spearmanr(dset[score].values, dset[c].values)
        res[(cname, c)] = dict(defv=defv, oth=oth, n_def=len(defv), n_oth=len(oth),
                               med_def=md, med_oth=mo, depl=depl, p=p, rho=rho, prho=prho)

# ── figure: rows = contrasts, cols = cohorts ──────────────────────────────────
rng = np.random.default_rng(11)
fig = plt.figure(figsize=(10.6, 8.4))
gs = gridspec.GridSpec(2, 3, figure=fig, hspace=0.50, wspace=0.32,
                       top=0.86, bottom=0.07, left=0.09, right=0.975)
fig.suptitle("Two binary classifications of 5S variants, each deficient vs all others",
             fontsize=12.5, fontweight="bold", y=0.965)
fig.text(0.5, 0.925, "per-cohort carrier frequency: deficient class vs all others",
         ha="center", fontsize=8.6, style="italic", color="#555")

ROW_LABEL = {0: "Incorporation-deficient\nvs all others",
             1: "Expression-deficient\nvs all others"}
for ri, (cname, flag, score, col, dset) in enumerate(CLASS):
    for ci, c in enumerate(COH):
        ax = fig.add_subplot(gs[ri, ci])
        r = res[(cname, c)]
        groups = [(f"{cname[:6].capitalize()}-\ndeficient", r["defv"], col),
                  ("All others", r["oth"], COL_OTH)]
        allyt = []
        for xc, (lab, arr, cc) in enumerate(groups):
            yt = logf(arr, c)
            tviolin(ax, xc, yt, cc)
            allyt.append(yt)
        # axis spans the shown decile ranges (P10..P90) across groups, plus label headroom
        ylo = min(np.percentile(v, 10) for v in allyt)
        yhi = max(np.percentile(v, 90) for v in allyt)
        rng_ = max(yhi - ylo, 0.5)
        ax.set_ylim(ylo - 0.10*rng_, yhi + 0.40*rng_)
        tb = mtransforms.blended_transform_factory(ax.transData, ax.transAxes)
        di = f"{'−' if r['depl']>=0 else '+'}{abs(r['depl']):.0f}%" if np.isfinite(r["depl"]) else "—"
        ax.text(0, 0.985, f"{stars(r['p'])}  {di}\n{pstr(r['p'])}  ρ={r['rho']:.2f}",
                transform=tb, ha="center", va="top", fontsize=7.6, fontweight="bold", color=col)
        ax.text(1, 0.985, "(ref)", transform=tb, ha="center", va="top", fontsize=7.5, color="#555")
        ax.set_xticks([0, 1]); ax.set_xticklabels([g[0] for g in groups], fontsize=8)
        ax.set_xlim(-.6, 1.6)
        if ri == 0:
            ax.set_title(c, color=COH_COL[c], fontweight="bold", fontsize=11, pad=6)
        if ci == 0:
            ax.set_ylabel(ROW_LABEL[ri] + f"\n\nper-individual freq (VAF≥{VAF_CONF*100:.1f}%, log)",
                          fontsize=8.5)
        ax.grid(axis="y", lw=.3, alpha=.4)

fig.savefig(OUTDIR / "98_binary_two_ways.pdf", bbox_inches="tight")
plt.close(fig)
print("wrote 98_binary_two_ways.pdf")

# ── tables ────────────────────────────────────────────────────────────────────
for cname, fname in [("incorporation", "98_incorporation_contrast.tsv"),
                     ("expression",    "98_expression_contrast.tsv")]:
    rows = []
    for c in COH:
        r = res[(cname, c)]
        rows.append(dict(cohort=c, contrast=f"{cname}-deficient vs all others",
                         n_deficient=r["n_def"], n_others=r["n_oth"],
                         median_freq_deficient=r["med_def"], median_freq_others=r["med_oth"],
                         depletion_pct=r["depl"], MW_p_one_sided=r["p"],
                         spearman_rho=r["rho"], spearman_p=r["prho"]))
    pd.DataFrame(rows).to_csv(TABLES / fname, sep="\t", index=False)
print("wrote 98 contrast tables")

print("\n=== summary (deficient vs all others) ===")
for cname, *_ in CLASS:
    print(f"-- {cname} --")
    for c in COH:
        r = res[(cname, c)]
        print(f"   {c:5s} depl={r['depl'] if np.isfinite(r['depl']) else float('nan'):6.1f}%  "
              f"MW p={r['p']:.2e} {stars(r['p'])}  rho={r['rho']:+.3f} (p={r['prho']:.2e})")
