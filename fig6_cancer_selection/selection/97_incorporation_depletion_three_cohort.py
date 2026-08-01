#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# 97_incorporation_depletion_three_cohort.py — per-variant carrier frequency by functional class (incorporation vs expression) across UKBB, HPRC, and GTEx.
# Author: Jonathan Bohlen
# (c) 2026 Jonathan Bohlen. Code accompanying Sengl et al. (2026),
# "Surveillance and selection of 5S ribosomal RNA genes in the human genome."
# Released under the MIT License; see LICENSE at the repository root.
# -----------------------------------------------------------------------------
"""
97_incorporation_depletion_three_cohort.py

Per-variant carrier-frequency analysis of 5S rRNA gene variants in three human
cohorts (UKBB, HPRC, GTEx), stratified by two functional axes (incorporation and
expression).

Two functional axes per assayed gene variant (functional_annotation):
  incorp_60s_mean  – relative incorporation into the 60S subunit (WT = 1)
  rna_expr_mean    – relative expression                          (WT = 1)
Defective cut for visualisation = < 0.50 (continuous scores used for all stats).

Population frequency (confident-carrier definition, VAF >= 0.30%) per cohort:
  UKBB  carriers (vaf_array >= 0.003)        / 490,075 samples
  HPRC  carrier array copies (array_member=1) / 43,724 copies
  GTEx  carrier samples (read_variant illumina, already confident) / 944 samples

Outputs (figures/06_population_genetics/):
  97_incorporation_depletion_MAIN.pdf   – B: 3-cohort depletion ; C: specificity vs expression
  97_incorporation_depletion_SUPP.pdf   – A: landscape ; D: dose-response ; robustness ; per-cohort detail
  97_master_table.tsv                   – per-variant merged data + per-cohort frequency
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
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D
from scipy.stats import mannwhitneyu, spearmanr
import statsmodels.api as sm

plt.rcParams.update({"pdf.fonttype": 42, "font.size": 9})

# ── paths / constants ─────────────────────────────────────────────────────────
DB     = Path(os.environ.get("FIVES_DB", "5S_rDNA.db"))
OUTDIR = Path(os.environ.get("FIVES_OUT", "output")) / "06_population_genetics"
OUTDIR.mkdir(parents=True, exist_ok=True)
TABLES = OUTDIR / "tables"
TABLES.mkdir(parents=True, exist_ok=True)

VAF_CONF   = 0.003          # confident-carrier VAF threshold, applied PER INDIVIDUAL in all
                            # three cohorts. HPRC: an individual carries if >=VAF_CONF of THEIR
                            # copies have the variant.
DEF_CUT    = 0.50           # defective cutoff for visualisation only
N_UKBB     = 490_075
GENE_LO, GENE_HI = 633, 746

COH = ["UKBB", "HPRC", "GTEx"]
COH_COL = {"UKBB": "#33608c", "HPRC": "#2ca25f", "GTEx": "#b8860b"}
COL_DEF, COL_NOR, COL_EXP = "#d62728", "#1f77b4", "#9467bd"

# ── 1. functional annotation (the two assay axes) ─────────────────────────────
con = sqlite3.connect(DB)
fa = pd.read_sql_query("""
    SELECT consensus_pos AS pos, ref_base AS ref, alt_base AS alt,
           rna_expr_mean AS expr, incorp_60s_mean AS incorp
    FROM functional_annotation
""", con)
print(f"functional_annotation: {len(fa)} variants "
      f"({fa.incorp.notna().sum()} with incorporation measured, "
      f"{fa.incorp.isna().sum()} incorp-NA = expression-defective extreme)")

# ── 2a. UKBB confident-carrier frequency ──────────────────────────────────────
uk = con.execute("""
    SELECT t2t_pos, t2t_ref, t2t_alt, vaf_array
    FROM ukbb_population_variants
    WHERE t2t_pos BETWEEN ? AND ?
""", (GENE_LO, GENE_HI)).fetchall()
uk_map = {}
for pos, ref, alt, blob in uk:
    vafs = np.frombuffer(blob, dtype=np.float32)
    uk_map[(pos, ref, alt)] = int((vafs >= VAF_CONF).sum())
fa["ukbb_carriers"] = [uk_map.get((p, r, a), 0) for p, r, a in zip(fa.pos, fa.ref, fa.alt)]
fa["UKBB"] = fa["ukbb_carriers"] / N_UKBB

# ── 2b. HPRC PER-INDIVIDUAL carrier frequency (matches UKBB/GTEx unit) ─────────
# Each HPRC assembly = one (diploid) individual. Within-individual VAF =
# carrier copies / that individual's total array copies; carrier if >= VAF_CONF.
asm_tot = pd.read_sql_query("""
    SELECT a.assembly_id AS asm, COUNT(*) AS tot
    FROM copy c
    JOIN haplotype h ON h.haplotype_id=c.haplotype_id
    JOIN assembly  a ON a.assembly_id=h.assembly_id
    WHERE a.cohort LIKE 'HPRC%' AND c.array_member=1
    GROUP BY a.assembly_id
""", con)
N_HPRC_IND = len(asm_tot)
tot = dict(zip(asm_tot.asm, asm_tot.tot))
hp = pd.read_sql_query(f"""
    SELECT v.consensus_pos AS pos, v.alt AS alt, a.assembly_id AS asm,
           COUNT(DISTINCT v.copy_id) AS carr
    FROM variant v
    JOIN copy c       ON c.copy_id=v.copy_id
    JOIN haplotype h  ON h.haplotype_id=c.haplotype_id
    JOIN assembly  a  ON a.assembly_id=h.assembly_id
    WHERE v.alignment_source='consensus_t2t' AND a.cohort LIKE 'HPRC%'
      AND c.array_member=1
      AND v.consensus_pos BETWEEN {GENE_LO} AND {GENE_HI}
    GROUP BY v.consensus_pos, v.alt, a.assembly_id
""", con)
hp["indiv_vaf"] = [carr / tot[a] for carr, a in zip(hp.carr, hp.asm)]
hp_carr = hp[hp.indiv_vaf >= VAF_CONF].groupby(["pos", "alt"]).size()
fa["hprc_carriers"] = [int(hp_carr.get((p, a), 0)) for p, a in zip(fa.pos, fa.alt)]
fa["HPRC"] = fa["hprc_carriers"] / N_HPRC_IND
HPRC_TOT = N_HPRC_IND

# ── 2c. GTEx carrier-sample frequency ─────────────────────────────────────────
GTEX_TOT = con.execute("SELECT COUNT(*) FROM assembly WHERE cohort='GTEx_v9_WGS'").fetchone()[0]
gt = pd.read_sql_query(f"""
    SELECT rv.consensus_pos AS pos, rv.alt AS alt,
           COUNT(DISTINCT rv.assembly_id) AS carriers
    FROM read_variant rv
    JOIN assembly a ON a.assembly_id=rv.assembly_id
    WHERE a.cohort='GTEx_v9_WGS' AND rv.modality='illumina'
      AND rv.vaf >= {VAF_CONF}
      AND rv.consensus_pos BETWEEN {GENE_LO} AND {GENE_HI}
    GROUP BY rv.consensus_pos, rv.alt
""", con)
gt_map = {(p, a): n for p, a, n in zip(gt.pos, gt.alt, gt.carriers)}
fa["gtex_carriers"] = [gt_map.get((p, a), 0) for p, a in zip(fa.pos, fa.alt)]
fa["GTEx"] = fa["gtex_carriers"] / GTEX_TOT
con.close()

print(f"per-individual denominators: UKBB={N_UKBB:,}  HPRC={HPRC_TOT}  GTEx={GTEX_TOT}  "
      f"(VAF_CONF={VAF_CONF*100:.1f}%)")

# ── 3. classes (viz only) ─────────────────────────────────────────────────────
fa["expressed"]   = fa.expr  >= DEF_CUT
fa["incorp_def"]  = fa.incorp < DEF_CUT          # NaN -> False (handled via expr-def)
fa["expr_def"]    = fa.expr  < DEF_CUT
def klass(r):
    if r.expr_def:                       return "expr_defective"
    if r.incorp_def:                     return "incorp_defective"   # expressed & low incorp
    return "competent"
fa["class"] = fa.apply(klass, axis=1)

# expressed subset used for the incorporation contrasts (isolates incorporation)
expr_ok = fa[fa.expressed & fa.incorp.notna()].copy()

fa.to_csv(OUTDIR / "97_master_table.tsv", sep="\t", index=False)
print("wrote 97_master_table.tsv")

# ── 4. per-cohort statistics ──────────────────────────────────────────────────
FLOOR = {c: max(fa[c][fa[c] > 0].min() / 2, 1e-7) for c in COH}   # log floor per cohort

def logf(x, c):
    return np.log10(np.where(x > 0, x, FLOOR[c]))

def tviolin(ax, x0, y, color, width=0.36):
    """Violin trimmed to the 10th-90th percentile (deciles); thin line = P10-P90,
    thick bar = IQR (Q25-Q75), white dot = median. KDE estimated on full data."""
    from scipy.stats import gaussian_kde
    y = np.asarray(y, float)
    p10, q25, med, q75, p90 = np.percentile(y, [10, 25, 50, 75, 90])
    if (p90 - p10) > 1e-6 and len(np.unique(y)) >= 3 and y.std() > 1e-6:
        kde = gaussian_kde(y)
        ys = np.linspace(p10, p90, 160)              # body trimmed to deciles
        d = kde(ys); d = d / d.max() * width
        ax.fill_betweenx(ys, x0 - d, x0 + d, color=color, alpha=0.55,
                         lw=0.6, edgecolor=color, zorder=2)
    else:
        ax.add_patch(Rectangle((x0 - width*0.35, p10), width*0.7,
                               max(p90 - p10, 0.03), fc=color, alpha=0.55, lw=0, zorder=2))
    ax.plot([x0, x0], [p10, p90], color="k", lw=0.8, zorder=3)   # 10th-90th decile range
    ax.plot([x0, x0], [q25, q75], color="k", lw=4, solid_capstyle="round", zorder=3)
    ax.scatter([x0], [med], s=16, color="w", edgecolor="k", lw=0.8, zorder=4)

stats = {}
for c in COH:
    comp = fa[fa["class"] == "competent"][c].values          # reference class
    inc  = fa[fa["class"] == "incorp_defective"][c].values   # incorporation-defective class
    exp  = fa[fa["class"] == "expr_defective"][c].values     # expression-defective class
    _, p_inc_mw = mannwhitneyu(inc, comp, alternative="less")
    _, p_exp_mw = mannwhitneyu(exp, comp, alternative="less")
    depl_inc = 100 * (np.median(comp) - np.median(inc)) / np.median(comp) if np.median(comp) > 0 else np.nan
    depl_exp = 100 * (np.median(comp) - np.median(exp)) / np.median(comp) if np.median(comp) > 0 else np.nan
    rho, p_rho = spearmanr(expr_ok.incorp.values, expr_ok[c].values)   # continuous, expressed
    # specificity regression on FULL assayed set (both axes), controlling expression
    reg = fa[fa.incorp.notna()].copy()
    X = sm.add_constant(reg[["incorp", "expr"]].values)
    y = logf(reg[c].values, c)
    m = sm.OLS(y, X).fit()
    stats[c] = dict(comp=comp, inc=inc, exp=exp,
                    p_inc_mw=p_inc_mw, p_exp_mw=p_exp_mw,
                    depl_inc=depl_inc, depl_exp=depl_exp, rho=rho, p_rho=p_rho,
                    b_inc=m.params[1], se_inc=m.bse[1], p_inc=m.pvalues[1],
                    b_exp=m.params[2], se_exp=m.bse[2])

def pstr(p):
    return "P<1e-4" if p < 1e-4 else f"P={p:.3g}"
def stars(p):
    return "***" if p < .001 else "**" if p < .01 else "*" if p < .05 else "n.s."

# ════════════════════════════ MAIN FIGURE ════════════════════════════════════
import matplotlib.transforms as mtransforms
rng = np.random.default_rng(7)
fig = plt.figure(figsize=(10.6, 8.8))
gs  = gridspec.GridSpec(2, 3, figure=fig, hspace=0.46, wspace=0.34,
                        top=0.83, bottom=0.085, left=0.085, right=0.975,
                        height_ratios=[1, 1])
fig.suptitle("Per-individual carrier frequency by functional variant class (three cohorts)",
             fontsize=12.5, fontweight="bold", x=0.5, y=0.975)
fig.text(0.5, 0.945,
         "functional classes: competent, incorporation-defective, expression-defective",
         ha="center", fontsize=8.6, style="italic", color="#555")

# ── Row 1 (Panel B): three-cohort, three functional classes ───────────────────
# competent (reference) / incorporation-defective / expression-defective
for j, c in enumerate(COH):
    ax = fig.add_subplot(gs[0, j])
    s = stats[c]
    groups = [("Competent", s["comp"], COL_NOR),
              ("Expr-\ndefective", s["exp"], COL_EXP),
              ("Incorp-\ndefective", s["inc"], COL_DEF)]
    allyt = []
    for xc, (lab, arr, col) in enumerate(groups):
        yt = logf(arr, c)
        tviolin(ax, xc, yt, col)
        allyt.append(yt)
    # axis spans the shown decile ranges (P10..P90) across groups, plus label headroom
    ylo = min(np.percentile(v, 10) for v in allyt)
    yhi = max(np.percentile(v, 90) for v in allyt)
    rng_ = max(yhi - ylo, 0.5)
    ax.set_ylim(ylo - 0.10*rng_, yhi + 0.42*rng_)
    tb = mtransforms.blended_transform_factory(ax.transData, ax.transAxes)
    di = f"−{s['depl_inc']:.0f}%" if np.isfinite(s["depl_inc"]) else "depleted"
    # one stats box per panel: continuous incorporation test (Spearman) + incorp-def depletion
    ax.text(0.03, 0.985,
            f"incorp. gradient\nSpearman ρ = {s['rho']:.2f}\n{pstr(s['p_rho'])}\nincorp-def: {di}",
            transform=ax.transAxes, ha="left", va="top", fontsize=7.4,
            bbox=dict(boxstyle="round,pad=.28", fc="#fff8d6", ec="goldenrod", alpha=.92))
    ax.set_xticks([0, 1, 2]); ax.set_xticklabels([g[0] for g in groups], fontsize=8)
    ax.set_xlim(-.6, 2.6)
    ax.set_title(f"{c}", color=COH_COL[c], fontweight="bold", fontsize=11, pad=6)
    if j == 0:
        ax.set_ylabel(f"per-individual carrier frequency\n(VAF ≥ {VAF_CONF*100:.1f}%, log scale)")
    ax.grid(axis="y", lw=.3, alpha=.4)
fig.text(0.04, 0.85, "B", fontsize=15, fontweight="bold")

# ── Row 2 left+mid (Panel C-left): 2D specificity heatmap (pooled, z-scored freq)
axh = fig.add_subplot(gs[1, 0:2])
reg = fa[fa.incorp.notna()].copy()
# pooled standardized frequency = mean of per-cohort z-scores of log freq
z = np.zeros(len(reg))
for c in COH:
    lf = logf(reg[c].values, c)
    z += (lf - lf.mean()) / lf.std()
reg["z"] = z / len(COH)
ex_edges = np.quantile(reg.expr,   np.linspace(0, 1, 5))
in_edges = np.quantile(reg.incorp, np.linspace(0, 1, 5))
H = np.full((4, 4), np.nan)
for ii in range(4):
    for jj in range(4):
        m = ((reg.incorp >= in_edges[ii]) & (reg.incorp <= in_edges[ii+1]) &
             (reg.expr   >= ex_edges[jj]) & (reg.expr   <= ex_edges[jj+1]))
        if m.sum() >= 2: H[ii, jj] = reg.z[m].mean()
im = axh.imshow(H, origin="lower", cmap="RdBu", vmin=-1, vmax=1, aspect="auto")
axh.set_xticks(range(4)); axh.set_yticks(range(4))
axh.set_xticklabels([f"{ex_edges[j]:.2f}–{ex_edges[j+1]:.2f}" for j in range(4)], fontsize=7)
axh.set_yticklabels([f"{in_edges[i]:.2f}–{in_edges[i+1]:.2f}" for i in range(4)], fontsize=7)
axh.set_xlabel("relative expression (quartile bins)")
axh.set_ylabel("relative incorporation (quartile bins)")
axh.set_title("Standardized frequency by incorporation and expression bins",
              fontsize=9.5)
cb = fig.colorbar(im, ax=axh, shrink=.8, pad=.02)
cb.set_label("pooled standardized\nfrequency (z)", fontsize=8)
fig.text(0.04, 0.45, "C", fontsize=15, fontweight="bold")

# ── Row 2 right (Panel C-right): forest plot of incorporation coefficient ──────
axf = fig.add_subplot(gs[1, 2])
ys = np.arange(len(COH))[::-1]
for y, c in zip(ys, COH):
    s = stats[c]
    axf.errorbar(s["b_inc"], y, xerr=1.96*s["se_inc"], fmt="o", color=COH_COL[c],
                 capsize=3, ms=7, lw=1.8)
    axf.text(s["b_inc"], y+.18, stars(s["p_inc"]), ha="center", fontsize=9, color=COH_COL[c])
axf.axvline(0, color="grey", ls="--", lw=1)
axf.set_yticks(ys); axf.set_yticklabels(COH)
axf.set_xlabel("incorporation effect on log₁₀ frequency\n(OLS, controlling for expression)", fontsize=8)
axf.set_title("Incorporation coefficient\n(expression-adjusted)", fontsize=9.5)
axf.grid(axis="x", lw=.3, alpha=.4)
axf.set_ylim(-.6, len(COH)-.4)

fig.savefig(OUTDIR / "97_incorporation_depletion_MAIN.pdf", bbox_inches="tight")
plt.close(fig)
print("wrote 97_incorporation_depletion_MAIN.pdf")

# ════════════════════════════ SUPPLEMENT ═════════════════════════════════════
figS = plt.figure(figsize=(12.5, 9))
gsS = gridspec.GridSpec(2, 3, figure=figS, hspace=0.40, wspace=0.34,
                        top=0.90, bottom=0.08, left=0.07, right=0.975)
figS.suptitle("Supplement — per-individual carrier frequency by functional class (three cohorts)",
              fontsize=12, fontweight="bold")

# A: functional landscape, coloured by pooled standardized freq
axA = figS.add_subplot(gsS[0, 0])
sc = axA.scatter(fa.expr, fa.incorp, c=[reg.set_index(["pos","alt"]).z.get((p,a), np.nan)
                                        for p,a in zip(fa.pos, fa.alt)],
                 cmap="RdBu", vmin=-1, vmax=1, s=26, edgecolors="k", linewidths=.3)
axA.axvline(DEF_CUT, color="grey", ls=":", lw=1); axA.axhline(DEF_CUT, color="grey", ls=":", lw=1)
axA.set_xlabel("relative expression"); axA.set_ylabel("relative incorporation")
axA.set_title("A  Functional landscape", loc="left", fontsize=10)
axA.set_xlim(0, min(3, fa.expr.max())); axA.set_ylim(0, min(3, fa.incorp.max()))
axA.annotate("incorp-defective", (0.05, 0.30), fontsize=7, color=COL_DEF)
figS.colorbar(sc, ax=axA, shrink=.8, label="pooled freq (z)")

# D: dose-response — freq vs incorporation (expressed), 3 cohorts overlaid
axD = figS.add_subplot(gsS[0, 1])
sub = expr_ok.sort_values("incorp")
edges = np.quantile(sub.incorp, np.linspace(0, 1, 7))
ctr = (edges[:-1] + edges[1:]) / 2
dose_tbl = {}
for c in COH:
    med = [np.median(sub[(sub.incorp >= lo) & (sub.incorp <= hi)][c])
           for lo, hi in zip(edges[:-1], edges[1:])]
    dose_tbl[c] = med
    axD.plot(ctr, med, "o-", color=COH_COL[c], lw=1.8, ms=5, label=c)
axD.axvline(DEF_CUT, color="grey", ls="--", lw=1)
axD.set_yscale("log"); axD.set_xlabel("relative incorporation"); axD.set_ylabel("median frequency")
axD.set_title("D  Dose–response (expressed variants)", loc="left", fontsize=10)
axD.legend(fontsize=8); axD.grid(lw=.3, alpha=.4)

# E: expression vs incorporation coefficients side by side (specificity)
axE = figS.add_subplot(gsS[0, 2])
w = .35
for k, (key, lab, col) in enumerate([("b_inc", "incorporation", "#c0392b"),
                                      ("b_exp", "expression", COL_EXP)]):
    xs = np.arange(len(COH)) + (k-0.5)*w
    vals = [stats[c][key] for c in COH]
    ses  = [stats[c]["se_inc" if key=="b_inc" else "se_exp"] for c in COH]
    axE.bar(xs, vals, w, yerr=[1.96*s for s in ses], color=col, capsize=3, label=lab)
axE.axhline(0, color="k", lw=.8)
axE.set_xticks(range(len(COH))); axE.set_xticklabels(COH)
axE.set_ylabel("effect on log₁₀ frequency"); axE.legend(fontsize=8)
axE.set_title("E  Incorporation vs expression coefficients", loc="left", fontsize=9.5)

# F: robustness — UKBB depletion across VAF thresholds
axF = figS.add_subplot(gsS[1, 0])
con = sqlite3.connect(DB)
uk2 = con.execute("""SELECT t2t_pos,t2t_ref,t2t_alt,vaf_array FROM ukbb_population_variants
                     WHERE t2t_pos BETWEEN ? AND ?""", (GENE_LO, GENE_HI)).fetchall()
con.close()
ukb_blob = {(p, r, a): np.frombuffer(b, np.float32) for p, r, a, b in uk2}
THS = [0.003, 0.005, 0.01, 0.02, 0.03]
depl_th, p_th = [], []
for T in THS:
    car = np.array([(ukb_blob.get((p, r, a), np.array([])) >= T).sum()
                    for p, r, a in zip(expr_ok.pos, expr_ok.ref, expr_ok.alt)], float)
    d = car[expr_ok.incorp_def.values]; n = car[~expr_ok.incorp_def.values]
    depl_th.append(100*(np.median(n)-np.median(d))/np.median(n) if np.median(n)>0 else np.nan)
    p_th.append(mannwhitneyu(d, n, alternative="less")[1])
axF.plot([t*100 for t in THS], depl_th, "o-", color=COH_COL["UKBB"], lw=1.8)
for t, dp, pv in zip(THS, depl_th, p_th):
    axF.annotate(stars(pv), (t*100, dp), fontsize=8, ha="center", va="bottom")
axF.set_xlabel("UKBB VAF threshold (%)"); axF.set_ylabel("median depletion (%)")
axF.set_title("F  UKBB threshold robustness", loc="left", fontsize=10)
axF.grid(lw=.3, alpha=.4)

# G: per-cohort detected fraction (zero-rate) defective vs competent
axG = figS.add_subplot(gsS[1, 1])
w = .35
det_tbl = {}
for k, (mask, lab, col) in enumerate([(expr_ok.incorp_def, "incorp-def", COL_DEF),
                                       (~expr_ok.incorp_def, "competent", COL_NOR)]):
    xs = np.arange(len(COH)) + (k-0.5)*w
    fracs = [(expr_ok[mask][c] > 0).mean() for c in COH]
    det_tbl[lab] = fracs
    axG.bar(xs, fracs, w, color=col, label=lab)
axG.set_xticks(range(len(COH))); axG.set_xticklabels(COH)
axG.set_ylabel("fraction of variants detected"); axG.legend(fontsize=8)
axG.set_title("G  Detection rate (zero-inflation)", loc="left", fontsize=10)

# H: stats summary text
axH = figS.add_subplot(gsS[1, 2]); axH.axis("off")
lines = ["Per-cohort statistics (expressed variants)", ""]
for c in COH:
    s = stats[c]
    lines += [f"{c}:",
              f"  incorp-def depletion {s['depl_inc']:.0f}%  MW {pstr(s['p_inc_mw'])} {stars(s['p_inc_mw'])}",
              f"  expr-def (control)    MW {pstr(s['p_exp_mw'])} {stars(s['p_exp_mw'])}",
              f"  Spearman ρ={s['rho']:.3f} ({pstr(s['p_rho'])})",
              f"  incorp β={s['b_inc']:.2f}±{s['se_inc']:.2f} {stars(s['p_inc'])}"
              f" | expr β={s['b_exp']:.2f}±{s['se_exp']:.2f}",
              ""]
lines += [f"n = {len(expr_ok)} expressed assayed variants",
          f"  ({int(expr_ok.incorp_def.sum())} incorp-defective, "
          f"{int((~expr_ok.incorp_def).sum())} competent)",
          f"frequency = confident carriers (VAF≥{VAF_CONF*100:.1f}%)"]
axH.text(0, 1, "\n".join(lines), va="top", family="monospace", fontsize=8)

figS.savefig(OUTDIR / "97_incorporation_depletion_SUPP.pdf", bbox_inches="tight")
plt.close(figS)
print("wrote 97_incorporation_depletion_SUPP.pdf")

# ════════════════════════ PER-PANEL DATA TABLES ══════════════════════════════
# master (per-variant, all panels A/B/D source) already at 97_master_table.tsv
fa.to_csv(TABLES / "97_master_per_variant.tsv", sep="\t", index=False)

# B — three-class depletion summary (cohort × class)
rowsB = []
for c in COH:
    s = stats[c]
    for cls, arr, pmw, depl in [("competent", s["comp"], np.nan, 0.0),
                                ("incorp_defective", s["inc"], s["p_inc_mw"], s["depl_inc"]),
                                ("expr_defective", s["exp"], s["p_exp_mw"], s["depl_exp"])]:
        rowsB.append(dict(cohort=c, klass=cls, n=len(arr),
                          median_freq=float(np.median(arr)), mean_freq=float(np.mean(arr)),
                          depletion_pct_vs_competent=depl, MW_p_vs_competent=pmw))
pd.DataFrame(rowsB).to_csv(TABLES / "97B_threeclass_summary.tsv", sep="\t", index=False)

# C-left — specificity heatmap matrix (incorp bin × expr bin -> pooled z)
Hrows = []
for ii in range(4):
    for jj in range(4):
        Hrows.append(dict(incorp_bin=f"{in_edges[ii]:.3f}-{in_edges[ii+1]:.3f}",
                          expr_bin=f"{ex_edges[jj]:.3f}-{ex_edges[jj+1]:.3f}",
                          pooled_z=H[ii, jj]))
pd.DataFrame(Hrows).to_csv(TABLES / "97C_heatmap_matrix.tsv", sep="\t", index=False)

# C-right / E — OLS coefficients (log10 freq ~ incorp + expr)
pd.DataFrame([dict(cohort=c, incorp_beta=stats[c]["b_inc"], incorp_se=stats[c]["se_inc"],
                   incorp_p=stats[c]["p_inc"], expr_beta=stats[c]["b_exp"],
                   expr_se=stats[c]["se_exp"], spearman_rho=stats[c]["rho"],
                   spearman_p=stats[c]["p_rho"]) for c in COH]
            ).to_csv(TABLES / "97CE_regression_coeffs.tsv", sep="\t", index=False)

# D — dose-response (incorp bin centre × per-cohort median frequency)
dfD = pd.DataFrame({"incorp_bin_centre": ctr, **{c: dose_tbl[c] for c in COH}})
dfD.to_csv(TABLES / "97D_dose_response.tsv", sep="\t", index=False)

# F — UKBB threshold robustness
pd.DataFrame({"vaf_threshold_pct": [t*100 for t in THS],
              "median_depletion_pct": depl_th, "MW_p": p_th}
            ).to_csv(TABLES / "97F_ukbb_threshold_robustness.tsv", sep="\t", index=False)

# G — detection rate (fraction with freq>0) per cohort × class
pd.DataFrame({"cohort": COH, **{lab: det_tbl[lab] for lab in det_tbl}}
            ).to_csv(TABLES / "97G_detection_rate.tsv", sep="\t", index=False)
print(f"wrote 7 per-panel tables → {TABLES}")

# ── console summary ───────────────────────────────────────────────────────────
print("\n=== summary ===")
for c in COH:
    s = stats[c]
    print(f"{c:5s} incorp-def depl={s['depl_inc']:5.1f}% (p={s['p_inc_mw']:.2e})  "
          f"expr-def(control) p={s['p_exp_mw']:.2e}  "
          f"rho={s['rho']:+.3f} (p={s['p_rho']:.2e})  "
          f"incorp_beta={s['b_inc']:+.3f} (p={s['p_inc']:.2e})")
