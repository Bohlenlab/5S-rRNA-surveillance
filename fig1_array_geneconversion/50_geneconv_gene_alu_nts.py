#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# 50_geneconv_gene_alu_nts.py — Compares the density and contiguity of large
# gene-conversion tracts across the 5S gene, the antisense Alu SINE, and the
# remaining non-transcribed spacer.
# Author: Jonathan Bohlen
# (c) 2026 Jonathan Bohlen. Code accompanying Sengl et al. (2026),
# "Surveillance and selection of 5S ribosomal RNA genes in the human genome."
# Released under the MIT License; see LICENSE at the repository root.
# -----------------------------------------------------------------------------
"""
50_geneconv_gene_alu_nts.py

Compares large gene-conversion tracts across three regions of the 5S repeat unit:
the 5S gene, the antisense Alu SINE, and the remaining non-transcribed spacer.
Computed per haplotype, enabling error bars and a paired test (each haplotype is
measured in all 3 regions).

Regions (consensus coords):
  5S gene            630-749   (120 bp)
  Alu SINE           787-1066  (280 bp)
  other NTS spacer   everything else (1768 bp)

Large conversion tract = a variant whose carrier copies form a contiguous run of
>= K copies (allowing single-copy skips, gap<=2).

Metrics per haplotype per region:
  (A) large-tract density  = #large-tract variants / region length (per kb)
  (B) conditional contiguity = among variants in >=K copies, fraction forming a
      >=K contiguous tract (controls for differing numbers of variants across
      regions); pooled with bootstrap-over-haplotypes 95% CI.

Outputs:
  figures/50_geneconv_gene_alu_nts.pdf
  exports/geneconv_gene_alu_nts_per_haplotype.csv
  exports/geneconv_gene_alu_nts_summary.csv
"""

import os
import sqlite3
from collections import defaultdict
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

DB = Path(os.environ.get("FIVES_DB", "5S_rDNA.db"))
FIG = Path(os.environ.get("FIVES_OUT", "output")) / "figures"
EXP = Path(os.environ.get("FIVES_DATA", "data")) / "exports"; EXP.mkdir(exist_ok=True)

GENE = (630, 749); ALU = (787, 1066)
LEN = {"gene": 120, "alu": 280, "other": 2168 - 120 - 280}
REGS = ["gene", "alu", "other"]
RLAB = {"gene": "5S gene", "alu": "Alu SINE", "other": "other NTS spacer"}
RCOL = {"gene": "#ff9900", "alu": "#7b68ee", "other": "#9aa0a6"}
K = 8; GAP_ALLOW = 2
rng = np.random.default_rng(0)


def region_of(pos):
    if GENE[0] <= pos <= GENE[1]: return "gene"
    if ALU[0] <= pos <= ALU[1]:  return "alu"
    return "other"


def max_tract(carriers):
    cs = sorted(set(carriers)); best = cur = 1
    for i in range(1, len(cs)):
        cur = cur+1 if cs[i]-cs[i-1] <= GAP_ALLOW else 1
        best = max(best, cur)
    return best


def main():
    con = sqlite3.connect(DB)
    rows = con.execute("""
        SELECT h.haplotype_id, a.sample_id, h.hap_label, c.copy_number,
               v.consensus_pos, v.alt
        FROM variant v JOIN copy c USING(copy_id)
        JOIN haplotype h USING(haplotype_id) JOIN assembly a USING(assembly_id)
        WHERE c.border_note='interior' AND a.cohort IN ('HPRC_Year1','HPRC_Release2')
          AND v.alignment_source='consensus_t2t' AND v.var_type='snp'
          AND v.masked=0
    """).fetchall()
    con.close()

    # group carriers per (hap,pos,alt)
    carr = defaultdict(lambda: defaultdict(list))
    meta = {}
    for hid, samp, hl, cn, pos, alt in rows:
        carr[hid][(pos, alt)].append(cn); meta[hid] = (samp, hl)

    # per-haplotype tallies
    per = []                    # density rows
    cond = defaultdict(lambda: defaultdict(lambda: [0, 0]))  # hid->reg->[contig>=K, eligible>=K]
    cond_rec = []               # per >=K-carrier variant: (hid, region, is_contiguous) for regression
    for hid, vd in carr.items():
        nlarge = {r: 0 for r in REGS}
        for (pos, alt), cps in vd.items():
            reg = region_of(pos)
            k = len(set(cps))
            if k < 2: continue
            mt = max_tract(cps)
            if mt >= K:
                nlarge[reg] += 1
            if k >= K:
                cond[hid][reg][1] += 1
                cond_rec.append((hid, reg, int(mt >= K)))
                if mt >= K:
                    cond[hid][reg][0] += 1
        samp, hl = meta[hid]
        row = {"sample_id": samp, "hap_label": hl, "haplotype_id": hid}
        for r in REGS:
            row[f"n_largetract_{r}"] = nlarge[r]
            row[f"density_{r}_per_kb"] = nlarge[r] / (LEN[r]/1000)
        per.append(row)
    perdf = pd.DataFrame(per)
    perdf.to_csv(EXP/"geneconv_gene_alu_nts_per_haplotype.csv", index=False)

    # ── stats: density mean + bootstrap 95% CI + paired Wilcoxon ──────────────
    # (NOT SEM: with n=466 SEM is ~1/21 of the SD and conveys precision-of-mean
    #  only; the data are zero-inflated/skewed. Use a non-parametric bootstrap CI
    #  of the mean, consistent with panel B.)
    dens = {r: perdf[f"density_{r}_per_kb"].values for r in REGS}
    mean = {r: dens[r].mean() for r in REGS}
    sd = {r: dens[r].std(ddof=1) for r in REGS}
    nH = len(perdf)
    bmean = {r: [] for r in REGS}
    idx = np.arange(nH)
    for _ in range(2000):
        bi = rng.choice(idx, size=nH, replace=True)
        for r in REGS:
            bmean[r].append(dens[r][bi].mean())
    ci = {r: (np.percentile(bmean[r], 2.5), np.percentile(bmean[r], 97.5)) for r in REGS}
    fried = stats.friedmanchisquare(dens["gene"], dens["alu"], dens["other"])
    pairs = {}
    for a, b in [("gene","alu"),("gene","other"),("alu","other")]:
        pairs[(a,b)] = stats.wilcoxon(dens[a], dens[b])[1]

    # ── conditional contiguity (pooled) + bootstrap-over-haplotypes CI ────────
    hids = list(cond.keys())
    def pooled_cond(sample_hids):
        num = {r: 0 for r in REGS}; den = {r: 0 for r in REGS}
        for hid in sample_hids:
            for r in REGS:
                c, e = cond[hid][r]; num[r] += c; den[r] += e
        return {r: (num[r]/den[r] if den[r] else np.nan) for r in REGS}, den
    obs_cond, obs_den = pooled_cond(hids)
    B = 1000
    boot = {r: [] for r in REGS}
    arr = np.array(hids)
    for _ in range(B):
        s = rng.choice(arr, size=len(arr), replace=True)
        pc, _ = pooled_cond(s)
        for r in REGS: boot[r].append(pc[r])
    cond_ci = {r: (np.nanpercentile(boot[r],2.5), np.nanpercentile(boot[r],97.5)) for r in REGS}

    # ── Panel B test: logistic regression of contiguity ~ region, cluster-robust
    # by haplotype (variants within an array are not independent). ORs vs reference.
    import statsmodels.api as sm
    import statsmodels.formula.api as smf
    cr = pd.DataFrame(cond_rec, columns=["hid", "region", "contig"])
    panelB = {}
    def fit(ref):
        m = smf.glm(f"contig ~ C(region, Treatment('{ref}'))", data=cr,
                    family=sm.families.Binomial())
        return m.fit(cov_type="cluster", cov_kwds={"groups": cr["hid"].values})
    rf = fit("other")
    for term in rf.params.index:
        if "gene" in term: panelB[("gene","other")] = (np.exp(rf.params[term]), rf.pvalues[term])
        if "alu" in term:  panelB[("alu","other")]  = (np.exp(rf.params[term]), rf.pvalues[term])
    rfa = fit("alu")
    for term in rfa.params.index:
        if "gene" in term: panelB[("gene","alu")] = (np.exp(rfa.params[term]), rfa.pvalues[term])
    print("Panel B logistic (contiguity ~ region, cluster-robust by haplotype):")
    for (a,b),(orr,p) in panelB.items():
        print(f"  {a} vs {b}: OR={orr:.3f}, p={p:.2e}")

    # ── summary export ────────────────────────────────────────────────────────
    summ = []
    for r in REGS:
        summ.append({"region": RLAB[r], "length_bp": LEN[r],
                     "density_per_kb_mean": round(mean[r],4), "density_per_kb_sd": round(sd[r],4),
                     "density_ci_low": round(ci[r][0],4), "density_ci_high": round(ci[r][1],4),
                     "cond_contiguous_frac": round(obs_cond[r],4),
                     "cond_ci_low": round(cond_ci[r][0],4), "cond_ci_high": round(cond_ci[r][1],4),
                     "n_eligible_ge_K": int(obs_den[r])})
    pd.DataFrame(summ).to_csv(EXP/"geneconv_gene_alu_nts_summary.csv", index=False)

    print(f"K={K}, gap_allow={GAP_ALLOW}, {len(perdf)} haplotypes")
    print("Density (large tracts /kb/hap, mean [bootstrap 95% CI]):")
    for r in REGS: print(f"  {RLAB[r]:18s}: {mean[r]:.3f} [{ci[r][0]:.3f}-{ci[r][1]:.3f}]")
    print(f"Friedman p={fried.pvalue:.2e}; paired Wilcoxon: "
          + ", ".join(f"{a}vs{b} p={p:.1e}" for (a,b),p in pairs.items()))
    print("Conditional contiguity (>=K carriers -> >=K contiguous tract):")
    for r in REGS:
        print(f"  {RLAB[r]:18s}: {obs_cond[r]*100:.1f}% [{cond_ci[r][0]*100:.1f}-{cond_ci[r][1]*100:.1f}]  (n={obs_den[r]})")

    # ── figure ─────────────────────────────────────────────────────────────────
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    x = np.arange(3)
    # A: density
    ax = axes[0]
    yerrA = np.array([[mean[r]-ci[r][0] for r in REGS], [ci[r][1]-mean[r] for r in REGS]])
    ax.bar(x, [mean[r] for r in REGS], yerr=yerrA,
           color=[RCOL[r] for r in REGS], edgecolor="black", linewidth=0.6, capsize=5)
    ax.set_xticks(x); ax.set_xticklabels([RLAB[r] for r in REGS])
    ax.set_ylabel(f"Large conversion tracts (≥{K} copies) per kb per haplotype")
    ax.set_title("A   Large-tract density by region (mean, bootstrap 95% CI)",
                 fontweight="bold", loc="left")
    ymax = max(ci[r][1] for r in REGS)
    def star(p): return "***" if p<1e-3 else "**" if p<1e-2 else "*" if p<0.05 else "ns"
    for i,((a,b),yy) in enumerate(zip([("gene","alu"),("gene","other"),("alu","other")],
                                      [ymax*1.08, ymax*1.22, ymax*1.36])):
        ia, ib = REGS.index(a), REGS.index(b)
        ax.plot([ia,ia,ib,ib],[yy,yy*1.01,yy*1.01,yy],lw=1,c="black")
        ax.text((ia+ib)/2, yy*1.015, star(pairs[(a,b)]), ha="center", fontsize=11)
    ax.set_ylim(0, ymax*1.5)
    ax.text(0.98,0.6,f"Friedman p={fried.pvalue:.0e}",transform=ax.transAxes,ha="right",fontsize=9)

    # B: conditional contiguity
    ax2 = axes[1]
    yv = [obs_cond[r]*100 for r in REGS]
    err = np.array([[obs_cond[r]*100-cond_ci[r][0]*100 for r in REGS],
                    [cond_ci[r][1]*100-obs_cond[r]*100 for r in REGS]])
    ax2.bar(x, yv, yerr=err, color=[RCOL[r] for r in REGS], edgecolor="black",
            linewidth=0.6, capsize=5)
    for i,r in enumerate(REGS):
        ax2.text(i, 1, f"n={obs_den[r]}", ha="center", va="bottom", fontsize=8)
    # significance brackets (logistic regression, cluster-robust by haplotype)
    ymaxB = max(cond_ci[r][1]*100 for r in REGS)
    for (a,b),yy in zip([("gene","alu"),("gene","other"),("alu","other")],
                        [ymaxB*1.06, ymaxB*1.18, ymaxB*1.30]):
        ia, ib = REGS.index(a), REGS.index(b)
        p = panelB.get((a,b), panelB.get((b,a),(np.nan,np.nan)))[1]
        ax2.plot([ia,ia,ib,ib],[yy,yy*1.01,yy*1.01,yy],lw=1,c="black")
        ax2.text((ia+ib)/2, yy*1.012, star(p), ha="center", fontsize=11)
    ax2.set_ylim(0, ymaxB*1.45)
    ax2.set_xticks(x); ax2.set_xticklabels([RLAB[r] for r in REGS])
    ax2.set_ylabel(f"% of ≥{K}-carrier variants forming a ≥{K} contiguous tract")
    ax2.set_title("B   Conditional contiguity (controls for variant abundance)\n"
                  "bars = bootstrap 95% CI; stars = logistic test (cluster-robust)",
                  fontweight="bold", loc="left", fontsize=10)

    fig.suptitle("Large gene-conversion tracts: 5S gene vs Alu SINE vs other spacer",
                 fontsize=13, fontweight="bold")
    fig.tight_layout(rect=[0,0,1,0.96])
    fig.savefig(FIG/"50_geneconv_gene_alu_nts.pdf", bbox_inches="tight", dpi=160)
    plt.close(fig)
    print(f"\nSaved figures/50_geneconv_gene_alu_nts.pdf + 2 CSVs")


if __name__ == "__main__":
    main()
