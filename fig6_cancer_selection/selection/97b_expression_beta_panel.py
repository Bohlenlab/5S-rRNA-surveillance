# -----------------------------------------------------------------------------
# 97b_expression_beta_panel.py — forest plots of the expression coefficient from the per-cohort log-frequency model, with substitution covariates.
# Author: Jonathan Bohlen
# (c) 2026 Jonathan Bohlen. Code accompanying Sengl et al. (2026),
# "Surveillance and selection of 5S ribosomal RNA genes in the human genome."
# Released under the MIT License; see LICENSE at the repository root.
# -----------------------------------------------------------------------------
"""Forest plot of the EXPRESSION coefficient from the model
log10(freq) ~ incorp + expr + is_transition + is_CpG, fit per cohort (UKBB, HPRC, GTEx) and formatted
identically to the incorporation-coefficient panel. Also emits a side-by-side incorporation vs expression
forest and a compact stacked version."""
import os
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import statsmodels.formula.api as smf
from pathlib import Path
plt.rcParams.update({"pdf.fonttype": 42, "font.size": 9})
OUT=Path(os.environ.get("FIVES_OUT","output"))/"06_population_genetics"
M=pd.read_csv(OUT/"97_master_table.tsv",sep="\t")
# substitution covariates from the consensus reference (1-based coords)
seq="".join(l.strip() for l in open(Path(os.environ.get("FIVES_REFS","refs"))/"5S_t2t_consensus.fa") if not l.startswith(">")).upper()
pur=set("AG")
M["is_transition"]=[(r in pur)==(a in pur) for r,a in zip(M.ref,M.alt)]
def is_cpg(p):
    b=seq[p-1]
    return (b=="C" and p<len(seq) and seq[p]=="G") or (b=="G" and p>1 and seq[p-2]=="C")
M["is_cpg"]=[is_cpg(p) for p in M.pos]
COH=["UKBB","HPRC","GTEx"]; COH_COL={"UKBB":"#33608c","HPRC":"#2ca25f","GTEx":"#b8860b"}
FLOOR={c:max(M[c][M[c]>0].min()/2,1e-7) for c in COH}
def stars(p): return "***" if p<.001 else "**" if p<.01 else "*" if p<.05 else "n.s."
reg=M[M.incorp.notna()].copy()
stats={}
for c in COH:
    reg["y"]=np.log10(np.where(reg[c]>0,reg[c],FLOOR[c]))
    m=smf.ols("y ~ incorp + expr + is_transition + is_cpg",data=reg).fit()
    stats[c]=dict(b_inc=m.params["incorp"],se_inc=m.bse["incorp"],p_inc=m.pvalues["incorp"],
                  b_exp=m.params["expr"], se_exp=m.bse["expr"], p_exp=m.pvalues["expr"])
    print(f"{c}: incorp b={m.params['incorp']:+.3f} p={m.pvalues['incorp']:.1e} | expr b={m.params['expr']:+.3f} p={m.pvalues['expr']:.2f}")
def forest(ax,kb,ks,kp,xlabel,title):
    ys=np.arange(len(COH))[::-1]
    for y,c in zip(ys,COH):
        s=stats[c]
        ax.errorbar(s[kb],y,xerr=1.96*s[ks],fmt="o",color=COH_COL[c],capsize=3,ms=7,lw=1.8)
        ax.text(s[kb],y+.18,stars(s[kp]),ha="center",fontsize=9,color=COH_COL[c])
    ax.axvline(0,color="grey",ls="--",lw=1); ax.set_yticks(ys); ax.set_yticklabels(COH)
    ax.set_xlabel(xlabel,fontsize=8); ax.set_title(title,fontsize=9.5)
    ax.grid(axis="x",lw=.3,alpha=.4); ax.set_ylim(-.6,len(COH)-.4)
xmax=max(abs(stats[c]["b_inc"])+2*stats[c]["se_inc"] for c in COH)*1.05
# standalone expression panel
fig,ax=plt.subplots(figsize=(3.2,3.3))
forest(ax,"b_exp","se_exp","p_exp","expression effect on log₁₀ frequency\n(OLS, controlling for incorporation)","Expression coefficient")
ax.set_xlim(-xmax,xmax)
fig.savefig(OUT/"97b_expression_beta_panel.pdf",bbox_inches="tight"); fig.savefig(OUT/"97b_expression_beta_panel.png",dpi=200,bbox_inches="tight"); plt.close(fig)
# side-by-side incorp | expression (shared x) for the direct specificity contrast
fig2,(a1,a2)=plt.subplots(1,2,figsize=(6.6,3.3),sharex=True,sharey=True)
forest(a1,"b_inc","se_inc","p_inc","incorporation effect on log₁₀ frequency\n(OLS, controlling for expression)","Incorporation coefficient")
forest(a2,"b_exp","se_exp","p_exp","expression effect on log₁₀ frequency\n(OLS, controlling for incorporation)","Expression coefficient")
a1.set_xlim(-xmax,xmax)
fig2.tight_layout(); fig2.savefig(OUT/"97b_specificity_incorp_vs_expr.pdf",bbox_inches="tight"); fig2.savefig(OUT/"97b_specificity_incorp_vs_expr.png",dpi=200,bbox_inches="tight"); plt.close(fig2)

# ── compact STACKED version, 4 cm x 4 cm (incorporation top, expression bottom, shared x) ──
CM=1/2.54
def forest_c(ax,kb,ks,kp,label):
    ys=np.arange(len(COH))[::-1]
    for y,c in zip(ys,COH):
        s=stats[c]
        ax.errorbar(s[kb],y,xerr=1.96*s[ks],fmt="o",color=COH_COL[c],capsize=1.6,ms=4,lw=1.1,elinewidth=1.1,mec="none")
        ax.text(s[kb],y+.30,stars(s[kp]),ha="center",va="bottom",fontsize=5.5,color=COH_COL[c])
    ax.axvline(0,color="grey",ls="--",lw=.7)
    ax.set_yticks(ys); ax.set_yticklabels(COH,fontsize=5.5); ax.set_ylim(-.7,len(COH)-.25)
    ax.text(0.02,0.90,label,transform=ax.transAxes,fontsize=6,fontweight="bold",va="top")
    ax.tick_params(length=2,pad=1.5)
    for sp in ["top","right"]: ax.spines[sp].set_visible(False)
figc,(c1,c2)=plt.subplots(2,1,figsize=(4*CM,4*CM),sharex=True,gridspec_kw=dict(hspace=0.22))
forest_c(c1,"b_inc","se_inc","p_inc","incorporation")
forest_c(c2,"b_exp","se_exp","p_exp","expression")
c2.set_xlim(-xmax,xmax); c2.set_xticks([-0.2,0,0.2]); c2.tick_params(axis="x",labelsize=5.5)
c2.set_xlabel("effect on log$_{10}$ frequency",fontsize=6,labelpad=2)
figc.savefig(OUT/"97b_specificity_stacked_4cm.pdf",bbox_inches="tight")
figc.savefig(OUT/"97b_specificity_stacked_4cm.png",dpi=600,bbox_inches="tight"); plt.close(figc)
print("wrote 97b_expression_beta_panel + 97b_specificity_incorp_vs_expr + 97b_specificity_stacked_4cm (pdf/png)")
