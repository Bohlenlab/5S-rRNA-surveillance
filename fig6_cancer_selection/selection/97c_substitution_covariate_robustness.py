#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# 97c_substitution_covariate_robustness.py — effect size and robustness of incorporation-defective 5S variant carrier frequency, controlling for substitution supply (Ts/Tv + CpG).
# Author: Jonathan Bohlen
# (c) 2026 Jonathan Bohlen. Code accompanying Sengl et al. (2026),
# "Surveillance and selection of 5S ribosomal RNA genes in the human genome."
# Released under the MIT License; see LICENSE at the repository root.
# -----------------------------------------------------------------------------
"""Effect-size and robustness analysis of incorporation-defective vs other 5S variant carrier frequency,
before and after controlling for mutational supply (transition/transversion + CpG context).
Continuous model: log10(freq) ~ incorp + expr (+ Ts/Tv + CpG).
Class model (for % effect): log10(freq) ~ is_incorp_def (+ Ts/Tv + CpG), expressed variants only.
% less frequent = (1 - 10^beta) * 100  (geometric-mean frequency)."""
import os
import pandas as pd, numpy as np, statsmodels.api as sm, matplotlib.pyplot as plt, glob
plt.rcParams.update({"font.family":"Arial","pdf.fonttype":42,"ps.fonttype":42,"axes.linewidth":0.8,"font.size":7})
CM=1/2.54
MT=glob.glob(f"{os.environ.get('FIVES_OUT','output')}/**/97_master_table.tsv",recursive=True)[0]
OUT=f"{os.environ.get('FIVES_OUT','output')}/Figure5"
fa=pd.read_csv(MT,sep="\t")
COH=["UKBB","HPRC","GTEx"]; COL={"UKBB":"#4a7fb5","HPRC":"#4c9e4c","GTEx":"#d0a13a"}

TS={("A","G"),("G","A"),("C","T"),("T","C")}
fa["is_transition"]=[1 if (r,a) in TS else 0 for r,a in zip(fa.ref,fa.alt)]
refseq=fa.drop_duplicates("pos").set_index("pos").ref.to_dict()
def is_cpg(pos,ref):
    if ref=="C": return int(refseq.get(pos+1)=="G")
    if ref=="G": return int(refseq.get(pos-1)=="C")
    return 0
fa["is_CpG"]=[is_cpg(p,r) for p,r in zip(fa.pos,fa.ref)]
reg=fa[fa.incorp.notna()].copy()
FLOOR={c:max(reg[c][reg[c]>0].min()/2,1e-7) for c in COH}
def logf(x,c): return np.log10(np.where(x>0,x,FLOOR[c]))

# ---- continuous incorporation beta (kept, for the table) ----
brow=[]
for c in COH:
    y=logf(reg[c].values,c)
    m0=sm.OLS(y,sm.add_constant(reg[["incorp","expr"]].values)).fit()
    m1=sm.OLS(y,sm.add_constant(reg[["incorp","expr","is_transition","is_CpG"]].values)).fit()
    brow+=[dict(cohort=c,model="incorp+expr",b_inc=m0.params[1],se=m0.bse[1],p=m0.pvalues[1]),
           dict(cohort=c,model="+Ts/Tv+CpG",b_inc=m1.params[1],se=m1.bse[1],p=m1.pvalues[1])]
pd.DataFrame(brow).to_csv(f"{OUT}/Figure5_incorporation_regression_robustness.tsv",sep="\t",index=False)

# ---- class model: % less frequent than neutral (incorp-def vs competent), expressed only ----
ex=reg[reg["class"].isin(["competent","incorp_defective"])].copy()
ex["is_incdef"]=(ex["class"]=="incorp_defective").astype(int)
def pct(b): return (1-10**b)*100
prow=[]
for c in COH:
    y=logf(ex[c].values,c)
    mr=sm.OLS(y,sm.add_constant(ex[["is_incdef"]].values)).fit()
    ma=sm.OLS(y,sm.add_constant(ex[["is_incdef","is_transition","is_CpG"]].values)).fit()
    for lab,m in [("raw",mr),("adjusted",ma)]:
        b,se,p=m.params[1],m.bse[1],m.pvalues[1]
        prow.append(dict(cohort=c,model=lab,pct_less=pct(b),
                         lo=pct(b+1.96*se),hi=pct(b-1.96*se),p=p))
# pooled across the 3 cohorts (cohort fixed effects, cluster-robust SE by variant)
lg=ex.melt(id_vars=["is_incdef","is_transition","is_CpG","pos","alt"],value_vars=COH,var_name="coh",value_name="freq")
lg["lf"]=[np.log10(f) if f>0 else np.log10(FLOOR[c]) for f,c in zip(lg.freq,lg.coh)]
Xd=pd.get_dummies(lg.coh,drop_first=True).astype(float)
for lab,cols in [("raw",["is_incdef"]),("adjusted",["is_incdef","is_transition","is_CpG"])]:
    X=sm.add_constant(pd.concat([lg[cols].astype(float),Xd],axis=1))
    mp=sm.OLS(lg.lf.astype(float),X).fit(cov_type="cluster",cov_kwds={"groups":(lg.pos.astype(str)+"_"+lg.alt)})
    b,se,p=mp.params["is_incdef"],mp.bse["is_incdef"],mp.pvalues["is_incdef"]
    prow.append(dict(cohort="Pooled",model=lab,pct_less=pct(b),lo=pct(b+1.96*se),hi=pct(b-1.96*se),p=p))
pe=pd.DataFrame(prow); pe.to_csv(f"{OUT}/Figure5_incdef_pct_depletion.tsv",sep="\t",index=False)
print(pe.round(1).to_string(index=False))
adj=pe[(pe.model=="adjusted")&(pe.cohort!="Pooled")]; raw=pe[(pe.model=="raw")&(pe.cohort!="Pooled")]
pa=pe[(pe.model=="adjusted")&(pe.cohort=="Pooled")].iloc[0]; pr=pe[(pe.model=="raw")&(pe.cohort=="Pooled")].iloc[0]
print(f"\nPer-cohort ADJUSTED: {adj.pct_less.min():.0f}-{adj.pct_less.max():.0f}% (raw {raw.pct_less.min():.0f}-{raw.pct_less.max():.0f}%)")
print(f"POOLED: raw {pr.pct_less:.0f}% (p={pr.p:.1e}) -> adjusted {pa.pct_less:.0f}% less frequent (95% CI {pa.lo:.0f}-{pa.hi:.0f}%, p={pa.p:.1e})")

# ---- figure: % less frequent than neutral, raw (hollow) vs Ts/Tv+CpG-adjusted (filled) ----
def star(p): return "***" if p<1e-3 else "**" if p<1e-2 else "*" if p<0.05 else "ns"
COL["Pooled"]="#333333"; ROWS=COH+["Pooled"]
fig,ax=plt.subplots(figsize=(8*CM,5.4*CM))
for i,c in enumerate(ROWS):
    r=pe[(pe.cohort==c)&(pe.model=="raw")].iloc[0]; a=pe[(pe.cohort==c)&(pe.model=="adjusted")].iloc[0]
    mk="D" if c=="Pooled" else "o"; ms=6 if c=="Pooled" else 5
    ax.errorbar(r.pct_less,i+0.16,xerr=[[r.pct_less-r.lo],[r.hi-r.pct_less]],fmt=mk,ms=ms,mfc="white",mec=COL[c],ecolor=COL[c],lw=1,capsize=2)
    ax.errorbar(a.pct_less,i-0.16,xerr=[[a.pct_less-a.lo],[a.hi-a.pct_less]],fmt=mk,ms=ms,color=COL[c],ecolor=COL[c],lw=1,capsize=2)
    ax.text(a.hi+1.5,i-0.16,star(a.p),va="center",fontsize=6,color=COL[c])
ax.axvline(0,color="grey",ls="--",lw=0.8)
ax.set_yticks(range(len(ROWS))); ax.set_yticklabels(ROWS,fontsize=7)
ax.set_xlabel("% less frequent than neutral variants\n(incorporation-defective, geometric mean)",fontsize=6.5)
ax.set_title("Incorporation-defective 5S variant carrier frequency\n(after controlling for substitution rate)",fontsize=6.5)
ax.tick_params(labelsize=6)
from matplotlib.lines import Line2D
ax.legend(handles=[Line2D([],[],marker="o",ls="",mfc="white",mec="k",label="raw"),
                   Line2D([],[],marker="o",ls="",color="k",label="+ Ts/Tv + CpG")],
          frameon=False,fontsize=5.5,loc="lower right")
for s in ["top","right"]: ax.spines[s].set_visible(False)
fig.tight_layout(); fig.savefig(f"{OUT}/Figure5_incdef_pct_depletion.pdf",bbox_inches="tight")
print("\nwrote Figure5_incdef_pct_depletion.pdf + .tsv")
