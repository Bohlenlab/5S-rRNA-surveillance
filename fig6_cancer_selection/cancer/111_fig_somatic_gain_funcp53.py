#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# 111_fig_somatic_gain_funcp53.py — bar panel of per-tumour somatic 5S gain prevalence by functional p53 status (competent vs inactive) in CPTAC and TCGA.
# Author: Jonathan Bohlen
# (c) 2026 Jonathan Bohlen. Code accompanying Sengl et al. (2026),
# "Surveillance and selection of 5S ribosomal RNA genes in the human genome."
# Released under the MIT License; see LICENSE at the repository root.
# -----------------------------------------------------------------------------
"""Somatic 5S GAIN per tumour, FUNCTIONAL-p53-competent vs -inactive, CPTAC + TCGA (gain only).
Functional p53-inactive = TP53 mut/del OR MDM2-amp OR MDM4-amp OR CDKN2A/ARF-del (+HPV for CPTAC).
CPTAC: surveillance_v2/tables/01_p53_status_expanded.tsv (p53_inactive). TCGA: built from
tp53/TCGA_tp53_status.tsv + tp53/tcga_cna_p53.tsv. Denominator = cases tested by the matched caller
(results_variants/somatic_calls_annotated.tsv). Bar panel with explicit p-values."""
import os, numpy as np, pandas as pd
from scipy import stats
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Patch
D=os.environ.get("FIVES_DATA","data"); OUT=os.environ.get("FIVES_OUT","output")
plt.rcParams.update({"pdf.fonttype":42,"ps.fonttype":42,"font.family":"Arial","font.size":7,
    "axes.linewidth":0.6,"xtick.major.width":0.6,"ytick.major.width":0.6,
    "xtick.major.size":2,"ytick.major.size":2})
CM=1/2.54; COMPc="#9b9b9b"; INACTc="#d62728"; THR=0.01

S=pd.read_csv(f"{D}/results_variants/somatic_calls_annotated.tsv",sep="\t")
S["cohortG"]=np.where(S.cohort.astype(str).str.startswith("TCGA"),"TCGA",np.where(S.cohort=="CPTAC","CPTAC","other"))
S=S[S.cohortG.isin(["CPTAC","TCGA"])].copy()

# ---- functional p53-inactive per case ----
cp=pd.read_csv(f"{D}/surveillance_v2/tables/01_p53_status_expanded.tsv",sep="\t")
cp_inact=dict(zip(cp.case, cp.p53_inactive.astype(int)))                      # CPTAC (incl. HPV)
tc=pd.read_csv(f"{D}/tp53/TCGA_tp53_status.tsv",sep="\t")[["case","tp53_status"]]
tcna=pd.read_csv(f"{D}/tp53/tcga_cna_p53.tsv",sep="\t")[["case","mdm2_amp","mdm4_amp","cdkn2a_del"]]
tf=tc.merge(tcna,on="case",how="outer").fillna(0)
tf["inact"]=(((tf.tp53_status=="deficient")|(tf.mdm2_amp>0)|(tf.mdm4_amp>0)|(tf.cdkn2a_del>0))).astype(int)
tc_inact=dict(zip(tf.case,tf.inact))
def func_p53(row):
    return cp_inact.get(row.case, np.nan) if row.cohortG=="CPTAC" else tc_inact.get(row.case, np.nan)
casetab=S.drop_duplicates("case")[["case","cohortG"]].copy()
casetab["inact"]=casetab.apply(func_p53,axis=1)

def wilson(k,n,z=1.96):
    # 95% CI on a proportion (tumour = Bernoulli replicate: has >=1 somatic 5S gain, or not)
    if n==0: return (np.nan,np.nan)
    p=k/n; den=1+z*z/n; ctr=(p+z*z/(2*n))/den; hw=z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/den
    return ctr-hw, ctr+hw
def agg(coh):
    ct=casetab[casetab.cohortG==coh].dropna(subset=["inact"])
    cov=len(ct); tot=(casetab.cohortG==coh).sum()
    sig=set(S[(S.cohortG==coh)&S.sig&(S.delta_vaf.abs()>=THR)&(S.direction=="gain")].case)
    ct=ct.assign(hasgain=ct.case.isin(sig).astype(int))
    out={}
    for v,lab in [(0,"comp"),(1,"inact")]:
        d=ct[ct.inact==v]; k=int(d.hasgain.sum()); N=len(d)
        out[lab]=dict(prev=k/N,k=k,N=N)
    # Fisher exact on 2x2 (group x has-gain) — tumours are the replicates
    tab=[[out["inact"]["k"],out["inact"]["N"]-out["inact"]["k"]],
         [out["comp"]["k"],out["comp"]["N"]-out["comp"]["k"]]]
    p=stats.fisher_exact(tab,alternative="two-sided")[1]
    rr=out["inact"]["prev"]/out["comp"]["prev"]
    return dict(cohort=coh,comp=out["comp"],inact=out["inact"],p=p,rr=rr,cov=cov,tot=tot)

res={c:agg(c) for c in ["CPTAC","TCGA"]}
for c in ["CPTAC","TCGA"]:
    r=res[c]; print(f"{c}: coverage {r['cov']}/{r['tot']} cases with functional-p53 | "
        f"competent {100*r['comp']['prev']:.1f}% ({r['comp']['k']}/{r['comp']['N']})  "
        f"inactive {100*r['inact']['prev']:.1f}% ({r['inact']['k']}/{r['inact']['N']})  "
        f"prev-ratio={r['rr']:.2f}  Fisher p={r['p']:.2e}")

# ---- figure: single compact panel ~4.5cm x 4cm, font 8, 1-pt lines ----
plt.rcParams.update({"font.size":8,"axes.linewidth":1.0,"xtick.major.width":1.0,
    "ytick.major.width":1.0,"xtick.major.size":2.5,"ytick.major.size":2.5})
def fmt_p(p): return (f"{p:.0e}".replace("e-0","e-") if p<1e-3 else f"{p:.2g}")
fig,ax=plt.subplots(figsize=(4.5*CM,4.0*CM),constrained_layout=True)
POS={"CPTAC":(0,1),"TCGA":(2.4,3.4)}; CTR={"CPTAC":0.5,"TCGA":2.9}; ymax=0
for coh in ["CPTAC","TCGA"]:
    r=res[coh]
    for x,col,key in [(POS[coh][0],COMPc,"comp"),(POS[coh][1],INACTc,"inact")]:
        val=100*r[key]["prev"]; lo,hi=[100*q for q in wilson(r[key]["k"],r[key]["N"])]
        ax.bar(x,val,0.85,color=col,edgecolor=col,linewidth=1.0,zorder=2)
        ax.errorbar(x,val,yerr=[[val-lo],[hi-val]],fmt="none",ecolor="0.2",elinewidth=1.0,capsize=2,capthick=1.0,zorder=3)
        ymax=max(ymax,hi)
for coh in ["CPTAC","TCGA"]:
    r=res[coh]; hi=100*wilson(r["inact"]["k"],r["inact"]["N"])[1]
    ax.text(CTR[coh],hi+0.5,f"p={fmt_p(r['p'])}",ha="center",va="bottom",fontsize=7)
ax.set_xticks([CTR["CPTAC"],CTR["TCGA"]]); ax.set_xticklabels(["CPTAC","TCGA"])
ax.set_ylabel("tumours with 5S gain (%)")
ax.set_ylim(0,ymax*1.30); ax.set_xlim(-0.7,4.1)
ax.spines[["top","right"]].set_visible(False)
ax.legend(handles=[Patch(fc=COMPc,ec=COMPc,label="p53 competent"),Patch(fc=INACTc,ec=INACTc,label="p53 inactive")],
    fontsize=6.5,frameon=False,loc="upper left",handlelength=0.9,handletextpad=0.4,labelspacing=0.2,borderpad=0.1)
for ext in ("pdf","png"):
    fig.savefig(f"{OUT}/surveillance_v2/figures/fig_somatic_5S_gain_funcp53.{ext}",dpi=300)
print("wrote fig_somatic_5S_gain_funcp53 (4.5x4cm, 300dpi)")
