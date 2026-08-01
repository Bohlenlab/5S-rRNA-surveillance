#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# 110_fig_somatic_gain_p53_panel.py — bar panel of per-tumour somatic 5S gain (and loss control) prevalence by TP53 status in CPTAC and TCGA.
# Author: Jonathan Bohlen
# (c) 2026 Jonathan Bohlen. Code accompanying Sengl et al. (2026),
# "Surveillance and selection of 5S ribosomal RNA genes in the human genome."
# Released under the MIT License; see LICENSE at the repository root.
# -----------------------------------------------------------------------------
"""Panel: per-tumour prevalence of somatic 5S GAIN in p53-mutant vs WT tumours (CPTAC + TCGA side by side),
with somatic LOSS as a p53-independent control.
Data: results_variants/pooled_flagexcl.tsv (per-tumour prevalence of >=1 significant somatic 5S event,
|dVAF|>=1%; genetic TP53-mut vs WT). filt='all' = all cancer types; 'excl' = flagged types dropped."""
import os
import pandas as pd, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Patch
D=os.environ.get("FIVES_DATA","data"); OUT=os.environ.get("FIVES_OUT","output")
plt.rcParams.update({"pdf.fonttype":42,"ps.fonttype":42,"font.family":"Arial","font.size":7,
    "axes.linewidth":0.6,"xtick.major.width":0.6,"ytick.major.width":0.6,
    "xtick.major.size":2,"ytick.major.size":2})
CM=1/2.54; WTc="#9b9b9b"; MUTc="#d62728"; FILT="all"
d=pd.read_csv(f"{D}/results_variants/pooled_flagexcl.tsv",sep="\t")
d=d[d.filt==FILT]
def wilson(k,n,z=1.96):
    if n==0: return (np.nan,np.nan)
    p=k/n; den=1+z*z/n; ctr=(p+z*z/(2*n))/den
    hw=z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/den
    return 100*(ctr-hw),100*(ctr+hw)
def stars(p): return "***" if p<1e-3 else "**" if p<1e-2 else "*" if p<0.05 else "†" if p<0.10 else "n.s."
COH=["CPTAC","TCGA"]; DIRS=["gain","loss"]
fig,axes=plt.subplots(1,2,figsize=(9.6*CM,5.4*CM),sharey=True,constrained_layout=True)
for ax,coh in zip(axes,COH):
    sub=d[d.cohort==coh]
    xg=[0,1]; centers={"gain":0,"loss":1.05}; w=0.38
    ymax=0
    for di,dr in enumerate(DIRS):
        row=sub[sub.direction==dr].iloc[0]
        cx=centers[dr]
        for j,(grp,col,ev,n) in enumerate([("WT",WTc,row.Ewt,row.nWT),("mut",MUTc,row.Emut,row.nMut)]):
            xpos=cx+(j-0.5)*w
            val=100*ev/n; lo,hi=wilson(ev,n)
            fc = col if dr=="gain" else "none"
            ec = col
            ax.bar(xpos,val,w,color=fc,edgecolor=ec,linewidth=0.8,zorder=2,
                   hatch=("" if dr=="gain" else ""))
            ax.errorbar(xpos,val,yerr=[[val-lo],[hi-val]],fmt="none",ecolor="0.25",
                        elinewidth=0.7,capsize=1.5,capthick=0.7,zorder=3)
            ymax=max(ymax,hi)
        rr=2**row.log2r
        ytop=max(wilson(row.Emut,row.nMut)[1],wilson(row.Ewt,row.nWT)[1])
        if dr=="gain":
            ax.text(cx,ytop+2.2,f"{stars(row.p)}",ha="center",va="bottom",fontsize=7,color=MUTc)
            ax.text(cx,ytop+5.6,f"{rr:.2f}×",ha="center",va="bottom",fontsize=5.6,color=MUTc)
        else:
            ax.text(cx,ytop+2.2,"n.s.",ha="center",va="bottom",fontsize=5.6,color="0.5")
    ax.set_xticks([centers["gain"],centers["loss"]]); ax.set_xticklabels(["Gain","Loss"],fontsize=7)
    ax.set_title(coh,fontsize=8,fontweight="bold")
    ax.set_xlim(-0.55,1.6); ax.spines[["top","right"]].set_visible(False)
    ax.tick_params(labelsize=6.5)
axes[0].set_ylabel("tumours with somatic 5S\nevent (%)",fontsize=7)
axes[0].set_ylim(0,48)
axes[0].legend(handles=[Patch(fc=WTc,ec=WTc,label="p53-WT"),
                        Patch(fc=MUTc,ec=MUTc,label="p53-mutant"),
                        Patch(fc="white",ec="0.4",label="open bars = loss (control)")],
               fontsize=5.6,frameon=False,loc="upper left",handlelength=1.1,handletextpad=0.5,
               labelspacing=0.35,borderpad=0.2)
fig.suptitle("Somatic 5S gain and loss prevalence by p53 status (CPTAC, TCGA)",fontsize=7.5)
for ext in ("pdf","png"):
    fig.savefig(f"{OUT}/surveillance_v2/figures/fig_somatic_5S_gain_p53_{FILT}.{ext}",
                dpi=300 if ext=="pdf" else 200,bbox_inches="tight")
print("wrote fig_somatic_5S_gain_p53_"+FILT)
for _,r in d.iterrows():
    print(f"{r.cohort:8s} {r.direction:4s}  mut {100*r.Emut/r.nMut:5.1f}% ({r.Emut}/{r.nMut})  "
          f"WT {100*r.Ewt/r.nWT:5.1f}% ({r.Ewt}/{r.nWT})  RR={2**r.log2r:.2f} p={r.p:.2g}")
