#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# fig4_coverage_panel.py — plot pooled RNA-seq read depth along the 5S repeat
# unit as a median line with quartile (q25-75) and decile (q10-90) bands.
# Author: Jonathan Bohlen
# (c) 2026 Jonathan Bohlen. Code accompanying Sengl et al. (2026),
# "Surveillance and selection of 5S ribosomal RNA genes in the human genome."
# Released under the MIT License; see LICENSE at the repository root.
# -----------------------------------------------------------------------------
"""RNA-seq coverage panel: per-donor pooled read depth along the 5S repeat unit,
median + quartile (q25-75) and decile (q10-90) shaded bands across donors."""
import os,numpy as np,pandas as pd,matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
plt.rcParams.update({"pdf.fonttype":42,"ps.fonttype":42,"font.family":"Arial",
    "axes.linewidth":0.6,"xtick.major.width":0.6,"ytick.major.width":0.6,
    "xtick.major.size":2,"ytick.major.size":2})
CM=1/2.54
DAT=str(Path(os.environ.get("FIVES_DATA","data"))/"figures"/"figure4_rnaseq"/"data")
OUT=str(Path(os.environ.get("FIVES_OUT","output"))/"Figure4")
GENE=(630,748); BLUE="#2171b5"; B1="#6baed6"; B2="#c6dbef"; ORange="#f4a300"
M=pd.read_csv(f"{DAT}/coverage_perdonor_matrix.tsv.gz",sep="\t",index_col="donor")
pos=M.columns.astype(int).values; A=M.values.astype(float)
q10,q25,q50,q75,q90=[np.percentile(A,p,axis=0) for p in (10,25,50,75,90)]
fig,ax=plt.subplots(figsize=(4*CM,4*CM))
ax.axvspan(GENE[0],GENE[1],color=ORange,alpha=.12,zorder=0,lw=0)
ax.fill_between(pos,np.clip(q10,1,None),np.clip(q90,1,None),color=B2,alpha=.7,lw=0,zorder=2,label="10–90%")
ax.fill_between(pos,np.clip(q25,1,None),np.clip(q75,1,None),color=B1,alpha=.7,lw=0,zorder=3,label="25–75%")
ax.plot(pos,np.clip(q50,1,None),color=BLUE,lw=1.0,zorder=4,label="median")
ax.set_yscale("log"); ax.set_ylim(5,max(q90)*1.3)
ax.set_xlim(pos.min(),pos.max()); ax.set_xticks([600,650,700,750])
ax.set_xlabel("Position on 5S rDNA repeat unit",fontsize=6)
ax.set_ylabel("Pooled RNA-seq depth\nper donor (reads)",fontsize=6)
ax.tick_params(labelsize=5.5,length=2,width=0.6)
ax.spines[["top","right"]].set_visible(False)
ax.text(np.mean(GENE),max(q90)*1.0,"5S rRNA gene",fontsize=4.8,ha="center",va="top",color="#8a5a00")
ax.legend(fontsize=4.4,loc="lower right",frameon=False,handlelength=1.0,handletextpad=0.4,
    labelspacing=0.25,borderpad=0.1)
ax.text(.02,.80,f"n={A.shape[0]} donors\nmedian {int(q50[690-pos.min()]):,}× at gene",
    transform=ax.transAxes,fontsize=4.4,va="top",ha="left",color="#333333")
fig.savefig(f"{OUT}/P0_rna_coverage_5S_profile.pdf",dpi=400,bbox_inches="tight")
fig.savefig(f"{OUT}/P0_rna_coverage_5S_profile.png",dpi=300,bbox_inches="tight"); plt.close(fig)
pd.DataFrame({"pos":pos,"q10":q10,"q25":q25,"q50_median":q50,"q75":q75,"q90":q90,
    "mean":A.mean(0),"n_donors":A.shape[0]}).to_csv(
    f"{OUT}/P0_rna_coverage_5S_profile.tsv",sep="\t",index=False)
print("wrote P0 coverage panel; gene-center median depth:",int(q50[690-pos.min()]))
