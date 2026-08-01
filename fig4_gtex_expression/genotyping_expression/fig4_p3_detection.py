#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# fig4_p3_detection.py — sliding-window profile of the fraction of tested 5S
# variants scored carrier-HIGH (expressed) along the gene, with a per-variant rug.
# Author: Jonathan Bohlen
# (c) 2026 Jonathan Bohlen. Code accompanying Sengl et al. (2026),
# "Surveillance and selection of 5S ribosomal RNA genes in the human genome."
# Released under the MIT License; see LICENSE at the repository root.
# -----------------------------------------------------------------------------
"""Sliding-window expression-detection profile along the 5S gene: fraction of tested
variants that are carrier-HIGH (expressed) in a +/-5 bp window. Rug = individual variants (det/not)."""
import os,numpy as np,pandas as pd,matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
plt.rcParams.update({"pdf.fonttype":42,"ps.fonttype":42,"font.family":"Arial",
    "axes.linewidth":0.6,"xtick.major.width":0.6,"ytick.major.width":0.6,
    "xtick.major.size":2,"ytick.major.size":2})
CM=1/2.54; RED="#d62728"; GREY="#9b9b9b"; ORange="#f4a300"; DGREY="#4d4d4d"
C=str(Path(os.environ.get("FIVES_DATA","data"))/"figures"/"wgs_rna_concordance")
OUT=str(Path(os.environ.get("FIVES_OUT","output"))/"Figure4")
RS=pd.read_csv(f"{C}/rank_skew_byvariant.tsv",sep="\t"); RS["variant"]=RS.pos.astype(str)+RS.alt
RS["hi"]=(RS.fdr<0.10)&(RS.dir.str.contains("HIGH"))
g=RS[(RS.pos>=630)&(RS.pos<=748)].copy()
GENE=(630,748); W=5; xs=np.arange(GENE[0],GENE[1]+1)
det=[];aucm=[];nt=[]
for x in xs:
    lo=min(max(x-W,GENE[0]),GENE[1]-2*W); hi=lo+2*W   # constant-width window, pinned inside the gene
    s=g[(g.pos>=lo)&(g.pos<=hi)]
    det.append(s.hi.mean() if len(s) else np.nan); aucm.append(s.auc.mean() if len(s) else np.nan); nt.append(len(s))
detp=np.array(det); aucp=np.array(aucm)

fig,ax=plt.subplots(figsize=(4*CM,4*CM))
ax.axvspan(*GENE,color=ORange,alpha=.10,zorder=0,lw=0)
ax.fill_between(xs,0,detp*100,color=RED,alpha=.30,lw=0,zorder=2)
ax.plot(xs,detp*100,color=RED,lw=1.0,zorder=3)
# rug of individual tested variants (red=expressed, grey=not)
for _,r in g.iterrows():
    ax.plot(r.pos,-7,marker="|",ms=3,mew=0.6,color=RED if r.hi else GREY,zorder=4)
ax.axhline(0,color="k",lw=0.5)
ax.set_ylim(-12,100); ax.set_xlim(626,752); ax.set_xticks([630,690,750])
ax.set_xlabel("Position on 5S rDNA repeat unit",fontsize=6)
ax.set_ylabel("Variants expressed in\n±5 bp window (%)",fontsize=6)
ax.tick_params(labelsize=5.5,length=2,width=0.6); ax.spines[["top","right"]].set_visible(False)
ax.text(689,99,"5S rRNA gene",fontsize=4.6,ha="center",va="top",color="#8a5a00")
fig.savefig(f"{OUT}/P3_detection_slidingwindow.pdf",dpi=400,bbox_inches="tight")
fig.savefig(f"{OUT}/P3_detection_slidingwindow.png",dpi=300,bbox_inches="tight"); plt.close(fig)
pd.DataFrame({"pos":xs,"n_tested_win":nt,"frac_expressed_win":det,"mean_auc_win":aucm}).to_csv(
    f"{OUT}/P3_detection_slidingwindow.tsv",sep="\t",index=False)
print("wrote P3; peak window det-frac:",np.nanmax(det).round(2),"at pos",xs[np.nanargmax(det)])
