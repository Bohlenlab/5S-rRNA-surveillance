#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# ex05j_p10_4x4_biogenesis.py — differential-expression volcano coloured by gene
# class: cytosolic ribosomal proteins, heat-shock/chaperones, and NPM1.
# Author: Jonathan Bohlen
# (c) 2026 Jonathan Bohlen. Code accompanying Sengl et al. (2026),
# "Surveillance and selection of 5S ribosomal RNA genes in the human genome."
# Released under the MIT License; see LICENSE at the repository root.
# -----------------------------------------------------------------------------
"""Differential-expression volcano: cytosolic ribosomal proteins in red, heat-shock/chaperone
genes in teal, and the ribosome-biogenesis factor NPM1 in orange. Writes P10_volcano_translation_4x4.pdf."""
import os,re,numpy as np,pandas as pd,matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
plt.rcParams.update({"pdf.fonttype":42,"ps.fonttype":42,"font.family":"Arial","axes.linewidth":1.0,
    "xtick.major.width":1.0,"ytick.major.width":1.0,"xtick.major.size":2.5,"ytick.major.size":2.5})
CM=1/2.54; RED="#d62728"; GREY="#cfcfcf"; DGREY="#4d4d4d"; ORANGE="#e6850e"; BLUE="#2171b5"
E=str(Path(os.environ.get("FIVES_DATA","data"))/"results"/"eqtl"/"extreme")
OUT=str(Path(os.environ.get("FIVES_OUT","output"))/"Figure4")
HSR={"HSPA1A","HSPA1B","HSPA6","HSPA8","HSPA4","HSPH1","HSP90AA1","HSP90AB1","DNAJB1","DNAJA1","DNAJA4",
     "DNAJB6","BAG3","HSPB1","HSF1","AHSA1","FKBP4","STIP1","CHORDC1","SERPINH1","HSPD1","HSPE1"}
FACT={"NPM1":"nucleolar (5S-RNP)"}          # only the ribosome-biogenesis factor
G=pd.read_csv(f"{E}/fullmodel_highexpr_DE.tsv",sep="\t")
cyto=lambda g: bool(re.match(r'^RP[LS]\d+[A-Z]?$',str(g))) or g in {"RPLP0","RPLP1","RPLP2","RPSA","FAU","UBA52"}
G["rp"]=G.gene.apply(cyto); G["hsr"]=G.gene.isin(HSR)
G["eff"]=G.eff*100   # ln-residual -> approx % expression change
# y-axis = -log10(BH-adjusted p / FDR), so the significance line sits at a round value (FDR 0.05 -> 1.3)
flo=G[G.fdr>0].fdr.min(); G["y"]=-np.log10(G.fdr.clip(lower=flo))
nsig=int((G.fdr<0.05).sum())
fig,ax=plt.subplots(figsize=(4*CM,4*CM))      # 4x4 cm
ns=G[~G.rp & ~G.hsr & ~G.gene.isin(FACT)]; rp=G[G.rp]; hs=G[G.hsr]; fa=G[G.gene.isin(FACT)]
ax.scatter(ns.eff,ns.y,s=2,c=GREY,linewidths=0,alpha=.5,rasterized=True,zorder=2)
ax.scatter(rp.eff,rp.y,s=7,c=RED,linewidths=0,alpha=.8,zorder=4)             # vector (not rasterized)
ax.scatter(hs.eff,hs.y,s=7,c=BLUE,linewidths=0,alpha=.8,zorder=5)            # vector (not rasterized)
ax.scatter(fa.eff,fa.y,s=14,c=ORANGE,edgecolor="k",linewidths=0.6,zorder=6)
ax.axvline(0,color="k",lw=1.0,ls=":"); ax.set_xlim(-7,7)
ax.axhline(-np.log10(0.05),color="grey",lw=1.0,ls="--")    # FDR 0.05 = 1.3
# NPM1 callout (ribosome-biogenesis factor). RP=red, HSR=teal, NPM1=orange.
for _,r in fa.iterrows():
    ax.annotate(r.gene,(r.eff,r.y),xytext=(4,2.4),fontsize=8,color="k",va="center",ha="left",
        fontweight="bold",arrowprops=dict(arrowstyle="-",lw=0.8,color="#888888",shrinkA=0,shrinkB=3))
ax.set_xlabel("Δ expression (%)",fontsize=8); ax.set_ylabel("$-$log$_{10}$ FDR",fontsize=8)
ax.tick_params(labelsize=8,length=2.5,width=1.0); ax.spines[["top","right"]].set_visible(False)
fig.savefig(f"{OUT}/P10_volcano_translation_4x4.pdf",dpi=300,bbox_inches="tight"); plt.close(fig)
print("wrote P10_volcano_translation_4x4.pdf (orange = NPM1 only)")
