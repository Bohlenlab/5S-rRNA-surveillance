#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# ex05g_dnarna_vaf_panel.py — plot the translation-module expression score against
# per-donor DNA VAF (%) and RNA VAF (%) of the gene-region 5S variants, on shared units.
# Author: Jonathan Bohlen
# (c) 2026 Jonathan Bohlen. Code accompanying Sengl et al. (2026),
# "Surveillance and selection of 5S ribosomal RNA genes in the human genome."
# Released under the MIT License; see LICENSE at the repository root.
# -----------------------------------------------------------------------------
"""Translation-module score vs variant VAF — DNA VAF (%) and RNA VAF (%) on the same units.
Produces (A) a single-panel overlay of both curves and (B) a two-panel version (DNA VAF | RNA VAF)."""
import os, sqlite3, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
from scipy import stats
plt.rcParams.update({"pdf.fonttype":42,"ps.fonttype":42,"font.family":"Arial",
    "axes.linewidth":0.6,"xtick.major.width":0.6,"ytick.major.width":0.6,"xtick.major.size":2,"ytick.major.size":2})
CM=1/2.54; BLUE="#2171b5"; RED="#d62728"
ROOT=str(Path(os.environ.get("FIVES_DATA","data"))); E=f"{ROOT}/results/eqtl/extreme"
OUT=str(Path(os.environ.get("FIVES_OUT","output"))/"Figure4")
DG=pd.read_pickle(f"{E}/audit_donor_resid.pkl")
res=pd.read_csv(f"{E}/fullmodel_highexpr_gsea.tsv",sep="\t")
TRG=[g for g in set().union(*[set(x.split(";")) for x in res[res.Term.str.contains("ranslat|ibosom|60S|40S|Peptide|rRNA",case=False)].head(12).Lead_genes]) if g in DG.columns]
score=DG[TRG].mean(1)
P=pd.read_csv(f"{E}/donor_variant_rnavaf.tsv",sep="\t"); GEN=["687G","701C","725C","726T","730G","733T","734T","743G"]
gg=P[P.variant.isin(GEN)]
rnav=gg.groupby("donor").rna_vaf.max()*100
dnav=gg.groupby("donor").wgs_vaf.max()*100                 # DNA VAF (%)
D=pd.DataFrame({"score":score}).join(rnav.rename("rna")).join(dnav.rename("dna")).dropna()
D["score"]=D.score*100   # ln-residual -> approx % expression deviation
print(f"carriers: {len(D)}  DNA VAF med {D.dna.median():.2f}% RNA VAF med {D.rna.median():.2f}%")

def binline(ax,x,y,nbin,color,label,lw=1.0,ms=3,elw=0.6,cap=1.5):
    qs=pd.qcut(x,nbin,labels=False,duplicates="drop")
    g=pd.DataFrame({"x":x,"y":y,"b":qs}).groupby("b")
    ax.errorbar(g.x.mean(),g.y.mean(),yerr=g.y.sem(),fmt="o-",color=color,ms=ms,lw=lw,
                capsize=cap,elinewidth=elw,mec=color,mfc=color,label=label)
    return stats.spearmanr(x,y)

# ---- (A) single-panel overlay (linear VAF %) ----
fig,ax=plt.subplots(figsize=(6.0*CM,4.6*CM))
rd,pd_=binline(ax,D.dna,D.score,5,BLUE,"DNA VAF")
rr,pr=binline(ax,D.rna,D.score,5,RED,"RNA VAF")
ax.axhline(0,color="grey",lw=0.4,ls=":")
ax.set_xlabel("Variant VAF (%)",fontsize=6); ax.set_ylabel("Translation\nprogram (%)",fontsize=6)
ax.tick_params(labelsize=5.5,length=2,width=0.6); ax.spines[["top","right"]].set_visible(False)
ax.text(.96,.96,f"DNA ρ={rd:+.2f}\nRNA ρ={rr:+.2f}",transform=ax.transAxes,fontsize=4.6,va="top",ha="right")
ax.legend(fontsize=5,frameon=False,loc="lower left")
fig.savefig(f"{OUT}/P13b_dnarna_vaf_overlay.pdf",dpi=400,bbox_inches="tight")
fig.savefig(f"{OUT}/P13b_dnarna_vaf_overlay.png",dpi=300,bbox_inches="tight"); plt.close(fig)

# ---- (B) two-panel homogenized (both VAF %) — 8x4 cm, dpi 300, font 8, 1 pt lines ----
plt.rcParams.update({"axes.linewidth":1.0,"xtick.major.width":1.0,"ytick.major.width":1.0})
fig,(a1,a2)=plt.subplots(1,2,figsize=(8*CM,4*CM),sharey=True)
binline(a1,D.dna,D.score,5,BLUE,"DNA",lw=1.0,ms=4,elw=1.0,cap=2)
binline(a2,D.rna,D.score,5,RED,"RNA",lw=1.0,ms=4,elw=1.0,cap=2)
for ax,(rho,p),col,lab in [(a1,stats.spearmanr(D.dna,D.score),BLUE,"Variant DNA-VAF (%)"),
                            (a2,stats.spearmanr(D.rna,D.score),RED,"Variant RNA-VAF (%)")]:
    ax.axhline(0,color="grey",lw=1.0,ls=":"); ax.set_xlabel(lab,fontsize=8)
    ax.tick_params(labelsize=8,length=3,width=1.0); ax.spines[["top","right"]].set_visible(False)
    ax.text(.96,.96,f"ρ={rho:+.2f}\np={p:.0e}",transform=ax.transAxes,fontsize=8,va="top",ha="right",color=col)
a1.set_ylabel("Translation\nprogram (%)",fontsize=8)
a1.set_title("DNA VAF",fontsize=8,color=BLUE); a2.set_title("RNA VAF",fontsize=8,color=RED)
fig.savefig(f"{OUT}/P13c_dnarna_vaf_2panel.pdf",dpi=300,bbox_inches="tight")
fig.savefig(f"{OUT}/P13c_dnarna_vaf_2panel.png",dpi=300,bbox_inches="tight"); plt.close(fig)
print("DNA VAF:",stats.spearmanr(D.dna,D.score),"\nRNA VAF:",stats.spearmanr(D.rna,D.score))
print("wrote P13b_overlay + P13c_2panel")
