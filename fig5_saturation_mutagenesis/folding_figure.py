#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# folding_figure.py — Render the 5S rRNA folding-energy figure: per-position ddG
# track plus ddG-vs-function and cross-algorithm validation panels.
# Author: Jonathan Bohlen
# (c) 2026 Jonathan Bohlen. Code accompanying Sengl et al. (2026),
# "Surveillance and selection of 5S ribosomal RNA genes in the human genome."
# Released under the MIT License; see LICENSE at the repository root.
# -----------------------------------------------------------------------------
"""5S rRNA folding-energy figure: ddG along the transcript + correlation panels.
(A) mean ddG per position with internal-promoter elements shaded;
(B) ddG vs 60S incorporation; (C) ddG vs expression; (D) algorithm validation RNAfold vs RNAstructure."""
import os, pandas as pd, numpy as np, matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from scipy.stats import spearmanr
plt.rcParams.update({"font.family":"Arial","pdf.fonttype":42,"ps.fonttype":42,"axes.linewidth":0.8,"font.size":7})
CM=1/2.54
OUT=os.environ.get("FIVES_OUT","output")
d=pd.read_csv(f"{OUT}/RNAfold_sense_per_variant_CORRECTED.csv")
m=pd.read_csv(f"{OUT}/folding_vs_function_CORRECTED.tsv",sep="\t")
ICR=[(50,60,"Box A"),(67,72,"IE"),(80,90,"Box C")]; ORANGE="#E8A33D"
def inicr(p): return any(s<=p<=e for s,e,_ in ICR)

fig=plt.figure(figsize=(17*CM,11*CM))
gs=GridSpec(2,3,height_ratios=[0.85,1.05],hspace=0.55,wspace=0.42,left=0.08,right=0.98,top=0.9,bottom=0.11)

# ---- A: per-substitution ddG along the sequence (colored by alt base) ----
axA=fig.add_subplot(gs[0,:])
ytop=d.ddG.max()
for s,e,lab in ICR:
    axA.axvspan(s-0.5,e+0.5,color=ORANGE,alpha=0.30,zorder=0,lw=0)
    axA.text((s+e)/2,ytop*1.06,lab,ha="center",va="bottom",fontsize=5.5,color="#9a6a12")
basecol={"A":"#e41a1c","C":"#000000","G":"#4daf4a","T":"#377eb8"}   # A red, C black, G green, T blue
basemk={"A":"o","C":"s","G":"^","T":"v"}
for b in "ACGT":
    sub=d[d.Alt==b]
    axA.scatter(sub.Pos,sub.ddG,s=9,c=basecol[b],marker=basemk[b],edgecolor="none",label=b,zorder=2)
axA.axhline(0,color="k",lw=0.7)
axA.set_xlim(0,120); axA.set_xlabel("position on 5S rRNA (sense strand)",fontsize=7)
axA.set_ylabel("ΔΔG (kcal/mol)",fontsize=6.5)
axA.set_title("Per-substitution ΔΔG along the 5S rRNA (positive = destabilizing)",fontsize=6.5)
axA.legend(title="substitution",fontsize=5.5,title_fontsize=5.5,ncol=4,frameon=False,loc="upper left",
           handletextpad=0.05,columnspacing=0.7,borderaxespad=0.2)
axA.tick_params(labelsize=6)
for sp in ["top","right"]: axA.spines[sp].set_visible(False)

def scatter(ax,x,y,xlab,ylab,title,logy=False,identity=False):
    s=m.dropna(subset=[x,y]) if x in m and y in m else d.dropna(subset=[x,y])
    src=s
    r=spearmanr(src[x],src[y])
    col=[ORANGE if inicr(p) else "#9aa0a6" for p in src["Pos" if "Pos" in src else src.index]]
    ax.scatter(src[x],src[y],s=6,c=col,edgecolor="none",alpha=0.8,zorder=2)
    # regression line (on the plotted scale)
    xx=np.linspace(src[x].min(),src[x].max(),50)
    b,a=np.polyfit(src[x],src[y],1); ax.plot(xx,a+b*xx,color="k",lw=0.8,zorder=3)
    if identity:
        lo,hi=min(src[x].min(),src[y].min()),max(src[x].max(),src[y].max()); ax.plot([lo,hi],[lo,hi],ls="--",color="#888",lw=0.7)
    ax.set_xlabel(xlab,fontsize=6.5); ax.set_ylabel(ylab,fontsize=6.5); ax.set_title(title,fontsize=6.5)
    if logy: ax.set_yscale("log")
    ax.tick_params(labelsize=6)
    for sp in ["top","right"]: ax.spines[sp].set_visible(False)
    ax.text(0.04,0.94,f"ρ={r[0]:+.2f}\np={r[1]:.0e}",transform=ax.transAxes,fontsize=6,va="top")

# ---- B: ddG vs incorporation ----
scatter(fig.add_subplot(gs[1,0]),"ddG","incorp_60s_mean","ΔΔG (kcal/mol)","60S incorporation\n(60S / total)","Folding vs incorporation",logy=True)
# ---- C: ddG vs expression ----
scatter(fig.add_subplot(gs[1,1]),"ddG","rna_expr_mean","ΔΔG (kcal/mol)","expression\n(expressed / plasmid)","Folding vs expression",logy=True)
# ---- D: validation ----
axD=fig.add_subplot(gs[1,2])
rv=spearmanr(d.ddG,d.ddG_RNAstructure)
axD.scatter(d.ddG,d.ddG_RNAstructure,s=6,c="#5b8c5a",edgecolor="none",alpha=0.7)
lo,hi=d.ddG.min(),d.ddG.max(); axD.plot([lo,hi],[lo,hi],ls="--",color="#888",lw=0.7)
axD.set_xlabel("RNAfold ΔΔG (kcal/mol)",fontsize=6.5); axD.set_ylabel("RNAstructure ΔΔG\n(kcal/mol)",fontsize=6.5)
axD.set_title("Independent-algorithm validation",fontsize=6.5); axD.tick_params(labelsize=6)
for sp in ["top","right"]: axD.spines[sp].set_visible(False)
axD.text(0.04,0.94,f"ρ={rv[0]:.2f}",transform=axD.transAxes,fontsize=6,va="top")

# legend for region colour
from matplotlib.lines import Line2D
fig.legend(handles=[Line2D([],[],marker='o',ls='',color=ORANGE,label='internal promoter (Box A/IE/C)',ms=4),
                    Line2D([],[],marker='o',ls='',color='#9aa0a6',label='other position',ms=4)],
           frameon=False,fontsize=5.5,loc="lower center",ncol=2,bbox_to_anchor=(0.5,0.0))
for fp in [f"{OUT}/Figure_folding_energy_CORRECTED.pdf",
           os.path.join(OUT, "manuscript_tables", "Figure5", "Figure_folding_energy.pdf")]:
    fig.savefig(fp,bbox_inches="tight")
print("wrote Figure_folding_energy_CORRECTED.pdf (+ manuscript_tables/Figure5 copy)")
print(f"panel B ρ={spearmanr(m.dropna(subset=['ddG','incorp_60s_mean']).ddG,m.dropna(subset=['ddG','incorp_60s_mean']).incorp_60s_mean)[0]:+.2f}")
