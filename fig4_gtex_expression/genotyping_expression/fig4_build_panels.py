#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# fig4_build_panels.py — build the GTEx RNA-seq expression panels: rank-skew
# volcano, per-locus carrier-vs-noncarrier VAF, and RNA-vs-DNA VAF concordance.
# Author: Jonathan Bohlen
# (c) 2026 Jonathan Bohlen. Code accompanying Sengl et al. (2026),
# "Surveillance and selection of 5S ribosomal RNA genes in the human genome."
# Released under the MIT License; see LICENSE at the repository root.
# -----------------------------------------------------------------------------
"""Build the GTEx RNA-seq publication panels.
Panels: P1 rank-skew volcano | P2 per-position carrier-vs-noncarrier dot plot |
P5 RNA-vs-DNA VAF concordance. Writes vector PDF + PNG per panel and a per-panel data TSV."""
import os,numpy as np,pandas as pd,matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
from matplotlib.lines import Line2D
from scipy.stats import spearmanr,pearsonr,mannwhitneyu
def bh(p):
    p=np.asarray(p,float); n=len(p); o=np.argsort(p)
    q=p[o]*n/np.arange(1,n+1); q=np.minimum.accumulate(q[::-1])[::-1]
    out=np.empty(n); out[o]=np.clip(q,0,1); return out
def fmtq(q):
    if q>=0.01: return f"{q:.3f}"
    if q>=1e-4: return f"{q:.1e}"
    return f"{q:.0e}"
plt.rcParams.update({"pdf.fonttype":42,"ps.fonttype":42,"font.family":"Arial",
    "axes.linewidth":0.6,"xtick.major.width":0.6,"ytick.major.width":0.6,
    "xtick.major.size":2,"ytick.major.size":2})
CM=1/2.54
ROOT=str(Path(os.environ.get("FIVES_DATA","data")))
C=f"{ROOT}/figures/wgs_rna_concordance"
DAT=f"{ROOT}/figures/figure4_rnaseq/data"
OUT=str(Path(os.environ.get("FIVES_OUT","output"))/"Figure4")
os.makedirs(OUT,exist_ok=True)
GEN=["687G","701C","725C","726T","730G","733T","734T","743G"]
RED="#d62728"; GREY="#9b9b9b"; DGREY="#4d4d4d"; ORange="#f4a300"
def style(ax):
    ax.spines[["top","right"]].set_visible(False)
    ax.tick_params(labelsize=5.5,length=2,width=0.6)
def save(fig,name):
    fig.savefig(f"{OUT}/{name}.pdf",dpi=400,bbox_inches="tight")
    fig.savefig(f"{OUT}/{name}.png",dpi=300,bbox_inches="tight"); plt.close(fig)

RS=pd.read_csv(f"{C}/rank_skew_byvariant.tsv",sep="\t"); RS["variant"]=RS.pos.astype(str)+RS.alt
RS["hi"]=(RS.fdr<0.10)&(RS.dir.str.contains("HIGH"))

# ───────────────────────── P1: rank-skew volcano (segmented y-axis) ─────────────────────────
d=RS.copy(); d["nlfdr"]=-np.log10(d.fdr.clip(lower=1e-300))
BRK=10; ymax=min(np.ceil(d.nlfdr.max()/10)*10,80); d["y"]=d.nlfdr.clip(upper=ymax)
ns=d[~d.hi]; hi=d[d.hi]
fig=plt.figure(figsize=(4.6*CM,4.4*CM))
gs=fig.add_gridspec(2,1,height_ratios=[1,4],hspace=0.10)   # top 20% (10..max), bottom 80% (0..10)
axt=fig.add_subplot(gs[0]); axb=fig.add_subplot(gs[1],sharex=axt)
for ax in (axt,axb):
    ax.scatter(ns.auc,ns.y,s=4,c=GREY,linewidths=0,alpha=.55,rasterized=True,zorder=2)
    ax.scatter(hi.auc,hi.y,s=5,c=RED,linewidths=0,alpha=.85,rasterized=True,zorder=3)
    ax.axvline(.5,color="k",lw=0.5,ls=":")
    ax.spines[["top","right"]].set_visible(False); ax.tick_params(labelsize=5.5,length=2,width=0.6)
axb.axhline(-np.log10(0.10),color="grey",lw=0.5,ls=":")
axt.set_ylim(11.2,ymax+2); axb.set_ylim(-0.6,10.6)   # break in the empty 10.6–11.2 gap -> no chopped dots
axb.set_xlim(0.28,0.85)
axt.set_yticks(np.arange(20,ymax+1,20)); axb.set_yticks([0,2,4,6,8,10])
axt.spines["bottom"].set_visible(False); axb.spines["top"].set_visible(False)
axt.tick_params(bottom=False,labelbottom=False)
dkw=dict(marker=[(-1,-0.5),(1,0.5)],markersize=3,linestyle="none",color="k",mec="k",mew=0.6,clip_on=False)
axt.plot([0,1],[0,0],transform=axt.transAxes,**dkw); axb.plot([0,1],[1,1],transform=axb.transAxes,**dkw)
axb.set_xlabel("Carrier RNA-VAF skew (AUC)",fontsize=6)
fig.text(-0.02,0.5,"$-$log$_{10}$ FDR",rotation=90,va="center",ha="center",fontsize=6)
axt.text(.02,.92,f"{int(d.hi.sum())} carrier-HIGH (0 carrier-LOW)",transform=axt.transAxes,
    fontsize=4.6,va="top",ha="left",color=RED)
save(fig,"P1_rankskew_volcano")
d[["variant","pos","alt","n_carr","n_nonc","auc","p","fdr","dir"]].assign(
    carrierHIGH_FDR10=d.hi).sort_values("fdr").to_csv(f"{OUT}/P1_rankskew_volcano.tsv",sep="\t",index=False)

# ───────────── P2: per-locus carrier vs non-carrier truncated violins (ranked test) ─────────────
pd2=pd.read_csv(f"{DAT}/panel4_perdonor_carrierHIGH_vaf.tsv",sep="\t")  # all 64 carrier-HIGH loci
VARS=["664C","687G","692A","726T","733T","734T","743G","690A","675G","648T"]  # curated set
PFLOOR=0.003  # display floor (%) for undetected (VAF==0)
def L(p): return np.log10(np.where(p>0,p,PFLOOR))   # p in percent -> log10
def vio(ax,vals_pct,xpos,color,width=0.34,fill=0.45):
    """truncated violin (density clipped to p10-p90) + IQR line + median tick"""
    Y=L(vals_pct); lo,hi=np.percentile(Y,[10,90]); Yc=Y[(Y>=lo)&(Y<=hi)]
    if hi-lo>0.05 and len(np.unique(Yc))>=3:
        for pc in ax.violinplot([Yc],positions=[xpos],showmedians=False,
                showextrema=False,widths=width)['bodies']:
            pc.set_facecolor(color); pc.set_alpha(fill); pc.set_edgecolor("none")
    q1,med,q3=np.percentile(Y,[25,50,75])
    ax.vlines(xpos,q1,q3,color="black",lw=0.7,zorder=5)
    ax.hlines(med,xpos-0.10,xpos+0.10,color="black",lw=0.9,zorder=6)
# Mann-Whitney carrier vs non-carrier per locus; BH-correct across the curated set; sort by carrier median
recs=[]
for v in VARS:
    sub=pd2[pd2.variant==v]
    ca=sub[sub.group=="carrier"].rna_vaf.values; nc=sub[sub.group=="noncarrier"].rna_vaf.values
    p=mannwhitneyu(ca,nc,alternative="two-sided").pvalue
    recs.append(dict(variant=v,n_carr=len(ca),n_nonc=len(nc),carr_median=np.median(ca),p_mwu=p))
st=pd.DataFrame(recs); st["q_BH"]=bh(st.p_mwu.values)
st=st.sort_values("carr_median",ascending=False).reset_index(drop=True)
loci=st.variant.tolist()
rng=np.random.default_rng(0)
fig,ax=plt.subplots(figsize=(10.5*CM,4.6*CM))
for i,v in enumerate(loci):
    sub=pd2[pd2.variant==v]
    nc=sub[sub.group=="noncarrier"].rna_vaf.values*100
    ca=sub[sub.group=="carrier"].rna_vaf.values*100
    vio(ax,nc,i-0.21,GREY)
    vio(ax,ca,i+0.21,RED)
    ax.scatter(i+0.21+rng.uniform(-0.10,0.10,len(ca)),L(ca),s=2.4,c=RED,
        linewidths=0,alpha=.75,rasterized=True,zorder=4)
    ax.text(i,1.02,fmtq(st.q_BH.iloc[i]),fontsize=3.6,ha="center",va="bottom",rotation=90,color="black")
    ax.text(i,np.log10(PFLOOR)-0.22,f"n={len(ca)}",fontsize=4,ha="center",va="top",color=DGREY)
ax.axhline(np.log10(PFLOOR),color="grey",lw=0.5,ls=":")
ax.text(len(loci)-0.5,np.log10(PFLOOR)+0.06,"0 (undetected)",fontsize=4.3,color="grey",va="bottom",ha="right")
ax.set_ylim(np.log10(PFLOOR)-0.45,1.55)
ax.set_yticks([-3,-2,-1,0,1]); ax.set_yticklabels(["0.001","0.01","0.1","1","10"])
ax.set_xticks(range(len(loci))); ax.set_xticklabels(loci,fontsize=5,rotation=90)
ax.set_xlim(-0.6,len(loci)-0.4)
ax.set_ylabel("Pooled RNA-VAF (%)",fontsize=6)
import matplotlib.patches as mp
ax.legend(handles=[mp.Patch(color=RED,alpha=.6,label="carrier"),
    mp.Patch(color=GREY,alpha=.6,label="non-carrier")],fontsize=4.6,loc="upper right",
    bbox_to_anchor=(1.0,1.16),ncol=2,handletextpad=0.4,columnspacing=0.8,
    frameon=False,borderpad=0.1,handlelength=1.0)
ax.text(0.0,1.04,"q: Mann–Whitney U, BH-corrected (×10)",transform=ax.transAxes,
    fontsize=4,va="bottom",ha="left",color=DGREY)
style(ax); save(fig,"P2_perposition_carrier_vs_noncarrier")
pd2[pd2.variant.isin(loci)].to_csv(f"{OUT}/P2_perposition_carrier_vs_noncarrier.tsv",sep="\t",index=False)
st.to_csv(f"{OUT}/P2_ranktest_stats.tsv",sep="\t",index=False)
print("P2 curated loci (by carrier median):",loci)
print(st.to_string(index=False))

# P3 (sliding-window expression-detection profile) is produced by a separate script.

# ───────────── P5: RNA vs DNA VAF concordance (coloured by variant) ─────────────
g=pd.read_csv(f"{DAT}/panelG_rna_vs_dna_vaf.tsv",sep="\t")
x=g.wgs_vaf.values*100; y=g.rna_vaf.values*100
rs,prs=spearmanr(x,y); rp,prp=pearsonr(np.log10(x),np.log10(y))
fold=np.median(x/y)
from scipy.stats import linregress
def ratio_med(v): s=g[g.variant==v]; return float((s.rna_vaf/s.wgs_vaf).median())
# Pearson correlations on raw VAF, computed across all variants and within each variant.
r_across=pearsonr(x,y)[0]
wr=[pearsonr(s.wgs_vaf,s.rna_vaf)[0] for v,s in g.groupby("variant") if len(s)>=15]
within_med=float(np.median(wr))
# three example variants highlighted
HIv,LOv,LO2v="687G","725C","733T"
HI_C,LO_C="#ff7f00","#377eb8"; LO2_CD="#6a51a3"; LO2_C="#b39ddb"   # 687G orange, 725C blue, 733T purple
fig,ax=plt.subplots(figsize=(4.4*CM,4.4*CM))
lo=min(x.min(),y.min())*0.6; hi=max(x.max(),y.max())*1.5
ax.plot([lo,hi],[lo,hi],color="k",lw=0.5,ls="--",zorder=1)                       # 1:1
ax.plot([lo,hi],[lo/6,hi/6],color="grey",lw=0.5,ls=":",zorder=1)                 # 1:6 (cohort median)
ax.text(hi*0.30,hi*0.30,"1:1",fontsize=3.8,color="k",ha="left",va="bottom",rotation=45)
ax.text(hi*0.30,hi*0.30/6,"÷6",fontsize=3.8,color="grey",ha="left",va="bottom",rotation=45)
def wfit(v,c,z,lw):                                # within-variant regression line (log-log)
    s=g[g.variant==v]; lr=linregress(np.log10(s.wgs_vaf*100),np.log10(s.rna_vaf*100))
    xx=np.array([s.wgs_vaf.min(),s.wgs_vaf.max()])*100
    ax.plot(xx,10**(lr.intercept+lr.slope*np.log10(xx)),color=c,lw=lw,zorder=z)
bg=g[~g.variant.isin([HIv,LOv,LO2v])]
ax.scatter(bg.wgs_vaf*100,bg.rna_vaf*100,s=3,c="#dddddd",linewidths=0,alpha=.45,rasterized=True,zorder=1)
# example variant LO2v plotted behind the other coloured dots
s2=g[g.variant==LO2v]
ax.scatter(s2.wgs_vaf*100,s2.rna_vaf*100,s=3.5,c=LO2_C,linewidths=0,alpha=.65,rasterized=True,zorder=2)
wfit(LO2v,LO2_CD,3,0.7)
for v,c in [(HIv,HI_C),(LOv,LO_C)]:
    s=g[g.variant==v]
    ax.scatter(s.wgs_vaf*100,s.rna_vaf*100,s=4,c=c,linewidths=0,alpha=.7,rasterized=True,zorder=4)
    if len(s)>=15: wfit(v,c,5,0.8)
for v,c,z in [(LO2v,LO2_CD,5),(HIv,HI_C,6),(LOv,LO_C,6)]:
    s=g[g.variant==v]
    ax.scatter(s.wgs_vaf.median()*100,s.rna_vaf.median()*100,s=24,c=c,marker="D",
        edgecolor="k",linewidths=0.5,zorder=z)
mx,my=g[g.variant==HIv].wgs_vaf.median()*100,g[g.variant==HIv].rna_vaf.median()*100
ax.annotate(f"{HIv}\nRNA≈DNA",(mx,my),fontsize=4.2,color=HI_C,va="bottom",ha="right",
    xytext=(-3,3),textcoords="offset points",fontweight="bold")
lx,ly=g[g.variant==LOv].wgs_vaf.median()*100,g[g.variant==LOv].rna_vaf.median()*100
ax.annotate(f"{LOv}\nRNA ÷{round(1/ratio_med(LOv))}",(lx,ly),fontsize=4.2,color=LO_C,va="top",ha="left",
    xytext=(4,-2),textcoords="offset points",fontweight="bold")
zx,zy=g[g.variant==LO2v].wgs_vaf.median()*100,g[g.variant==LO2v].rna_vaf.median()*100
ax.annotate(f"{LO2v}\nRNA ÷{round(1/ratio_med(LO2v))}",(zx,zy),fontsize=4.2,color=LO2_CD,va="top",ha="left",
    xytext=(4,-2),textcoords="offset points")
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlim(lo,hi); ax.set_ylim(lo,hi); ax.set_aspect("equal")
ax.set_xlabel("DNA VAF (%)",fontsize=6); ax.set_ylabel("RNA VAF (%)",fontsize=6)
ax.text(.03,.97,f"Pearson r: within {within_med:.2f}, across {r_across:.2f}\nRNA ≈ {fold:.0f}× < DNA",
    transform=ax.transAxes,fontsize=4.2,va="top",ha="left")
style(ax); save(fig,"P5_rna_vs_dna_vaf")
g.assign(rna_dna_ratio=g.rna_vaf/g.wgs_vaf).to_csv(f"{OUT}/P5_rna_vs_dna_vaf.tsv",sep="\t",index=False)
print(f"P5: n={len(g)} Pearson across {r_across:.3f} | within-variant median Pearson {within_med:.3f} | fold {fold:.1f}")
print("panels written ->",OUT); print(sorted(os.listdir(OUT)))
