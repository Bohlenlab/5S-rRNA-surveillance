# -----------------------------------------------------------------------------
# make_fig5_region.py — region-resolved box-whisker of 5S variant effect on
# expression and 60S incorporation, split by internal-promoter element.
# Author: Jonathan Bohlen
# (c) 2026 Jonathan Bohlen. Code accompanying Sengl et al. (2026),
# "Surveillance and selection of 5S ribosomal RNA genes in the human genome."
# Released under the MIT License; see LICENSE at the repository root.
# -----------------------------------------------------------------------------
"""Region-resolved box-whisker: variant effect on EXPRESSION and 60S INCORPORATION,
split by internal-promoter element (Box A 50-60, IE 67-72, Box C 80-90) vs other positions.
Box = Q1-Q3, whiskers = deciles (P10-P90), no outliers, log-y, neutral line at 1. Rank test (MWU) vs 'other'."""
import sqlite3, os, pandas as pd, numpy as np, matplotlib.pyplot as plt
from pathlib import Path
from scipy.stats import mannwhitneyu
plt.rcParams.update({"font.family":"Arial","pdf.fonttype":42,"ps.fonttype":42,"axes.linewidth":0.8,"font.size":7})
CM=1/2.54
OUT=str(Path(os.environ.get("FIVES_OUT","output"))/"Figure5"); os.makedirs(OUT,exist_ok=True)
DB=str(Path(os.environ.get("FIVES_DB","5S_rDNA.db")))
d=pd.read_sql("SELECT gene_pos,ref_base,alt_base,rna_expr_mean,incorp_60s_mean FROM functional_annotation",sqlite3.connect(DB))
def region(p): return "Box A" if 50<=p<=60 else "IE" if 67<=p<=72 else "Box C" if 80<=p<=90 else "other"
d["region"]=d.gene_pos.apply(region)
d.to_csv(f"{OUT}/Figure5_region_source_data.tsv",sep="\t",index=False)
order=["Box A","IE","Box C","other"]; ORANGE,GREY="#E8A33D","#b5b5b5"; cols=[ORANGE,ORANGE,ORANGE,GREY]
oth=d[d.region=="other"]
def stars(p): return "***" if p<1e-3 else "**" if p<1e-2 else "*" if p<0.05 else "ns"
# stats table (quartiles + deciles + fraction censored + rank p vs other)
rows=[]
for r in order:
    for col,nm in [("rna_expr_mean","expression"),("incorp_60s_mean","incorporation")]:
        s=d[d.region==r][col].dropna()
        p=mannwhitneyu(s,oth[col].dropna())[1] if r!="other" else np.nan
        rows.append(dict(region=r,readout=nm,n=len(s),p10=round(s.quantile(.1),4),q1=round(s.quantile(.25),4),
            median=round(s.median(),4),q3=round(s.quantile(.75),4),p90=round(s.quantile(.9),4),
            frac_below_0p05=round((s<0.05).mean(),3),p_vs_other=(None if r=="other" else float(f"{p:.2e}"))))
stat=pd.DataFrame(rows); stat.to_csv(f"{OUT}/Figure5_region_stats.tsv",sep="\t",index=False)
# figure
fig,axes=plt.subplots(1,2,figsize=(8*CM,5*CM))
for ax,col,ylab,ylim,clip in [(axes[0],"rna_expr_mean","Expression (expressed / plasmid)",(9e-3,12),0.01),
                          (axes[1],"incorp_60s_mean","Incorporation (60S / total)",(0.15,7),None)]:
    data=[d[d.region==r][col].dropna().values for r in order]
    if clip: data=[np.clip(v,clip,None) for v in data]     # saturate display floor (censored 0.001 -> 0.01)
    bp=ax.boxplot(data,whis=(10,90),showfliers=False,patch_artist=True,widths=0.62,
        medianprops=dict(color="k",lw=1),boxprops=dict(lw=0.8,edgecolor="k"),
        whiskerprops=dict(lw=0.8,color="k"),capprops=dict(lw=0.8,color="k"))
    for patch,c in zip(bp["boxes"],cols): patch.set_facecolor(c)
    ax.set_yscale("log"); ax.axhline(1,color="grey",ls="--",lw=0.7,zorder=0)
    ax.set_xticks(range(1,5)); ax.set_xticklabels(order,fontsize=6.5)
    ax.set_ylabel(ylab,fontsize=6.5); ax.set_ylim(*ylim); ax.tick_params(labelsize=6,length=2)
    for sp in ["top","right"]: ax.spines[sp].set_visible(False)
    for i,r in enumerate(order[:3]):
        s=d[d.region==r][col].dropna(); p=mannwhitneyu(s,oth[col].dropna())[1]
        ax.text(i+1,ylim[1]*0.42,stars(p),ha="center",va="bottom",fontsize=6.5,fontweight="bold")
fig.suptitle("Internal-promoter (Box A / IE / Box C) variants: expression vs 60S incorporation",fontsize=6.5,y=1.02)
fig.tight_layout(); fig.savefig(f"{OUT}/Figure5_region_expr_incorp.pdf",bbox_inches="tight")
print("wrote:",OUT); print(stat.to_string(index=False))
