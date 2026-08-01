#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# 24g_violin_raw_vs_neighbor.py — Two-panel violins of raw versus neighbour-normalised methylation by 5S-gene variant class.
# Author: Jonathan Bohlen
# (c) 2026 Jonathan Bohlen. Code accompanying Sengl et al. (2026),
# "Surveillance and selection of 5S ribosomal RNA genes in the human genome."
# Released under the MIT License; see LICENSE at the repository root.
# -----------------------------------------------------------------------------
"""Two-panel violin plots of whole-copy CpG methylation by 5S-gene variant consequence:
(A) raw whole-copy CpG methylation;
(B) neighbour-normalised methylation (copy minus the mean of its +/-K flanking copies).
Groups are compared by Kruskal-Wallis across classes and by Mann-Whitney U against the
no-gene-SNV baseline.

Input : 5S_rDNA.db (tables variant, copy_methylation, copy);
        97_master_table.tsv (saturation-mutagenesis functional classes).
Output: <FIVES_OUT>/Figure3/FGH_violin_raw_vs_neighbor_normalized.pdf

Paths are read from environment variables: FIVES_DB, FIVES_DATA, FIVES_OUT."""
import os
import sqlite3, numpy as np, pandas as pd, matplotlib.pyplot as plt
from scipy.stats import kruskal, mannwhitneyu, gaussian_kde
plt.rcParams.update({"font.family":"Arial","pdf.fonttype":42,"ps.fonttype":42,"axes.linewidth":1.0,"font.size":8})
CM=1/2.54; K=3
DB=os.environ.get("FIVES_DB","5S_rDNA.db")
MT=os.environ.get("FIVES_DATA","data")+"/97_master_table.tsv"
OUT=os.environ.get("FIVES_OUT","output")+"/Figure3"
con=sqlite3.connect(DB); mt=pd.read_csv(MT,sep="\t")[["pos","alt","class"]]
gv=pd.read_sql("SELECT copy_id,consensus_pos,alt FROM variant WHERE region='gene' AND alt IN ('A','C','G','T')",con).merge(mt,left_on=["consensus_pos","alt"],right_on=["pos","alt"])
cnt=gv.groupby(["copy_id","class"]).size().unstack(fill_value=0).reset_index()
for c in ["competent","expr_defective","incorp_defective"]:
    if c not in cnt: cnt[c]=0
cnt=cnt.rename(columns={"competent":"n_neutral","expr_defective":"n_exprdef","incorp_defective":"n_incdef"})
cm=pd.read_sql("SELECT copy_id,mean_meth FROM copy_methylation WHERE n_conf_calls>=10",con)
cp=pd.read_sql("SELECT copy_id,haplotype_id,unit_start_local FROM copy WHERE array_member=1",con)
d=cm.merge(cp,on="copy_id").merge(cnt[["copy_id","n_neutral","n_exprdef","n_incdef"]],on="copy_id",how="left").fillna(0)
d["n_total"]=d.n_neutral+d.n_exprdef+d.n_incdef
d=d.sort_values(["haplotype_id","unit_start_local"]).reset_index(drop=True)
nb=np.full(len(d),np.nan); vals=d.mean_meth.values
for h,g in d.groupby("haplotype_id"):
    idx=g.index.values; m=vals[idx]; n=len(m)
    for i in range(n):
        j=[k for k in range(max(0,i-K),min(n,i+K+1)) if k!=i]
        if j: nb[idx[i]]=m[j].mean()
d["neighbor_meth"]=nb; d=d.dropna(subset=["neighbor_meth"])
d["raw"]=d.mean_meth*100; d["norm"]=(d.mean_meth-d.neighbor_meth)*100

groups=[("no gene\nSNV",d[d.n_total==0],"#999999"),("neutral",d[d.n_neutral>=1],"#4c9e4c"),
        ("expr-def",d[d.n_exprdef>=1],"#c0392b"),("incorp-def",d[d.n_incdef>=1],"#d0a13a")]
def violin(ax,col,ylab,ref,reflab,ylim,title):
    data=[g[col].values for _,g,_ in groups]; base=data[0]; kw=kruskal(*data)[1]
    for i,((lab,g,c),y) in enumerate(zip(groups,data)):
        p10,q25,med,q75,p90=np.percentile(y,[10,25,50,75,90])
        kde=gaussian_kde(y); ys=np.linspace(ylim[0],ylim[1],200); dd=kde(ys); dd=dd/dd.max()*0.36
        ax.fill_betweenx(ys,i-dd,i+dd,color=c,alpha=0.6,lw=1.0,edgecolor=c,zorder=2)
        ax.plot([i,i],[p10,p90],color="k",lw=1.0,zorder=3); ax.plot([i,i],[q25,q75],color="k",lw=3.2,solid_capstyle="round",zorder=3)
        ax.scatter([i],[med],s=16,color="w",edgecolor="k",lw=1.0,zorder=4)
        if i>0:
            p=mannwhitneyu(y,base)[1]; s="***" if p<1e-3 else "**" if p<1e-2 else "*" if p<0.05 else "ns"
            ax.text(i,ylim[1]-(ylim[1]-ylim[0])*0.04,s,ha="center",fontsize=8)
    ax.axhline(ref,color="#c0392b",ls="--",lw=1.0,zorder=1); ax.text(3.45,ref+(ylim[1]-ylim[0])*0.01,reflab,color="#c0392b",fontsize=8,va="bottom",ha="right")
    ax.set_xticks(range(4)); ax.set_xticklabels(["no SNV","neutral","expr-def","inc-def"],fontsize=8,rotation=35,ha="right"); ax.set_ylim(*ylim)
    ax.set_ylabel(ylab,fontsize=8); ax.set_title(f"{title}\nKW p={kw:.1e}",fontsize=8); ax.tick_params(labelsize=8)
    for s in ["top","right"]: ax.spines[s].set_visible(False)

# explicit 4x4 cm PLOTTING AREA per panel (axes size fixed; labels/title live in the margins)
FIGW,FIGH=13*CM,7*CM; fig=plt.figure(figsize=(FIGW,FIGH))
aw,ah=4*CM/FIGW,4*CM/FIGH; lm=1.9*CM/FIGW; bm=2.0*CM/FIGH; gap=2.3*CM/FIGW
axA=fig.add_axes([lm,bm,aw,ah]); axB=fig.add_axes([lm+aw+gap,bm,aw,ah])
violin(axA,"raw","CpG methylation (%)",65,"<65%",(0,100),"Raw copy methylation")
lo,hi=np.percentile(d.norm,[0.5,99.5])
violin(axB,"norm","methylation − neighbours (Δ%)",0,"local mean",(min(lo,-hi),max(hi,-lo)),"Neighbour-normalised")
fig.savefig(f"{OUT}/FGH_violin_raw_vs_neighbor_normalized.pdf",bbox_inches="tight",dpi=300)
print("raw medians:",{l:round(np.median(g.raw),1) for l,g,_ in groups})
print("neighbour-normalised medians:",{l:round(np.median(g.norm),2) for l,g,_ in groups})
print("wrote FGH_violin_raw_vs_neighbor_normalized.pdf")
