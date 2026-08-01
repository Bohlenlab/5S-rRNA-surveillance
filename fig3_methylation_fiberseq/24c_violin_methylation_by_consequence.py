#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# 24c_violin_methylation_by_consequence.py — Violin plots of whole-copy CpG methylation by functional class of carried 5S-gene variant.
# Author: Jonathan Bohlen
# (c) 2026 Jonathan Bohlen. Code accompanying Sengl et al. (2026),
# "Surveillance and selection of 5S ribosomal RNA genes in the human genome."
# Released under the MIT License; see LICENSE at the repository root.
# -----------------------------------------------------------------------------
"""Violin plots of whole-copy CpG methylation for the four copy classes:
no gene SNV / carries neutral / carries expr-def / carries incorp-def gene variant.
Functional class from saturation mutagenesis, DEF_CUT=0.50 (<50% of WT = defective).

Input : 5S_rDNA.db (tables variant, copy_methylation, copy);
        97_master_table.tsv (saturation-mutagenesis functional classes).
Output: <FIVES_OUT>/Figure3/FGH_violin_methylation_by_consequence.pdf

Paths are read from environment variables: FIVES_DB, FIVES_DATA, FIVES_OUT."""
import os
import sqlite3, numpy as np, pandas as pd, matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from scipy.stats import kruskal, mannwhitneyu, gaussian_kde
plt.rcParams.update({"font.family":"Arial","pdf.fonttype":42,"ps.fonttype":42,"axes.linewidth":0.8,"font.size":7})
CM=1/2.54
DB=os.environ.get("FIVES_DB","5S_rDNA.db")
MT=os.environ.get("FIVES_DATA","data")+"/97_master_table.tsv"
OUT=os.environ.get("FIVES_OUT","output")+"/Figure3"
con=sqlite3.connect(DB); mt=pd.read_csv(MT,sep="\t")[["pos","alt","class"]]
gv=pd.read_sql("SELECT copy_id,consensus_pos,alt FROM variant WHERE region='gene' AND alt IN ('A','C','G','T')",con)
gv=gv.merge(mt,left_on=["consensus_pos","alt"],right_on=["pos","alt"],how="inner")
cnt=gv.groupby(["copy_id","class"]).size().unstack(fill_value=0).reset_index()
for c in ["competent","expr_defective","incorp_defective"]:
    if c not in cnt: cnt[c]=0
cnt=cnt.rename(columns={"competent":"n_neutral","expr_defective":"n_exprdef","incorp_defective":"n_incdef"})
cm=pd.read_sql("SELECT copy_id,mean_meth FROM copy_methylation WHERE n_conf_calls>=10",con)
cp=pd.read_sql("SELECT copy_id FROM copy WHERE array_member=1",con)
d=cm.merge(cp,on="copy_id").merge(cnt[["copy_id","n_neutral","n_exprdef","n_incdef"]],on="copy_id",how="left").fillna(0)
d["n_total"]=d.n_neutral+d.n_exprdef+d.n_incdef; d["meth"]=d.mean_meth*100

groups=[("no gene\nSNV",d[d.n_total==0],"#999999"),("neutral",d[d.n_neutral>=1],"#4c9e4c"),
        ("expr-def",d[d.n_exprdef>=1],"#c0392b"),("incorp-def",d[d.n_incdef>=1],"#d0a13a")]
data=[g.meth.values for _,g,_ in groups]; base=data[0]
kw=kruskal(*data)[1]

fig,ax=plt.subplots(figsize=(8*CM,6*CM))
for i,((lab,g,col),y) in enumerate(zip(groups,data)):
    p10,q25,med,q75,p90=np.percentile(y,[10,25,50,75,90])
    kde=gaussian_kde(y); ys=np.linspace(0,100,200); dd=kde(ys); dd=dd/dd.max()*0.36
    ax.fill_betweenx(ys,i-dd,i+dd,color=col,alpha=0.6,lw=0.6,edgecolor=col,zorder=2)
    ax.plot([i,i],[p10,p90],color="k",lw=0.8,zorder=3)
    ax.plot([i,i],[q25,q75],color="k",lw=3.5,solid_capstyle="round",zorder=3)
    ax.scatter([i],[med],s=16,color="w",edgecolor="k",lw=0.8,zorder=4)
    ax.text(i,-8,f"n={len(y):,}",ha="center",fontsize=5,color="#444")
    if i>0:
        p=mannwhitneyu(y,base)[1]; s="***" if p<1e-3 else "**" if p<1e-2 else "*" if p<0.05 else "ns"
        ax.text(i,103,s,ha="center",fontsize=7)
ax.axhline(65,color="#c0392b",ls="--",lw=0.8,zorder=1); ax.text(3.45,66,"<65%",color="#c0392b",fontsize=5.5,va="bottom",ha="right")
ax.set_xticks(range(4)); ax.set_xticklabels([g[0] for g in groups],fontsize=6.5)
ax.set_ylim(-12,112); ax.set_ylabel("CpG methylation (%)",fontsize=7); ax.set_yticks([0,20,40,60,80,100])
ax.set_title(f"Whole-copy methylation by 5S-gene variant consequence\nKruskal-Wallis p={kw:.1e}",fontsize=6.5)
ax.tick_params(labelsize=6)
for s in ["top","right"]: ax.spines[s].set_visible(False)
fig.tight_layout(); fig.savefig(f"{OUT}/FGH_violin_methylation_by_consequence.pdf",bbox_inches="tight")
print("median methylation %:", {lab:round(np.median(g.meth),1) for lab,g,_ in groups})
print(f"KW p={kw:.2e}; wrote FGH_violin_methylation_by_consequence.pdf")
