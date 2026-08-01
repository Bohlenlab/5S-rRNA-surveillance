#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# 24d_array_positional_bias.py — Tests whether 5S-gene variants are biased toward array edges (positional confounder check).
# Author: Jonathan Bohlen
# (c) 2026 Jonathan Bohlen. Code accompanying Sengl et al. (2026),
# "Surveillance and selection of 5S ribosomal RNA genes in the human genome."
# Released under the MIT License; see LICENSE at the repository root.
# -----------------------------------------------------------------------------
"""Tests whether 5S-gene variants (total, and each functional subclass) are biased toward array
edges, a positional confounder for the variant-methylation association.
edge_score per copy = 2*|pos_frac-0.5| in [0,1]  (0 = array centre, 1 = array edge).
Restricted to arrays with n_copies>=10 so array position is resolvable; groups are compared by
Mann-Whitney U against the no-gene-SNV baseline.

Input : 5S_rDNA.db (tables variant, copy, haplotype);
        97_master_table.tsv (saturation-mutagenesis functional classes).
Output: <FIVES_OUT>/Figure3/FGH_array_positional_bias.{tsv,pdf}

Paths are read from environment variables: FIVES_DB, FIVES_DATA, FIVES_OUT."""
import os
import sqlite3, numpy as np, pandas as pd, matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu
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
cp=pd.read_sql("SELECT copy_id,haplotype_id,copy_number FROM copy WHERE array_member=1",con)
hap=pd.read_sql("SELECT haplotype_id,n_copies FROM haplotype",con)
d=cp.merge(hap,on="haplotype_id"); d=d[d.n_copies>=10].copy()          # need enough copies for position
d["pos_frac"]=(d.copy_number-0.5)/d.n_copies
d["edge_score"]=2*(d.pos_frac-0.5).abs()                                # 0=centre, 1=edge
d=d.merge(cnt[["copy_id","n_neutral","n_exprdef","n_incdef"]],on="copy_id",how="left").fillna(0)
d["n_total"]=d.n_neutral+d.n_exprdef+d.n_incdef
print(f"copies with resolved position (n_copies>=10): {len(d)}")

base=d[d.n_total==0].edge_score.values
groups=[("no gene SNV",d[d.n_total==0]),("any gene SNV",d[d.n_total>=1]),
        ("neutral",d[d.n_neutral>=1]),("expr-def",d[d.n_exprdef>=1]),("incorp-def",d[d.n_incdef>=1])]
print(f"\n{'group':14s}{'n':>7s}{'mean edge_score':>16s}{'% in outer 40%':>16s}{'MWU vs baseline':>18s}")
rows=[]
for lab,g in groups:
    es=g.edge_score.values; pe=100*(es>0.6).mean()
    p=np.nan if lab=="no gene SNV" else mannwhitneyu(es,base)[1]
    print(f"{lab:14s}{len(g):>7d}{es.mean():>16.3f}{pe:>15.1f}%{('' if lab=='no gene SNV' else f'{p:.2e}'):>18s}")
    rows.append(dict(group=lab,n=len(g),mean_edge=es.mean(),pct_outer40=pe,p_vs_baseline=p))
pd.DataFrame(rows).to_csv(f"{OUT}/FGH_array_positional_bias.tsv",sep="\t",index=False)

# figure: edge_score distribution by group (violin), edge=1 at top
fig,ax=plt.subplots(figsize=(8*CM,5.5*CM))
from scipy.stats import gaussian_kde
cols=["#999999","#3b6fb0","#4c9e4c","#c0392b","#d0a13a"]
for i,((lab,g),col) in enumerate(zip(groups,cols)):
    y=g.edge_score.values; kde=gaussian_kde(y); ys=np.linspace(0,1,150); dd=kde(ys); dd=dd/dd.max()*0.38
    ax.fill_betweenx(ys,i-dd,i+dd,color=col,alpha=0.6,lw=0.5,edgecolor=col)
    ax.scatter([i],[np.median(y)],s=14,color="w",edgecolor="k",lw=0.8,zorder=4)
    if lab!="no gene SNV":
        p=mannwhitneyu(y,base)[1]; s="***" if p<1e-3 else "**" if p<1e-2 else "*" if p<0.05 else "ns"
        ax.text(i,1.04,s,ha="center",fontsize=6.5)
ax.set_xticks(range(len(groups))); ax.set_xticklabels([g[0] for g in groups],fontsize=6,rotation=20,ha="right")
ax.set_ylabel("array position\n(0 = centre, 1 = edge)",fontsize=6.5); ax.set_ylim(-0.05,1.12)
ax.set_title("5S-gene variants vs array position",fontsize=6.5); ax.tick_params(labelsize=6)
for s in ["top","right"]: ax.spines[s].set_visible(False)
fig.tight_layout(); fig.savefig(f"{OUT}/FGH_array_positional_bias.pdf",bbox_inches="tight")
print("\nwrote FGH_array_positional_bias.pdf + .tsv")
