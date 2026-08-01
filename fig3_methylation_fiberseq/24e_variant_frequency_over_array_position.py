#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# 24e_variant_frequency_over_array_position.py — 5S-gene variant frequency as a function of fractional array position.
# Author: Jonathan Bohlen
# (c) 2026 Jonathan Bohlen. Code accompanying Sengl et al. (2026),
# "Surveillance and selection of 5S ribosomal RNA genes in the human genome."
# Released under the MIT License; see LICENSE at the repository root.
# -----------------------------------------------------------------------------
"""5S-gene variant frequency as a function of fractional array position.
x = fractional array position (0..1, both ends = edges, 0.5 = centre).
(A) mean gene SNVs per copy (total).  (B) relative enrichment per functional class.
Restricted to arrays with n_copies>=10 so array position is resolvable.

Input : 5S_rDNA.db (tables variant, copy, haplotype);
        97_master_table.tsv (saturation-mutagenesis functional classes).
Output: <FIVES_OUT>/Figure3/FGH_variant_frequency_over_array_position.pdf

Paths are read from environment variables: FIVES_DB, FIVES_DATA, FIVES_OUT."""
import os
import sqlite3, numpy as np, pandas as pd, matplotlib.pyplot as plt
plt.rcParams.update({"font.family":"Arial","pdf.fonttype":42,"ps.fonttype":42,"axes.linewidth":1.0,"font.size":8})
CM=1/2.54
DB=os.environ.get("FIVES_DB","5S_rDNA.db")
MT=os.environ.get("FIVES_DATA","data")+"/97_master_table.tsv"
OUT=os.environ.get("FIVES_OUT","output")+"/Figure3"
con=sqlite3.connect(DB); mt=pd.read_csv(MT,sep="\t")[["pos","alt","class"]]
gv=pd.read_sql("SELECT copy_id,consensus_pos,alt FROM variant WHERE region='gene' AND alt IN ('A','C','G','T')",con).merge(mt,left_on=["consensus_pos","alt"],right_on=["pos","alt"])
cnt=gv.groupby(["copy_id","class"]).size().unstack(fill_value=0).reset_index()
for c in ["competent","expr_defective","incorp_defective"]:
    if c not in cnt: cnt[c]=0
cnt=cnt.rename(columns={"competent":"n_neutral","expr_defective":"n_exprdef","incorp_defective":"n_incdef"})
cp=pd.read_sql("SELECT copy_id,haplotype_id,copy_number FROM copy WHERE array_member=1",con)
hap=pd.read_sql("SELECT haplotype_id,n_copies FROM haplotype",con)
d=cp.merge(hap,on="haplotype_id"); d=d[d.n_copies>=10].copy()
d["pos"]=(d.copy_number-0.5)/d.n_copies                                  # 0..1 along array
d=d.merge(cnt[["copy_id","n_neutral","n_exprdef","n_incdef"]],on="copy_id",how="left").fillna(0)
d["n_total"]=d.n_neutral+d.n_exprdef+d.n_incdef
NB=20; d["bin"]=np.clip((d.pos*NB).astype(int),0,NB-1); xc=(np.arange(NB)+0.5)/NB

fig,(axA,axB)=plt.subplots(1,2,figsize=(14*CM,4.3*CM))
# A: total variants/copy vs position
g=d.groupby("bin").n_total; m=g.mean().reindex(range(NB)); se=g.sem().reindex(range(NB))
axA.fill_between(xc,m-se,m+se,color="#3b6fb0",alpha=0.25,lw=0)
axA.plot(xc,m,color="#3b6fb0",lw=1.0,marker="o",ms=3)
axA.axhline(d.n_total.mean(),color="grey",ls="--",lw=1.0)
axA.set_ylim(0,(m+se).max()*1.1)
axA.set_xlabel("array position (0,1 = edges · 0.5 = centre)",fontsize=8)
axA.set_ylabel("5S-gene SNVs per copy",fontsize=8); axA.set_title("Total gene variant frequency",fontsize=8)
axA.tick_params(labelsize=8)
# B: relative enrichment (bin mean / overall mean) per class
for col,lab,c in [("n_total","total","#3b6fb0"),("n_neutral","neutral","#4c9e4c"),("n_exprdef","expr-def","#c0392b"),("n_incdef","incorp-def","#d0a13a")]:
    mm=d.groupby("bin")[col].mean().reindex(range(NB))/d[col].mean()
    axB.plot(xc,mm,lw=1.0,marker="o",ms=2.5,color=c,label=lab)
axB.axhline(1,color="grey",ls="--",lw=1.0)
axB.set_xlabel("array position (0,1 = edges · 0.5 = centre)",fontsize=8)
axB.set_ylabel("relative enrichment\n(bin mean / array mean)",fontsize=8); axB.set_title("By functional consequence",fontsize=8)
axB.legend(frameon=False,fontsize=7,loc="upper center",ncol=2); axB.tick_params(labelsize=8)
for ax in (axA,axB):
    for s in ["top","right"]: ax.spines[s].set_visible(False)
fig.tight_layout(); fig.savefig(f"{OUT}/FGH_variant_frequency_over_array_position.pdf",bbox_inches="tight",dpi=300)
# quantify edge vs middle
edge=d[(d.pos<0.15)|(d.pos>0.85)].n_total.mean(); mid=d[(d.pos>0.40)&(d.pos<0.60)].n_total.mean()
print(f"mean gene SNVs/copy: EDGE (outer 15%)={edge:.3f}  MIDDLE (40-60%)={mid:.3f}  -> {edge/mid:.1f}x higher at edges")
for col,lab in [("n_neutral","neutral"),("n_exprdef","expr-def"),("n_incdef","incorp-def")]:
    e=d[(d.pos<0.15)|(d.pos>0.85)][col].mean(); mm=d[(d.pos>0.40)&(d.pos<0.60)][col].mean()
    print(f"  {lab:10s}: edge={e:.3f} middle={mm:.3f} -> {e/mm:.1f}x")
print("wrote FGH_variant_frequency_over_array_position.pdf")
