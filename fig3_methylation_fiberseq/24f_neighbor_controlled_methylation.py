#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# 24f_neighbor_controlled_methylation.py — Neighbor-controlled test of methylation association with 5S-gene variant functional class.
# Author: Jonathan Bohlen
# (c) 2026 Jonathan Bohlen. Code accompanying Sengl et al. (2026),
# "Surveillance and selection of 5S ribosomal RNA genes in the human genome."
# Released under the MIT License; see LICENSE at the repository root.
# -----------------------------------------------------------------------------
"""Neighbour-controlled test of methylation association with 5S-gene variant functional class.
Each copy is controlled against the local methylation of its neighbouring copies in the same array
(capturing the per-array pattern, including telomere/centromere asymmetry and array-size
differences, non-parametrically).
neighbor_meth = mean methylation of the +/-K copies flanking each copy (by genomic order, excl. self).
Tests whether a functional-consequence class predicts methylation beyond the local baseline, via a
donor-clustered logistic model fitted with and without the neighbour term.

Input : 5S_rDNA.db (tables variant, copy_methylation, copy, haplotype, assembly);
        97_master_table.tsv (saturation-mutagenesis functional classes).
Output: console summary.

Paths are read from environment variables: FIVES_DB, FIVES_DATA, FIVES_OUT."""
import os
import sqlite3, numpy as np, pandas as pd, statsmodels.api as sm
from scipy.stats import mannwhitneyu
K=3
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
hap=pd.read_sql("SELECT haplotype_id,assembly_id FROM haplotype",con)
asm=pd.read_sql("SELECT assembly_id,sample_id AS donor FROM assembly",con)
d=cm.merge(cp,on="copy_id").merge(hap,on="haplotype_id").merge(asm,on="assembly_id")
d=d.merge(cnt[["copy_id","n_neutral","n_exprdef","n_incdef"]],on="copy_id",how="left").fillna(0)
d["n_total"]=d.n_neutral+d.n_exprdef+d.n_incdef

# --- neighbour methylation: mean of +/-K flanking copies in the same array (by genomic order), excl. self ---
d=d.sort_values(["haplotype_id","unit_start_local"]).reset_index(drop=True)
nb=np.full(len(d),np.nan); vals=d.mean_meth.values; hid=d.haplotype_id.values
for h,g in d.groupby("haplotype_id"):
    idx=g.index.values; m=vals[idx]; n=len(m)
    for i in range(n):
        lo=max(0,i-K); hi=min(n,i+K+1); j=[k for k in range(lo,hi) if k!=i]
        if j: nb[idx[i]]=m[j].mean()
d["neighbor_meth"]=nb
d=d.dropna(subset=["neighbor_meth"])
d["resid"]=d.mean_meth-d.neighbor_meth          # copy meth relative to its local array environment
d["hypo"]=(d.mean_meth<0.65).astype(int)
for c in ["n_neutral","n_exprdef","n_incdef"]: d[c+"_f"]=(d[c]>=1).astype(int)
print(f"copies with neighbour baseline: {len(d)}  (K={K})")

# --- descriptive: methylation RELATIVE TO NEIGHBOURS by group (residual, in % points) ---
print("\nmean methylation MINUS local neighbour mean (negative = hypomethylated vs local array context):")
base=d[d.n_total==0].resid.values
for lab,mask in [("no gene SNV",d.n_total==0),("neutral",d.n_neutral>=1),("expr-def",d.n_exprdef>=1),("incorp-def",d.n_incdef>=1)]:
    r=d[mask].resid.values*100
    p=np.nan if lab=="no gene SNV" else mannwhitneyu(d[mask].resid,base)[1]
    print(f"   {lab:12s} n={len(r):6d}  resid={r.mean():+.2f}%  (MWU vs baseline p={'' if lab=='no gene SNV' else f'{p:.2e}'})")

# --- model: hypomethylation controlling for LOCAL neighbour methylation (donor-clustered) ---
def fit(cols):
    m=sm.GLM(d.hypo,sm.add_constant(d[cols]),family=sm.families.Binomial()).fit(cov_type="cluster",cov_kwds={"groups":d.donor})
    return {c:(np.exp(m.params[c]),m.pvalues[c]) for c in cols}
naive=fit(["n_neutral_f","n_exprdef_f","n_incdef_f"])
ctrl =fit(["n_neutral_f","n_exprdef_f","n_incdef_f","neighbor_meth"])
print("\nhypomethylation OR (OR<1 = more methylated), naive vs neighbour-controlled:")
print(f"{'class':10s}{'OR naive':>20s}{'OR + neighbour ctrl':>24s}")
for c in ["n_neutral_f","n_exprdef_f","n_incdef_f"]:
    print(f"{c.replace('_f',''):10s}{naive[c][0]:>10.3f} (p={naive[c][1]:.1e}){ctrl[c][0]:>13.3f} (p={ctrl[c][1]:.1e})")
print(f"{'neighbor':10s}{'':>20s}{ctrl['neighbor_meth'][0]:>13.3g} (p={ctrl['neighbor_meth'][1]:.1e})")
