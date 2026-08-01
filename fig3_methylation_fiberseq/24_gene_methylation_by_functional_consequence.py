#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# 24_gene_methylation_by_functional_consequence.py — Stratifies whole-copy CpG methylation by the functional consequence of 5S-gene SNVs.
# Author: Jonathan Bohlen
# (c) 2026 Jonathan Bohlen. Code accompanying Sengl et al. (2026),
# "Surveillance and selection of 5S ribosomal RNA genes in the human genome."
# Released under the MIT License; see LICENSE at the repository root.
# -----------------------------------------------------------------------------
"""Stratifies whole-copy CpG methylation by the functional consequence of 5S-gene SNVs.

Per copy, gene SNVs are counted in each functional class (from the saturation-mutagenesis
master table: competent / expr_defective / incorp_defective) and tested against methylation,
controlling for donor pseudoreplication:
  - joint donor-clustered logistic:  hypomethylated ~ z(n_competent)+z(n_exprdef)+z(n_incdef)
  - per-class per-donor Spearman(meth, count) -> one-sample Wilcoxon (donor = unit)
Outcome = whole-copy CpG methylation; hypomethylated = mean_meth < 0.65.

Input : 5S_rDNA.db (tables variant, copy_methylation, copy, haplotype, assembly);
        97_master_table.tsv (saturation-mutagenesis functional classes).
Output: <FIVES_OUT>/Figure3/FGH_methylation_by_functional_consequence.tsv

Paths are read from environment variables: FIVES_DB, FIVES_DATA, FIVES_OUT."""
import os
import sqlite3, numpy as np, pandas as pd, statsmodels.api as sm
from scipy.stats import spearmanr, wilcoxon
DB=os.environ.get("FIVES_DB","5S_rDNA.db")
MT=os.environ.get("FIVES_DATA","data")+"/97_master_table.tsv"
OUT=os.environ.get("FIVES_OUT","output")+"/Figure3"
con=sqlite3.connect(DB)
mt=pd.read_csv(MT,sep="\t")[["pos","alt","class"]]

# --- per-copy gene-SNV counts by functional class ---
gv=pd.read_sql("SELECT copy_id,consensus_pos,alt FROM variant WHERE region='gene' AND alt IN ('A','C','G','T')",con)
gv=gv.merge(mt,left_on=["consensus_pos","alt"],right_on=["pos","alt"],how="inner")
cnt=gv.groupby(["copy_id","class"]).size().unstack(fill_value=0).reset_index()
for c in ["competent","expr_defective","incorp_defective"]:
    if c not in cnt: cnt[c]=0
cnt=cnt.rename(columns={"competent":"n_neutral","expr_defective":"n_exprdef","incorp_defective":"n_incdef"})

# --- copies with methylation + donor, filtered (array_member=1, n_conf_calls>=10) ---
cm=pd.read_sql("SELECT copy_id,mean_meth,n_conf_calls FROM copy_methylation WHERE n_conf_calls>=10",con)
cp=pd.read_sql("SELECT copy_id,haplotype_id,array_member FROM copy WHERE array_member=1",con)
hap=pd.read_sql("SELECT haplotype_id,assembly_id FROM haplotype",con)
asm=pd.read_sql("SELECT assembly_id,sample_id AS donor FROM assembly",con)
d=cm.merge(cp,on="copy_id").merge(hap,on="haplotype_id").merge(asm,on="assembly_id")
d=d.merge(cnt[["copy_id","n_neutral","n_exprdef","n_incdef"]],on="copy_id",how="left").fillna({"n_neutral":0,"n_exprdef":0,"n_incdef":0})
d["n_total"]=d.n_neutral+d.n_exprdef+d.n_incdef
d["hypo"]=(d.mean_meth<0.65).astype(int)
print(f"copies={len(d)}  donors={d.donor.nunique()}  hypometh={d.hypo.mean():.2%}")
print(f"per-copy gene SNV counts (mean): neutral={d.n_neutral.mean():.3f}  exprdef={d.n_exprdef.mean():.3f}  incdef={d.n_incdef.mean():.3f}")

def z(x): return (x-x.mean())/x.std() if x.std()>0 else x*0
# --- total gene-SNV count slope (donor-clustered logistic) ---
m0=sm.GLM(d.hypo,sm.add_constant(z(d.n_total)),family=sm.families.Binomial()).fit(cov_type="cluster",cov_kwds={"groups":d.donor})
print(f"\n[sanity] hypo ~ z(total gene SNV): OR/SD={np.exp(m0.params[1]):.3f} p={m0.pvalues[1]:.2e}  (OR<1 = more variants->more methylation)")

# --- joint donor-clustered logistic by class ---
X=pd.DataFrame({"neutral":z(d.n_neutral),"exprdef":z(d.n_exprdef),"incdef":z(d.n_incdef)})
mj=sm.GLM(d.hypo,sm.add_constant(X),family=sm.families.Binomial()).fit(cov_type="cluster",cov_kwds={"groups":d.donor})
print("\n[JOINT logistic] hypomethylation odds per +1 SD of each class (OR<1 = MORE methylation):")
for c in ["neutral","exprdef","incdef"]:
    print(f"   {c:8s}: OR/SD={np.exp(mj.params[c]):.3f}  p={mj.pvalues[c]:.2e}")

# --- per-class per-donor Spearman -> Wilcoxon (donor = unit) ---
print("\n[per-donor] Spearman(mean_meth, class count) then Wilcoxon vs 0 (positive rho = more variants->more meth):")
res=[]
for c,col in [("neutral","n_neutral"),("exprdef","n_exprdef"),("incdef","n_incdef")]:
    rhos=[]
    for don,g in d.groupby("donor"):
        if g[col].nunique()<2 or len(g)<5: continue
        r=spearmanr(g.mean_meth,g[col])[0]
        if np.isfinite(r): rhos.append(r)
    rhos=np.array(rhos); W,p=wilcoxon(rhos)
    print(f"   {c:8s}: median rho={np.median(rhos):+.3f}  n_donor={len(rhos)}  Wilcoxon p={p:.2e}")
    res.append(dict(cls=c,median_rho=np.median(rhos),n_donor=len(rhos),wilcoxon_p=p,
                    joint_OR=np.exp(mj.params[c]),joint_p=mj.pvalues[c]))
pd.DataFrame(res).to_csv(f"{OUT}/FGH_methylation_by_functional_consequence.tsv",sep="\t",index=False)
print(f"\nwrote {OUT}/FGH_methylation_by_functional_consequence.tsv")
