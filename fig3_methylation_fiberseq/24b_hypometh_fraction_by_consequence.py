#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# 24b_hypometh_fraction_by_consequence.py — Fraction of lowly-methylated copies by functional class of carried 5S-gene variant.
# Author: Jonathan Bohlen
# (c) 2026 Jonathan Bohlen. Code accompanying Sengl et al. (2026),
# "Surveillance and selection of 5S ribosomal RNA genes in the human genome."
# Released under the MIT License; see LICENSE at the repository root.
# -----------------------------------------------------------------------------
"""Fraction of lowly-methylated copies (whole-copy CpG methylation < 65%) among copies carrying
each functional class of 5S-gene variant, versus copies with no gene SNV.
Reports raw % (pooled copies), donor-averaged % (pseudoreplication control), and a
donor-clustered logistic test of each class flag against the no-gene-SNV baseline.

Input : 5S_rDNA.db (tables variant, copy_methylation, copy, haplotype, assembly);
        97_master_table.tsv (saturation-mutagenesis functional classes).
Output: <FIVES_OUT>/Figure3/FGH_hypometh_fraction_by_consequence.{tsv,pdf}

Paths are read from environment variables: FIVES_DB, FIVES_DATA, FIVES_OUT."""
import os
import sqlite3, numpy as np, pandas as pd, statsmodels.api as sm, matplotlib.pyplot as plt
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
cp=pd.read_sql("SELECT copy_id,haplotype_id FROM copy WHERE array_member=1",con)
hap=pd.read_sql("SELECT haplotype_id,assembly_id FROM haplotype",con)
asm=pd.read_sql("SELECT assembly_id,sample_id AS donor FROM assembly",con)
d=cm.merge(cp,on="copy_id").merge(hap,on="haplotype_id").merge(asm,on="assembly_id")
d=d.merge(cnt[["copy_id","n_neutral","n_exprdef","n_incdef"]],on="copy_id",how="left").fillna(0)
d["n_total"]=d.n_neutral+d.n_exprdef+d.n_incdef
d["hypo"]=(d.mean_meth<0.65).astype(int)

# groups: no gene SNV (baseline) + carries >=1 of each class (non-exclusive)
groups={"no gene SNV":d.n_total==0,"carries neutral":d.n_neutral>=1,
        "carries expr-def":d.n_exprdef>=1,"carries incorp-def":d.n_incdef>=1}
base=d[d.n_total==0]
rows=[]
for name,mask in groups.items():
    g=d[mask]
    perdon=g.groupby("donor").hypo.mean()          # per-donor % hypo, then average (pseudorep control)
    rows.append(dict(group=name,n_copies=len(g),n_donor=g.donor.nunique(),
                     pct_hypo_raw=100*g.hypo.mean(),
                     pct_hypo_donoravg=100*perdon.mean(),donoravg_se=100*perdon.sem()))
res=pd.DataFrame(rows)
# donor-clustered logistic test of each class flag vs baseline (no gene SNV)
for cls,col in [("carries neutral","n_neutral"),("carries expr-def","n_exprdef"),("carries incorp-def","n_incdef")]:
    sub=d[(d.n_total==0)|(d[col]>=1)].copy(); sub["flag"]=(sub[col]>=1).astype(int)
    m=sm.GLM(sub.hypo,sm.add_constant(sub.flag),family=sm.families.Binomial()).fit(cov_type="cluster",cov_kwds={"groups":sub.donor})
    res.loc[res.group==cls,"OR_vs_baseline"]=np.exp(m.params[1]); res.loc[res.group==cls,"p_vs_baseline"]=m.pvalues[1]
res.to_csv(f"{OUT}/FGH_hypometh_fraction_by_consequence.tsv",sep="\t",index=False)
print(res.round(3).to_string(index=False))

# bar figure: % hypomethylated copies by group (donor-averaged, SE)
fig,ax=plt.subplots(figsize=(8*CM,5*CM))
cols=["#999999","#4c9e4c","#c0392b","#d0a13a"]
y=res.pct_hypo_donoravg.values; e=res.donoravg_se.values
ax.bar(range(len(res)),y,yerr=e,color=cols,width=0.66,error_kw=dict(lw=0.8,capsize=2))
base_pct=res.iloc[0].pct_hypo_donoravg; ax.axhline(base_pct,color="grey",ls="--",lw=0.7)
for i,r in res.iterrows():
    if pd.notna(r.get("p_vs_baseline")):
        s="***" if r.p_vs_baseline<1e-3 else "**" if r.p_vs_baseline<1e-2 else "*" if r.p_vs_baseline<0.05 else "ns"
        ax.text(i,y[i]+e[i]+0.6,s,ha="center",fontsize=6)
ax.set_xticks(range(len(res))); ax.set_xticklabels([g.replace("carries ","") for g in res.group],fontsize=6.5)
ax.set_ylabel("% lowly-methylated copies\n(whole-copy CpG < 65%)",fontsize=6.5)
ax.set_title("Hypomethylated-copy fraction by 5S-gene variant consequence",fontsize=6.5)
ax.tick_params(labelsize=6)
for s in ["top","right"]: ax.spines[s].set_visible(False)
fig.tight_layout(); fig.savefig(f"{OUT}/FGH_hypometh_fraction_by_consequence.pdf",bbox_inches="tight")
print("\nwrote FGH_hypometh_fraction_by_consequence.pdf + .tsv")
