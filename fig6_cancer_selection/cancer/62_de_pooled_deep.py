#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# 62_de_pooled_deep.py — pooled deep-coverage TCGA differential expression of the rna_excess predictor with a TP53 interaction and cytosolic-RP readout.
# Author: Jonathan Bohlen
# (c) 2026 Jonathan Bohlen. Code accompanying Sengl et al. (2026),
# "Surveillance and selection of 5S ribosomal RNA genes in the human genome."
# Released under the MIT License; see LICENSE at the repository root.
# -----------------------------------------------------------------------------
"""Pooled deep-coverage TCGA differential expression: single DESeq2 over all samples with RNA gene-region
depth >= DE_MIN_DEPTH, cancer type as a covariate (not stratifier). Design
~ ctype + age_z + SEX + SV1..k(RUVr) + rna_excess + tp53 + rna_excess:tp53.
Extracts the rna_excess MAIN coefficient + the rna_excess:tp53 INTERACTION; reads out the cytosolic-RP module."""
import os, warnings, numpy as np, pandas as pd
warnings.simplefilter("ignore")
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
R=os.environ.get("FIVES_DATA","data"); RV=f"{R}/results_variants"
MIN_DEPTH=float(os.environ.get("DE_MIN_DEPTH","300")); N_CPUS=6
M=pd.read_pickle(f"{R}/tcga_expr_counts.pkl")
dm=pd.read_csv(f"{RV}/tcga_donor_metrics.tsv",sep="\t").set_index("donor")
cov=pd.read_csv(f"{RV}/tcga_covariates.tsv",sep="\t").set_index("case")
dep=pd.read_csv(f"{RV}/tcga_rna_depth.tsv",sep="\t").set_index("case").depth
lib=M.sum(0).replace(0,1); cpm=lambda e: np.log1p((M.loc[e] if e in M.index else pd.Series(0,index=M.columns))/lib*1e6)
sex=pd.Series(np.where(cpm("ENSG00000129824")+cpm("ENSG00000067048")>cpm("ENSG00000229807"),"M","F"),index=M.columns)
cases=[c for c in M.columns if c in cov.index and c in dm.index and cov.loc[c,"ctype"] is not np.nan
       and cov.loc[c,"tp53"] in("WT","deficient") and pd.notna(dm.loc[c,"rna_excess_z"]) and pd.notna(cov.loc[c,"age_yr"])
       and dep.get(c,0)>=MIN_DEPTH]
print(f"pooled deep-tail: {len(cases)} tumours at >= {MIN_DEPTH}x | mut {sum(cov.loc[c,'tp53']=='deficient' for c in cases)}")
print("per ctype:", dict(pd.Series([cov.loc[c,'ctype'] for c in cases]).value_counts()))
cnt=M[cases].T.values.astype(int); keep=((cnt>=10).mean(0)>=0.5); cnt=cnt[:,keep]; gn=M.index[keep]
ob=pd.DataFrame(index=cases)
ob["age_z"]=(cov.loc[cases,"age_yr"].astype(float)-cov.loc[cases,"age_yr"].astype(float).mean())/cov.loc[cases,"age_yr"].astype(float).std()
ob["SEX"]=sex.reindex(cases).values
ct=cov.loc[cases,"ctype"].astype(str); ct=ct.where(ct.isin(ct.value_counts()[ct.value_counts()>=10].index),"other"); ob["ctype"]=ct.values
ob["tp53"]=pd.Categorical(cov.loc[cases,"tp53"].values,categories=["WT","deficient"]); ob["metric"]=dm.loc[cases,"rna_excess_z"].astype(float).values
# RUVr SVs
X=cnt.astype(float); L=X.sum(1,keepdims=True); L[L==0]=1; Y=np.log1p(X/L*1e4); Y=Y-Y.mean(0)
Dc=[np.ones(len(cases)),ob.age_z.values,ob.metric.values,(ob.tp53=="deficient").astype(float).values,ob.metric.values*(ob.tp53=="deficient").astype(float).values]
for cat in ["SEX","ctype"]:
    for lv in pd.get_dummies(ob[cat],drop_first=True).T.values: Dc.append(lv.astype(float))
D=np.column_stack(Dc); beta,*_=np.linalg.lstsq(D,Y,rcond=None); Rr=Y-D@beta
U,S,_=np.linalg.svd(Rr,full_matrices=False); k=8; SV=(U[:,:k]*S[:k]); SV=(SV-SV.mean(0))/(SV.std(0)+1e-9)
md=ob[["age_z","SEX","ctype","tp53","metric"]].copy()
md["SEX"]=md.SEX.astype("category"); md["ctype"]=md.ctype.astype("category")
for j in range(k): md[f"SV{j+1}"]=SV[:,j]
sv=" + ".join(f"SV{j+1}" for j in range(k))
dds=DeseqDataSet(counts=pd.DataFrame(cnt,index=cases,columns=gn),metadata=md,
                 design=f"~ ctype + age_z + SEX + {sv} + metric + tp53 + metric:tp53",quiet=True,n_cpus=N_CPUS)
dds.deseq2(); cols=list(dds.obsm["design_matrix"].columns)
def coef(sel_fn):
    sel=[c for c in cols if sel_fn(c)]
    if not sel: return None
    vec=np.array([1.0 if c==sel[0] else 0.0 for c in cols]); st=DeseqStats(dds,contrast=vec,quiet=True); st.summary(); return st.results_df
rm=coef(lambda c:c=="metric"); ri=coef(lambda c:"metric" in c and c!="metric" and (":" in c or "tp53" in c))
mod=pd.read_csv(f"{RV}/translation_dosage_module.tsv",sep="\t"); rp=set(mod[mod.is_RP==True].ensg)
for lab,r in [("MAIN 5S->RP",rm),("INTERACTION xTP53",ri)]:
    if r is None: continue
    r=r.copy(); r["e"]=[e.split(".")[0] for e in r.index]; z=r[r.e.isin(rp)]
    print(f"  {lab:18}: RP-module mean z={z.stat.mean():+.2f} ({(z.stat<0).sum()}/{len(z)} down)  median LFC={z.log2FoldChange.median():+.3f}")
