#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# 66_expanded_p53_classification.py — classify CPTAC tumours as p53-inactive by TP53 mutation or non-mutational p53-pathway lesions (MDM2/MDM4 amp, CDKN2A/ARF del, HPV).
# Author: Jonathan Bohlen
# (c) 2026 Jonathan Bohlen. Code accompanying Sengl et al. (2026),
# "Surveillance and selection of 5S ribosomal RNA genes in the human genome."
# Released under the MIT License; see LICENSE at the repository root.
# -----------------------------------------------------------------------------
"""Expanded p53-inactive classification for CPTAC. p53_inactive = TP53-mutant OR MDM2-amp OR
MDM4-amp OR CDKN2A/ARF-deletion OR HPV+ (genetic/molecular). Reclassifies genetically-WT tumours carrying a
non-mutational p53-pathway lesion into the 'inactive' group. HPV not available for CPTAC via cBioPortal ->
flagged NA. Inputs: covariates (tp53) + cptac_cna_p53.tsv (MDM2/MDM4/CDKN2A GISTIC).
Output: surveillance_v2/tables/01_p53_status_expanded.tsv"""
import os, pandas as pd
ROOT=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT=os.path.join(ROOT,"surveillance_v2","tables"); os.makedirs(OUT,exist_ok=True)
cov=pd.read_csv(os.path.join(ROOT,"results_variants","cptac_covariates.tsv"),sep="\t")
cna=pd.read_csv(os.path.join(ROOT,"tp53","cptac_cna_p53.tsv"),sep="\t")[["case","mdm2_amp","mdm4_amp","cdkn2a_del"]]
m=cov[["case","ctype","tp53"]].merge(cna,on="case",how="left")
m=m[m.tp53.isin(["WT","deficient"])].copy()
m["tp53_mut"]=(m.tp53=="deficient").astype(int)
for c in ["mdm2_amp","mdm4_amp","cdkn2a_del"]: m[c]=m[c].fillna(0).astype(int)
m["hpv_pos"]=pd.NA   # not accessible via cBioPortal for CPTAC HNSC/CESC
m["cna_available"]=cov.set_index("case").reindex(m.case).index.isin(cna.case).astype(int) if False else m.mdm2_amp.notna().astype(int)
m["p53_inactive"]=((m.tp53_mut==1)|(m.mdm2_amp==1)|(m.mdm4_amp==1)|(m.cdkn2a_del==1)).astype(int)
def mech(r):
    ms=[k for k,v in [("TP53mut",r.tp53_mut),("MDM2amp",r.mdm2_amp),("MDM4amp",r.mdm4_amp),("CDKN2Adel",r.cdkn2a_del)] if v]
    return ";".join(ms) if ms else "p53_competent"
m["mechanism"]=m.apply(mech,axis=1)
# flag squamous/HPV-driven types where HPV status is unavailable
m["hpv_relevant"]=m.ctype.isin(["HNSC","LSCC","CESC","HNSCC","CSCC"]).astype(int)
m.to_csv(os.path.join(OUT,"01_p53_status_expanded.tsv"),sep="\t",index=False)
# report
tot=len(m)
print(f"CPTAC cases classified: {tot} | p53-inactive (expanded): {int(m.p53_inactive.sum())} ({100*m.p53_inactive.mean():.0f}%) | mut-only baseline: {int(m.tp53_mut.sum())}")
recl=m[(m.tp53_mut==0)&(m.p53_inactive==1)]
print(f"WT tumours RECLASSIFIED to inactive (non-mutational lesion): {len(recl)} of {int((m.tp53_mut==0).sum())} WT")
print("\nper type: WT_orig -> competent_new (reclassified out) | inactive_total")
rows=[]
for t,g in m.groupby("ctype"):
    wt=int((g.tp53_mut==0).sum()); comp=int(((g.tp53_mut==0)&(g.p53_inactive==0)).sum()); reclt=wt-comp
    rows.append((t,len(g),int(g.tp53_mut.sum()),wt,reclt,comp,int(g.p53_inactive.sum())))
rep=pd.DataFrame(rows,columns=["ctype","n","tp53mut","WT_orig","reclassified","competent_new","inactive_total"]).sort_values("reclassified",ascending=False)
print(rep.to_string(index=False))
rep.to_csv(os.path.join(OUT,"01_reclassification_by_type.tsv"),sep="\t",index=False)
