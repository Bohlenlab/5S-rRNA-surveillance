# -----------------------------------------------------------------------------
# 102_somatic_calls_incorporation.py — incorporation-defective fraction of significant somatic 5S gains vs germline background and by p53 status (CPTAC, TCGA).
# Author: Jonathan Bohlen
# (c) 2026 Jonathan Bohlen. Code accompanying Sengl et al. (2026),
# "Surveillance and selection of 5S ribosomal RNA genes in the human genome."
# Released under the MIT License; see LICENSE at the repository root.
# -----------------------------------------------------------------------------
"""Incorporation analysis on matched somatic calls (matched Fisher + beta-binomial, tumour vs matched
blood-normal). Events analysed = significant somatic GAINS of variant 5S. Two tests:
(A) Fisher test of whether somatic gains differ from the germline background in incorporation-defective
    (incorp<0.5) fraction.
(B) p53-STRATIFIED: compares the incorporation-defective fraction of somatic gains between p53-deficient and
    p53-WT tumours; per-tumour Poisson model with a cancer-type covariate and an offset for total gains.
    CPTAC + TCGA, pooled + per-cohort."""
import os
import pandas as pd, numpy as np, warnings; warnings.simplefilter("ignore")
from scipy.stats import fisher_exact, mannwhitneyu
import statsmodels.api as sm, statsmodels.formula.api as smf
T=os.environ.get("FIVES_DATA","data"); DEFCUT=0.50
fa=pd.read_csv(f"{T}/results_variants/funcannot_ei.tsv",sep="\t"); fa["key"]=fa.pos.astype(str)+fa.alt
incorp={k:v for k,v in zip(fa.key,fa.incorp_60s_mean) if pd.notna(v)}
S=pd.read_csv(f"{T}/results_variants/somatic_calls_annotated.tsv",sep="\t")
S["key"]=S.pos.astype(str)+S.alt; S["incorp"]=S.key.map(incorp); S["idef"]=S.incorp<DEFCUT
sg=S[(S.sig==True)&(S.direction=="gain")&S.incorp.notna()].copy()   # significant somatic GAINS, annotated
print(f"significant somatic gains with incorp annotation: {len(sg)} (CPTAC {int((sg.cohort=='CPTAC').sum())}, TCGA {int((sg.cohort=='TCGA').sum())})")
# germline background incorp-def rate (matched CPTAC blood-normal carriers)
gm=pd.read_csv(f"{T}/cptac_germline_carriers.tsv",sep="\t"); gm["key"]=gm.pos.astype(str)+gm.alt; gm["incorp"]=gm.key.map(incorp); gm=gm[gm.incorp.notna()]
germ_def=100*(gm.incorp<DEFCUT).mean()
print(f"\n=== (A) somatic GAINS vs germline background: incorp-defective enrichment ===")
print(f"  germline carriers: {germ_def:.1f}% incorp-defective (mean incorp {gm.incorp.mean():.2f}, n={len(gm)})")
for coh,d in [("ALL",sg),("CPTAC",sg[sg.cohort=='CPTAC']),("TCGA",sg[sg.cohort=='TCGA'])]:
    a=(d.incorp<DEFCUT).sum(); b=(d.incorp>=DEFCUT).sum(); c=(gm.incorp<DEFCUT).sum(); e=(gm.incorp>=DEFCUT).sum()
    orr,pf=fisher_exact([[a,b],[c,e]])
    print(f"  {coh:<6} gains: {100*a/(a+b):.1f}% incorp-def (mean {d.incorp.mean():.2f}), OR vs germline={orr:.2f} p={pf:.3f}  {'ENRICHED' if orr>1 else 'depleted'}")
print("  (OR>1 = gains higher incorp-defective fraction than germline; OR<1 = lower)\n")

print("=== (B) p53-stratified: do p53-def tumours GAIN incorp-defective variants more than WT? ===")
for coh in ["ALL","CPTAC","TCGA"]:
    d=sg if coh=="ALL" else sg[sg.cohort==coh]; d=d[d.tp53.isin(["WT","deficient"])]
    # event-level: fraction of somatic gains that are incorp-def, WT-tumour vs def-tumour
    wt=d[d.tp53=="WT"]; de=d[d.tp53=="deficient"]
    a=(de.incorp<DEFCUT).sum(); b=(de.incorp>=DEFCUT).sum(); c=(wt.incorp<DEFCUT).sum(); e=(wt.incorp>=DEFCUT).sum()
    orr,pf=fisher_exact([[a,b],[c,e]],alternative="greater")
    print(f"  {coh:<6}: gains in def-tumours {100*a/(a+b):.1f}% incorp-def (n={len(de)}) vs WT-tumours {100*c/(c+e):.1f}% (n={len(wt)}) | Fisher OR(def/WT)={orr:.2f} p(def>WT)={pf:.3f}")
    # per-tumour: n incorp-def gains ~ p53 + offset(log total gains) + ctype
    g=d.groupby("case").agg(n_tot=("idef","size"),n_def=("idef","sum"),tp53=("tp53","first"),ctype=("ctype","first")).reset_index()
    g["is_def"]=(g.tp53=="deficient").astype(int)
    if coh=="ALL" and g.ctype.nunique()>1:
        try:
            m=smf.glm("n_def ~ is_def + C(ctype)",data=g[g.n_tot>0],family=sm.families.Poisson(),offset=np.log(g[g.n_tot>0].n_tot)).fit()
            print(f"         per-tumour Poisson n_def~p53def+ctype+offset(log gains): beta={m.params['is_def']:+.3f} p={m.pvalues['is_def']:.3f} RR={np.exp(m.params['is_def']):.2f}")
        except Exception as ex: print("   poisson failed",ex)
print("\n  (OR/RR>1 = p53-deficient tumours gain incorp-defective variants at higher rate)")
sg.to_csv(f"{T}/results_variants/somatic_gains_incorp_annotated.tsv",sep="\t",index=False)
print("wrote somatic_gains_incorp_annotated.tsv")
