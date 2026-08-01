# -----------------------------------------------------------------------------
# 101_somatic_p53_incorporation.py — somatic-vs-germline and p53-stratified incorporation-defective 5S variant burden in CPTAC and TCGA tumours.
# Author: Jonathan Bohlen
# (c) 2026 Jonathan Bohlen. Code accompanying Sengl et al. (2026),
# "Surveillance and selection of 5S ribosomal RNA genes in the human genome."
# Released under the MIT License; see LICENSE at the repository root.
# -----------------------------------------------------------------------------
"""SOMATIC and p53-stratified incorporation analysis.
(A) SOMATIC-ONLY (CPTAC, matched blood-normal): somatic variant = tumour carrier (VAF>=0.3%) that is NOT a
    germline carrier in the SAME donor. Tests whether somatic (tumour-acquired) variants differ from germline
    variants in incorporation-defective (incorp<0.5) fraction (Mann-Whitney + Fisher).
(B) p53-STRATIFIED: compares the incorporation-defective somatic burden between p53-deficient and p53-WT
    tumours. CPTAC (somatic) + TCGA (tumour carriers). Per-donor Poisson model with a cancer-type covariate
    and an offset for total burden."""
import os
import pandas as pd, numpy as np, warnings; warnings.simplefilter("ignore")
from scipy.stats import mannwhitneyu, fisher_exact
import statsmodels.api as sm, statsmodels.formula.api as smf
T=os.environ.get("FIVES_DATA","data"); GENE=(633,746); DEFCUT=0.50; VC=0.003
fa=pd.read_csv(f"{T}/results_variants/funcannot_ei.tsv",sep="\t"); fa["key"]=fa.pos.astype(str)+fa.alt
incorp={k:v for k,v in zip(fa.key,fa.incorp_60s_mean) if pd.notna(v)}
def ann(df): df=df[(df.pos>=GENE[0])&(df.pos<=GENE[1])].copy(); df["key"]=df.pos.astype(str)+df.alt; df["incorp"]=df.key.map(incorp); return df
tum=ann(pd.read_csv(f"{T}/cptac_wgs_carriers.tsv",sep="\t"))
germ=ann(pd.read_csv(f"{T}/cptac_germline_carriers.tsv",sep="\t"))
cov=pd.read_csv(f"{T}/results_variants/cptac_covariates_expanded.tsv",sep="\t").set_index("case")
matched=sorted(set(tum.donor)&set(germ.donor))
gset={(d,k) for d,k in zip(germ.donor,germ.key)}          # germline carrier (donor,variant)
tm=tum[tum.donor.isin(matched) & (tum.dna_vaf>=VC)].copy()
tm["somatic"]=[(d,k) not in gset for d,k in zip(tm.donor,tm.key)]
tm=tm[tm.incorp.notna()]
som=tm[tm.somatic]; ger=tm[~tm.somatic]  # tumour events split into somatic (de novo) vs germline-shared
print("=== (A) SOMATIC vs GERMLINE incorporation (CPTAC, matched, event-level) ===")
for lab,s in [("somatic (de novo)",som),("germline-shared",ger)]:
    print(f"  {lab:<18} n={len(s):>5}  %incorp-def(<0.5)={100*(s.incorp<DEFCUT).mean():.1f}  mean incorp={s.incorp.mean():.2f}")
u,p=mannwhitneyu(som.incorp,ger.incorp,alternative="less")  # somatic incorp LOWER (more defective)?
a=(som.incorp<DEFCUT).sum(); b=(som.incorp>=DEFCUT).sum(); c=(ger.incorp<DEFCUT).sum(); dd=(ger.incorp>=DEFCUT).sum()
orr,pf=fisher_exact([[a,b],[c,dd]])
print(f"  MannWhitney somatic<germline incorp p={p:.3f} | Fisher incorp-def OR(somatic/germline)={orr:.2f} p={pf:.3f}")
print("  (OR>1 = somatic higher incorp-defective fraction; OR<1 = lower)\n")

def p53burden(df,cohortlabel,tp53map,covdf):
    """per-donor counts of incorp-def vs well-incorp variants, then p53-def vs WT."""
    df=df[df.incorp.notna()].copy(); df["idef"]=df.incorp<DEFCUT
    g=df.groupby("donor").agg(n_tot=("idef","size"),n_def=("idef","sum")).reset_index()
    g["tp53"]=g.donor.map(tp53map); g=g[g.tp53.isin(["WT","deficient"])]
    g["frac_def"]=g.n_def/g.n_tot; g["ctype"]=g.donor.map(lambda d: covdf.loc[d,"ctype"] if d in covdf.index else np.nan)
    wt=g[g.tp53=="WT"]; de=g[g.tp53=="deficient"]
    print(f"=== (B) p53-stratified incorp-def burden — {cohortlabel} (n WT={len(wt)} def={len(de)}) ===")
    print(f"  mean n_incorp-def/tumour: WT={wt.n_def.mean():.2f} def={de.n_def.mean():.2f} | frac-def: WT={wt.frac_def.mean():.3f} def={de.frac_def.mean():.3f}")
    _,pn=mannwhitneyu(de.n_def,wt.n_def,alternative="greater")   # def MORE incorp-def variants?
    _,pfr=mannwhitneyu(de.frac_def.dropna(),wt.frac_def.dropna(),alternative="greater")
    print(f"  MannWhitney def>WT: n_incorp-def p={pn:.3f} | fraction-incorp-def p={pfr:.3f}")
    # negative-binomial: n_def ~ p53 + offset(log n_tot) + ctype  (enrichment of incorp-def among the load)
    gg=g[g.n_tot>0].dropna(subset=["ctype"]).copy(); gg["is_def"]=(gg.tp53=="deficient").astype(int)
    try:
        m=smf.glm("n_def ~ is_def + C(ctype)",data=gg,family=sm.families.Poisson(),offset=np.log(gg.n_tot)).fit()
        print(f"  Poisson n_def ~ p53def + offset(log n_tot) + ctype:  p53def beta={m.params['is_def']:+.3f} p={m.pvalues['is_def']:.3f} (RR={np.exp(m.params['is_def']):.2f})\n")
    except Exception as e: print("  poisson failed",e,"\n")
    return g

tp53_cptac=cov.tp53.to_dict()
p53burden(som,"CPTAC SOMATIC",tp53_cptac,cov)           # somatic-only
p53burden(tm,"CPTAC all-tumour",tp53_cptac,cov)          # all tumour (germline balanced by p53)
# TCGA (tumour carriers; no matched normal -> all-tumour, germline balanced across p53)
tct=ann(pd.read_csv(f"{T}/tcga_tumour_carriers.tsv",sep="\t")) if os.path.exists(f"{T}/tcga_tumour_carriers.tsv") else None
covT=pd.read_csv(f"{T}/results_variants/tcga_covariates.tsv",sep="\t").set_index("case")
if tct is not None:
    if "dna_vaf" in tct.columns: tct=tct[tct.dna_vaf>=VC]
    p53burden(tct,"TCGA all-tumour",covT.tp53.to_dict(),covT)
print("done")
