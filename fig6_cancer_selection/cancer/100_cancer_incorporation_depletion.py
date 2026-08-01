# -----------------------------------------------------------------------------
# 100_cancer_incorporation_depletion.py — compare incorporation-related 5S variant carrier-frequency structure across germline (UKBB/HPRC/GTEx) and cancer (CPTAC/TCGA) cohorts.
# Author: Jonathan Bohlen
# (c) 2026 Jonathan Bohlen. Code accompanying Sengl et al. (2026),
# "Surveillance and selection of 5S ribosomal RNA genes in the human genome."
# Released under the MIT License; see LICENSE at the repository root.
# -----------------------------------------------------------------------------
"""Compare incorporation-related carrier-frequency structure of 5S variants between cancer cohorts (TCGA,
CPTAC) and germline cohorts (UKBB, HPRC, GTEx). Within each cohort, test whether variant carrier frequency
scales with incorporation among EXPRESSED variants using denominator-invariant statistics (Spearman + OLS
slope) plus a Mann-Whitney contrast of incorporation-defective vs competent variants. GENE region 633-746,
DEF_CUT 0.50."""
import os
import pandas as pd, numpy as np, warnings; warnings.simplefilter("ignore")
from scipy.stats import spearmanr, mannwhitneyu
import statsmodels.api as sm
T=os.environ.get("FIVES_DATA","data"); OUT=os.environ.get("FIVES_OUT","output"); GENE=(633,746); DEF=0.50
M=pd.read_csv(f"{OUT}/06_population_genetics/97_master_table.tsv",sep="\t")  # germline baseline
M["key"]=M.pos.astype(str)+M.alt
# ---- cancer cohort per-variant carrier frequency (VAF>=0.3%, per individual) ----
cp=pd.read_csv(f"{T}/cptac_wgs_carriers.tsv",sep="\t"); cp=cp[(cp.pos>=GENE[0])&(cp.pos<=GENE[1])]
N_CPTAC=cp.donor.nunique()  # every donor in file is a carrier; denominator-invariant statistic
cpc=cp[cp.dna_vaf>=0.003].groupby(["pos","alt"]).donor.nunique();
tc=pd.read_csv(f"{T}/results_variants/tcga_pervariant.tsv",sep="\t")
tc=tc[(tc.pos>=GENE[0])&(tc.pos<=GENE[1])]; N_TCGA=tc.donor.nunique()
tcc=tc[(tc.carrier==1)|(tc.dna_vaf>=0.003)].groupby(["pos","alt"]).donor.nunique()
print(f"denominators: CPTAC={N_CPTAC} TCGA={N_TCGA} | UKBB=490075 (from master)",flush=True)
M["CPTAC"]=[cpc.get((p,a),0)/N_CPTAC for p,a in zip(M.pos,M.alt)]
M["TCGA"] =[tcc.get((p,a),0)/N_TCGA for p,a in zip(M.pos,M.alt)]
COH=["UKBB","HPRC","GTEx","CPTAC","TCGA"]; GERM=["UKBB","HPRC","GTEx"]; CANC=["CPTAC","TCGA"]
# ---- per-cohort incorporation-depletion stat (among EXPRESSED variants) ----
print(f"\n{'cohort':<8}{'n_expr':>7}{'rho(incorp,freq)':>18}{'p':>10}{'OLS incorp beta':>16}{'p':>10}{'MW p(inc-def<)':>15}")
res={}
for c in COH:
    e=M[M.expressed & M.incorp.notna() & (M[c]>0 if False else True)].copy()  # expressed + incorp measured
    e=e[e.incorp.notna() & M.expressed]
    rho,prho=spearmanr(e.incorp, e[c])
    fl=max(e[c][e[c]>0].min()/2,1e-8) if (e[c]>0).any() else 1e-8
    y=np.log10(e[c]+fl); X=sm.add_constant(e[["incorp","expr"]].astype(float))
    ols=sm.OLS(y,X).fit(); beta,pb=ols.params["incorp"],ols.pvalues["incorp"]
    incdef=e[e.incorp<DEF][c]; comp=e[e.incorp>=DEF][c]
    _,pmw=mannwhitneyu(incdef,comp,alternative="less") if len(incdef)>0 and len(comp)>0 else (np.nan,np.nan)
    res[c]=dict(n=len(e),rho=rho,prho=prho,beta=beta,pb=pb,pmw=pmw,med_incdef=incdef.median(),med_comp=comp.median())
    print(f"{c:<8}{len(e):>7}{rho:>18.3f}{prho:>10.2e}{beta:>16.3f}{pb:>10.2e}{pmw:>15.2e}",flush=True)
print("\n(positive rho / positive OLS beta = incorporation-defective variants at lower frequency)")
print("(compare germline vs cancer cohorts)")
R=pd.DataFrame(res).T; R.to_csv(f"{T}/results_variants/incorp_depletion_5cohort.tsv",sep="\t")
germ_beta=np.mean([res[c]["beta"] for c in GERM]); canc_beta=np.mean([res[c]["beta"] for c in CANC])
print(f"\nMEAN OLS incorp-beta: germline={germ_beta:+.3f} | cancer={canc_beta:+.3f}")
print("wrote incorp_depletion_5cohort.tsv")
