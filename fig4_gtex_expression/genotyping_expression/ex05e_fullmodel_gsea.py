#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# ex05e_fullmodel_gsea.py — full covariate model for high-expresser donors:
# regress out covariates + hidden factors, donor-level Welch t test, then GSEA.
# Author: Jonathan Bohlen
# (c) 2026 Jonathan Bohlen. Code accompanying Sengl et al. (2026),
# "Surveillance and selection of 5S ribosomal RNA genes in the human genome."
# Released under the MIT License; see LICENSE at the repository root.
# -----------------------------------------------------------------------------
"""Full covariate model comparing high-expresser vs remaining donors.
Sample-level regress-out of: tissue (one-hot) + RIN + ischemic + sex + genotype PC1-5
+ K hidden factors (PCs of the [covariates+high_status]-residualized matrix, so they
cannot absorb the high-expresser effect; PEER-analogous). Then donor-level high-vs-rest
Welch t test (no pseudo-replication) -> GSEA (Hallmark + Reactome)."""
import re,gzip,numpy as np,pandas as pd,scanpy as sc,anndata as ad,gseapy as gp
from scipy import stats
from sklearn.decomposition import TruncatedSVD
import os
from pathlib import Path
ROOT=str(Path(os.environ.get("FIVES_DATA","data"))); OUT=f"{ROOT}/results/eqtl/extreme"
GTF=str(Path(os.environ.get("FIVES_REFS","refs"))/"gencode.v49.annotation.gtf.gz")
import sqlite3,os
K_HF=30; THR=0.015
EXCL=os.environ.get("EXCL_CULTURED"); ADJ=os.environ.get("ADJ_COMP")
SUF="_compadj" if ADJ else ("_tissueonly" if EXCL else "")
CULT={"Cells - Cultured fibroblasts","Cells - EBV-transformed lymphocytes","Cells - Leukemia cell line (CML)"}
EPI=["CDH1","GRHL2","ESRP1","ESRP2","EHF","ELF3","DSP","PKP3","PKP2","TJP3","RAB25","ST14","SPINT1","SPINT2","PROM2","LLGL2","AP1M2","CD2AP","TACSTD2","EPCAM","KRT8","KRT18","KRT19","CLDN3","CLDN4","CLDN7","LAD1","PIGR","MAL2","OCLN","CGN","MARVELD2"]
MES=["MYL9","MYH11","TAGLN","CNN1","ACTN1","DES","CSRP1","PPP1R14A","DSTN","ACTA2","MYLK","TPM1","TPM2","LMOD1","PALLD","SYNM","VIM","FN1","COL1A1","COL1A2","COL3A1","COL6A1","COL6A2","COL6A3","SPARC","LUM","DCN"]
meta=pd.read_csv(f"{ROOT}/metadata/GTEx_Analysis_2025-08-22_v11_Annotations_SampleAttributesDS.txt",
   sep="\t",usecols=["SAMPID","SMTSD","SMRIN","SMTSISCH"],low_memory=False).set_index("SAMPID")
e2s={}
for ln in gzip.open(GTF,"rt"):
    if ln[0]=="#" or "\tgene\t" not in ln: continue
    e2s[re.search(r'gene_id "([^"]+)"',ln).group(1).split(".")[0]]=re.search(r'gene_name "([^"]+)"',ln).group(1)
db=sqlite3.connect(str(Path(os.environ.get("FIVES_DB","5S_rDNA.db"))))
sexd=pd.read_sql("SELECT sample_id,sex FROM assembly WHERE cohort='GTEx_v9_WGS'",db).set_index("sample_id").sex
pc=pd.read_csv(f"{ROOT}/results/eqtl_inputs/genotype_20PCs.eigenvec.txt",sep="\t")
pc["donor"]=pc.IID.str.extract(r'(GTEX-[A-Z0-9]+)'); pc=pc.drop_duplicates("donor").set_index("donor")[[f"PC{i}" for i in range(1,6)]]
P=pd.read_csv(f"{OUT}/donor_variant_rnavaf.tsv",sep="\t"); GEN=["687G","701C","725C","726T","730G","733T","734T","743G"]
mv=P[P.variant.isin(GEN)].groupby("donor").rna_vaf.max(); high=set(mv[mv>=THR].index)

A=ad.read_h5ad(os.environ.get("FIVES_GTEX_COUNTS","GTEx_v11_bulk_gene_counts.h5ad"))
A.var["ensg"]=[v.split(".")[0] for v in A.var_names]; A=A[:,~A.var.ensg.duplicated()].copy(); A.var_names=A.var.ensg.values
A.obs=A.obs.join(meta,how="left"); A=A[A.obs.SMTSD.notna()].copy()
if EXCL: A=A[~A.obs.SMTSD.isin(CULT)].copy(); print(f"EXCLUDED cultured cells -> {A.n_obs} samples")
sc.pp.normalize_total(A,target_sum=1e4); sc.pp.log1p(A)
X=A.X.tocsc(); mean=np.asarray(X.mean(0)).ravel(); det=np.asarray((X>0).mean(0)).ravel()
keep=(mean>0.1)&(det>0.2)&np.array([v in e2s for v in A.var_names]); A=A[:,keep].copy()
ob=A.obs.copy(); ob["donor"]=ob.donor.astype(str)
ob["sex"]=ob.donor.map(sexd).map({"M":0,"F":1}); ob=ob.join(pc,on="donor")
ok=ob[["SMRIN","SMTSISCH","sex"]+[f"PC{i}" for i in range(1,6)]].notna().all(1).values
A=A[ok].copy(); ob=ob[ok]; print(f"samples after covariate filter: {A.n_obs}; genes {A.n_vars}")
Y=A.X.toarray().astype(np.float32)
# known covariate design C (sample level)
tis=pd.get_dummies(ob.SMTSD,drop_first=True).values.astype(np.float32)
base=ob[["SMRIN","SMTSISCH","sex","PC1","PC2","PC3","PC4","PC5"]].values.astype(np.float32)
if ADJ:   # add epithelial + mesenchymal composition proxies (from log-norm expression) as covariates
    _symk=np.array([e2s.get(v,"") for v in A.var_names])
    _ms=lambda gl: np.asarray(Y[:,np.isin(_symk,gl)].mean(1)).ravel()
    base=np.c_[base,_ms(EPI)[:,None],_ms(MES)[:,None]]; print(f"ADJ: added epi+mes composition covariates")
num=(base-base.mean(0))/(base.std(0)+1e-9)
C=np.c_[np.ones(len(ob),np.float32),tis,num]
h=ob.donor.isin(high).values.astype(np.float32)[:,None]
def resid(D,Y):
    B,*_=np.linalg.lstsq(D,Y,rcond=None); return Y-D@B
# protected hidden factors: PCs of [C,h]-residualized Y
Rh=resid(np.c_[C,h],Y)
svd=TruncatedSVD(n_components=K_HF,random_state=0); HF=svd.fit_transform(Rh).astype(np.float32); del Rh
print(f"hidden factors: {K_HF} (protected, orthogonal to high-status)")
# regress out known covariates + hidden factors (NOT h); donor-level test
Yr=resid(np.c_[C,HF],Y); del Y
DG=pd.DataFrame(Yr,columns=[e2s[v] for v in A.var_names]); DG["donor"]=ob.donor.values
# sample-level module residual scores (full-model residuals incl. hidden factors) for per-tissue analysis
_cyto=lambda s:bool(re.match(r'^RP[LS]\d+[A-Z]?$',str(s))) or s in {"RPLP0","RPLP1","RPLP2","RPSA","FAU","UBA52"}
_rpc=[c for c in DG.columns if c!="donor" and _cyto(c)]
EPI=["CDH1","GRHL2","ESRP1","ESRP2","EHF","ELF3","DSP","PKP3","PKP2","TJP3","RAB25","ST14","SPINT1","SPINT2","PROM2","LLGL2","AP1M2","CD2AP","TACSTD2","EPCAM","KRT8","KRT18","KRT19","CLDN3","CLDN4","CLDN7","LAD1","PIGR","MAL2","OCLN","CGN","MARVELD2"]
MES=["MYL9","MYH11","TAGLN","CNN1","ACTN1","DES","CSRP1","PPP1R14A","DSTN","ACTA2","MYLK","TPM1","TPM2","LMOD1","PALLD","SYNM","VIM","FN1","COL1A1","COL1A2","COL3A1","COL6A1","COL6A2","COL6A3","SPARC","LUM","DCN"]
_epi=[g for g in EPI if g in DG.columns]; _mes=[g for g in MES if g in DG.columns]
pd.DataFrame({"donor":ob.donor.values,"tissue":ob.SMTSD.values,"high":ob.donor.isin(high).values,
    "rp_resid":DG[_rpc].mean(1).values,"epi_resid":DG[_epi].mean(1).values,"mes_resid":DG[_mes].mean(1).values
    }).to_csv(f"{OUT}/fullmodel_sample_rpscore{SUF}.tsv",sep="\t",index=False)
print(f"module sizes: RP={len(_rpc)} EPI={len(_epi)} MES={len(_mes)}")
# comprehensive sample x pathway module-score table (Hallmark + custom) for the per-tissue pathway volcano
_HSR=["HSPA1A","HSPA1B","HSPA6","HSPA8","HSPA4","HSPH1","HSP90AA1","HSP90AB1","DNAJB1","DNAJA1","DNAJA4","DNAJB6","BAG3","HSPB1","HSF1","AHSA1","FKBP4","STIP1","CHORDC1","SERPINH1","HSPD1","HSPE1"]
try: _lib=gp.get_library("MSigDB_Hallmark_2020")
except Exception as _ex: _lib={}; print("hallmark fetch failed:",_ex)
_lib={"HALLMARK_"+k.replace(" ","_"):v for k,v in _lib.items()}
_lib.update({"_TRANSLATION_RP":_rpc,"_EPITHELIAL":EPI,"_MESENCHYMAL":MES,"_HEAT_SHOCK":_HSR})
_ps=pd.DataFrame({"donor":ob.donor.values,"tissue":ob.SMTSD.values,"high":ob.donor.isin(high).values})
for _nm,_gl in _lib.items():
    _g=[x for x in _gl if x in DG.columns]
    if len(_g)>=8: _ps[_nm]=DG[_g].mean(1).values
_ps.to_csv(f"{OUT}/fullmodel_sample_pathscores{SUF}.tsv",sep="\t",index=False)
print(f"saved {_ps.shape[1]-3} pathway scores x {len(_ps)} samples")
DG=DG.groupby("donor").mean(); DG=DG.loc[:,~DG.columns.duplicated()]
isx=DG.index.isin(high); print(f"donor-level: {isx.sum()} high vs {(~isx).sum()} rest")
t,p=stats.ttest_ind(DG[isx].values,DG[~isx].values,axis=0,equal_var=False,nan_policy="omit")
G=pd.DataFrame({"gene":DG.columns,"t":t,"p":p,"eff":DG[isx].mean()-DG[~isx].mean()}).dropna().sort_values("t",ascending=False)
o=np.argsort(G.p.values); q=np.empty(len(G)); q[o]=np.minimum.accumulate((G.p.values[o]*len(G)/np.arange(1,len(G)+1))[::-1])[::-1]
G["fdr"]=np.clip(q,0,1); G.to_csv(f"{OUT}/fullmodel_highexpr_DE{SUF}.tsv",sep="\t",index=False)
print(f"genes FDR<0.1: {(G.fdr<0.1).sum()}")
rnk=G[["gene","t"]].drop_duplicates("gene")
r=gp.prerank(rnk=rnk,gene_sets=["MSigDB_Hallmark_2020","Reactome_2022"],min_size=10,max_size=500,
    permutation_num=1000,seed=7,threads=8,no_plot=True,outdir=None,verbose=False)
res=r.res2d.copy(); fc=[c for c in res.columns if "FDR" in c][0]
res["NES"]=res.NES.astype(float); res[fc]=res[fc].astype(float); res.to_csv(f"{OUT}/fullmodel_highexpr_gsea{SUF}.tsv",sep="\t",index=False)
print("\n=== FULL MODEL: top 12 DOWN pathways ===")
for _,x in res.sort_values("NES").head(12).iterrows(): print(f"  NES={x.NES:+.2f} FDR={x[fc]:.3f}  {x.Term.split('__')[-1]}")
dn=res.sort_values("NES").reset_index(drop=True)
key=dn[dn.Term.str.contains('ranslat|ibosom|60S|40S|Peptide chain',case=False)]
print(f"\ntranslation/ribosome best rank: {key.index.min()+1 if len(key) else 'NONE'} of {len(dn)}; "
      f"top transl NES={key.NES.min() if len(key) else float('nan'):+.2f} FDR={key[fc].min() if len(key) else float('nan'):.3f}")
