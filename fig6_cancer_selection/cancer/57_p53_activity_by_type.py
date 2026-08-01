# -----------------------------------------------------------------------------
# 57_p53_activity_by_type.py — per-cancer-type p53 activation/p21/proliferation scores and their correlation with the RP-module DE signal (CPTAC).
# Author: Jonathan Bohlen
# (c) 2026 Jonathan Bohlen. Code accompanying Sengl et al. (2026),
# "Surveillance and selection of 5S ribosomal RNA genes in the human genome."
# Released under the MIT License; see LICENSE at the repository root.
# -----------------------------------------------------------------------------
import os
import pandas as pd,numpy as np,glob
from scipy import stats
R=os.environ.get("FIVES_DATA","data")
# MDM2-excluded p53 activation target genes (MDM2 dropped to avoid amplicon confound)
P53={"CDKN1A":"ENSG00000124762","GADD45A":"ENSG00000116717","BBC3":"ENSG00000105327","BAX":"ENSG00000087088",
"SESN1":"ENSG00000080546","SESN2":"ENSG00000130766","RRM2B":"ENSG00000048392","ZMAT3":"ENSG00000172667",
"BTG2":"ENSG00000159388","FAS":"ENSG00000026103","TNFRSF10B":"ENSG00000120889","FDXR":"ENSG00000161513",
"TP53I3":"ENSG00000115129","DDB2":"ENSG00000134574","AEN":"ENSG00000181026","PHLDA3":"ENSG00000174307",
"GDF15":"ENSG00000130513","EDA2R":"ENSG00000131080"}
CDKN1A="ENSG00000124762"
PROLIF={"MKI67":"ENSG00000148773","PCNA":"ENSG00000132646","TOP2A":"ENSG00000131747","MCM2":"ENSG00000073111",
"E2F1":"ENSG00000101412","BUB1":"ENSG00000169679","CCNB1":"ENSG00000134057","MYC":"ENSG00000136997"}
M=pd.read_pickle(f"{R}/cptac_expr_counts.pkl"); M.index=[e.split(".")[0] for e in M.index]
cov=pd.read_csv(f"{R}/results_variants/cptac_covariates.tsv",sep="\t").set_index("case")
dm=set(pd.read_csv(f"{R}/results_variants/cptac_donor_metrics.tsv",sep="\t").donor)
lib=M.sum(0).replace(0,1); logcpm=np.log1p(M/lib*1e6)
mod=pd.read_csv(f"{R}/results_variants/translation_dosage_module.tsv",sep="\t"); rp=set(mod[mod.is_RP==True].ensg)
rpz={}
for f in glob.glob(f"{R}/out_de/de_rna_*.tsv"):
    t=f.split("de_rna_")[1].replace(".tsv","")
    if t=="allctype": continue
    d=pd.read_csv(f,sep="\t"); d["e"]=[e.split(".")[0] for e in d.ensg]; rpz[t]=d[d.e.isin(rp)].z_main.mean()
def score(genes,cases):
    g=[e for e in genes if e in M.index]; s=logcpm.loc[g,cases]
    return s.sub(s.mean(1),axis=0).div(s.std(1)+1e-9,axis=0).mean(0)
rows=[]
for t in rpz:
    cases=[c for c in M.columns if c in cov.index and cov.loc[c,"ctype"]==t and cov.loc[c,"tp53"] in("WT","deficient") and c in dm]
    if len(cases)<25: continue
    st=cov.loc[cases,"tp53"].values
    p53=score(P53.values(),cases); p21=score([CDKN1A],cases); prol=score(PROLIF.values(),cases)
    rows.append(dict(type=t,rp_z=rpz[t],p53_WT=p53[st=="WT"].mean(),p53gap=p53[st=="WT"].mean()-p53[st=="deficient"].mean(),
                     p21_WT=p21[st=="WT"].mean(),p21gap=p21[st=="WT"].mean()-p21[st=="deficient"].mean(),
                     prolif=prol.mean()))
D=pd.DataFrame(rows).sort_values("rp_z")
print("=== MDM2-clean p53 activation + p21-alone + proliferation, per type ===")
print(D.round(3).to_string(index=False))
D2=D[D.type!="CCRCC"]
for lab,col in [("p53gap(clean)","p53gap"),("p53_WT(clean)","p53_WT"),("p21gap","p21gap"),("p21_WT","p21_WT"),("prolif","prolif")]:
    r,p=stats.spearmanr(D2.rp_z,D2[col]); print(f"Spearman(RP z, {lab:16}) [excl CCRCC, n={len(D2)}]: rho={r:+.2f} p={p:.3f}")
