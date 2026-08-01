# -----------------------------------------------------------------------------
# 58_comprehensive_p53_table.py — per-cancer-type table of TP53 status, non-genetic p53-pathway lesions, p53-activity gap, and RP-module DE statistics (CPTAC).
# Author: Jonathan Bohlen
# (c) 2026 Jonathan Bohlen. Code accompanying Sengl et al. (2026),
# "Surveillance and selection of 5S ribosomal RNA genes in the human genome."
# Released under the MIT License; see LICENSE at the repository root.
# -----------------------------------------------------------------------------
import os
import pandas as pd,numpy as np,glob
from scipy import stats
R=os.environ.get("FIVES_DATA","data")
P53={"CDKN1A":"ENSG00000124762","GADD45A":"ENSG00000116717","BBC3":"ENSG00000105327","BAX":"ENSG00000087088",
"SESN1":"ENSG00000080546","SESN2":"ENSG00000130766","RRM2B":"ENSG00000048392","ZMAT3":"ENSG00000172667",
"BTG2":"ENSG00000159388","FAS":"ENSG00000026103","TNFRSF10B":"ENSG00000120889","FDXR":"ENSG00000161513",
"TP53I3":"ENSG00000115129","DDB2":"ENSG00000134574","AEN":"ENSG00000181026","PHLDA3":"ENSG00000174307",
"GDF15":"ENSG00000130513","EDA2R":"ENSG00000131080"}
M=pd.read_pickle(f"{R}/cptac_expr_counts.pkl"); M.index=[e.split(".")[0] for e in M.index]
cov=pd.read_csv(f"{R}/results_variants/cptac_covariates.tsv",sep="\t").set_index("case")
dm=set(pd.read_csv(f"{R}/results_variants/cptac_donor_metrics.tsv",sep="\t").donor)
cna=pd.read_csv(f"{R}/tp53/cptac_cna_p53.tsv",sep="\t").set_index("case")
lib=M.sum(0).replace(0,1); logcpm=np.log1p(M/lib*1e6)
mod=pd.read_csv(f"{R}/results_variants/translation_dosage_module.tsv",sep="\t"); rp=set(mod[mod.is_RP==True].ensg)
rpz={};iz={}
for f in glob.glob(f"{R}/out_de/de_rna_*.tsv"):
    t=f.split("de_rna_")[1].replace(".tsv","")
    if t=="allctype": continue
    d=pd.read_csv(f,sep="\t"); d["e"]=[e.split(".")[0] for e in d.ensg]; r=d[d.e.isin(rp)]
    rpz[t]=r.z_main.mean(); iz[t]=r.z_inter.mean() if "z_inter" in r else np.nan
tg=[e for e in P53.values() if e in M.index]
rows=[]
for t in rpz:
    cs=[c for c in M.columns if c in cov.index and cov.loc[c,"ctype"]==t and cov.loc[c,"tp53"] in("WT","deficient") and c in dm]
    if len(cs)<25: continue
    st=cov.loc[cs,"tp53"]; wtc=[c for c in cs if st[c]=="WT"]; mutc=[c for c in cs if st[c]=="deficient"]
    sub=logcpm.loc[tg,cs]; z=sub.sub(sub.mean(1),axis=0).div(sub.std(1)+1e-9,axis=0).mean(0)
    gap=z[wtc].mean()-z[mutc].mean(); gp=stats.mannwhitneyu(z[wtc],z[mutc]).pvalue if len(mutc)>=5 else np.nan
    # among WT: molecular non-genetic p53 inactivation
    wcna=cna.reindex(wtc)
    m2=wcna.mdm2_amp.fillna(0); m4=wcna.mdm4_amp.fillna(0); cd=wcna.cdkn2a_del.fillna(0)
    anyinact=((m2+m4+cd)>0).mean()
    rows.append(dict(type=t,n=len(cs),pct_TP53mut=round(100*len(mutc)/len(cs)),
        WT_MDM2amp=round(100*m2.mean()),WT_MDM4amp=round(100*m4.mean()),WT_CDKN2Adel=round(100*cd.mean()),
        WT_anyNonGen=round(100*anyinact),p53act_gap=round(gap,2),gap_p=gp,
        RP_main_z=round(rpz[t],2),TP53_inter_z=round(iz[t],2)))
D=pd.DataFrame(rows).sort_values("RP_main_z")
D.to_csv(f"{R}/tp53/cptac_p53_comprehensive.tsv",sep="\t",index=False)
pd.set_option("display.width",200,"display.max_columns",20)
print(D.to_string(index=False))
print("\nWT_anyNonGen = % of TP53-WT tumours with MDM2amp/MDM4amp/CDKN2Adel (non-genetic p53 inactivation; HPV not included)")
print("p53act_gap = WT minus mut p53-activity")
