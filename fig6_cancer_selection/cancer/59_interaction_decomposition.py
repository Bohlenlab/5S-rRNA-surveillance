# -----------------------------------------------------------------------------
# 59_interaction_decomposition.py — decompose the per-cancer-type 5S-dose to RP-module slope into TP53-WT and TP53-mutant components (CPTAC).
# Author: Jonathan Bohlen
# (c) 2026 Jonathan Bohlen. Code accompanying Sengl et al. (2026),
# "Surveillance and selection of 5S ribosomal RNA genes in the human genome."
# Released under the MIT License; see LICENSE at the repository root.
# -----------------------------------------------------------------------------
import os
import pandas as pd,numpy as np,glob
R=os.environ.get("FIVES_DATA","data")
mod=pd.read_csv(f"{R}/results_variants/translation_dosage_module.tsv",sep="\t"); rp=set(mod[mod.is_RP==True].ensg)
rows=[]
for f in glob.glob(f"{R}/out_de/de_rna_*.tsv"):
    t=f.split("de_rna_")[1].replace(".tsv","")
    if t=="allctype": continue
    d=pd.read_csv(f,sep="\t"); d["e"]=[e.split(".")[0] for e in d.ensg]; r=d[d.e.isin(rp)].copy()
    r["mut_lfc"]=r.lfc_main+r.lfc_inter               # slope in TP53-mut = WT slope + interaction
    wt=r.lfc_main.mean(); mut=r["mut_lfc"].mean(); inter=r.lfc_inter.mean()
    blunt=(r["mut_lfc"]>r.lfc_main).mean()            # frac RP genes where mut slope is higher (less repression)
    rows.append(dict(type=t,n=int(d.n.iloc[0]),n_mut=int(d.n_mut.iloc[0]),
        WT_slope=round(wt,3),MUT_slope=round(mut,3),interaction=round(inter,3),
        WT_z=round(r.z_main.mean(),2),inter_z=round(r.z_inter.mean(),2),
        pct_genes_blunted=round(100*blunt),
        direction="BLUNTED (mut less RP-down)" if inter>0 else "mut MORE RP-down"))
D=pd.DataFrame(rows).sort_values("WT_slope")
pd.set_option("display.width",220,"display.max_columns",20)
print("=== per-type 5S->RP slope decomposition: WT(p53-active) vs MUT subset ===")
print(D[["type","n","n_mut","WT_slope","MUT_slope","interaction","WT_z","inter_z","pct_genes_blunted","direction"]].to_string(index=False))
print("\nWT_slope = RP-module mean log2FC per SD rna_excess in TP53-WT")
print("MUT_slope = same in TP53-mut (= WT_slope + interaction).  interaction>0 = blunting.")
D.to_csv(f"{R}/tp53/cptac_interaction_decomp.tsv",sep="\t",index=False)
