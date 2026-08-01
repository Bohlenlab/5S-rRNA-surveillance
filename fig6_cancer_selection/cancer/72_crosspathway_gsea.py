# -----------------------------------------------------------------------------
# 72_crosspathway_gsea.py — preranked Hallmark GSEA of the 5S-dose transcriptional program in p53-competent vs p53-inactive tumours.
# Author: Jonathan Bohlen
# (c) 2026 Jonathan Bohlen. Code accompanying Sengl et al. (2026),
# "Surveillance and selection of 5S ribosomal RNA genes in the human genome."
# Released under the MIT License; see LICENSE at the repository root.
# -----------------------------------------------------------------------------
"""Cross-pathway preranked GSEA (MSigDB Hallmark) of the 5S-dose transcriptional program in p53-competent
tumours and in p53-inactive tumours. Ranks genes by the competent meta slope (and the inactive slope) and
runs preranked GSEA vs MSigDB Hallmark. Paired-bar figure: top pathways' NES in competent vs inactive."""
import pandas as pd, numpy as np, glob, os, warnings; warnings.simplefilter("ignore")
import gseapy as gp
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy import stats
R=os.environ.get("FIVES_DATA","data"); FIG=f'{os.environ.get("FIVES_OUT","output")}/surveillance_v2/figures'; TAB=f'{os.environ.get("FIVES_OUT","output")}/surveillance_v2/tables'
mod=pd.read_csv(f"{R}/results_variants/translation_dosage_module.tsv",sep="\t"); RP=set(mod[mod.is_RP==True].ensg)
e2s=pd.read_csv(f"{R}/results_variants/ensg2symbol.tsv",sep="\t").drop_duplicates("ensg").set_index("ensg").symbol
d={}
for f in glob.glob(f"{R}/out_de_cptacexp/de_rna_*.tsv"):
    t=f.split("de_rna_")[1].replace(".tsv","")
    if t=="allctype": continue
    df=pd.read_csv(f,sep="\t"); df["e"]=[e.split(".")[0] for e in df.ensg]; d[t]=df
COMP=[t for t,df in d.items() if df[df.e.isin(RP)].z_main.mean()<0]
def meta(col_l,col_s):
    A=pd.concat([d[t][["e",col_l,col_s]].rename(columns={col_l:"l",col_s:"s"}) for t in COMP]).dropna(); A=A[A.s>0]
    cnt=A.groupby("e").size(); A=A[A.e.isin(set(cnt[cnt==len(COMP)].index))]
    g=A.groupby("e").apply(lambda x:pd.Series({"z":np.sum(x.l/x.s**2)/np.sum(1/x.s**2)/np.sqrt(1/np.sum(1/x.s**2))}))
    return g.z
def rnk(z):
    r=pd.DataFrame({"g":[e2s.get(e,None) for e in z.index],"s":z.values}).dropna()
    r=r.groupby("g").s.apply(lambda v:v.loc[v.abs().idxmax()]).reset_index(); return r.sort_values("s")
def gsea(z):
    r=rnk(z)
    pre=gp.prerank(rnk=r,gene_sets="MSigDB_Hallmark_2020",min_size=8,max_size=500,permutation_num=1000,seed=0,outdir=None,no_plot=True,threads=4)
    res=pre.res2d.copy(); res["NES"]=pd.to_numeric(res["NES"],errors="coerce"); res["FDR q-val"]=pd.to_numeric(res["FDR q-val"],errors="coerce")
    res["Term"]=res.Term.str.replace("MSigDB_Hallmark_2020__","",regex=False); return res.set_index("Term")
print(f"competent types: {COMP}\nrunning GSEA (competent)..."); gc=gsea(meta("lfc_main","se_main"))
print("running GSEA (inactive)..."); gx=gsea(meta("lfc_inact","se_inact"))
# top pathways by |competent NES|
top=gc.reindex(gc.NES.abs().sort_values(ascending=False).index).head(15)
comp_nes=top.NES; inact_nes=gx.NES.reindex(top.index)
merged=pd.DataFrame({"competent_NES":comp_nes,"competent_FDR":top["FDR q-val"],"inactive_NES":inact_nes}); merged.to_csv(f"{TAB}/crosspathway_gsea.tsv",sep="\t")
print("\nTop pathways (competent NES | inactive NES | FDR):"); print(merged.round(2).to_string())
# figure: paired horizontal bars
y=np.arange(len(top))[::-1]; fig,ax=plt.subplots(figsize=(10,8))
ax.barh(y+0.2,comp_nes.values,height=0.4,color=np.where(comp_nes<0,"#c0392b","#2471a3"),label="p53-COMPETENT")
ax.barh(y-0.2,inact_nes.values,height=0.4,color="#bbbbbb",label="p53-INACTIVE")
for yi,term in zip(y,top.index):
    fdr=top.loc[term,"FDR q-val"]; ax.text(0.02,yi+0.2,("*" if fdr<0.05 else ""),va="center",fontsize=12,color="k")
ax.set_yticks(y); ax.set_yticklabels([t.replace("_"," ").title()[:34] for t in top.index],fontsize=9)
ax.axvline(0,c="k",lw=.6); ax.set_xlabel("GSEA NES (5S-dose effect)"); ax.legend(loc="lower right")
ax.set_title("Cross-pathway program of variant-5S dose: p53-competent vs inactive (Hallmark)\n"
             "* FDR<0.05 (competent). Negative NES = down with 5S dose.",fontsize=10)
fig.tight_layout(); fig.savefig(f"{FIG}/fig_crosspathway_gsea.png",dpi=150,bbox_inches="tight"); plt.close(fig)
print("\nwrote fig_crosspathway_gsea.png + crosspathway_gsea.tsv")
