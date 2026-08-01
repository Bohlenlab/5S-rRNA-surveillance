#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# 56_de_cancer.py — per-cancer-type differential expression of the 5S-variant dosage predictor with a TP53 interaction (CPTAC).
# Author: Jonathan Bohlen
# (c) 2026 Jonathan Bohlen. Code accompanying Sengl et al. (2026),
# "Surveillance and selection of 5S ribosomal RNA genes in the human genome."
# Released under the MIT License; see LICENSE at the repository root.
# -----------------------------------------------------------------------------
"""CPTAC 5S-variant -> gene-expression differential expression, per cancer type, with a TP53 interaction.
Design = ~ age_z + SEX + race + SV1..k(RUVr) + <predictor> + tp53 + <predictor>:tp53 , per cancer type.
Extracts the <predictor> MAIN coefficient AND the <predictor>:tp53_deficient INTERACTION.
Predictor via DE_METRIC (rna_excess_z / mut_z / cn_z).
Outputs out_de/de_<metric>_<ctype>.tsv (per-gene main + interaction).
Then meta-analysis + cytosolic-RP module readout when run with --summarize."""
import sys, os, re, glob, warnings, numpy as np, pandas as pd
warnings.simplefilter("ignore")
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
R=os.environ.get("FIVES_DATA","data"); RV=f"{R}/results_variants"; OUT=f"{R}/out_de{os.environ.get('DE_SUFFIX','')}"; os.makedirs(OUT,exist_ok=True)
MIN_CASE=25; N_CPUS=int(os.environ.get("DE_CPUS","4")); METRIC=os.environ.get("DE_METRIC","rna_excess_z")
SHORT={"rna_excess_z":"rna","mut_z":"mut","cn_z":"cn"}

def ruvr_svs(raw, D, k_max=15, seed=0):
    X=raw.astype(np.float64); lib=X.sum(1,keepdims=True); lib[lib==0]=1
    Y=np.log1p(X/lib*1e4); Y=Y-Y.mean(0)
    beta,*_=np.linalg.lstsq(D,Y,rcond=None); Rr=Y-D@beta
    U,S,_=np.linalg.svd(Rr,full_matrices=False); ev=S**2
    rng=np.random.default_rng(seed); null=np.zeros((3,len(ev)))
    for p in range(3):
        cols=rng.choice(Rr.shape[1],min(2000,Rr.shape[1]),replace=False)
        Rp=np.column_stack([rng.permutation(Rr[:,j]) for j in cols])
        sp=np.linalg.svd(Rp-Rp.mean(0),compute_uv=False)**2; sp=sp*(ev.sum()/(Rr.shape[1]/Rp.shape[1])/sp.sum())
        null[p,:len(sp)]=sp[:len(ev)]
    thr=null.mean(0); k=int(min(k_max,max(1,np.sum(ev[:k_max]>thr[:k_max]))))
    SV=U[:,:k]*S[:k]; SV=(SV-SV.mean(0))/(SV.std(0)+1e-9); return SV,k

def load_data():
    M=pd.read_pickle(os.environ.get("DE_MATRIX",f"{R}/cptac_expr_counts.pkl"))   # genes x samples (int)
    dm=pd.read_csv(os.environ.get("DE_METRICS",f"{RV}/cptac_donor_metrics.tsv"),sep="\t").set_index("donor")
    cov=pd.read_csv(os.environ.get("DE_COV",f"{RV}/cptac_covariates.tsv"),sep="\t").set_index("case")
    # derive SEX from expression: XIST vs RPS4Y1/DDX3Y
    def gg(ensg): return M.loc[ensg] if ensg in M.index else pd.Series(0,index=M.columns)
    lib=M.sum(0).replace(0,1); cpm=lambda e: np.log1p(gg(e)/lib*1e6)
    xist=cpm("ENSG00000229807"); ychr=cpm("ENSG00000129824")+cpm("ENSG00000067048")  # RPS4Y1+DDX3Y
    sex=pd.Series(np.where(ychr>xist,"M","F"),index=M.columns)
    df=os.environ.get("DE_DEPTH_FILE")
    if df:
        dep=pd.read_csv(df,sep="\t").set_index("case").depth; cov["_depth"]=cov.index.map(dep)
    else: cov["_depth"]=1e9
    return M,dm,cov,sex

def run_ctype(M,dm,cov,sex,ctype):
    tag=f"{SHORT[METRIC]}_{re.sub(r'[^A-Za-z0-9]+','_',ctype)}"
    if os.environ.get("DE_RESUME") and os.path.exists(f"{OUT}/de_{tag}.tsv"): return
    cases=[c for c in M.columns if cov.get("ctype",{}).get(c)==ctype] if False else \
          [c for c in M.columns if (c in cov.index and cov.loc[c,"ctype"]==ctype and c in dm.index
           and cov.loc[c,"tp53"] in ("WT","deficient") and pd.notna(dm.loc[c,METRIC]) and pd.notna(cov.loc[c,"age_yr"])
           and cov.loc[c,"_depth"]>=float(os.environ.get("DE_MIN_DEPTH","0")))]
    if len(cases)<MIN_CASE: print(f"  {ctype}: only {len(cases)} cases, skip"); return
    cnt=M[cases].T.values.astype(int)                        # samples x genes
    keep=((cnt>=10).mean(0)>=0.5); cnt=cnt[:,keep]; gn=M.index[keep]
    ob=pd.DataFrame(index=cases)
    ob["age_z"]=(cov.loc[cases,"age_yr"].astype(float)-cov.loc[cases,"age_yr"].astype(float).mean())/cov.loc[cases,"age_yr"].astype(float).std()
    ob["SEX"]=sex.reindex(cases).values
    race=cov.loc[cases,"race"].astype(str); race=race.where(race.isin(race.value_counts()[race.value_counts()>=5].index),"other")
    ob["race"]=race.values
    ob["tp53"]=pd.Categorical(cov.loc[cases,"tp53"].values,categories=["WT","deficient"])
    ob["metric"]=dm.loc[cases,METRIC].astype(float).values
    cnterm=""
    if os.environ.get("DE_ADD_CN") and "cn_z" in dm.columns and METRIC!="cn_z":
        ob["cn_z"]=dm.loc[cases,"cn_z"].astype(float).fillna(0).values; cnterm=" + cn_z"
    # RUVr design matrix (intercept + covars + predictor + tp53 + interaction) for residualization
    Dcols=[np.ones(len(cases)), ob.age_z.values, ob.metric.values, (ob.tp53=="deficient").astype(float).values,
           ob.metric.values*(ob.tp53=="deficient").astype(float).values]
    if cnterm: Dcols.append(ob.cn_z.values)
    for cat in ["SEX","race"]:
        for lv in pd.get_dummies(ob[cat],drop_first=True).T.values: Dcols.append(lv.astype(float))
    D=np.column_stack(Dcols)
    SV,k=ruvr_svs(cnt,D)
    keepcols=["age_z","SEX","race","tp53","metric"]+(["cn_z"] if cnterm else [])
    md=ob[keepcols].copy()
    md["SEX"]=md.SEX.astype("category"); md["race"]=md.race.astype("category")
    for j in range(k): md[f"SV{j+1}"]=SV[:,j]
    svterm=" + ".join(f"SV{j+1}" for j in range(k))
    cdf=pd.DataFrame(cnt,index=cases,columns=gn)
    design=f"~ age_z + SEX + race{cnterm} + {svterm} + metric + tp53 + metric:tp53"
    try:
        dds=DeseqDataSet(counts=cdf,metadata=md,design=design,quiet=True,n_cpus=N_CPUS); dds.deseq2()
        cols=list(dds.obsm["design_matrix"].columns)
        ivar=[c for c in cols if ("metric" in c and c!="metric" and (":" in c or "tp53" in c))]  # interaction col
        def coef(pred):  # 'main'=competent slope, 'inter'=interaction, 'inactive'=competent+interaction slope
            if pred=="main": vec=np.array([1.0 if c=="metric" else 0.0 for c in cols])
            elif pred=="inter":
                if not ivar: return None
                vec=np.array([1.0 if c==ivar[0] else 0.0 for c in cols])
            elif pred=="inactive":  # main + interaction = slope in the p53-inactive group (full-model SE)
                vec=np.array([1.0 if (c=="metric" or (ivar and c==ivar[0])) else 0.0 for c in cols])
            st=DeseqStats(dds,contrast=vec,quiet=True); st.summary(); return st.results_df
        rm=coef("main"); ri=coef("inter"); rx=coef("inactive") if ivar else None
        if rm is None: print(f"  {ctype}: no metric coef"); return
        out=pd.DataFrame(index=rm.index)
        out["lfc_main"]=rm.log2FoldChange; out["se_main"]=rm.lfcSE; out["z_main"]=rm.stat; out["p_main"]=rm.pvalue
        if ri is not None:
            out["lfc_inter"]=ri.log2FoldChange; out["se_inter"]=ri.lfcSE; out["z_inter"]=ri.stat; out["p_inter"]=ri.pvalue
        if rx is not None:
            out["lfc_inact"]=rx.log2FoldChange; out["se_inact"]=rx.lfcSE; out["z_inact"]=rx.stat; out["p_inact"]=rx.pvalue
        out["ctype"]=ctype; out["n"]=len(cases); out["k_sv"]=k; out["n_mut"]=int((ob.tp53=="deficient").sum())
        out.index.name="ensg"; out.to_csv(f"{OUT}/de_{tag}.tsv",sep="\t")
        print(f"  {ctype:8} n={len(cases)} (mut {int((ob.tp53=='deficient').sum())}) genes={len(out)} k_sv={k} | main z(RP-ish) done, inter={'yes' if ri is not None else 'NO'}",flush=True)
    except Exception as e:
        print(f"  FAIL {ctype}: {type(e).__name__}: {str(e)[:120]}",flush=True)

def summarize():
    # RP module (cytosolic RP = is_RP True)
    mod=pd.read_csv(f"{RV}/translation_dosage_module.tsv",sep="\t"); rp=set(mod[mod.is_RP==True].ensg)
    files=glob.glob(f"{OUT}/de_{SHORT[METRIC]}_*.tsv"); allr=[]
    for f in files:
        d=pd.read_csv(f,sep="\t"); d["ensg_"]=[e.split(".")[0] for e in d.ensg]; allr.append(d)
    if not allr: print("no results"); return
    A=pd.concat(allr)
    # fixed-effect inverse-variance meta per gene, for main + interaction
    def meta(col_lfc,col_se):
        sub=A.dropna(subset=[col_lfc,col_se]); w=1/sub[col_se]**2
        g=sub.groupby("ensg_").apply(lambda x:pd.Series({"beta":np.sum(x[col_lfc]/x[col_se]**2)/np.sum(1/x[col_se]**2),
                                                          "se":np.sqrt(1/np.sum(1/x[col_se]**2)),"nt":len(x)}))
        g["z"]=g.beta/g.se; return g
    print(f"\n=== METRIC={METRIC} | cytosolic-RP module readout ({len(rp)} genes) ===")
    for lab,cl,cs in [("MAIN (5S->RP)","lfc_main","se_main"),("INTERACTION (x TP53-mut)","lfc_inter","se_inter")]:
        if cl not in A.columns: continue
        g=meta(cl,cs); rpz=g[g.index.isin(rp)]
        print(f"  {lab:26}: RP-module mean meta-z = {rpz.z.mean():+.2f}  (n_RP={len(rpz)}, {(rpz.z<0).sum()} down / {(rpz.z>0).sum()} up); all-gene median z={g.z.median():+.2f}")
    A.to_csv(f"{OUT}/de_{SHORT[METRIC]}_allctype.tsv",sep="\t",index=False)
    print(f"saved {OUT}/de_{SHORT[METRIC]}_allctype.tsv")

if __name__=="__main__":
    if "--summarize" in sys.argv: summarize(); sys.exit()
    M,dm,cov,sex=load_data()
    types=cov.ctype.value_counts(); types=[t for t in types.index if pd.notna(t) and types[t]>=MIN_CASE]
    if "--ctype" in sys.argv: types=[sys.argv[sys.argv.index("--ctype")+1]]
    print(f"[DE_METRIC={METRIC}] cancer types: {types}",flush=True)
    for t in types: run_ctype(M,dm,cov,sex,t)
    summarize()
