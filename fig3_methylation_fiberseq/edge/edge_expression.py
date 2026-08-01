# -----------------------------------------------------------------------------
# edge_expression.py — Association between nascent-RNA expression of 5S-gene variants and array-edge copy placement.
# Author: Jonathan Bohlen
# (c) 2026 Jonathan Bohlen. Code accompanying Sengl et al. (2026),
# "Surveillance and selection of 5S ribosomal RNA genes in the human genome."
# Released under the MIT License; see LICENSE at the repository root.
# -----------------------------------------------------------------------------
import os, pandas as pd, numpy as np, re
from scipy import stats
import statsmodels.api as sm

HPRC=os.environ.get("FIVES_DATA","data")
SAMPLES=["NA18505","NA18508","NA18522","NA18879"]; GS,GE=629,748; TEL,CEN=30,15

# continuous RNA-VAF per (sample,pos,alt) from gold pileup (0 if allele absent)
def pileup_vaf(f):
    v={}
    for ln in open(f):
        p=ln.rstrip("\n").split("\t")
        if len(p)<5: continue
        pos,ref,alt,dp,ad=p
        try: dp=int(dp); ads=[int(x) for x in ad.split(",")]
        except: continue
        for i,al in enumerate(alt.split(",")):
            if al in("<*>",".",""): continue
            a=ads[i+1] if i+1<len(ads) else 0
            v[(int(pos),al)]=a/dp if dp>0 else 0.0
    return v
rvaf={s:pileup_vaf(f'{os.environ.get("FIVES_DATA","data")}/nascent_edge/{s}_gold_pileup.tsv') for s in SAMPLES}

meth=pd.read_csv(f"{HPRC}/hifi_array_percopy_meth.tsv",sep="\t")
def tel_end(s,hap):
    m=meth[(meth['sample']==s)&(meth.hap==hap)].sort_values('copy_number')
    if len(m)<8: return None
    k=max(3,len(m)//4); return 'low' if m.head(k).mean_meth.mean()<m.tail(k).mean_meth.mean() else 'high'
def parse(s):
    if pd.isna(s) or str(s).strip() in ("none",""): return []
    return [(int(m.group(1))-GS+1,m.group(3)) for m in (re.match(r"\s*(\d+):([ACGT])>([ACGT])",t) for t in str(s).split(";")) if m and GS<=int(m.group(1))<GE]

rows=[]; tot_copies={}
for s in SAMPLES:
    tc=0
    for hap in ("hap1","hap2"):
        db=pd.read_csv(f"{HPRC}/databases/{s}_{hap}.tsv",sep="\t"); n=len(db); tc+=n; te=tel_end(s,hap)
        for _,r in db.iterrows():
            cid=int(r.copy_id)
            if te=='low': rt,rc=cid-1,n-cid
            elif te=='high': rt,rc=n-cid,cid-1
            else: rt=rc=99999
            cls='tel_edge' if rt<TEL else 'cen_edge' if rc<CEN else 'interior'
            for gp,alt in parse(r.gene_variants):
                rows.append(dict(sample=s,gene_pos=gp,alt=alt,cls=cls))
    tot_copies[s]=tc
A=pd.DataFrame(rows)
V=A.groupby(["sample","gene_pos","alt"]).agg(
    n_copies=("cls","size"), n_tel=("cls",lambda x:(x=='tel_edge').sum()),
    n_cen=("cls",lambda x:(x=='cen_edge').sum()), n_int=("cls",lambda x:(x=='interior').sum())).reset_index()
V["n_edge"]=V.n_tel+V.n_cen
V["rna_vaf"]=[rvaf[s].get((p,a),0.0) for s,p,a in zip(V["sample"],V.gene_pos,V.alt)]
V["dna_vaf"]=[nc/tot_copies[s] for s,nc in zip(V["sample"],V.n_copies)]
V["log_rna"]=np.log10(V.rna_vaf+1e-4)
V.to_csv(f'{os.environ.get("FIVES_DATA","data")}/nascent_edge/variant_expression_table.tsv',sep="\t",index=False)
print(f"=== {len(V)} variants, RNA-VAF continuous (median {V.rna_vaf.median():.4f}, {(V.rna_vaf>0).sum()} with >0 RNA reads) ===\n")

print("--- does RNA expression track EDGE copies more than INTERIOR copies? (Spearman) ---")
for col in ["n_edge","n_tel","n_cen","n_int","n_copies","dna_vaf"]:
    r,p=stats.spearmanr(V[col],V.rna_vaf); print(f"  RNA-VAF vs {col:9s}: rho={r:+.2f} p={p:.3f}")

print("\n--- regression: RNA-VAF ~ n_edge + n_interior (partial: which copies carry expression?) ---")
X=sm.add_constant(V[["n_edge","n_int"]].astype(float)); m=sm.OLS(V.rna_vaf,X).fit()
for k in ["n_edge","n_int"]: print(f"  {k:8s}: beta={m.params[k]:+.5f}  p={m.pvalues[k]:.3f}")
print("\n--- telomere vs centromere vs interior partials ---")
X=sm.add_constant(V[["n_tel","n_cen","n_int"]].astype(float)); m=sm.OLS(V.rna_vaf,X).fit()
for k in ["n_tel","n_cen","n_int"]: print(f"  {k:8s}: beta={m.params[k]:+.5f}  p={m.pvalues[k]:.3f}")

print("\n--- expression EFFICIENCY: over/under-expression vs edge fraction ---")
# residual of RNA-VAF given DNA copy fraction; is it higher when copies are edge?
V["frac_edge"]=V.n_edge/V.n_copies
res=sm.OLS(V.rna_vaf,sm.add_constant(V.dna_vaf.astype(float))).fit().resid
r,p=stats.spearmanr(V.frac_edge,res); print(f"  Spearman(edge-fraction, RNA-VAF residual over DNA): rho={r:+.2f} p={p:.3f}")
r2,p2=stats.spearmanr(V.n_tel/V.n_copies,res); print(f"  Spearman(telomere-fraction, residual):            rho={r2:+.2f} p={p2:.3f}")
