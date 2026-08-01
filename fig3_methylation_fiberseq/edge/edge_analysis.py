# -----------------------------------------------------------------------------
# edge_analysis.py — Association between nascent-RNA detection of 5S-gene variants and array-edge copy placement.
# Author: Jonathan Bohlen
# (c) 2026 Jonathan Bohlen. Code accompanying Sengl et al. (2026),
# "Surveillance and selection of 5S ribosomal RNA genes in the human genome."
# Released under the MIT License; see LICENSE at the repository root.
# -----------------------------------------------------------------------------
import os, pandas as pd, numpy as np, re
from scipy import stats
import statsmodels.api as sm

HPRC=os.environ.get("FIVES_DATA","data")
SAMPLES=["NA18505","NA18508","NA18522","NA18879"]
GS,GE=629,748; TEL,CEN=30,15

# 1. nascent detected variants (AD>=5, VAF>=0.5%) per sample -> set of (gene_pos, alt)
det={}
for s in SAMPLES:
    d=pd.read_csv(f'{os.environ.get("FIVES_DATA","data")}/nascent_edge/{s}_rna_vs_lr.tsv',sep="\t",comment="#")
    d=d[(d.ad>=5)&(d.vaf>=0.005)]; det[s]=set(zip(d.pos.astype(int),d.alt))

# 2. telomere end per (sample,hap) from methylation gradient (telomere = hypomethylated)
meth=pd.read_csv(f"{HPRC}/hifi_array_percopy_meth.tsv",sep="\t")
def tel_end(s,hap):
    m=meth[(meth['sample']==s)&(meth.hap==hap)].sort_values('copy_number')
    if len(m)<8: return None,None,None
    k=max(3,len(m)//4); lo=m.head(k).mean_meth.mean(); hi=m.tail(k).mean_meth.mean()
    return ('low' if lo<hi else 'high'),lo,hi

# 3. assembly per-copy 5S gene variants
def parse(s):
    if pd.isna(s) or str(s).strip() in ("none",""): return []
    out=[]
    for tok in str(s).split(";"):
        mm=re.match(r"\s*(\d+):([ACGT])>([ACGT])",tok)
        if mm and GS<=int(mm.group(1))<GE: out.append((int(mm.group(1))-GS+1,mm.group(3)))
    return out

print("=== orientation (telomere end) per hap, via methylation ===")
rows=[]
for s in SAMPLES:
    for hap in ("hap1","hap2"):
        db=pd.read_csv(f"{HPRC}/databases/{s}_{hap}.tsv",sep="\t"); n=len(db)
        te,lo,hi=tel_end(s,hap)
        print(f"  {s} {hap}: {n} copies | meth low-end={lo}, high-end={hi} -> telomere={te}")
        for _,r in db.iterrows():
            cid=int(r.copy_id)
            if te=='low': rt,rc=cid-1,n-cid
            elif te=='high': rt,rc=n-cid,cid-1
            else: rt=rc=None
            cls=('unknown' if rt is None else 'tel_edge' if rt<TEL else 'cen_edge' if rc<CEN else 'interior')
            for gp,alt in parse(r.gene_variants):
                rows.append(dict(sample=s,hap=hap,copy_id=cid,n=n,gene_pos=gp,alt=alt,cls=cls,
                                 detected=(gp,alt) in det[s]))
A=pd.DataFrame(rows)
# aggregate to unique donor variant
def agg(g):
    return pd.Series(dict(n_copies=len(g), n_tel=(g.cls=='tel_edge').sum(), n_cen=(g.cls=='cen_edge').sum(),
        n_int=(g.cls=='interior').sum(), has_tel=(g.cls=='tel_edge').any(), has_edge=g.cls.isin(['tel_edge','cen_edge']).any(),
        frac_edge=g.cls.isin(['tel_edge','cen_edge']).mean(), detected=g.detected.iloc[0]))
V=A.groupby(["sample","gene_pos","alt"]).apply(agg).reset_index()
V["detected"]=V.detected.astype(bool); V["has_edge"]=V.has_edge.astype(bool); V["has_tel"]=V.has_tel.astype(bool)
print(f"\n=== {len(V)} unique assembly 5S variants (donor-level); {V.detected.sum()} detected in nascent RNA ===")

print("\n--- detection rate by edge presence (RAW) ---")
for col,lab in [("has_edge","≥1 edge copy (tel30|cen15)"),("has_tel","≥1 telomere-edge copy")]:
    t=pd.crosstab(V[col],V.detected)
    a=V[V[col]].detected.mean(); b=V[~V[col]].detected.mean()
    OR,p=stats.fisher_exact([[t.loc[True,True],t.loc[True,False]],[t.loc[False,True],t.loc[False,False]]]) if t.shape==(2,2) else (np.nan,np.nan)
    print(f"  {lab}: detected {a:.0%} (with) vs {b:.0%} (without)  OR={OR:.2f} Fisher p={p:.3f}")

print("\n--- copy-number is a confounder (detection tracks copy count) ---")
print(f"  detected variants: median {V[V.detected].n_copies.median():.0f} copies; non-detected: {V[~V.detected].n_copies.median():.0f}")
print(f"  MWU frac_edge detected vs not: p={stats.mannwhitneyu(V[V.detected].frac_edge,V[~V.detected].frac_edge).pvalue:.3f}")

print("\n--- logistic: detected ~ log(n_copies) + edge  (edge effect ON TOP of copy number) ---")
V["lognc"]=np.log10(V.n_copies)
for edgecol in ["has_edge","frac_edge","has_tel"]:
    X=sm.add_constant(V[["lognc",edgecol]].astype(float));
    try:
        m=sm.Logit(V.detected.astype(int),X).fit(disp=0)
        print(f"  {edgecol:10s}: beta={m.params[edgecol]:+.2f}  p={m.pvalues[edgecol]:.3f}   (copies beta={m.params['lognc']:+.2f} p={m.pvalues['lognc']:.3f})")
    except Exception as e: print(f"  {edgecol}: {e}")
V.to_csv(f'{os.environ.get("FIVES_DATA","data")}/nascent_edge/variant_edge_table.tsv',sep="\t",index=False)
print("\nwrote variant_edge_table.tsv")
