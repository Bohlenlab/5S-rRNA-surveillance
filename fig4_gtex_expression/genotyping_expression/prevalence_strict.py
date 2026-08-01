#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# prevalence_strict.py — Per-donor prevalence of expressed 5S variants under a
# strict, background-exceeding expression call, with a permutation null.
# Author: Jonathan Bohlen
# (c) 2026 Jonathan Bohlen. Code accompanying Sengl et al. (2026),
# "Surveillance and selection of 5S ribosomal RNA genes in the human genome."
# Released under the MIT License; see LICENSE at the repository root.
# -----------------------------------------------------------------------------
"""Per-donor prevalence of expressed 5S variants under a strict call: a carrier is
counted as expressing a variant when its pooled RNA-VAF exceeds the maximum RNA-VAF
seen in any non-carrier at that position (with alt-depth >= 5) — a level never
observed as background. Prevalence is the fraction of donors expressing at least one
carried variant, referenced to a permutation null and stratified by gene-region read
depth. Also writes a per-variant carrier-vs-non-carrier VAF distribution table.

Input : 5S_rDNA.db (read_variant, assembly) + per-donor RNA pileups under FIVES_DATA.
Output: <FIVES_OUT>/wgs_rna_concordance/proband_expression_strict_bydepth.tsv,
        pervariant_vaf_dist.tsv
"""
import sqlite3,glob,gzip,os,numpy as np,pandas as pd
from collections import defaultdict,Counter
DB=os.environ.get("FIVES_DB","5S_rDNA.db")
RNADIR=os.path.join(os.environ.get("FIVES_DATA","data"),"gtex_rna_pileups_RNAONLY")
OUT=os.path.join(os.environ.get("FIVES_OUT","output"),"wgs_rna_concordance")
os.makedirs(OUT,exist_ok=True)
GENE=(630,748); MIN_DP=20; MIN_AD=5; rng=np.random.default_rng(7); NPERM=500
wgs=defaultdict(dict)
for sid,pos,ref,alt,v in sqlite3.connect(DB).execute(
  "SELECT a.sample_id,rv.consensus_pos,rv.ref,rv.alt,rv.vaf FROM read_variant rv JOIN assembly a "
  "ON a.assembly_id=rv.assembly_id WHERE a.cohort='GTEx_v9_WGS' AND rv.modality='illumina' AND rv.consensus_pos BETWEEN ? AND ?",GENE):
    for al in str(alt).split(','):
        if al not in ('<*>','','.',ref): wgs[sid][(pos,al)]=float(v)
rdon=sorted(set(os.path.basename(f).split('.')[0] for f in glob.glob(RNADIR+"/*.pileup.tsv.gz")) & set(wgs))
pDP=defaultdict(lambda:defaultdict(int)); pAD=defaultdict(lambda:defaultdict(int))
for d in rdon:
    for fn in glob.glob(f"{RNADIR}/{d}.*.pileup.tsv.gz"):
        with gzip.open(fn,'rt') as f:
            for ln in f:
                p=ln.rstrip("\n").split('\t')
                if len(p)<5: continue
                try: pos=int(p[0]); dp=int(p[3])
                except: continue
                if not(GENE[0]<=pos<=GENE[1]): continue
                pDP[d][pos]+=dp
                for i,a in enumerate(p[2].split(',')):
                    if a in('<*>','','.',p[1]): continue
                    try: adi=int(p[4].split(',')[i+1])
                    except: continue
                    if adi>0: pAD[d][(pos,a)]+=adi
def vaf(d,pa): return pAD[d].get(pa,0)/pDP[d][pa[0]] if pDP[d][pa[0]] else 0.0
N=len(rdon)
gene_depth={d:np.mean([pDP[d].get(p,0) for p in range(GENE[0],GENE[1]+1)]) for d in rdon}
alleles=sorted(set(a for d in rdon for a in wgs[d]))
assd={pa:[d for d in rdon if pDP[d][pa[0]]>=MIN_DP] for pa in alleles}
nonc={pa:[d for d in assd[pa] if pa not in wgs[d]] for pa in alleles}
alleles=[pa for pa in alleles if len(nonc[pa])>=10]
thrmax={pa:max([vaf(d,pa) for d in nonc[pa]]+[0]) for pa in alleles}
def call(d,pa): return pAD[d].get(pa,0)>=MIN_AD and vaf(d,pa)>thrmax[pa]
carrier_of={pa:[d for d in assd[pa] if pa in wgs[d]] for pa in alleles}
obs=set(d for pa in alleles for d in carrier_of[pa] if call(d,pa))
# permutation null
null=[]
for _ in range(NPERM):
    ex=set()
    for pa in alleles:
        k=len(carrier_of[pa])
        for d in rng.choice(assd[pa],size=min(k,len(assd[pa])),replace=False):
            if call(d,pa): ex.add(d)
    null.append(len(ex))
null=np.array(null); nullpct=100*null.mean()/N
print(f"STRICT call (VAF > max non-carrier, AD>=5), all loci | {N} donors")
print(f"observed {len(obs)} ({100*len(obs)/N:.1f}%) | null {null.mean():.1f}+/-{null.std():.1f} ({nullpct:.1f}%) | "
      f"excess {len(obs)-null.mean():.0f} ({100*(len(obs)-null.mean())/N:.1f}%)  perm p={(np.sum(null>=len(obs))+1)/(NPERM+1):.2g}")
# depth strata
rows=[]
for lab,td in [("all",0),("depth>=500x",500),("depth>=1000x",1000),("depth>=2000x",2000)]:
    sub=[d for d in rdon if gene_depth[d]>=td]; ex=[d for d in obs if d in set(sub)]
    rows.append(dict(stratum=lab,n_probands=len(sub),n_express=len(ex),pct_express=round(100*len(ex)/len(sub),1),
                     null_pct=round(nullpct,1),excess_pct=round(100*len(ex)/len(sub)-nullpct,1)))
sdf=pd.DataFrame(rows); sdf.to_csv(f"{OUT}/proband_expression_strict_bydepth.tsv",sep="\t",index=False)
print(sdf.to_string(index=False))
drive=Counter()
for pa in alleles:
    for d in carrier_of[pa]:
        if call(d,pa): drive[f"{pa[0]}{pa[1]}"]+=1
print("\ntop drivers:",", ".join(f"{v}:{c}" for v,c in drive.most_common(10)))

# per-variant carrier vs non-carrier VAF distribution
vrows=[]
for pa in alleles:
    cv=np.array([vaf(d,pa) for d in carrier_of[pa]]); nv=np.array([vaf(d,pa) for d in nonc[pa]])
    if len(cv)<5: continue
    vrows.append(dict(variant=f"{pa[0]}{pa[1]}",n_carr=len(cv),
        carr_mean=cv.mean(),carr_p90=np.percentile(cv,90),carr_p95=np.percentile(cv,95),carr_max=cv.max(),
        nonc_p90=np.percentile(nv,90),nonc_p95=np.percentile(nv,95),nonc_max=nv.max(),nonc_mean=nv.mean()))
pd.DataFrame(vrows).to_csv(f"{OUT}/pervariant_vaf_dist.tsv",sep="\t",index=False)
print(f"\nsaved per-variant VAF distribution ({len(vrows)} variants) -> pervariant_vaf_dist.tsv")
