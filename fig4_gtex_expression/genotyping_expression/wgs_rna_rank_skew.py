#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# wgs_rna_rank_skew.py — threshold-free rank test of whether DNA (WGS) carriers of
# each gene-region 5S variant have skewed pooled RNA-VAF, via Mann-Whitney U / AUC.
# Author: Jonathan Bohlen
# (c) 2026 Jonathan Bohlen. Code accompanying Sengl et al. (2026),
# "Surveillance and selection of 5S ribosomal RNA genes in the human genome."
# Released under the MIT License; see LICENSE at the repository root.
# -----------------------------------------------------------------------------
"""Threshold-free carrier rank test. For each gene (pos,alt) WGS variant: rank all
assessable donors by their pooled RNA-VAF at that position, then test whether the
(WGS) carriers are ranked differently from non-carriers. Mann-Whitney U (carrier vs
non-carrier RNA-VAF) with AUC = rank-biserial effect size (fraction of carrier/non-carrier
pairs in which the carrier has the higher RNA-VAF). Uses the continuous VAF with no
detection threshold. Writes a per-variant statistics table and an AUC vs FDR volcano."""
import sqlite3, glob, gzip, os, numpy as np
from pathlib import Path
from collections import defaultdict
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
DB=str(Path(os.environ.get("FIVES_DB","5S_rDNA.db")))
RNADIR=str(Path(os.environ.get("FIVES_DATA","data"))/"results"/"gtex_rna_pileups_RNAONLY")
OUT=str(Path(os.environ.get("FIVES_DATA","data"))/"figures"/"wgs_rna_concordance")
GENE=(630,748); MIN_DP=20; MIN_CARR=3; MIN_NONC=10
wgs=defaultdict(dict)
for sid,pos,ref,alt,vaf in sqlite3.connect(DB).execute(
  "SELECT a.sample_id,rv.consensus_pos,rv.ref,rv.alt,rv.vaf FROM read_variant rv JOIN assembly a "
  "ON a.assembly_id=rv.assembly_id WHERE a.cohort='GTEx_v9_WGS' AND rv.modality='illumina' "
  "AND rv.consensus_pos BETWEEN ? AND ?",GENE):
    for al in str(alt).split(','):
        if al not in ('<*>','','.',ref): wgs[sid][(pos,al)]=float(vaf)
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
def vaf(d,pos,alt):
    return pAD[d].get((pos,alt),0)/pDP[d][pos] if pDP[d][pos] else 0.0
alleles=sorted(set(a for d in rdon for a in wgs[d]))
rows=[]
for (pos,alt) in alleles:
    assd=[d for d in rdon if pDP[d][pos]>=MIN_DP]
    carr=[d for d in assd if (pos,alt) in wgs[d]]; nonc=[d for d in assd if (pos,alt) not in wgs[d]]
    if len(carr)<MIN_CARR or len(nonc)<MIN_NONC: continue
    cv=np.array([vaf(d,pos,alt) for d in carr]); nv=np.array([vaf(d,pos,alt) for d in nonc])
    U,p=mannwhitneyu(cv,nv,alternative='two-sided')
    auc=U/(len(cv)*len(nv))                                  # P(carrier VAF > non-carrier VAF)
    nc_rate=np.mean(nv>0)                                    # fraction of non-carriers with nonzero RNA-VAF (background)
    rows.append(dict(pos=pos,alt=alt,n_carr=len(carr),n_nonc=len(nonc),auc=auc,p=p,
                     med_carr=np.median(cv),med_nonc=np.median(nv),
                     carr_frac_pos=np.mean(cv>0),nonc_frac_pos=nc_rate))
import pandas as pd
D=pd.DataFrame(rows)
D["fdr"]=multipletests(D.p,method='fdr_bh')[1]
D["dir"]=np.where(D.auc>0.5,"carrier-HIGH (expression)","carrier-LOW (silencing)")
D=D.sort_values("p")
D.to_csv(f"{OUT}/rank_skew_byvariant.tsv",sep="\t",index=False)
print(f"{len(rdon)} donors | {len(D)} gene variants tested (>= {MIN_CARR} carriers, {MIN_NONC} non-carriers)\n")
sig=D[D.fdr<0.10]
print(f"FDR<0.10: {len(sig)}  ({(sig.auc>0.5).sum()} carrier-HIGH / {(sig.auc<0.5).sum()} carrier-LOW)\n")
pd.set_option('display.width',200)
print("=== top 25 by p (AUC = fraction of (carrier,non-carrier) pairs where carrier VAF is higher) ===")
show=D.head(25).copy()
show["var"]=show.pos.astype(str)+show.alt
print(show[["var","n_carr","n_nonc","auc","p","fdr","carr_frac_pos","nonc_frac_pos","dir"]].to_string(index=False,
      formatters={"auc":"{:.3f}".format,"p":"{:.1e}".format,"fdr":"{:.1e}".format,
                  "carr_frac_pos":"{:.2f}".format,"nonc_frac_pos":"{:.2f}".format}))

# Volcano: AUC vs -log10 FDR
fig,ax=plt.subplots(figsize=(8,6))
col=np.where(D.fdr<0.10,np.where(D.auc>0.5,"#1b7837","#b2182b"),"#bbbbbb")
ax.scatter(D.auc,-np.log10(D.fdr),c=col,s=20+ D.n_carr.clip(upper=60))
for _,r in D[D.fdr<0.10].iterrows(): ax.annotate(f"{r.pos}{r.alt}",(r.auc,-np.log10(r.fdr)),fontsize=6)
ax.axvline(0.5,ls='--',color='#999'); ax.axhline(-np.log10(0.10),ls=':',color='#999')
ax.set_xlabel("AUC  (<0.5 carrier RNA-VAF ranks lower   |   >0.5 carrier RNA-VAF ranks higher)")
ax.set_ylabel("-log10 FDR"); ax.set_title(f"Carrier RNA-VAF rank skew per gene variant (n={len(rdon)} donors)")
plt.tight_layout(); plt.savefig(f"{OUT}/fig15_rank_skew_volcano.png",dpi=150); plt.close()
print(f"\n-> rank_skew_byvariant.tsv + fig15_rank_skew_volcano.png in {OUT}")
