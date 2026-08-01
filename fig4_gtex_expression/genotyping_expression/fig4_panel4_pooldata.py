#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# fig4_panel4_pooldata.py — pool per-donor RNA-VAF at the AUC carrier-HIGH loci
# (rank-skew FDR<0.10) for all assessable donors, labelled carrier/non-carrier.
# Author: Jonathan Bohlen
# (c) 2026 Jonathan Bohlen. Code accompanying Sengl et al. (2026),
# "Surveillance and selection of 5S ribosomal RNA genes in the human genome."
# Released under the MIT License; see LICENSE at the repository root.
# -----------------------------------------------------------------------------
"""Pool per-donor RNA-VAF at the AUC carrier-HIGH loci (FDR<0.10) for all assessable
donors, labelled carrier/non-carrier, for the per-donor prevalence panel."""
import sqlite3,glob,gzip,os,numpy as np,pandas as pd
from pathlib import Path
from collections import defaultdict
DB=str(Path(os.environ.get("FIVES_DB","5S_rDNA.db")))
RNADIR=str(Path(os.environ.get("FIVES_DATA","data"))/"results"/"gtex_rna_pileups_RNAONLY")
C=str(Path(os.environ.get("FIVES_DATA","data"))/"figures"/"wgs_rna_concordance")
OUT=str(Path(os.environ.get("FIVES_DATA","data"))/"figures"/"figure4_rnaseq"/"data")
RS=pd.read_csv(f"{C}/rank_skew_byvariant.tsv",sep="\t")
hi=RS[(RS.fdr<0.10)&(RS.dir.str.contains("HIGH"))]
GEN=set((int(p),a) for p,a in zip(hi.pos,hi.alt))
POSSET={p for p,_ in GEN}; MIN_DP=20
wgs=defaultdict(set)
for sid,pos,ref,alt in sqlite3.connect(DB).execute(
  "SELECT a.sample_id,rv.consensus_pos,rv.ref,rv.alt FROM read_variant rv JOIN assembly a "
  "ON a.assembly_id=rv.assembly_id WHERE a.cohort='GTEx_v9_WGS' AND rv.modality='illumina' "
  "AND rv.consensus_pos IN (%s)"%",".join(map(str,sorted(POSSET)))):
    for al in str(alt).split(','):
        if (pos,al) in GEN: wgs[sid].add((pos,al))
rdon=sorted(set(os.path.basename(f).split('.')[0] for f in glob.glob(RNADIR+"/*.pileup.tsv.gz")))
pDP=defaultdict(lambda:defaultdict(int)); pAD=defaultdict(lambda:defaultdict(int))
for i,d in enumerate(rdon):
    for fn in glob.glob(f"{RNADIR}/{d}.*.pileup.tsv.gz"):
        with gzip.open(fn,'rt') as f:
            for ln in f:
                p=ln.rstrip("\n").split('\t')
                if len(p)<5: continue
                try: pos=int(p[0]); dp=int(p[3])
                except: continue
                if pos not in POSSET: continue
                pDP[d][pos]+=dp
                for j,a in enumerate(p[2].split(',')):
                    if (pos,a) not in GEN: continue
                    try: adi=int(p[4].split(',')[j+1])
                    except: continue
                    if adi>0: pAD[d][(pos,a)]+=adi
    if (i+1)%150==0: print(f"  pooled {i+1}/{len(rdon)}")
rows=[]
for d in rdon:
    for (pos,al) in sorted(GEN):
        dp=pDP[d].get(pos,0)
        if dp<MIN_DP: continue
        ad=pAD[d].get((pos,al),0)
        rows.append(dict(donor=d,variant=f"{pos}{al}",pos=pos,alt=al,
            group="carrier" if (pos,al) in wgs.get(d,()) else "noncarrier",
            pooled_DP=dp,pooled_AD=ad,rna_vaf=ad/dp))
df=pd.DataFrame(rows)
df.to_csv(f"{OUT}/panel4_perdonor_carrierHIGH_vaf.tsv",sep="\t",index=False)
print("wrote",len(df),"rows |",df.donor.nunique(),"donors |",len(GEN),"loci")
