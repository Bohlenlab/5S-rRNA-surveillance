#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# fig4_coverage_pooldata.py — pool per-donor RNA-seq read depth along the 5S
# repeat unit (summed across a donor's tissue pileups) into a coverage matrix.
# Author: Jonathan Bohlen
# (c) 2026 Jonathan Bohlen. Code accompanying Sengl et al. (2026),
# "Surveillance and selection of 5S ribosomal RNA genes in the human genome."
# Released under the MIT License; see LICENSE at the repository root.
# -----------------------------------------------------------------------------
"""Pool per-donor RNA-seq read depth along the 5S repeat unit (sum across all of a donor's
tissue pileups), for the coverage-distribution panel. Window spans the 5S gene + flanks."""
import glob,gzip,os,numpy as np,pandas as pd
from pathlib import Path
from collections import defaultdict
RNADIR=str(Path(os.environ.get("FIVES_DATA","data"))/"results"/"gtex_rna_pileups_RNAONLY")
OUT=str(Path(os.environ.get("FIVES_DATA","data"))/"figures"/"figure4_rnaseq"/"data")
LO,HI=600,780                     # window (gene region 630-748 + flanks)
positions=list(range(LO,HI+1))
files=glob.glob(RNADIR+"/*.pileup.tsv.gz")
donors=sorted(set(os.path.basename(f).split('.')[0] for f in files))
didx={d:i for i,d in enumerate(donors)}
M=np.zeros((len(donors),len(positions)),dtype=np.int64)
for k,f in enumerate(files):
    d=os.path.basename(f).split('.')[0]; r=didx[d]
    with gzip.open(f,'rt') as fh:
        for ln in fh:
            p=ln.rstrip("\n").split('\t')
            if len(p)<4: continue
            try: pos=int(p[0]); dp=int(p[3])
            except: continue
            if LO<=pos<=HI: M[r,pos-LO]+=dp
    if (k+1)%3000==0: print(f"  {k+1}/{len(files)} files")
df=pd.DataFrame(M,index=donors,columns=positions); df.index.name="donor"
df.to_csv(f"{OUT}/coverage_perdonor_matrix.tsv.gz",sep="\t",compression="gzip")
print("wrote matrix",M.shape,"->",f"{OUT}/coverage_perdonor_matrix.tsv.gz")
print("median pooled depth at gene center (pos 690):",int(np.median(M[:,690-LO])))
