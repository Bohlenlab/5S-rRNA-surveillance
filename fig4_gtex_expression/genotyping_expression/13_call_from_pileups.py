#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# 13_call_from_pileups.py — re-derive WGS 5S variant calls and filter-dependent
# QC from saved per-donor pileups, without re-streaming.
# Author: Jonathan Bohlen
# (c) 2026 Jonathan Bohlen. Code accompanying Sengl et al. (2026),
# "Surveillance and selection of 5S ribosomal RNA genes in the human genome."
# Released under the MIT License; see LICENSE at the repository root.
# -----------------------------------------------------------------------------
"""Post-streaming caller: re-derive variant calls + filter-dependent QC from the saved
per-donor pileups ({donor}.pileup.tsv.gz). No re-streaming. Change EXCL or the
thresholds below and re-run to re-filter the whole cohort.
Reads control_depth/readlen/slice_reads from the existing {donor}.qc.tsv (filter-independent,
written by the streaming step)."""
import os, gzip, glob, csv, statistics, sys
from pathlib import Path
CD=os.environ.get("COHORT_DIR",str(Path(os.environ.get("FIVES_DATA","data"))/"results"/"wgs"/"cohort"))
MIN_AD     = int(os.environ.get("MIN_AD","5"))
MIN_VAF    = float(os.environ.get("MIN_VAF","0.003"))
MIN_CALLABLE_DP = 50
EXCL       = set(range(1,90))|set(range(2034,2169))   # WGS: edges only
GENE=(630,748)
def region(p): return "gene" if GENE[0]<=p<=GENE[1] else ("nts_pre" if p<GENE[0] else "nts_post")

pgs=sorted(glob.glob(CD+"/*.pileup.tsv.gz"))
print(f"calling from {len(pgs)} pileups | EXCL={sorted(EXCL)[:2]}..{sorted(EXCL)[-2:]} | AD>={MIN_AD} VAF>={MIN_VAF}")
for pg in pgs:
    d=os.path.basename(pg)[:-len(".pileup.tsv.gz")]
    qf=f"{CD}/{d}.qc.tsv"; control=readlen=slice_reads=None
    if os.path.exists(qf):
        r=list(csv.DictReader(open(qf),delimiter='\t'))
        if r:
            control=r[0].get('control_depth'); readlen=r[0].get('readlen'); slice_reads=r[0].get('slice_reads')
    depths={}; calls=[]
    with gzip.open(pg,'rt') as fh:
        for ln in fh:
            p=ln.rstrip("\n").split("\t")
            if len(p)<5: continue
            try: pos=int(p[0]); dp=int(p[3])
            except: continue
            depths[pos]=dp
            if pos in EXCL: continue
            ad=[int(x) for x in p[4].split(",") if x!="."]
            if not ad: continue
            for i,a in enumerate(p[2].split(",")):
                if a in ("<*>",".",""): continue
                adi=ad[i+1] if i+1<len(ad) else 0
                vaf=adi/dp if dp else 0
                if adi>=MIN_AD and vaf>=MIN_VAF:
                    calls.append((pos,p[1],a,dp,adi,round(vaf,5),region(pos)))
    nonx=[v for k,v in depths.items() if k not in EXCL]
    med=statistics.median(nonx) if nonx else 0
    try: est=round(med/float(control),2) if control and float(control)>0 else ""
    except: est=""
    callf=round(sum(1 for v in nonx if v>=MIN_CALLABLE_DP)/len(nonx),4) if nonx else 0
    with open(f"{CD}/{d}.variants.tsv","w") as f:
        f.write("donor_id\tconsensus_pos\tref\talt\tdepth\talt_depth\tvaf\tregion\n")
        for c in sorted(calls): f.write(d+"\t"+"\t".join(map(str,c))+"\n")
    with open(qf,"w") as f:
        f.write("donor_id\tmedian_depth\tcontrol_depth\test_copies\tn_variants\tcallable_fraction\tslice_reads\treadlen\tstatus\n")
        f.write(f"{d}\t{med}\t{control or ''}\t{est}\t{len(calls)}\t{callf}\t{slice_reads or ''}\t{readlen or ''}\tok\n")
print("CALL_FROM_PILEUPS_DONE")
