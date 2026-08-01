#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# 23_call_from_rna_pileups.py — call 5S variants from per-(donor,tissue,sample)
# RNA pileups, producing per-tissue calls and pooled-per-donor calls.
# Author: Jonathan Bohlen
# (c) 2026 Jonathan Bohlen. Code accompanying Sengl et al. (2026),
# "Surveillance and selection of 5S ribosomal RNA genes in the human genome."
# Released under the MIT License; see LICENSE at the repository root.
# -----------------------------------------------------------------------------
"""Post-streaming RNA caller. From per-(donor,tissue,sample) pileups, produce
(a) per-tissue calls (tissue provenance) and (b) pooled-per-donor calls (sum DP/AD across a donor's
tissues, for depth). RNA mask = full exclusion set (positions maskable to pseudogene contamination in
RNA). Re-runnable. Env: DONORS, MIN_AD, MIN_VAF."""
import os, gzip, glob, csv, collections
from pathlib import Path
CD=str(Path(os.environ.get("FIVES_DATA","data"))/"results"/"rnaseq"/"cohort")
MIN_AD=int(os.environ.get("MIN_AD","5")); MIN_VAF=float(os.environ.get("MIN_VAF","0.005"))
EXCL=set(range(1,90))|set(range(790,933))|set(range(974,1058))|set(range(2034,2169))
GENE=(630,748); DONORS=[d for d in os.environ.get("DONORS","").split(",") if d]
def region(p): return "gene" if GENE[0]<=p<=GENE[1] else ("nts_pre" if p<GENE[0] else "nts_post")
def parse(pg):
    d={}
    with gzip.open(pg,'rt') as fh:
        for ln in fh:
            p=ln.rstrip("\n").split("\t")
            if len(p)<5: continue
            try: pos=int(p[0]); dp=int(p[3])
            except: continue
            ad=[int(x) for x in p[4].split(",") if x!="."]
            if not ad: continue
            altd={}
            for i,a in enumerate(p[2].split(",")):
                if a not in ("<*>",".",""): altd[a]=ad[i+1] if i+1<len(ad) else 0
            d[pos]=(p[1],dp,altd)
    return d
bydonor=collections.defaultdict(list)
for pg in sorted(glob.glob(CD+"/*.pileup.tsv.gz")):
    name=os.path.basename(pg)[:-len(".pileup.tsv.gz")]; donor=name.split(".")[0]
    if DONORS and donor not in DONORS: continue
    parts=name.split("."); bydonor[donor].append((parts[1] if len(parts)>=3 else "NA",pg))
pt=open(CD+"/_rna_pertissue.variants.tsv","w"); pt.write("donor_id\ttissue\tconsensus_pos\tref\talt\tdepth\talt_depth\tvaf\tregion\n")
pl=open(CD+"/_rna_pooled.variants.tsv","w"); pl.write("donor_id\tconsensus_pos\tref\talt\tdepth\talt_depth\tvaf\tregion\tn_tissues\n")
nd=0
for donor,lst in bydonor.items():
    nd+=1; pool=collections.defaultdict(lambda:["",0,collections.Counter(),collections.Counter()])
    for tissue,pg in lst:
        for pos,(ref,dp,altd) in parse(pg).items():
            if pos in EXCL: continue
            e=pool[pos]; e[0]=ref; e[1]+=dp
            for a,ad in altd.items():
                e[2][a]+=ad
                if dp and ad>=MIN_AD and ad/dp>=MIN_VAF:
                    e[3][a]+=1
                    pt.write(f"{donor}\t{tissue}\t{pos}\t{ref}\t{a}\t{dp}\t{ad}\t{round(ad/dp,5)}\t{region(pos)}\n")
    for pos,(ref,dps,ads,tc) in pool.items():
        for a,adsum in ads.items():
            vaf=adsum/dps if dps else 0
            if adsum>=MIN_AD and vaf>=MIN_VAF:
                pl.write(f"{donor}\t{pos}\t{ref}\t{a}\t{dps}\t{adsum}\t{round(vaf,5)}\t{region(pos)}\t{tc[a]}\n")
pt.close(); pl.close(); print(f"RNA_CALL_DONE donors={nd}")
