# -----------------------------------------------------------------------------
# fig2_rescue_error_vs_mutation.py — position-stratified genotype-linkage and
# assembly-support tests distinguishing recurrent error from shared mutation in
# short-read/HiFi false-positive 5S calls.
# Author: Jonathan Bohlen
# (c) 2026 Jonathan Bohlen. Code accompanying Sengl et al. (2026),
# "Surveillance and selection of 5S ribosomal RNA genes in the human genome."
# Released under the MIT License; see LICENSE at the repository root.
# -----------------------------------------------------------------------------
import sqlite3, os, json, numpy as np
from collections import defaultdict, Counter
DB=os.environ.get("FIVES_DB","5S_rDNA.db"); BAM=os.environ.get("FIVES_DATA","data")+"/bam"
BOUNDARY=set(range(1,90))|set(range(2034,2169))
SR_MIN_AD=3; LR_AD_MIN=5; LR_VAF=0.003; COHORT="HPRC_Year1"
EXCL_ALL={"HG02486"}; EXCL_HIFI={"HG02818"}
con=sqlite3.connect(DB); gt={}
for sid,pos,ref,alt in con.execute("SELECT a.sample_id,v.consensus_pos,v.ref,v.alt FROM variant v JOIN copy c USING(copy_id) JOIN haplotype h USING(haplotype_id) JOIN assembly a USING(assembly_id) WHERE a.cohort=? AND v.alignment_source='gene_unit_t2t'",(COHORT,)):
    if int(pos) in BOUNDARY: continue
    gt.setdefault(sid,set()).add((int(pos),alt))
con.close()
pooled_gt=set().union(*gt.values())   # any (pos,alt) real in >=1 assembly
def parse(path,min_ad):
    out={}
    for line in open(path):
        p=line.rstrip().split("\t")
        if len(p)<5: continue
        try: pos=int(p[0]); dp=int(p[3])
        except: continue
        if dp==0 or pos in BOUNDARY: continue
        alts=[a for a in p[2].split(",") if a and a!="<*>"]; ads=p[4].split(",")
        for i,alt in enumerate(alts):
            try: ad=int(ads[i+1])
            except: continue
            if ad>=min_ad: out[(pos,alt)]=(ad,ad/dp)
    return out
sr_fp={}; hifi_fp={}
for sid in gt:
    if sid in EXCL_ALL: continue
    srp=f"{BAM}/{sid}/{sid}_illumina.tsv"; hfp=f"{BAM}/{sid}/{sid}_hifi_variants.tsv"
    if not os.path.exists(srp): continue
    src=set(parse(srp,SR_MIN_AD)); hfc=set()
    if os.path.exists(hfp) and sid not in EXCL_HIFI:
        hfc={k for k,(ad,v) in parse(hfp,1).items() if ad>=LR_AD_MIN and v>=LR_VAF}
    sr_fp[sid]=src-gt[sid]; hifi_fp[sid]=hfc-gt[sid]
donors=[s for s in sr_fp if hifi_fp.get(s)]; N=len(donors)
# per-position donor counts
ns=Counter(); nh=Counter(); nb=Counter()
for s in donors:
    for k in sr_fp[s]: ns[k]+=1
    for k in hifi_fp[s]: nh[k]+=1
    for k in (sr_fp[s]&hifi_fp[s]): nb[k]+=1
obs_double=sum(nb.values())
exp_double=sum(ns[k]*nh[k]/N for k in set(ns)&set(nh))
print(f"donors={N}")
print(f"(1) POSITION-STRATIFIED genotype-linkage test:")
print(f"    observed double-positive donor-incidences: {obs_double}")
print(f"    expected if SR/HiFi independent within each site: {exp_double:.0f}")
print(f"    excess (donor-specific co-occurrence) = {obs_double-exp_double:.0f}  ({obs_double/exp_double:.2f}x)")
# (2) cross-donor assembly support
cocall_sites=set(nb)   # (pos,alt) double-positive in >=1 donor
sr_only_sites=set(ns)-set(nh)
def frac_real(s): return np.mean([k in pooled_gt for k in s]) if s else float("nan")
print(f"(2) CROSS-DONOR ASSEMBLY SUPPORT (fraction of sites that are real GT in >=1 donor):")
print(f"    SR∩HiFi co-call sites:        {frac_real(cocall_sites)*100:.0f}%  (n={len(cocall_sites)})")
print(f"    SR-only FP sites (no HiFi):   {frac_real(sr_only_sites)*100:.0f}%  (n={len(sr_only_sites)})")
print(f"    all HiFi-FP sites:            {frac_real(set(nh))*100:.0f}%  (n={len(nh)})")
# among co-call sites: split by how many donors share (frequently-shared vs rare)
freq_shared=[k for k in nb if nb[k]>=5]; rare_shared=[k for k in nb if nb[k]==1]
print(f"(3) at FREQUENTLY-shared co-call sites (>=5 donors, n={len(freq_shared)}): real-in-assembly = {frac_real(set(freq_shared))*100:.0f}%")
print(f"    at singleton co-call sites (1 donor, n={len(rare_shared)}): real-in-assembly = {frac_real(set(rare_shared))*100:.0f}%")
