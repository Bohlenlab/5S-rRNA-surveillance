# -----------------------------------------------------------------------------
# fig2_rescue_permutation_null.py — permutation null test for enrichment of
# short-read/HiFi shared false-positive 5S calls (uniform and propensity-preserving nulls).
# Author: Jonathan Bohlen
# (c) 2026 Jonathan Bohlen. Code accompanying Sengl et al. (2026),
# "Surveillance and selection of 5S ribosomal RNA genes in the human genome."
# Released under the MIT License; see LICENSE at the repository root.
# -----------------------------------------------------------------------------
import sqlite3, os, json, numpy as np
DB=os.environ.get("FIVES_DB","5S_rDNA.db"); BAM=os.environ.get("FIVES_DATA","data")+"/bam"
BOUNDARY=set(range(1,90))|set(range(2034,2169))
SR_MIN_AD=3; LR_AD_MIN=5; LR_VAF=0.003; COHORT="HPRC_Year1"
EXCL_ALL={"HG02486"}; EXCL_HIFI={"HG02818"}
# consensus ref per position -> build eligible (pos,alt) universe (3 non-ref alts per non-boundary pos)
cons="".join(json.load(open(os.environ.get("FIVES_DATA","data")+"/consensus_reference.json"))["consensus"]).upper().replace("-","")
nonbound=[p for p in range(1,2169) if p not in BOUNDARY]
universe=[(p,a) for p in nonbound for a in "ACGT" if a!=cons[p-1]]   # 5832
uni_idx=np.arange(len(universe)); key2i={k:i for i,k in enumerate(universe)}
# GT
con=sqlite3.connect(DB); gt={}
for sid,pos,ref,alt in con.execute("SELECT a.sample_id,v.consensus_pos,v.ref,v.alt FROM variant v JOIN copy c USING(copy_id) JOIN haplotype h USING(haplotype_id) JOIN assembly a USING(assembly_id) WHERE a.cohort=? AND v.alignment_source='gene_unit_t2t'",(COHORT,)):
    if int(pos) in BOUNDARY: continue
    gt.setdefault(sid,set()).add((int(pos),alt))
con.close()
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
            if ad>=min_ad and (pos,alt) in key2i: out[(pos,alt)]=(ad,ad/dp)
    return out
sr_fp={}; hifi_fp={}
for sid in gt:
    if sid in EXCL_ALL: continue
    srp=f"{BAM}/{sid}/{sid}_illumina.tsv"; hfp=f"{BAM}/{sid}/{sid}_hifi_variants.tsv"
    if not os.path.exists(srp): continue
    src=set(parse(srp,SR_MIN_AD))
    hfc=set()
    if os.path.exists(hfp) and sid not in EXCL_HIFI:
        hfc={k for k,(ad,v) in parse(hfp,1).items() if ad>=LR_AD_MIN and v>=LR_VAF}
    sr_fp[sid]=src-gt[sid]; hifi_fp[sid]=hfc-gt[sid]
donors=[s for s in sr_fp if hifi_fp.get(s)]
obs=sum(len(sr_fp[s]&hifi_fp[s]) for s in donors)
print(f"donors={len(donors)}  observed SR-FP∩HiFi-FP overlap={obs}")
# realized HiFi-FP universe (positions HiFi ever FPs) for marginal null
pooled=[k for s in donors for k in hifi_fp[s]]
from collections import Counter; freq=Counter(pooled)

real_idx=np.array([key2i[k] for k in freq]); real_wt=np.array([freq[k] for k in freq],dtype=float); real_wt/=real_wt.sum()
rng=np.random.default_rng(0); NP=2000
def run(sample_fn,label):
    null=np.empty(NP)
    sr_i={s:np.array([key2i[k] for k in sr_fp[s]]) for s in donors}
    sr_set={s:set(sr_i[s].tolist()) for s in donors}
    for j in range(NP):
        tot=0
        for s in donors:
            n=len(hifi_fp[s])
            if n==0: continue
            draw=sample_fn(n)
            tot+=len(set(draw.tolist())&sr_set[s])
        null[j]=tot
    enr=obs/max(null.mean(),1e-9); p=(np.sum(null>=obs)+1)/(NP+1)
    print(f"  [{label}] null mean={null.mean():.1f} sd={null.std():.1f}  enrichment={enr:.1f}x  empirical p={p:.2e} (>= obs in {int((null>=obs).sum())}/{NP})")
# uniform null: draw n from full 5832 universe
run(lambda n: rng.choice(uni_idx,size=n,replace=False), "UNIFORM (all non-boundary sites)")
# marginal null: draw n from realized HiFi-FP positions weighted by frequency
run(lambda n: rng.choice(real_idx,size=min(n,len(real_idx)),replace=False,p=real_wt), "MARGINAL (HiFi-FP propensity-preserving)")
