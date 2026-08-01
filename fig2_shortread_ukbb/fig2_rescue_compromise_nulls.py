# -----------------------------------------------------------------------------
# fig2_rescue_compromise_nulls.py — two compromise null models (error-matched and
# genotype-conditioned) for enrichment of short-read/HiFi shared false-positive 5S calls.
# Author: Jonathan Bohlen
# (c) 2026 Jonathan Bohlen. Code accompanying Sengl et al. (2026),
# "Surveillance and selection of 5S ribosomal RNA genes in the human genome."
# Released under the MIT License; see LICENSE at the repository root.
# -----------------------------------------------------------------------------
"""HiFi-rescue enrichment — two compromise nulls between the over-lenient UNIFORM and the
over-conservative MARGINAL null.
 Option 1 (ERROR-MATCHED positional null): like the marginal null, but positional weights are
   the pure-ERROR propensity. Positions that are real GT in >=1 donor are de-weighted to the
   baseline (mean) error rate instead of their signal-inflated HiFi-FP frequency, so the null no
   longer launders real recurrent variants into itself.
 Option 2 (GENOTYPE-CONDITIONED null): fix each position's SR-FP and HiFi-FP donor-set sizes,
   randomize WHICH donors carry the HiFi-FP (breaks genotype linkage). Excess co-incidence in the
   SAME donors = real. Analytic E[o_k]=s_k*h_k/N + permutation p; stratified by HiFi recurrence.
Variant and false-positive loading follow the standard procedure."""
import sqlite3, os, json, numpy as np
from collections import Counter, defaultdict
DB=os.environ.get("FIVES_DB","5S_rDNA.db"); BAM=os.environ.get("FIVES_DATA","data")+"/bam"
BOUNDARY=set(range(1,90))|set(range(2034,2169))
SR_MIN_AD=3; LR_AD_MIN=5; LR_VAF=0.003; COHORT="HPRC_Year1"
EXCL_ALL={"HG02486"}; EXCL_HIFI={"HG02818"}
cons="".join(json.load(open(os.environ.get("FIVES_DATA","data")+"/consensus_reference.json"))["consensus"]).upper().replace("-","")
nonbound=[p for p in range(1,2169) if p not in BOUNDARY]
universe=[(p,a) for p in nonbound for a in "ACGT" if a!=cons[p-1]]
uni_idx=np.arange(len(universe)); key2i={k:i for i,k in enumerate(universe)}
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
REAL=set().union(*gt.values())                       # (pos,alt) real GT in >=1 donor
obs=sum(len(sr_fp[s]&hifi_fp[s]) for s in donors)
obs_real=sum(len((sr_fp[s]&hifi_fp[s])&REAL) for s in donors)
print(f"donors={len(donors)}  observed overlap={obs}  (at real-GT positions: {obs_real}={100*obs_real/obs:.0f}%)")

pooled=[k for s in donors for k in hifi_fp[s]]; freq=Counter(pooled)
rng=np.random.default_rng(0); NP=2000
sr_i={s:np.array([key2i[k] for k in sr_fp[s]]) for s in donors}
sr_set={s:set(sr_i[s].tolist()) for s in donors}
def run(sample_fn,label):
    null=np.empty(NP)
    for j in range(NP):
        tot=0
        for s in donors:
            n=len(hifi_fp[s])
            if n: tot+=len(set(sample_fn(n).tolist())&sr_set[s])
        null[j]=tot
    enr=obs/max(null.mean(),1e-9); p=(np.sum(null>=obs)+1)/(NP+1)
    print(f"  [{label}] null mean={null.mean():.1f} sd={null.std():.1f}  enrichment={enr:.2f}x  p={p:.1e}")
    return enr

# context: the two bracketing nulls
real_idx=np.array([key2i[k] for k in freq]); real_wt=np.array([freq[k] for k in freq],float); real_wt/=real_wt.sum()
run(lambda n: rng.choice(uni_idx,size=n,replace=False),"UNIFORM (bracket: lenient)")
run(lambda n: rng.choice(real_idx,size=min(n,len(real_idx)),replace=False,p=real_wt),"MARGINAL (bracket: conservative)")

# ---- OPTION 1: error-matched positional null ----
nonreal_freqs=[freq[k] for k in freq if k not in REAL]
baseline=float(np.mean(nonreal_freqs)) if nonreal_freqs else 1.0
ek=list(freq.keys()); e_idx=np.array([key2i[k] for k in ek])
e_wt=np.array([(baseline if k in REAL else freq[k]) for k in ek],float); e_wt/=e_wt.sum()
print(f"  (error-matched: {sum(k in REAL for k in ek)}/{len(ek)} HiFi-FP positions are real-GT -> de-weighted to baseline {baseline:.2f})")
run(lambda n: rng.choice(e_idx,size=min(n,len(e_idx)),replace=False,p=e_wt),"OPTION 1: ERROR-MATCHED")

# ---- OPTION 2: genotype-conditioned null ----
didx={s:i for i,s in enumerate(donors)}; N=len(donors)
sr_don=defaultdict(set); hf_don=defaultdict(set)
for s in donors:
    for k in sr_fp[s]: sr_don[k].add(didx[s])
    for k in hifi_fp[s]: hf_don[k].add(didx[s])
both=[k for k in hf_don if k in sr_don]
o_k={k:len(sr_don[k]&hf_don[k]) for k in both}
obs2=sum(o_k.values()); exp2=sum(len(sr_don[k])*len(hf_don[k])/N for k in both)
# permutation: randomize which donors carry each position's HiFi-FP
NP2=2000; null2=np.empty(NP2)
sr_arr={k:np.fromiter(sr_don[k],int) for k in both}; hk={k:len(hf_don[k]) for k in both}
srmask={k:np.zeros(N,bool) for k in both}
for k in both: srmask[k][sr_arr[k]]=True
for j in range(NP2):
    tot=0
    for k in both:
        d=rng.choice(N,size=hk[k],replace=False)
        tot+=int(srmask[k][d].sum())
    null2[j]=tot
p2=(np.sum(null2>=obs2)+1)/(NP2+1)
print(f"  [OPTION 2: GENOTYPE-CONDITIONED] obs={obs2} analytic-exp={exp2:.1f} -> {obs2/exp2:.2f}x ; "
      f"perm null mean={null2.mean():.1f} sd={null2.std():.1f} -> {obs2/null2.mean():.2f}x  p={p2:.1e}")
# recurrence stratification by HiFi-FP incidence h_k
print("  recurrence-stratified (by HiFi-FP donor incidence h_k): obs/exp enrichment")
for lab,lo,hi in [("h=1",1,1),("h=2-4",2,4),("h>=5",5,10**9)]:
    ks=[k for k in both if lo<=hk[k]<=hi]
    o=sum(o_k[k] for k in ks); e=sum(len(sr_don[k])*len(hf_don[k])/N for k in ks)
    print(f"     {lab:6s} n_pos={len(ks):4d}  obs={o:4d} exp={e:6.1f}  enrichment={o/e if e else float('nan'):.2f}x")
