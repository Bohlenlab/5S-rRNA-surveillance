#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# fig4_panel4_build.py — per-donor prevalence panel over the AUC carrier-HIGH 5S
# loci: histogram of expressed-variants-per-carrier vs a per-locus chance null.
# Author: Jonathan Bohlen
# (c) 2026 Jonathan Bohlen. Code accompanying Sengl et al. (2026),
# "Surveillance and selection of 5S ribosomal RNA genes in the human genome."
# Released under the MIT License; see LICENSE at the repository root.
# -----------------------------------------------------------------------------
"""AUC-based per-donor prevalence panel (threshold-free framing).
Expressed set = AUC carrier-HIGH loci (FDR<0.10). Per-donor call calibrated against
the non-carrier null at each locus (RNA-VAF > non-carrier Pth percentile, AD>=3); default P=q99
(~1% FPR). Set P4_PCTL env to test other cutoffs. Histogram of #expressed variants per carrier
(observed vs chance-null), + fraction expressing >=1 with bootstrap 95% CI and chance correction."""
import os,numpy as np,pandas as pd,matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
plt.rcParams.update({"pdf.fonttype":42,"ps.fonttype":42,"font.family":"Arial",
    "axes.linewidth":0.6,"xtick.major.width":0.6,"ytick.major.width":0.6,
    "xtick.major.size":2,"ytick.major.size":2})
CM=1/2.54; RED="#d62728"; GREY="#9b9b9b"; DGREY="#4d4d4d"
DAT=str(Path(os.environ.get("FIVES_DATA","data"))/"figures"/"figure4_rnaseq"/"data")
OUT=str(Path(os.environ.get("FIVES_OUT","output"))/"Figure4")
df=pd.read_csv(f"{DAT}/panel4_perdonor_carrierHIGH_vaf.tsv",sep="\t")
PCTL=float(os.environ.get("P4_PCTL","99")); SUF="" if PCTL==99 else f"_q{PCTL:g}".replace(".","p")  # q99 = canonical
MIN_AD=3; rng=np.random.default_rng(7); NSIM=2000; NBOOT=2000; MAXBIN=4

# per-locus non-carrier calibration: threshold = nonc PCTL-th pctile, FPR = realized nonc rate above it
loci=sorted(df.variant.unique()); thr={}; fpr={}
for v in loci:
    nc=df[(df.variant==v)&(df.group=="noncarrier")].rna_vaf.values
    t=np.quantile(nc,PCTL/100) if len(nc) else np.inf
    thr[v]=t; fpr[v]=float(np.mean(nc>t)) if len(nc) else 0.0
# carriers: per donor list of carried loci + detection
car=df[df.group=="carrier"].copy()
car["det"]=(car.rna_vaf>car.variant.map(thr))&(car.pooled_AD>=MIN_AD)
g=car.groupby("donor")
donors=sorted(df.donor.unique())   # ALL assessable donors (non-carriers of every locus -> 0 expressed)
obs=g.det.sum().reindex(donors).fillna(0).astype(int)          # # expressed per donor
carried_fpr={d:[fpr[v] for v in car[car.donor==d].variant] for d in donors}
Ncar=len(donors)
def histvec(counts):
    h=np.zeros(MAXBIN+1)
    for c in counts:
        h[min(int(c),MAXBIN)]+=1
    return 100*h/len(counts)
obs_h=histvec(obs.values)
obs_ge1=float(np.mean(obs.values>=1))
# chance-null: per carrier draw Bernoulli(fpr) over carried loci
null_ge1=[]; null_hs=[]
for _ in range(NSIM):
    cc=[int(np.sum(rng.random(len(f))<f)) for f in carried_fpr.values()]
    null_hs.append(histvec(cc)); null_ge1.append(np.mean(np.array(cc)>=1))
null_h=np.mean(null_hs,axis=0); null_ge1=np.array(null_ge1)
null_ge1_mean=null_ge1.mean(); null_ge1_lo,null_ge1_hi=np.percentile(null_ge1,[2.5,97.5])
# bootstrap CI on observed >=1
ov=obs.values
boot=[np.mean(rng.choice(ov,Ncar,replace=True)>=1) for _ in range(NBOOT)]
obs_lo,obs_hi=np.percentile(boot,[2.5,97.5])
attr=obs_ge1-null_ge1_mean
attr_lo=obs_lo-null_ge1_hi; attr_hi=obs_hi-null_ge1_lo
# all-assessable-donor denominator (donors assessable at >=1 carrier-HIGH locus)
n_assess=df.donor.nunique()
print(f"assessable donors (denominator): {Ncar} | of which carry >=1 carrier-HIGH: {(obs.index.isin(car.donor)).sum()}")
print(f"observed express >=1: {100*obs_ge1:.1f}% (95% CI {100*obs_lo:.1f}-{100*obs_hi:.1f})")
print(f"chance-null >=1: {100*null_ge1_mean:.1f}% (CI {100*null_ge1_lo:.1f}-{100*null_ge1_hi:.1f})")
print(f"genotype-attributable: {100*attr:.1f}% (CI {100*attr_lo:.1f}-{100*attr_hi:.1f})")
print(f"as % of all assessable donors: {100*obs.values.sum()/0+0 if False else 100*np.sum(obs.values>=1)/n_assess:.1f}% express >=1")

# ---- plot ----
fig,ax=plt.subplots(figsize=(4.8*CM,4.4*CM))
x=np.arange(MAXBIN+1); w=0.4
ax.bar(x-w/2,obs_h,w,color=RED,label="observed",zorder=3)
ax.bar(x+w/2,null_h,w,color="none",edgecolor=GREY,lw=0.8,label="chance",zorder=3)
ax.set_xticks(x); ax.set_xticklabels([str(i) for i in range(MAXBIN)]+[f"{MAXBIN}+"],fontsize=5.5)
ax.set_xlabel("Expressed 5S variants per donor",fontsize=6)
ax.set_ylabel("Donors (%)",fontsize=6)
ax.tick_params(labelsize=5.5,length=2,width=0.6)
ax.spines[["top","right"]].set_visible(False)
ax.legend(fontsize=4.8,frameon=False,loc="upper right",handletextpad=0.4,borderpad=0.1)
ax.text(.42,.80,f"call: VAF > nonc P{int(PCTL)}\nexpress ≥1: {100*obs_ge1:.0f}% (CI {100*obs_lo:.0f}–{100*obs_hi:.0f})\n"
    f"chance {100*null_ge1_mean:.0f}%\nattributable {100*attr:.0f}% (CI {100*attr_lo:.0f}–{100*attr_hi:.0f})",
    transform=ax.transAxes,fontsize=4.6,va="top",ha="left")
fig.savefig(f"{OUT}/P4_auc_prevalence_histogram{SUF}.pdf",dpi=400,bbox_inches="tight")
fig.savefig(f"{OUT}/P4_auc_prevalence_histogram{SUF}.png",dpi=300,bbox_inches="tight"); plt.close(fig)

# ---- data tables ----
pd.DataFrame({"n_expressed_per_carrier":[str(i) for i in range(MAXBIN)]+[f"{MAXBIN}+"],
    "observed_pct":obs_h,"chance_null_pct":null_h}).to_csv(
    f"{OUT}/P4_auc_prevalence_histogram{SUF}.tsv",sep="\t",index=False)
pd.DataFrame([dict(variant=v,nonc_thr=thr[v],nonc_FPR=fpr[v],
    n_carrier=int((car.variant==v).sum()),n_carrier_detected=int(car[(car.variant==v)].det.sum()))
    for v in loci]).to_csv(f"{OUT}/P4_perlocus_calibration{SUF}.tsv",sep="\t",index=False)
with open(f"{OUT}/P4_estimate{SUF}.txt","w") as fh:
    fh.write(f"Figure 4 Panel 4 — AUC-based prevalence (threshold-free framing), call = nonc P{int(PCTL)}\n")
    fh.write(f"Set: {len(loci)} AUC carrier-HIGH loci (rank-skew FDR<0.10). Per-donor call: pooled RNA-VAF\n")
    fh.write(f"> non-carrier {int(PCTL)}th percentile at that locus (AD>={MIN_AD}); ~{100-PCTL:.0f}% FPR by construction.\n\n")
    fh.write(f"Denominator = ALL assessable donors: {Ncar}  (of which {int((obs.index.isin(car.donor)).sum())} carry >=1 carrier-HIGH variant)\n")
    fh.write(f"Express >=1 variant : {100*obs_ge1:.1f}%  (bootstrap 95% CI {100*obs_lo:.1f}-{100*obs_hi:.1f})\n")
    fh.write(f"Chance-null         : {100*null_ge1_mean:.1f}%  (95% CI {100*null_ge1_lo:.1f}-{100*null_ge1_hi:.1f})\n")
    fh.write(f"Genotype-attributable: {100*attr:.1f}%  (95% CI {100*attr_lo:.1f}-{100*attr_hi:.1f})\n")
print("wrote P4 ->",OUT)
