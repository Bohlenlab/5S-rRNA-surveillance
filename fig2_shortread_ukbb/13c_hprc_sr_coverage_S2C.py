#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# 13c_hprc_sr_coverage_S2C.py — per-position short-read depth across the 5S rDNA
# repeat unit (HPRC Illumina WGS): median with inter-sample IQR ribbon.
# Author: Jonathan Bohlen
# (c) 2026 Jonathan Bohlen. Code accompanying Sengl et al. (2026),
# "Surveillance and selection of 5S ribosomal RNA genes in the human genome."
# Released under the MIT License; see LICENSE at the repository root.
# -----------------------------------------------------------------------------
"""Short-read depth across the 5S rDNA repeat unit (HPRC Year 1).

Illumina WGS reads extracted from the chr1q42 5S array are re-aligned to a single
2168 bp population-consensus repeat unit; depth is computed per position for every
HPRC Year-1 sample (samtools depth -a). Shows per-sample profiles + median with the
inter-sample IQR ribbon."""
import os, subprocess, glob, numpy as np
from pathlib import Path
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

ST=os.environ.get("SAMTOOLS","samtools")
HPRC=Path(os.environ.get("FIVES_DATA","data"))
OUT=Path(os.environ.get("FIVES_OUT","output"))/"SupplementaryFigure2"
UNIT=2168; GENE=(630,748); EXCL=(89,2034)   # boundary-exclusion edges kept out of the array caller

bams=sorted(glob.glob(str(HPRC/"bam/*/*_illumina.sorted.bam")))
mat=[]
for b in bams:
    if not Path(b+".bai").exists(): subprocess.run([ST,"index",b])
    r=subprocess.run([ST,"depth","-a","-d","0",b],capture_output=True,text=True)
    dp={}
    for l in r.stdout.splitlines():
        p=l.split("\t")
        if len(p)>=3:
            try: dp[int(p[1])]=int(p[2])
            except: pass
    mat.append(np.array([dp.get(i,np.nan) for i in range(1,UNIT+1)],float))
mat=np.stack(mat); n=len(bams)
pos=np.arange(1,UNIT+1)
med=np.nanmedian(mat,axis=0)
p25=np.nanpercentile(mat,25,axis=0); p75=np.nanpercentile(mat,75,axis=0)
p5 =np.nanpercentile(mat, 5,axis=0); p95=np.nanpercentile(mat,95,axis=0)
core=slice(EXCL[0],EXCL[1])
core_med=np.nanmedian(med[core])
print(f"n={n}  core median depth {core_med:.0f}x  (CV {np.nanstd(med[core])/np.nanmean(med[core]):.3f})")

# ---- profile table ----
with open(OUT/"S2C_read_depth_profile.tsv","w") as o:
    o.write("unit_pos\tmedian_dp\tp5_dp\tp25_dp\tp75_dp\tp95_dp\n")
    for i in range(UNIT):
        o.write(f"{i+1}\t{med[i]:.0f}\t{p5[i]:.0f}\t{p25[i]:.0f}\t{p75[i]:.0f}\t{p95[i]:.0f}\n")

# ---- figure ----
plt.rcParams.update({"font.family":"Arial","pdf.fonttype":42,"font.size":8,"axes.linewidth":0.8})
CM=1/2.54
fig,ax=plt.subplots(figsize=(11*CM,4.3*CM))   # 30% flatter than 6.2cm
BLUE="#2166ac"
ax.fill_between(pos,p25,p75,color=BLUE,alpha=0.22,lw=0,zorder=2,label="inter-sample IQR")
ax.plot(pos,med,color=BLUE,lw=1.1,zorder=3,label=f"median (n={n})")
ax.axhline(core_med,color=BLUE,ls="--",lw=0.6,alpha=0.5)
# gene + boundary-exclusion shading
ax.axvspan(*GENE,color="#74c476",alpha=0.18,lw=0,zorder=0)
ax.axvspan(1,EXCL[0],color="grey",alpha=0.12,lw=0,zorder=0)
ax.axvspan(EXCL[1],UNIT,color="grey",alpha=0.12,lw=0,zorder=0)
ax.text(sum(GENE)/2,ax.get_ylim()[1] if False else 7600,"5S gene",fontsize=6,ha="center",color="#2d8f4e")
ax.text(45,300,"excl.",fontsize=5.5,ha="center",color="grey")
ax.text(2101,300,"excl.",fontsize=5.5,ha="center",color="grey")
ax.set_xlim(1,UNIT); ax.set_ylim(0,8000)
ax.text(0.985,0.06,f"core median {core_med:.0f}x",transform=ax.transAxes,ha="right",fontsize=6.5,color=BLUE)
ax.set_xlabel("position in 5S rDNA repeat unit (bp)")
ax.set_ylabel("short-read depth")
ax.set_title("Even coverage of the 5S consensus (HPRC Illumina WGS)",fontsize=8,loc="left",fontweight="bold")
ax.legend(fontsize=6,frameon=False,loc="upper left")
for s in ["top","right"]: ax.spines[s].set_visible(False)
fig.savefig(OUT/"S2C_read_depth_profile_HPRC.pdf",bbox_inches="tight",dpi=300)
fig.savefig(OUT/"S2C_read_depth_profile_HPRC.png",bbox_inches="tight",dpi=200)
print("wrote S2C_read_depth_profile_HPRC.pdf/png + S2C_read_depth_profile.tsv")
