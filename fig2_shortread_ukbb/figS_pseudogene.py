#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# figS_pseudogene.py — supplementary figure assessing dispersed 5S pseudogene
# read contamination of short-read 5S variant calling.
# Author: Jonathan Bohlen
# (c) 2026 Jonathan Bohlen. Code accompanying Sengl et al. (2026),
# "Surveillance and selection of 5S ribosomal RNA genes in the human genome."
# Released under the MIT License; see LICENSE at the repository root.
# -----------------------------------------------------------------------------
"""Supplementary figure assessing dispersed 5S pseudogene read contamination of short-read 5S variant calling.
(A) identity of the 339 dispersed loci; (B) T2T ground-truth read-level confusion matrix (reads simulated
from CHM13's real uncollapsed array + real pseudogene loci, labelled by origin, aligned to GRCh38);
(C) 150-mer alignability track across the array + flank; (D) bleedthrough competitive vs naive;
(E) bidirectional leakage is robust to read substitutions (0-12 mismatches/read)."""
import os, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

SCRATCH=os.environ.get("FIVES_DATA","data")+"/pseudo"
OUT=os.environ.get("FIVES_OUT","output")+"/Supplementary_Tables"
ident=np.array([float(l) for l in open(SCRATCH+"/disp_identities.txt")])

plt.rcParams.update({"font.family":"Arial","pdf.fonttype":42,"font.size":8,"axes.linewidth":0.8})
CM=1/2.54
fig=plt.figure(figsize=(17*CM,19*CM))
gs=GridSpec(3,2,figure=fig,hspace=0.62,wspace=0.42,height_ratios=[1,0.9,1])
axA=fig.add_subplot(gs[0,0]); axB=fig.add_subplot(gs[0,1])
axC=fig.add_subplot(gs[1,:])
axD=fig.add_subplot(gs[2,0]); axE=fig.add_subplot(gs[2,1])

# ---- A: identity distribution of dispersed loci ----
ax=axA
ax.hist(ident,bins=np.arange(75,101,1),color="#4393c3",edgecolor="white",lw=0.4)
ax.axvline(95,color="#b2182b",ls="--",lw=0.9); ax.axvline(99,color="#7a0177",ls=":",lw=0.9)
ax.text(94.5,ax.get_ylim()[1]*0.95,f"{int((ident>=95).sum())} loci >=95%",color="#b2182b",fontsize=6.5,va="top",ha="right")
ax.text(99.3,ax.get_ylim()[1]*0.6,f"{int((ident>=99).sum())} >=99%",color="#7a0177",fontsize=6.5,va="top")
ax.set_xlabel("% identity of dispersed 5S locus to array gene"); ax.set_ylabel("number of dispersed loci")
ax.set_title(f"A   {len(ident)} dispersed 5S loci genome-wide",fontsize=8,loc="left",fontweight="bold")
for s in ["top","right"]: ax.spines[s].set_visible(False)

# ---- B: T2T ground-truth confusion matrix (reads from the extracted 5S gene region) ----
# rows = read origin (CHM13 truth), cols = where they align in GRCh38
ax=axB
dests=["5S array","dispersed 5S\npseudogene","other\nlocus","un-\nmapped"]
M=np.array([[399956,      0,     0, 0],       # ARRAY gene-ROI origin
            [     1, 383711, 16382, 0]],dtype=float)  # PSEUDO origin
rowtot=M.sum(1,keepdims=True); P=M/rowtot*100
im=ax.imshow(P,cmap="Blues",vmin=0,vmax=100,aspect="auto")
ax.set_xticks(range(4)); ax.set_xticklabels(dests,fontsize=5.8)
ax.set_yticks([0,1]); ax.set_yticklabels(["5S array\n(gene ROI,\nn=400k)","pseudogene\n(n=400k)"],fontsize=6.3)
for i in range(2):
    for j in range(4):
        v=P[i,j]
        if M[i,j]==0: txt="0"
        elif v<0.05:  txt=f"{v:.4f}%"      # tiny nonzero -> percent, not a raw count
        else:         txt=f"{v:.1f}%"
        ax.text(j,i,txt,ha="center",va="center",fontsize=6,
                color="white" if v>55 else "#222")
# outline the two decision-critical cells
import matplotlib.patches as mp
ax.add_patch(mp.Rectangle((-0.5,-0.5),1,1,fill=False,ec="#238b45",lw=1.6))   # array->array (correct)
ax.add_patch(mp.Rectangle((-0.5, 0.5),1,1,fill=False,ec="#b2182b",lw=1.6))   # pseudo->array (contamination)
ax.set_xlabel("aligns in GRCh38 to  ->",fontsize=6.5); ax.set_ylabel("read origin (CHM13 T2T truth)",fontsize=6.5)
ax.set_title("B   Ground-truth confusion matrix",fontsize=8,loc="left",fontweight="bold")
ax.text(1.5,2.12,"0 array reads -> pseudogene;\n0.0003% of pseudogene reads -> array",
        fontsize=5.4,ha="center",va="top",color="#444")

# ---- C: 150-mer alignability track ----
ax=axC
d=np.genfromtxt(SCRATCH+"/mappability.tsv",names=True,dtype=None,encoding=None)
WLO,WHI=228610069,228646359
x=(d["pos"]-WLO)/1000.0    # kb relative to window start
ax.fill_between(x,0,d["mappability"],step="mid",color="#9ecae1",lw=0)
ax.plot(x,d["mappability"],lw=0.5,color="#2166ac",drawstyle="steps-mid")
ax.axvspan((WLO-WLO)/1000,(WHI-WLO)/1000,color="#fde0dd",lw=0,zorder=0,alpha=0.6)
ax.text((WHI-WLO)/2000,1.06,"5S array extraction window (collapsed repeat)",fontsize=6,ha="center",color="#b2182b")
ax.set_ylim(0,1.15); ax.set_xlim(x.min(),x.max())
ax.set_xlabel("position along chr1 relative to array-window start (kb)")
ax.set_ylabel("150-mer alignability\n(1 / # best genomic hits)")
ax.axhline(1.0,ls=":",lw=0.6,color="#888")
ax.text((WHI-WLO)/4000,0.30,"array body:\nmedian 1/14\n(multi-mapping)",fontsize=6,color="#08519c",ha="center")
ax.text(x.max()*0.72,0.55,"unique flank:\n99% single-copy",fontsize=6,color="#238b45",ha="center")
ax.set_title("C   Array is multi-mapping, flanks are unique (GRCh38 150-mer alignability)",fontsize=8,loc="left",fontweight="bold")
for s in ["top","right"]: ax.spines[s].set_visible(False)

# ---- D: bleedthrough competitive vs naive ----
ax=axD
ax.bar([0,1],[0.0,24.6],width=0.55,color=["#238b45","#b2182b"],edgecolor="white")
ax.text(0,1.3,"0%",ha="center",fontsize=8,color="#238b45",fontweight="bold")
ax.text(1,25.4,"24.6%",ha="center",fontsize=8,color="#b2182b",fontweight="bold")
ax.set_xticks([0,1]); ax.set_xticklabels(["competitive\nalignment\n(this study)","naive\nsequence-based"],fontsize=6.5)
ax.set_ylabel("% dispersed reads reaching array pool"); ax.set_ylim(0,29)
ax.set_title("D   Bleedthrough into variant calling",fontsize=8,loc="left",fontweight="bold")
for s in ["top","right"]: ax.spines[s].set_visible(False)

# ---- E: bidirectional leakage vs substitutions ----
ax=axE
subs=[0,1,2,3,5,8,12]
a2p=[0,0,0,0,0,0.004,0.091]
p2a=[0,0,0,0,0,0,0.00025]
ax.axvspan(0.5,2.5,color="#fff2cc",lw=0,zorder=0); ax.text(1.5,4e-4,"typical real\nvariant load",fontsize=5.8,color="#b58900",ha="center",va="bottom")
ax.plot(subs,a2p,"-o",color="#b2182b",ms=3,lw=1.3,label="array -> pseudogene")
ax.plot(subs,p2a,"-s",color="#2166ac",ms=3,lw=1.3,label="pseudogene -> array")
ax.axhline(1e-4,ls=":",lw=0.6,color="#888"); ax.text(12.4,1.15e-4,"~single-read limit",fontsize=5,color="#888",ha="right",va="bottom")
ax.set_yscale("symlog",linthresh=1e-4)
ax.set_yticks([0,1e-4,1e-3,1e-2,1e-1])
ax.set_xlabel("substitutions per 150 bp read"); ax.set_ylabel("% reads mis-mapped to other 5S locus")
ax.set_ylim(-1.3e-4,0.2); ax.set_xlim(-0.3,12.5)   # bottom below 0 so the zero points sit clear of the axis
ax.legend(fontsize=6.5,frameon=False,loc="upper left")
ax.set_title("E   Leakage stays ~0 with mismatches",fontsize=8,loc="left",fontweight="bold")
for s in ["top","right"]: ax.spines[s].set_visible(False)

fig.savefig(os.path.join(OUT,"FigS_pseudogene_contamination.pdf"),bbox_inches="tight",dpi=300)
fig.savefig(os.path.join(OUT,"FigS_pseudogene_contamination.png"),bbox_inches="tight",dpi=200)
print("wrote FigS_pseudogene_contamination.pdf (5 panels)")
