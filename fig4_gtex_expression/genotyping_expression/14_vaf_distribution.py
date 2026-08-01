#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# 14_vaf_distribution.py — plot the VAF distribution of called WGS 5S variants
# across the cohort, overall and stratified by repeat-unit region.
# Author: Jonathan Bohlen
# (c) 2026 Jonathan Bohlen. Code accompanying Sengl et al. (2026),
# "Surveillance and selection of 5S ribosomal RNA genes in the human genome."
# Released under the MIT License; see LICENSE at the repository root.
# -----------------------------------------------------------------------------
"""VAF distribution of called 5S variants across the WGS cohort, overall and by region.
Reads the per-donor *.variants.tsv variant call tables."""
import os
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt, numpy as np, glob, csv
from pathlib import Path
CD=str(Path(os.environ.get("FIVES_DATA","data"))/"results"/"wgs"/"cohort"); FIG=str(Path(os.environ.get("FIVES_OUT","output"))/"figures")
vafs={'gene':[],'nts_pre':[],'nts_post':[]}
for f in glob.glob(CD+"/*.variants.tsv"):
    for v in csv.DictReader(open(f),delimiter='\t'): vafs[v['region']].append(float(v['vaf']))
allv=np.array(vafs['gene']+vafs['nts_pre']+vafs['nts_post'])
print(f"total calls: {len(allv)}")
for r in vafs:
    a=np.array(vafs[r]); print(f"  {r}: n={len(a)} medianVAF={np.median(a)*100:.2f}% <1%:{(a<0.01).mean()*100:.0f}% >35%:{(a>0.35).mean()*100:.0f}%")
bins=np.linspace(0,1,51); fig,ax=plt.subplots(1,2,figsize=(13,4.8))
ax[0].hist(allv,bins=bins,color='#4477aa',edgecolor='w'); ax[0].set_yscale('log')
ax[0].axvline(0.003,color='g',ls=':',label='VAF floor 0.3%')
ax[0].set_xlabel('VAF'); ax[0].set_ylabel('# called variants (log)'); ax[0].set_title(f'All called variants (n={len(allv)})'); ax[0].legend()
cols={'gene':'#cc3311','nts_pre':'#4477aa','nts_post':'#999999'}
ax[1].hist([vafs['nts_post'],vafs['nts_pre'],vafs['gene']],bins=bins,stacked=True,
           color=[cols['nts_post'],cols['nts_pre'],cols['gene']],label=['nts_post','nts_pre','gene'])
ax[1].set_yscale('log'); ax[1].set_xlabel('VAF'); ax[1].set_ylabel('# called variants (log)'); ax[1].set_title('By region (stacked)'); ax[1].legend()
fig.suptitle('GTEx WGS cohort — VAF distribution of called 5S variants (AD>=5, VAF>=0.3%)',fontweight='bold')
plt.tight_layout(); plt.savefig(FIG+'/05_vaf_distribution.png',dpi=150); plt.savefig(FIG+'/05_vaf_distribution.pdf')
print("saved 05_vaf_distribution")
