#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# 20_gtex_rnaseq_stream.py — stream GTEx RNA-seq per sample, realign the RNA5S
# loci (+optional unmapped reads) to the 5S consensus, and write one pileup + QC
# per (donor, tissue, sample).
# Author: Jonathan Bohlen
# (c) 2026 Jonathan Bohlen. Code accompanying Sengl et al. (2026),
# "Surveillance and selection of 5S ribosomal RNA genes in the human genome."
# Released under the MIT License; see LICENSE at the repository root.
# -----------------------------------------------------------------------------
"""Stream GTEx RNA-seq: per-SAMPLE slice (RNA5S loci [+ unmapped if DEEP]) -> realign to 5S
consensus -> one pileup per (donor,tissue,sample) + QC. Tissue provenance preserved (the atomic
unit is the sample = one tissue); pooling/calling is post-streaming.
Resumable. Env: DONORS=comma-list (default all RNASEQ donors), DEEP=1 (loci+unmapped) / 0 (loci-only),
MAX_PER_DONOR=N (0=all)."""
import os, subprocess, csv, re
from pathlib import Path
PROJ=os.environ.get("GTEX_GCP_PROJECT","")                       # requester-pays GCP project id
RNAB=os.environ.get("GTEX_RNA_BUCKET","gs://BUCKET/GTEx_Analysis_2025-08-22_v11_RNAseq_BAM_files")
CONS=str(Path(os.environ.get("FIVES_REFS","refs"))/"5S_t2t_consensus.fa")
BED=str(Path(os.environ.get("FIVES_REFS","refs"))/"rna5s_all_loci.bed")
META=str(Path(os.environ.get("FIVES_DATA","data"))/"metadata"/"GTEx_Analysis_2025-08-22_v11_Annotations_SampleAttributesDS.txt")
OUT=str(Path(os.environ.get("FIVES_DATA","data"))/"results"/"rnaseq"/"cohort")
DEEP=os.environ.get("DEEP","1")=="1"; MAXP=int(os.environ.get("MAX_PER_DONOR","0"))
DONORS=[d for d in os.environ.get("DONORS","").split(",") if d]
THREADS=6; GENE=(630,748); os.makedirs(OUT,exist_ok=True)
CN=open(CONS).readline().strip().lstrip(">").split()[0]
def sh(c,env=None): subprocess.run(c,shell=True,executable="/bin/bash",env=env)
def o(c,env=None): return subprocess.run(c,shell=True,executable="/bin/bash",capture_output=True,text=True,env=env).stdout
def token(): return o("gcloud auth print-access-token").strip()
def slug(s): return re.sub(r'[^A-Za-z0-9]+','_',s).strip('_')
def gb(bam,env): return o(f"samtools depth -a -r {CN}:{GENE[0]}-{GENE[1]} {bam} 2>/dev/null | awk '{{s+=$3}}END{{if(NR)printf \"%.1f\",s/NR; else print 0}}'",env).strip()

perdonor={}
for r in csv.DictReader(open(META),delimiter='\t'):
    if r['SMAFRZE']!='RNASEQ': continue
    d='-'.join(r['SAMPID'].split('-')[:2])
    if DONORS and d not in DONORS: continue
    perdonor.setdefault(d,[]).append((r['SAMPID'],r['SMTSD']))
samples=[]
for d,lst in perdonor.items():
    lst=sorted(lst); lst=lst[:MAXP] if MAXP else lst
    samples+=[(d,s,t) for s,t in lst]
print(f"{len(samples)} samples / {len(perdonor)} donors; DEEP={DEEP}",flush=True)

for i,(d,sampid,tissue) in enumerate(samples,1):
    base=f"{OUT}/{d}.{slug(tissue)}.{sampid}"; pg=base+".pileup.tsv.gz"; qc=base+".qc.tsv"
    if os.path.exists(pg) and os.path.exists(qc): print(f"[{i}/{len(samples)}] {sampid} skip",flush=True); continue
    env=dict(os.environ,GCS_OAUTH_TOKEN=token(),GCS_REQUESTER_PAYS_PROJECT=PROJ)
    BAM=f"{RNAB}/{sampid}.v11.Aligned.sortedByCoord.out.patched.md.bam"
    wd=f"{OUT}/_tmp_{sampid}"; os.makedirs(wd,exist_ok=True)
    try:
        sh(f'samtools view -b -L {BED} "{BAM}" 2>/dev/null | samtools collate -u -O - 2>/dev/null | samtools fastq -n 2>/dev/null > {wd}/loci.fq',env)
        sh(f'bwa mem -B 6 -O 8 -L 5,5 -T 30 -t {THREADS} {CONS} {wd}/loci.fq 2>/dev/null | samtools view -b -q 30 -F 0x904 | samtools sort -o {wd}/loci.bam 2>/dev/null',env)
        sh(f"samtools index {wd}/loci.bam",env)
        loci_reads=int((o(f'wc -l < {wd}/loci.fq').strip() or 0))//4; gb_loci=gb(f"{wd}/loci.bam",env)
        unmap_reads=0; allbam=f"{wd}/loci.bam"
        if DEEP:
            sh(f"""samtools view -b "{BAM}" '*' 2>/dev/null | samtools collate -u -O - 2>/dev/null | samtools fastq -n 2>/dev/null > {wd}/un.fq""",env)
            unmap_reads=int((o(f'wc -l < {wd}/un.fq').strip() or 0))//4
            sh(f'bwa mem -B 6 -O 8 -L 5,5 -T 30 -t {THREADS} {CONS} {wd}/un.fq 2>/dev/null | samtools view -b -q 30 -F 0x904 | samtools sort -o {wd}/un.bam 2>/dev/null',env)
            sh(f"samtools merge -f {wd}/all.bam {wd}/loci.bam {wd}/un.bam 2>/dev/null; samtools index {wd}/all.bam",env); allbam=f"{wd}/all.bam"
        gb_all=gb(allbam,env)
        tot=int((o(f'samtools idxstats "{BAM}" 2>/dev/null | awk \'{{m+=$3}}END{{print m}}\'',env).strip() or 0))
        cons_reads=int((o(f"samtools view -c {allbam}").strip() or 0))
        rpm=round(cons_reads/tot*1e6,2) if tot else ""
        sh(f'bcftools mpileup -Ou -d 200000 -q 30 -Q 30 -a FORMAT/AD,FORMAT/DP -f {CONS} {allbam} 2>/dev/null | bcftools query -f %POS"\\t"%REF"\\t"%ALT"\\t["%DP"]\\t["%AD"]\\n" | gzip -c > {pg}',env)
        with open(qc,"w") as f:
            f.write("donor_id\ttissue\tsample_id\ttotal_reads\tloci_reads\tunmap_reads\tcons_reads\trpm\tgb_depth_loci\tgb_depth_all\tdeep\n")
            f.write(f"{d}\t{tissue}\t{sampid}\t{tot}\t{loci_reads}\t{unmap_reads}\t{cons_reads}\t{rpm}\t{gb_loci}\t{gb_all}\t{int(DEEP)}\n")
        print(f"[{i}/{len(samples)}] {sampid} {slug(tissue)}: gb_loci={gb_loci} gb_all={gb_all} rpm={rpm}",flush=True)
    except Exception as e:
        print(f"[{i}/{len(samples)}] {sampid} ERROR {e}",flush=True)
    finally:
        sh(f"rm -rf {wd}")
print("RNA_STREAM_DONE",flush=True)
