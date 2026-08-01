#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# 10_gtex_wgs_cohort.py — stream GTEx WGS CRAMs, realign the 5S cluster slice to
# the 5S consensus, and write per-donor variant + copy-number QC TSVs.
# Author: Jonathan Bohlen
# (c) 2026 Jonathan Bohlen. Code accompanying Sengl et al. (2026),
# "Surveillance and selection of 5S ribosomal RNA genes in the human genome."
# Released under the MIT License; see LICENSE at the repository root.
# -----------------------------------------------------------------------------
"""Stream GTEx WGS CRAMs, call 5S variants from the hg38 RNA5S cluster slice,
write per-donor variant + QC TSVs. Resumable."""
import os, sys, subprocess, csv, statistics
from pathlib import Path
PROJ=os.environ.get("GTEX_GCP_PROJECT","")                       # requester-pays GCP project id
WGSB=os.environ.get("GTEX_WGS_BUCKET","gs://BUCKET/GTEx_Analysis_2021-02-11_v9_WGS_CRAM_files")
REF38=str(Path(os.environ.get("FIVES_REFS","refs"))/"Homo_sapiens_assembly38.fasta")
CONS=str(Path(os.environ.get("FIVES_REFS","refs"))/"5S_t2t_consensus.fa")
META=str(Path(os.environ.get("FIVES_DATA","data"))/"metadata"/"GTEx_Analysis_2025-08-22_v11_Annotations_SampleAttributesDS.txt")
CLUSTER="chr1:228605000-228650000"
CTRL="chr12:50000000-50100000"          # single-copy autosomal control window for CN
OUT=str(Path(os.environ.get("FIVES_DATA","data"))/"results"/"wgs"/"cohort")
N_DONORS=int(os.environ.get("N_DONORS","100"))
MIN_AD, MIN_VAF, MIN_CALLABLE_DP = 5, 0.003, 50
EXCL=set(range(1,90))|set(range(2034,2169))  # WGS: edges only (790-932/974-1057 are RNA/pseudogene artefacts, real in WGS cluster slice)
GENE=(630,748); THREADS=6
os.makedirs(OUT, exist_ok=True)
CN=open(CONS).readline().strip().lstrip(">").split()[0]

def sh(cmd, env=None): subprocess.run(cmd, shell=True, executable="/bin/bash", env=env)
def out(cmd, env=None): return subprocess.run(cmd, shell=True, executable="/bin/bash",
                                              capture_output=True, text=True, env=env).stdout
def token(): return out("gcloud auth print-access-token").strip()
def region(p): return "gene" if GENE[0]<=p<=GENE[1] else ("nts_pre" if p<GENE[0] else "nts_post")

def wgs_donors():
    s=set()
    with open(META) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if row["SMAFRZE"]=="WGS":
                s.add("-".join(row["SAMPID"].split("-")[:2]))
    return sorted(s)[:N_DONORS]

def process(d, env):
    vtsv=f"{OUT}/{d}.variants.tsv"; qtsv=f"{OUT}/{d}.qc.tsv"
    if os.path.exists(vtsv) and os.path.exists(qtsv): return "skip"
    cram=out(f'gsutil -u {PROJ} ls "{WGSB}/{d}-*.cram" 2>/dev/null | head -1', env).strip()
    if not cram: 
        open(qtsv,"w").write("donor_id\tstatus\n%s\tno_cram\n"%d); return "no_cram"
    wd=f"{OUT}/_tmp_{d}"; os.makedirs(wd, exist_ok=True)
    sh(f'samtools view -b -T {REF38} "{cram}" {CLUSTER} 2>/dev/null | samtools collate -u -O - 2>/dev/null '
       f'| samtools fastq -1 {wd}/r1.fq -2 {wd}/r2.fq -s {wd}/s.fq -0 /dev/null -n 2>/dev/null', env)
    rl=out(f"head -2 {wd}/r1.fq | tail -1 | tr -d '\\n' | wc -c").strip()
    nreads=(int(out(f"wc -l < {wd}/r1.fq") or 0)+int(out(f"wc -l < {wd}/s.fq") or 0))//4
    sh(f'( bwa mem -B 6 -O 8 -L 5,5 -T 30 -t {THREADS} {CONS} {wd}/r1.fq {wd}/r2.fq 2>/dev/null; '
       f'bwa mem -B 6 -O 8 -L 5,5 -T 30 -t {THREADS} {CONS} {wd}/s.fq 2>/dev/null | grep -v "^@" ) '
       f'| samtools view -b -q 30 -F 0x904 - | samtools sort -o {wd}/d.bam 2>/dev/null', env)
    sh(f"samtools index {wd}/d.bam", env)
    sh(f'bcftools mpileup -Ou -d 200000 -q 30 -Q 30 -a FORMAT/AD,FORMAT/DP -f {CONS} {wd}/d.bam 2>/dev/null '
       f'| bcftools query -f %POS"\\t"%REF"\\t"%ALT"\\t["%DP"]\\t["%AD"]\\n" > {wd}/pile.tsv', env)
    depths={}; calls=[]
    for ln in open(f"{wd}/pile.tsv"):
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
    creads=out(f'samtools view -c -T {REF38} "{cram}" {CTRL} 2>/dev/null', env).strip()
    rli=int(rl) if rl.isdigit() else 150
    ctrl=(int(creads or 0)*rli)/100000.0   # CTRL window = 100kb
    nonx=[v for k,v in depths.items() if k not in EXCL and 90<=k<=2033]
    med=statistics.median(nonx) if nonx else 0
    est_cn=round(med/ctrl,2) if ctrl else ""
    callable_frac=round(sum(1 for v in nonx if v>=MIN_CALLABLE_DP)/len(nonx),4) if nonx else 0
    with open(vtsv,"w") as f:
        f.write("donor_id\tconsensus_pos\tref\talt\tdepth\talt_depth\tvaf\tregion\n")
        for c in sorted(calls): f.write(d+"\t"+"\t".join(map(str,c))+"\n")
    with open(qtsv,"w") as f:
        f.write("donor_id\tmedian_depth\tcontrol_depth\test_copies\tn_variants\tcallable_fraction\tslice_reads\treadlen\tstatus\n")
        f.write(f"{d}\t{med}\t{round(ctrl,1)}\t{est_cn}\t{len(calls)}\t{callable_frac}\t{nreads}\t{rl}\tok\n")
    sh(f"gzip -c {wd}/pile.tsv > {OUT}/{d}.pileup.tsv.gz")
    sh(f"rm -rf {wd}")
    return f"ok med={med} cn={est_cn} nvar={len(calls)}"

def main():
    donors=wgs_donors()
    print(f"{len(donors)} donors; out={OUT}", flush=True)
    for i,d in enumerate(donors,1):
        env=dict(os.environ, GCS_OAUTH_TOKEN=token(), GCS_REQUESTER_PAYS_PROJECT=PROJ)
        try:
            r=process(d, env)
        except Exception as e:
            r=f"ERROR {e}"
        print(f"[{i}/{len(donors)}] {d}: {r}", flush=True)
    print("COHORT_RUN_DONE", flush=True)

if __name__=="__main__": main()
