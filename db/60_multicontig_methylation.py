#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# 60_multicontig_methylation.py — regenerate per-copy/per-position methylation exports for multi-contig haplotypes.
# Author: Jonathan Bohlen
# (c) 2026 Jonathan Bohlen. Code accompanying Sengl et al. (2026),
# "Surveillance and selection of 5S ribosomal RNA genes in the human genome."
# Released under the MIT License; see LICENSE at the repository root.
# -----------------------------------------------------------------------------
"""
60 — Methylation re-derivation for multi-contig-migrated haplotypes.

Regenerates per-copy and per-position methylation for migrated multi-contig
haplotypes from the cached 5S reads.

Per-hap reference:
  multi-contig hap : concatenation of its per-copy full units (mc_gene.fa,
                     copy_number order, consensus-framed gene~630) -> copies sit
                     at known cumulative offsets.
  single-contig hap: existing genomic array FASTA {sample}_{hap}_5S_array.fa
                     (used only as the competition partner; not re-output).

Competitive alignment (minimap2 -y -N1) tags reads by hap; modkit extract on the
hap reference; CpGs assigned to copy by offset; wpos -> region (gene 630-748,
ALU 787-1066), REPEAT_LEN 2168.

Outputs: methylation/copy_meth_export_mc[_hifi].tsv (per-copy) and
         methylation/copy_meth_pos_export_mc[_hifi].tsv (per-position),
keyed by (sample, hap_label, copy_number) for loading into SQLite.

Paths and the tool-binary directory are read from environment variables:
  FIVES_DATA  input/derived-data directory (sequences/, databases/, methylation/)
  FIVES_BIN   directory containing external tools (minimap2, samtools, modkit)

Usage: 60_multicontig_methylation.py --samples S1,S2 [--hifi] [--twostep]
"""
import argparse, os, subprocess, sys, tempfile, gzip
from pathlib import Path
import numpy as np, pandas as pd

T2T=Path(os.environ.get("FIVES_DATA", "data")); HPRC=T2T
SEQS=HPRC/"sequences"; DBMC=HPRC/"databases_mc"; DB=HPRC/"databases"; METH=HPRC/"methylation"
ENV=os.environ.get("FIVES_BIN", "")
MINIMAP2,SAMTOOLS,MODKIT=f"{ENV}/minimap2",f"{ENV}/samtools",f"{ENV}/modkit"
REPEAT_LEN=2168; GENE_S,GENE_E=630,748; ALU_S,ALU_E=787,1066; MOD_HI,MOD_LO=0.8,0.2

def run(cmd): subprocess.run(cmd,shell=True,check=True,executable="/bin/bash")
def region_of(p): return "gene" if GENE_S<=p<=GENE_E else ("nts_pre" if p<GENE_S else "nts_post")

def read_fasta(path):
    op=gzip.open if str(path).endswith(".gz") else open
    seqs={}; h=None; buf=[]
    with op(path,"rt") as fh:
        for ln in fh:
            if ln.startswith(">"):
                if h: seqs[h]="".join(buf)
                h=ln[1:].split()[0]; buf=[]
            else: buf.append(ln.strip())
    if h: seqs[h]="".join(buf)
    return seqs

def build_ref(sample,hap):
    """Return (ref_fa_path, header, copymap) ; copymap=list of (copy_number,start,end). None if unavailable."""
    tag=f"{sample}__{hap}__arr"
    mc=SEQS/sample/f"{sample}_{hap}_mc_gene.fa"
    out=METH/sample/f"{sample}_{hap}_mcmethref.fa"
    out.parent.mkdir(parents=True,exist_ok=True)
    if mc.exists():                       # multi-contig: concat units in copy order
        seqs=read_fasta(mc)               # headers copyNNN
        items=sorted(seqs.items(), key=lambda kv:int(kv[0].replace("copy","")))
        cur=0; cmap=[]; cat=[]
        for cid,s in items:
            n=int(cid.replace("copy","")); cmap.append((n,cur,cur+len(s))); cat.append(s); cur+=len(s)
        with open(out,"w") as f: f.write(f">{tag}\n{''.join(cat)}\n")
        run(f"'{SAMTOOLS}' faidx '{out}'")
        return out,tag,cmap
    arr=SEQS/f"{sample}_{hap}_5S_array.fa"  # single-contig: existing genomic array
    db=DB/f"{sample}_{hap}.tsv"
    if not arr.exists() or not db.exists(): return None
    # offset from header {contig}:{start}-{end}
    off=0
    with open(arr) as fh:
        hdr=fh.readline().strip().split()
        try: off=int(hdr[1].split(":")[1].split("-")[0])
        except Exception: off=0
    d=pd.read_csv(db,sep="\t"); d=d[~d["border_note"].astype(str).str.contains("border",na=False)]
    cmap=[(int(r.copy_id),int(r.unit_start_local-off),int(r.unit_end_local-off)) for _,r in d.iterrows()]
    # continuous genomic array reference
    s=read_fasta(arr); seq=list(s.values())[0]
    with open(out,"w") as f: f.write(f">{tag}\n{seq}\n")
    run(f"'{SAMTOOLS}' faidx '{out}'")
    return out,tag,cmap

TAGGER=r"""
import sys
def hap_of(rn):
    # header tag: {sample}__{hap}__arr
    try: return rn.split('__')[1]
    except: return 'NA'
def emit(buf):
    if not buf: return
    prim=None; ph=None; sec=set()
    for line in buf:
        f=line.split('\t'); fl=int(f[1])
        if fl&0x800: continue
        if fl&0x100: sec.add(hap_of(f[2]))
        else: prim=line; ph=hap_of(f[2])
    if prim is None: return
    yh='ambiguous' if (sec-{ph}) else ph
    sys.stdout.write(prim.rstrip('\n')+f'\tYH:Z:{yh}\n')
cur=None;buf=[]
for line in sys.stdin:
    if line[0]=='@': sys.stdout.write(line); continue
    q=line.split('\t',1)[0]
    if q!=cur: emit(buf); buf=[]; cur=q
    buf.append(line)
emit(buf)
"""

def build_compref(sample,hap):
    """Continuous competition reference (FAIR hap assignment): multi-contig hap ->
    saved genomic 5S contigs; single-contig -> existing _5S_array.fa. Headers
    retagged {sample}__{hap}__N so the tagger reads hap from split('__')[1]."""
    tag=f"{sample}__{hap}"; out=METH/sample/f"{sample}_{hap}_compref.fa"
    out.parent.mkdir(parents=True,exist_ok=True)
    contigs=SEQS/sample/f"{sample}_{hap}_5S_contigs.fa.gz"
    if contigs.exists():
        seqs=read_fasta(contigs)
        with open(out,"w") as f:
            for i,(nm,s) in enumerate(seqs.items()): f.write(f">{tag}__{i}\n{s}\n")
    else:
        arr=SEQS/f"{sample}_{hap}_5S_array.fa"
        if not arr.exists(): return None
        s=list(read_fasta(arr).values())[0]
        with open(out,"w") as f: f.write(f">{tag}__arr\n{s}\n")
    run(f"'{SAMTOOLS}' faidx '{out}'"); return out

def _assign(m,cmap,sample,h):
    """Vectorized copy/region assignment from a modkit df -> (copy_df, pos_df)."""
    cmap=sorted(cmap,key=lambda c:c[1])
    starts=np.array([c[1] for c in cmap],dtype=np.int64); ends=np.array([c[2] for c in cmap],dtype=np.int64); nums=np.array([c[0] for c in cmap],dtype=np.int64)
    pos=m.ref_position.values.astype(np.int64); mq=m.mod_qual.values
    idx=np.searchsorted(ends,pos); ic=np.clip(idx,0,len(nums)-1)
    valid=(idx<len(nums))&(starts[ic]<=pos)&(pos<ends[ic])
    keep=valid&((mq<=MOD_LO)|(mq>=MOD_HI))
    if not keep.any(): return None,None
    cn=nums[ic][keep]; wpos=pos[keep]-starts[ic][keep]; mqk=mq[keep]
    inr=(wpos>=0)&(wpos<REPEAT_LEN); cn=cn[inr]; wpos=wpos[inr]; mqk=mqk[inr]
    if len(cn)==0: return None,None
    is_meth=(mqk>=MOD_HI).astype(np.int64)
    region=np.where((wpos>=GENE_S)&(wpos<=GENE_E),"gene",np.where(wpos<GENE_S,"nts_pre","nts_post"))
    is_alu=((wpos>=ALU_S)&(wpos<=ALU_E)).astype(np.int64)
    mm=pd.DataFrame({"cn":cn,"wpos":wpos,"is_meth":is_meth,"region":region,"is_alu":is_alu})
    ov=mm.groupby("cn")["is_meth"].agg(n_conf_calls="count",n_meth="sum")
    reg=mm.groupby(["cn","region"])["is_meth"].agg(n="count",s="sum").unstack("region")
    alu=mm[mm.is_alu==1].groupby("cn")["is_meth"].agg(alu_n="count",alu_meth="sum")
    cc=ov.copy()
    for rn in ("nts_pre","gene","nts_post"):
        cc[f"{rn}_n"]=reg[("n",rn)] if ("n",rn) in reg.columns else 0
        cc[f"{rn}_meth"]=reg[("s",rn)] if ("s",rn) in reg.columns else 0
    cc=cc.join(alu).fillna(0); cc["mean_meth"]=cc.n_meth/cc.n_conf_calls
    cc=cc.reset_index().rename(columns={"cn":"copy_number"}); cc.insert(0,"hap",h); cc.insert(0,"sample",sample)
    pp=mm.groupby(["cn","wpos"])["is_meth"].agg(n_conf="count",n_meth="sum").reset_index()
    pp=pp.rename(columns={"cn":"copy_number","wpos":"wpos_bin"}); pp.insert(0,"hap",h); pp.insert(0,"sample",sample)
    print(f"  {sample} {h}: {len(mm)} CpG calls, {int(cc.copy_number.nunique())} copies",flush=True)
    return cc,pp

def process_sample(sample, haps_all, haps_out, hifi, twostep=False):
    md=METH/sample; tag="hifi" if hifi else "ont"
    reads=md/(f"{sample}_hifi_meth_reads.fq.gz" if hifi else f"{sample}_meth_reads.fq.gz")
    if not reads.exists(): return [],[]
    preset="map-hifi" if hifi else "map-ont"; sfx="2s" if twostep else "mc"
    refs={h:build_ref(sample,h) for h in haps_all}; refs={h:r for h,r in refs.items() if r}
    if twostep:
        comprefs={h:build_compref(sample,h) for h in haps_all}; comprefs={h:p for h,p in comprefs.items() if p}
    else:
        comprefs={h:refs[h][0] for h in refs}
    if not any(h in comprefs for h in haps_out): return [],[]
    comp=md/f"{sample}_{tag}_{sfx}_comp.fa"
    with open(comp,"w") as o:
        for p in comprefs.values(): o.write(Path(p).read_text())
    run(f"'{SAMTOOLS}' faidx '{comp}'")
    bam=md/f"{sample}_{tag}_{sfx}_aln.bam"
    with tempfile.NamedTemporaryFile("w",suffix=".py",delete=False) as t: t.write(TAGGER); tg=t.name
    try:
        run(f"'{MINIMAP2}' -ax {preset} -y -N1 -t 8 '{comp}' '{reads}' 2>/dev/null | '{ENV}/python3' {tg} | '{SAMTOOLS}' view -b -F 0x4 | '{SAMTOOLS}' sort -o '{bam}'")
        run(f"'{SAMTOOLS}' index '{bam}'")
    finally: os.unlink(tg)
    copy_rows=[]; pos_rows=[]
    for h in haps_out:
        if h not in refs: continue
        refp,htag,cmap=refs[h]; mod=md/f"{sample}_{h}_{tag}_{sfx}_modkit.tsv"
        if mod.exists(): mod.unlink()
        if twostep:
            # extract hap-assigned reads (with MM/ML), re-align to mc_gene concat-units
            fq=md/f"{sample}_{h}_{tag}_2s.fq"; abam=md/f"{sample}_{h}_{tag}_2s_assign.bam"
            run(f"'{SAMTOOLS}' view -bh -d 'YH:{h}' '{bam}' | '{SAMTOOLS}' fastq -T MM,ML - > '{fq}' 2>/dev/null")
            if os.path.getsize(fq)<10: continue
            run(f"'{MINIMAP2}' -ax {preset} -y -t 4 '{refp}' '{fq}' 2>/dev/null | '{SAMTOOLS}' view -b -F 0x904 | '{SAMTOOLS}' sort -o '{abam}'"); run(f"'{SAMTOOLS}' index '{abam}'")
            run(f"'{MODKIT}' extract full --cpg --mapped-only --reference '{refp}' -t 4 '{abam}' '{mod}' 2>/dev/null")
        else:
            fbam=md/f"{sample}_{h}_{tag}_mc_only.bam"
            run(f"'{SAMTOOLS}' view -bh -d 'YH:{h}' '{bam}' | '{SAMTOOLS}' sort -o '{fbam}'"); run(f"'{SAMTOOLS}' index '{fbam}'")
            run(f"'{MODKIT}' extract full --cpg --mapped-only --reference '{refp}' -t 4 '{fbam}' '{mod}' 2>/dev/null")
        m=pd.read_csv(mod,sep="\t"); m=m[(m.mod_code=="m")&(m.chrom==htag)]
        if m.empty: continue
        cc,pp=_assign(m,cmap,sample,h)
        if cc is not None: copy_rows.append(cc); pos_rows.append(pp)
    return copy_rows,pos_rows

def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--samples",required=True)
    ap.add_argument("--hifi",action="store_true"); ap.add_argument("--out",default="mc")
    ap.add_argument("--outhaps",default="",help="validation: force these hap labels as output (comma)")
    ap.add_argument("--twostep",action="store_true",help="fair competition (continuous refs) then re-align to unit refs")
    a=ap.parse_args()
    import sqlite3
    # determine, per sample, all haps + which are multi-contig (have mc_gene.fa)
    samples=a.samples.split(",")
    allc=[]; allp=[]
    for s in samples:
        haps=[p.name.split("_")[-2]+"" for p in []]  # placeholder
        # haps from databases_mc (multi) + databases (all)
        multi=[f.name[len(s)+1:-4] for f in DBMC.glob(f"{s}_*.tsv")]
        allh=sorted({f.name[len(s)+1:-4] for f in DB.glob(f"{s}_*.tsv")} | set(multi))
        out_haps=a.outhaps.split(",") if a.outhaps else multi  # validation override / multi-contig
        cr,pr=process_sample(s,allh,out_haps,a.hifi,a.twostep)
        allc+=cr; allp+=pr
    suf="_hifi" if a.hifi else ""
    src="HPRC_HiFi" if a.hifi else "HPRC_Year1_ONT"
    cdf=pd.concat(allc,ignore_index=True) if allc else pd.DataFrame()
    if len(cdf): cdf["source"]=src
    pdf=pd.concat(allp,ignore_index=True) if allp else pd.DataFrame()
    cdf.to_csv(METH/f"copy_meth_export_{a.out}{suf}.tsv",sep="\t",index=False)
    pdf.to_csv(METH/f"copy_meth_pos_export_{a.out}{suf}.tsv",sep="\t",index=False)
    print(f"wrote {len(cdf)} copy rows, {len(pdf)} pos rows")

if __name__=="__main__": main()
