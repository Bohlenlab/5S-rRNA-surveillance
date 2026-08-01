#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# 59c_rederive_gene_unit_t2t.py — additively derive T2T-polarized gene_unit_t2t per-copy variants.
# Author: Jonathan Bohlen
# (c) 2026 Jonathan Bohlen. Code accompanying Sengl et al. (2026),
# "Surveillance and selection of 5S ribosomal RNA genes in the human genome."
# Released under the MIT License; see LICENSE at the repository root.
# -----------------------------------------------------------------------------
"""
59c — Additively derive gene_unit_t2t (T2T-polarized per-copy variants) for every
haplotype that has 'gene_unit' but lacks 'gene_unit_t2t'. Each copy's gene-region
sequence is aligned to the T2T consensus (MAFFT --add, orientation auto-detected)
and its deviations are recorded. Inserts ADDITIVELY (preserves existing
consensus_t2t / gene_unit / nts variants).

Paths are read from environment variables:
  FIVES_DB    path to 5S_rDNA.db
  FIVES_DATA  input derived-data directory (sequences/)
  FIVES_REFS  reference directory (consensus FASTA)

Run: python3 59c_rederive_gene_unit_t2t.py
"""
import os, re, subprocess, tempfile, sqlite3
from pathlib import Path
T2T=Path(os.environ.get("FIVES_DATA", "data"))
DB=Path(os.environ.get("FIVES_DB", "5S_rDNA.db")); SEQ=T2T/"sequences"; T2T_FA=Path(os.environ.get("FIVES_REFS", "refs"))/"5S_t2t_consensus.fa"
MAFFT="mafft"
COMP=str.maketrans("ACGTacgt","TGCAtgca"); ORIENT_SNP_THRESHOLD=100

def parse_fasta(path):
    seqs,h,buf={},None,[]
    for line in open(path):
        line=line.rstrip()
        if line.startswith(">"):
            if h: seqs[h]="".join(buf)
            h=line[1:].split()[0]; buf=[]
        elif h: buf.append(line.upper())
    if h: seqs[h]="".join(buf)
    return seqs
def read_fasta_string(text):
    seqs,h,buf={},None,[]
    for line in text.splitlines():
        line=line.rstrip()
        if line.startswith(">"):
            if h: seqs[h]="".join(buf)
            h=line[1:].split()[0]; buf=[]
        elif h: buf.append(line.upper())
    if h: seqs[h]="".join(buf)
    return seqs
def _mafft_add(aln_path,ref_seq,ref_id):
    with tempfile.NamedTemporaryFile(mode="w",suffix=".fa",delete=False) as tmp:
        tmp.write(f">{ref_id}\n{ref_seq}\n"); rf=tmp.name
    try:
        r=subprocess.run([MAFFT,"--add",rf,"--quiet",str(aln_path)],capture_output=True,text=True,check=True)
    finally: os.unlink(rf)
    return read_fasta_string(r.stdout)
def _snp_count(a,b): return sum(1 for x,y in zip(a,b) if x!="-" and y!="-" and x!=y)
def _region(p): return "nts_pre" if p<630 else "gene" if p<=748 else "nts_post"
def derive_t2t_variants(aln,t2t_seq):
    L=len(t2t_seq); aug=_mafft_add(aln,t2t_seq,"T2T_ref"); trow=None; cs={}
    for rid,seq in aug.items():
        if rid.startswith("T2T"): trow=seq
        else:
            m=re.match(r"copy(\d+)",rid)
            if m: cs[int(m.group(1))]=seq
    if trow is None or not cs: return {},"fwd"
    samp=list(cs.values())[:10]; avg=sum(_snp_count(trow,s) for s in samp)/len(samp)
    if avg>ORIENT_SNP_THRESHOLD:
        rc=t2t_seq.translate(COMP)[::-1]; aug=_mafft_add(aln,rc,"T2T_rc"); trow=None; cs={}
        for rid,seq in aug.items():
            if rid.startswith("T2T"): trow=seq
            else:
                m=re.match(r"copy(\d+)",rid)
                if m: cs[int(m.group(1))]=seq
        if trow is None: return {},"rc"
        c2r={}; cc=0
        for col,b in enumerate(trow):
            if b!="-": cc+=1; c2r[col]=cc
        res={}
        for cn,seq in cs.items():
            v=[]
            for col,prc in c2r.items():
                if seq[col]=="-": continue
                if seq[col]!=trow[col]:
                    pos=L-prc+1; v.append((pos,t2t_seq[pos-1],seq[col].translate(COMP),_region(pos)))
            res[cn]=v
        return res,"rc"
    c2t={}; cc=0
    for col,b in enumerate(trow):
        if b!="-": cc+=1; c2t[col]=cc
    res={}
    for cn,seq in cs.items():
        v=[]
        for col,tp in c2t.items():
            if seq[col]!="-" and seq[col]!=trow[col]: v.append((tp,trow[col],seq[col],_region(tp)))
        res[cn]=v
    return res,"fwd"

def main():
    t2t_seq=list(parse_fasta(T2T_FA).values())[0]
    con=sqlite3.connect(DB)
    need=con.execute("""
      SELECT a.sample_id,h.hap_label,h.haplotype_id FROM haplotype h JOIN assembly a USING(assembly_id)
      WHERE (SELECT COUNT(*) FROM variant v JOIN copy c2 USING(copy_id)
             WHERE c2.haplotype_id=h.haplotype_id AND v.alignment_source='gene_unit_t2t')=0
        AND (SELECT COUNT(*) FROM variant v JOIN copy c2 USING(copy_id)
             WHERE c2.haplotype_id=h.haplotype_id AND v.alignment_source='gene_unit')>0
      ORDER BY a.sample_id,h.hap_label""").fetchall()
    print(f"{len(need)} haplotypes need gene_unit_t2t",flush=True)
    ok=fail=tot=0
    for sid,hl,hid in need:
        aln=SEQ/sid/f"{sid}_{hl}_gene_aln.fa"
        if not aln.exists(): print(f"  SKIP {sid} {hl}: no gene_aln"); fail+=1; continue
        cidmap={cn:cid for cn,cid in con.execute(
            "SELECT copy_number,copy_id FROM copy WHERE haplotype_id=? AND array_member=1",(hid,))}
        try: cv,orient=derive_t2t_variants(aln,t2t_seq)
        except Exception as e: print(f"  FAIL {sid} {hl}: {e}"); fail+=1; continue
        cids=tuple(cidmap.values())
        if cids: con.execute("DELETE FROM variant WHERE alignment_source='gene_unit_t2t' AND copy_id IN (%s)"%",".join("?"*len(cids)),cids)
        n=0
        for cn,vl in cv.items():
            cid=cidmap.get(cn)
            if cid is None: continue
            for pos,ref,alt,region in vl:
                con.execute("INSERT INTO variant (copy_id,alignment_source,consensus_pos,ref,alt,region) VALUES (?,'gene_unit_t2t',?,?,?,?)",(cid,int(pos),ref,alt,region)); n+=1
        con.commit(); ok+=1; tot+=n
        if ok%25==0: print(f"  ... {ok}/{len(need)}",flush=True)
    print(f"done: {ok} haps, {tot} gene_unit_t2t variants ({fail} failed)")
    n0=con.execute("""SELECT COUNT(*) FROM haplotype h WHERE
      (SELECT COUNT(*) FROM variant v JOIN copy c2 USING(copy_id) WHERE c2.haplotype_id=h.haplotype_id AND v.alignment_source='gene_unit_t2t')=0
      AND (SELECT COUNT(*) FROM copy c WHERE c.haplotype_id=h.haplotype_id AND c.array_member=1)>0""").fetchone()[0]
    print(f"array haplotypes still lacking gene_unit_t2t: {n0}")
    con.close()

if __name__=="__main__": main()
