#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# 59_migrate_multicontig.py — replace dominant-contig copies with the complete multi-contig array in the database.
# Author: Jonathan Bohlen
# (c) 2026 Jonathan Bohlen. Code accompanying Sengl et al. (2026),
# "Surveillance and selection of 5S ribosomal RNA genes in the human genome."
# Released under the MIT License; see LICENSE at the repository root.
# -----------------------------------------------------------------------------
"""
59_migrate_multicontig.py — DB-wide multi-contig array recovery.

For multi-contig haplotypes, replace the dominant-contig-only copies with the
complete flank-ordered array (from databases_mc/{sample}_{hap}.tsv). Delete-reinsert
is safe because the only copy_id-keyed data on these haplotypes is methylation
(regenerated separately) and variants (re-derived separately).

Adds schema:
  copy:      source_contig TEXT, contig_rank INTEGER, array_member INTEGER DEFAULT 1
  haplotype: n_array_contigs INT, array_fragmented INT, array_order_resolved TEXT,
             n_copies_dominant_legacy INT, contigs_fasta TEXT, array_contigs TEXT
Single-contig haps get backfilled defaults (n_array_contigs=1, order='single', …).

Orphan/dispersed 5S (array_member=0) are stored but excluded from array analyses.

Paths are read from environment variables:
  FIVES_DB    path to 5S_rDNA.db
  FIVES_DATA  input derived-data directory (databases_mc/)

Usage: 59_migrate_multicontig.py [--dry-run]
"""
import os, re, sqlite3, sys
from pathlib import Path
import pandas as pd

T2T  = Path(os.environ.get("FIVES_DATA", "data"))
HPRC = T2T
DB   = Path(os.environ.get("FIVES_DB", "5S_rDNA.db"))
MCDIR = HPRC / "databases_mc"
DRY = "--dry-run" in sys.argv
NTS_POST_START_1 = 749

def parse_variants_str(s):
    if not s or str(s).strip().lower() in ("none",".","","nan"): return []
    out=[]
    for tok in str(s).split(";"):
        m=re.match(r"\s*(\d+):([A-Za-z]+)>([A-Za-z]+)",tok.strip())
        if m: out.append((int(m.group(1)),m.group(2),m.group(3)))
    return out
def region_of(p): return "nts_pre" if p<630 else "gene" if p<=748 else "nts_post"
def si(v):
    try: f=float(v); return int(f) if pd.notna(f) else None
    except: return None
def sf(v):
    try: f=float(v); return f if pd.notna(f) else None
    except: return None
def ss(v):
    s=str(v).strip() if v is not None else ""
    return None if s.lower() in ("nan","none","",".") else s

def ensure_schema(con):
    cc={r[1] for r in con.execute("PRAGMA table_info(copy)")}
    for col,typ in [("source_contig","TEXT"),("contig_rank","INTEGER"),("array_member","INTEGER DEFAULT 1")]:
        if col not in cc: con.execute(f"ALTER TABLE copy ADD COLUMN {col} {typ}")
    hc={r[1] for r in con.execute("PRAGMA table_info(haplotype)")}
    for col,typ in [("n_array_contigs","INTEGER"),("array_fragmented","INTEGER"),
                    ("array_order_resolved","TEXT"),("n_copies_dominant_legacy","INTEGER"),
                    ("contigs_fasta","TEXT"),("array_contigs","TEXT")]:
        if col not in hc: con.execute(f"ALTER TABLE haplotype ADD COLUMN {col} {typ}")

def main():
    con=sqlite3.connect(DB); con.execute("PRAGMA foreign_keys=OFF")
    ensure_schema(con)
    tsvs=sorted(MCDIR.glob("*.tsv"))
    print(f"{len(tsvs)} multi-contig TSVs; DRY-RUN={DRY}")
    # hap lookup
    haps={}
    for sid,hl,hid in con.execute("""SELECT a.sample_id,h.hap_label,h.haplotype_id
                                     FROM haplotype h JOIN assembly a USING(assembly_id)"""):
        haps[(sid,hl)]=hid
    tot_before=tot_after=tot_disp=0
    for f in tsvs:
        df=pd.read_csv(f,sep="\t")
        sid=str(df.sample_id.iloc[0]); hl=str(df.haplotype.iloc[0])
        hid=haps.get((sid,hl))
        if hid is None: print(f"  [skip] {sid} {hl} not in DB"); continue
        before=con.execute("SELECT COUNT(*) FROM copy WHERE haplotype_id=?",(hid,)).fetchone()[0]
        arr=df[df.array_member==1].copy(); disp=df[df.array_member==0].copy()
        order_res=ss(df.order_resolved.iloc[0])
        arr_contigs=sorted(set(arr.source_contig.dropna()))       # array-fragment contigs only
        all_contigs=sorted(set(df.source_contig.dropna()))        # all 5S-bearing (in saved FASTA)
        n_arr_ctg=len(arr_contigs); contigs=arr_contigs
        contig_fa=f"sequences/{sid}/{sid}_{hl}_5S_contigs.fa.gz"
        tot_before+=before; tot_after+=len(arr); tot_disp+=len(disp)
        if DRY:
            print(f"  {sid} {hl}: {before} -> {len(arr)} array (+{len(disp)} disp), "
                  f"{len(contigs)} contigs [{order_res}]")
            continue
        # delete methylation + variants + copies for this hap
        cids=[r[0] for r in con.execute("SELECT copy_id FROM copy WHERE haplotype_id=?",(hid,))]
        if cids:
            ph=",".join("?"*len(cids))
            for tbl in ("copy_methylation","copy_methylation_hifi","copy_meth_pos","copy_meth_pos_hifi","variant"):
                con.execute(f"DELETE FROM {tbl} WHERE copy_id IN ({ph})",cids)
            con.execute(f"DELETE FROM copy WHERE copy_id IN ({ph})",cids)
        # insert array copies + variants
        def insert_copy(r,member,cat,copy_num=None):
            cn = copy_num if copy_num is not None else si(r.get("copy_id"))
            cid=con.execute("""INSERT INTO copy
              (haplotype_id,copy_number,unit_start_local,unit_end_local,unit_length_bp,
               spacing_to_next_bp,gene_lo_local,gene_hi_local,gene_pct_identity,
               gene_mismatches,gene_gaps,n_snv_gene,n_snv_5s_gene,n_snv_nts_pre,
               n_snv_nts_post,category,border_note,source_contig,contig_rank,array_member)
              VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)""",
              (hid,cn,si(r.get("unit_start_local")),si(r.get("unit_end_local")),
               si(r.get("unit_length_bp")),sf(r.get("spacing_to_next_bp")),
               si(r.get("gene_lo_local")),si(r.get("gene_hi_local")),sf(r.get("gene_pct_identity")),
               si(r.get("gene_mismatches")),si(r.get("gene_gaps")),si(r.get("n_snv_gene")),
               si(r.get("n_snv_5s_gene")),si(r.get("n_snv_nts_pre")),si(r.get("n_snv_nts_post")),
               ss(r.get("category2")) or cat,ss(r.get("border_note")),
               ss(r.get("source_contig")),si(r.get("contig_rank")),member)).lastrowid
            return cid
        for _,r in arr.iterrows():
            cid=insert_copy(r,1,"array")
            vr=[]
            for p0,ref,alt in parse_variants_str(r.get("gene_variants","")):
                vr.append((cid,"gene_unit",p0+1,ref,alt,region_of(p0+1)))
            for p0,ref,alt in parse_variants_str(r.get("nts_pre_variants","")):
                vr.append((cid,"nts_pre_aln",p0+1,ref,alt,"nts_pre"))
            for p0,ref,alt in parse_variants_str(r.get("nts_post_variants","")):
                vr.append((cid,"nts_post_aln",NTS_POST_START_1+p0,ref,alt,"nts_post"))
            if vr: con.executemany("""INSERT INTO variant
                (copy_id,alignment_source,consensus_pos,ref,alt,region) VALUES (?,?,?,?,?,?)""",vr)
        for i,(_,r) in enumerate(disp.iterrows(),1):
            insert_copy(r,0,"dispersed_5S",copy_num=len(arr)+i)
        # update haplotype metadata
        con.execute("""UPDATE haplotype SET n_copies=?, n_array_contigs=?, array_fragmented=?,
            array_order_resolved=?, n_copies_dominant_legacy=?, contigs_fasta=?, array_contigs=?,
            array_chrom=? WHERE haplotype_id=?""",
            (len(arr),n_arr_ctg,1 if n_arr_ctg>=2 else 0,order_res,before,contig_fa,",".join(all_contigs),
             ss(arr.sort_values("copy_id").source_contig.iloc[0]) if len(arr) else None,hid))
        print(f"  {sid} {hl}: {before} -> {len(arr)} array on {n_arr_ctg} contig(s) (+{len(disp)} disp) [{order_res}]")
    # backfill single-contig haps (unchanged ones)
    if not DRY:
        con.execute("""UPDATE haplotype SET n_array_contigs=1, array_fragmented=0,
            array_order_resolved='single', n_copies_dominant_legacy=n_copies
            WHERE n_array_contigs IS NULL""")
        con.execute("""UPDATE copy SET array_member=1 WHERE array_member IS NULL""")
        # source_contig=array_chrom, contig_rank=1 for single-contig copies missing it
        con.execute("""UPDATE copy SET source_contig=(SELECT array_chrom FROM haplotype h
            WHERE h.haplotype_id=copy.haplotype_id), contig_rank=1
            WHERE source_contig IS NULL""")
        con.commit()
    print(f"\nTOTAL: before={tot_before} -> array={tot_after} (+{tot_after-tot_before}); dispersed stored={tot_disp}")
    con.close()

if __name__=="__main__":
    main()
