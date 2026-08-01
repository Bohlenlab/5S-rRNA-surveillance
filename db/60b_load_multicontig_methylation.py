#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# 60b_load_multicontig_methylation.py — load multi-contig methylation exports into the database.
# Author: Jonathan Bohlen
# (c) 2026 Jonathan Bohlen. Code accompanying Sengl et al. (2026),
# "Surveillance and selection of 5S ribosomal RNA genes in the human genome."
# Released under the MIT License; see LICENSE at the repository root.
# -----------------------------------------------------------------------------
"""
60b — Load multi-contig methylation exports into SQLite.

Maps (sample, hap_label, copy_number) -> copy_id (array_member=1 copies) and
inserts into copy_methylation / copy_methylation_hifi (per-copy) and
copy_meth_pos / copy_meth_pos_hifi (per-position). Idempotent: deletes any
existing rows for the loaded copy_ids first (the migrated haps' methylation was
already deleted, so normally a no-op).

Paths are read from environment variables:
  FIVES_DB    path to 5S_rDNA.db
  FIVES_DATA  input derived-data directory (methylation/)

Usage: 60b_load_multicontig_methylation.py            # ONT
       60b_load_multicontig_methylation.py --hifi     # HiFi
"""
import argparse, os, sqlite3
from pathlib import Path
import pandas as pd

T2T=Path(os.environ.get("FIVES_DATA", "data"))
DB=Path(os.environ.get("FIVES_DB", "5S_rDNA.db")); MCDIR=T2T/"methylation"

def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--hifi",action="store_true"); a=ap.parse_args()
    suf="_hifi" if a.hifi else ""
    ctab="copy_methylation_hifi" if a.hifi else "copy_methylation"
    ptab="copy_meth_pos_hifi" if a.hifi else "copy_meth_pos"
    cexp=MCDIR/f"copy_meth_export_mc{suf}.tsv"; pexp=MCDIR/f"copy_meth_pos_export_mc{suf}.tsv"
    cdf=pd.read_csv(cexp,sep="\t"); pdf=pd.read_csv(pexp,sep="\t")
    con=sqlite3.connect(DB)
    lut={(s,hl,cn):cid for s,hl,cn,cid in con.execute(
        """SELECT a.sample_id,h.hap_label,c.copy_number,c.copy_id FROM copy c
           JOIN haplotype h USING(haplotype_id) JOIN assembly a USING(assembly_id)
           WHERE c.array_member=1""").fetchall()}
    def cid(r): return lut.get((r["sample"],r["hap"],int(r["copy_number"])))
    cdf["copy_id"]=cdf.apply(cid,axis=1); pdf["copy_id"]=pdf.apply(cid,axis=1)
    miss=cdf.copy_id.isna().sum()
    print(f"{ctab}: {len(cdf)} export rows, {miss} unmapped (dropped)")
    cdf=cdf.dropna(subset=["copy_id"]); pdf=pdf.dropna(subset=["copy_id"])
    cdf["copy_id"]=cdf.copy_id.astype(int); pdf["copy_id"]=pdf.copy_id.astype(int)
    ids=tuple(set(cdf.copy_id))
    if ids:
        con.execute(f"DELETE FROM {ctab} WHERE copy_id IN (%s)"%",".join("?"*len(ids)),ids)
        con.execute(f"DELETE FROM {ptab} WHERE copy_id IN (%s)"%",".join("?"*len(ids)),ids)
    con.executemany(f"""INSERT INTO {ctab}
        (copy_id,source,n_conf_calls,n_meth,mean_meth,nts_pre_n,nts_pre_meth,gene_n,gene_meth,
         nts_post_n,nts_post_meth,alu_n,alu_meth) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)""",
        [(int(r.copy_id),r.source,int(r.n_conf_calls),int(r.n_meth),float(r.mean_meth),
          int(r.nts_pre_n),int(r.nts_pre_meth),int(r.gene_n),int(r.gene_meth),
          int(r.nts_post_n),int(r.nts_post_meth),int(r.alu_n),int(r.alu_meth)) for r in cdf.itertuples()])
    con.executemany(f"INSERT INTO {ptab} (copy_id,wpos_bin,n_conf,n_meth) VALUES (?,?,?,?)",
        [(int(r.copy_id),int(r.wpos_bin),int(r.n_conf),int(r.n_meth)) for r in pdf.itertuples()])
    con.commit()
    print(f"loaded {len(cdf)} {ctab} rows, {len(pdf)} {ptab} rows; copies={len(ids)}")
    print("now total:", con.execute(f"SELECT COUNT(*) FROM {ctab}").fetchone()[0], ctab)
    con.close()

if __name__=="__main__": main()
