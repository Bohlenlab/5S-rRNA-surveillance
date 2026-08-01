#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# 59b_rederive_consensus_migrated.py — re-derive consensus_t2t variants for the multi-contig-migrated haplotypes.
# Author: Jonathan Bohlen
# (c) 2026 Jonathan Bohlen. Code accompanying Sengl et al. (2026),
# "Surveillance and selection of 5S ribosomal RNA genes in the human genome."
# Released under the MIT License; see LICENSE at the repository root.
# -----------------------------------------------------------------------------
"""
59b — Re-derive consensus_t2t variants for the multi-contig-migrated haplotypes.

After migration, the migrated haplotypes' consensus_t2t variants were deleted with
their old copies. This reuses the pass2_derive routine from the consensus-rederivation
module (12_consensus_rederive.py) against the SAVED population consensus, for every
haplotype with a multi-contig array (contigs_fasta set). Array copies only.

Usage: 59b_rederive_consensus_migrated.py
"""
import importlib.util, json, sqlite3
from pathlib import Path
SCRIPTS = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("cr", SCRIPTS / "12_consensus_rederive.py")
cr = importlib.util.module_from_spec(spec); spec.loader.exec_module(cr)

def main():
    j = json.load(open(cr.CONS_OUT))
    consensus = list(j["consensus"]); t2t_seq = j["chm13"]; gap_pos = set(j["consensus_gap_positions"])
    cache_dir = cr.HPRC / "cache"; cache_dir.mkdir(exist_ok=True)
    con = cr.db_connect(); cr.ensure_columns(con)
    haps = con.execute("""SELECT a.sample_id, h.hap_label, h.haplotype_id
                          FROM haplotype h JOIN assembly a USING(assembly_id)
                          WHERE h.contigs_fasta IS NOT NULL ORDER BY a.sample_id,h.hap_label""").fetchall()
    print(f"{len(haps)} migrated haplotypes to re-derive consensus_t2t")
    ok=snp=indel=fail=0
    for sid, hap, hid in haps:
        cmap = {cn: cid for cn, cid in con.execute(
            "SELECT copy_number, copy_id FROM copy WHERE haplotype_id=? AND array_member=1", (hid,))}
        task = dict(sample_id=sid, hap_label=hap, haplotype_id=hid, copy_id_map=cmap,
                    t2t_seq=t2t_seq, consensus=consensus, gap_positions=gap_pos, cache_dir=cache_dir)
        res = cr.pass2_derive(task)
        if res["status"] != "ok":
            print(f"  FAIL {sid} {hap}: {res.get('error')}"); fail+=1; continue
        cids = tuple(cmap.values())
        if cids:
            con.execute("DELETE FROM variant WHERE alignment_source='consensus_t2t' AND copy_id IN (%s)"
                        % ",".join("?"*len(cids)), cids)
        con.executemany("""INSERT INTO variant (copy_id,alignment_source,consensus_pos,ref,alt,region,
                           chm13_base,var_type,masked) VALUES (?,?,?,?,?,?,?,?,?)""", res["rows"])
        con.commit(); ok+=1; snp+=res["n_snp"]; indel+=res["n_indel"]
    print(f"done: {ok} ok, {fail} fail; {snp} SNP, {indel} indel")
    print("consensus_t2t by cohort:")
    for c,n in con.execute("""SELECT a.cohort,COUNT(*) FROM variant v JOIN copy c2 USING(copy_id)
        JOIN haplotype h USING(haplotype_id) JOIN assembly a USING(assembly_id)
        WHERE v.alignment_source='consensus_t2t' GROUP BY a.cohort ORDER BY a.cohort"""):
        print(f"  {c or 'NULL':16s}: {n:,}")
    con.close()

if __name__=="__main__":
    main()
