#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# 2E_fair_perdonor_vaf.py — per-diploid-individual UK Biobank vs assembly VAF comparison across assembly scopes.
# Author: Jonathan Bohlen
# (c) 2026 Jonathan Bohlen. Code accompanying Sengl et al. (2026),
# "Surveillance and selection of 5S ribosomal RNA genes in the human genome."
# Released under the MIT License; see LICENSE at the repository root.
# -----------------------------------------------------------------------------
"""
2E_fair_perdonor_vaf.py
UK Biobank vs assembly per-variant VAF comparison with BOTH sides on a per-diploid-individual
denominator (assembly VAF pooled across both haplotypes per donor). Run across assembly scopes
(EUR-matched Release2, all Release2, all HPRC).

Filters: gene_unit_t2t SNPs, array_member=1, positions 467-967, confident carriers VAF>=0.3%.
UKBB per-individual VAF is the mean within-person VAF (confident carriers) from the input table.

Outputs (<FIVES_OUT>/Figure2/):
  2E_UKBB_vs_assembly_VAF_perdonor_<scope>.tsv   per-variant table per scope
  2E_UKBB_vs_assembly_VAF_perdonor.stats.txt     scope comparison
"""
import os
import sqlite3
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, pearsonr

DB   = Path(os.environ.get("FIVES_DB", "5S_rDNA.db"))
FIG2 = Path(os.environ.get("FIVES_OUT", "output")) / "Figure2"
VAF_CONF = 0.003; LO, HI = 467, 967

SCOPES = {
    "EUR_Release2": "a.cohort='HPRC_Release2' AND a.superpopulation='EUR'",
    "ALL_Release2": "a.cohort='HPRC_Release2'",
    "ALL_HPRC":     "a.cohort IN ('HPRC_Year1','HPRC_Release2')",
}

uk = pd.read_csv(FIG2 / "2E_UKBB_vs_assembly_VAF.tsv", sep="\t")[["pos","alt","ukbb_conf_mean_vaf_pct"]]

def compute(where):
    con = sqlite3.connect(DB)
    tot_donor = dict(con.execute(f"""
        SELECT a.assembly_id, COUNT(*) FROM copy c
        JOIN haplotype h ON h.haplotype_id=c.haplotype_id JOIN assembly a ON a.assembly_id=h.assembly_id
        WHERE {where} AND c.array_member=1 GROUP BY a.assembly_id""").fetchall())
    tot_hap = dict(con.execute(f"""
        SELECT h.haplotype_id, COUNT(*) FROM copy c
        JOIN haplotype h ON h.haplotype_id=c.haplotype_id JOIN assembly a ON a.assembly_id=h.assembly_id
        WHERE {where} AND c.array_member=1 GROUP BY h.haplotype_id""").fetchall())
    rows = con.execute(f"""
        SELECT v.consensus_pos, v.alt, a.assembly_id, h.haplotype_id, COUNT(DISTINCT v.copy_id)
        FROM variant v JOIN copy c ON c.copy_id=v.copy_id
        JOIN haplotype h ON h.haplotype_id=c.haplotype_id JOIN assembly a ON a.assembly_id=h.assembly_id
        WHERE v.alignment_source='gene_unit_t2t' AND v.var_type='snp' AND c.array_member=1
          AND {where} AND v.consensus_pos BETWEEN {LO} AND {HI}
        GROUP BY v.consensus_pos, v.alt, a.assembly_id, h.haplotype_id""").fetchall()
    con.close()
    df = pd.DataFrame(rows, columns=["pos","alt","donor","hap","carr"])
    df["hap_vaf"] = df.carr / df.hap.map(tot_hap)
    don = df.groupby(["pos","alt","donor"])["carr"].sum().reset_index()
    don["donor_vaf"] = don.carr / don.donor.map(tot_donor)
    def pv(d, col, name):
        d = d[d[col] >= VAF_CONF]
        return (d.groupby(["pos","alt"])[col].mean() * 100).rename(name)
    m = (pd.concat([pv(df,"hap_vaf","asm_perhap_vaf_pct"),
                    pv(don,"donor_vaf","asm_perdonor_vaf_pct")], axis=1).reset_index()
            .merge(uk, on=["pos","alt"], how="inner").dropna())
    return len(tot_donor), m

def st(a, u):
    rs,_ = spearmanr(a, u); rp,_ = pearsonr(a, u)
    return rs, rp, np.median(a), np.median(u), np.median(u)/np.median(a)

out = [f"Per-donor UKBB-vs-assembly VAF: EUR vs ALL HPRC",
       f"gene_unit_t2t SNPs, array_member=1, pos {LO}-{HI}, confident carriers VAF>={VAF_CONF*100:.1f}%",
       f"UKBB = per-diploid within-person VAF (fixed); assembly = per-donor diploid VAF", ""]
for name, where in SCOPES.items():
    nd, m = compute(where)
    m.to_csv(FIG2 / f"2E_UKBB_vs_assembly_VAF_perdonor_{name}.tsv", sep="\t", index=False)
    rs_h, rp_h, ma_h, mu, ratio_h = st(m.asm_perhap_vaf_pct.values, m.ukbb_conf_mean_vaf_pct.values)
    rs_d, rp_d, ma_d, _,  ratio_d = st(m.asm_perdonor_vaf_pct.values, m.ukbb_conf_mean_vaf_pct.values)
    out += [f"[{name}]  donors={nd}  n_variants_matched={len(m)}  median UKBB VAF={mu:.3f}%",
            f"   per-haplotype: Spearman ρ={rs_h:.3f}  Pearson r={rp_h:.3f}  "
            f"med asm={ma_h:.3f}%  UKBB/asm ratio={ratio_h:.2f}",
            f"   per-donor    : Spearman ρ={rs_d:.3f}  Pearson r={rp_d:.3f}  "
            f"med asm={ma_d:.3f}%  UKBB/asm ratio={ratio_d:.2f}", ""]
txt = "\n".join(out)
(FIG2 / "2E_UKBB_vs_assembly_VAF_perdonor.stats.txt").write_text(txt)
print(txt)
