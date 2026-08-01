#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# 11_load_gtex_wgs_to_db.py — load GTEx WGS per-donor variant/QC TSVs into the
# 5S database as cohort GTEx_v9_WGS (idempotent, DB backed up first).
# Author: Jonathan Bohlen
# (c) 2026 Jonathan Bohlen. Code accompanying Sengl et al. (2026),
# "Surveillance and selection of 5S ribosomal RNA genes in the human genome."
# Released under the MIT License; see LICENSE at the repository root.
# -----------------------------------------------------------------------------
"""Idempotently load GTEx WGS per-donor variant/QC TSVs into 5S_rDNA.db as cohort GTEx_v9_WGS.
Backs up the DB first. assembly upsert (ON CONFLICT, preserves assembly_id); read_variant +
wgs_cohort_qc are delete-then-insert scoped to this cohort."""
import os, sys, sqlite3, csv, glob, shutil, datetime, argparse
from pathlib import Path
SUBJ=str(Path(os.environ.get("FIVES_DATA","data"))/"metadata"/"GTEx_Analysis_2025-08-22_v11_Annotations_SubjectPhenotypesDS.txt")
SAMP=str(Path(os.environ.get("FIVES_DATA","data"))/"metadata"/"GTEx_Analysis_2025-08-22_v11_Annotations_SampleAttributesDS.txt")
COHORT="GTEx_v9_WGS"; MOD="illumina"; SOURCE="GTEx_AnVIL_phs000424_v9WGS"
RACE={'1':('Asian','EAS'),'2':('Black','AFR'),'3':('White','EUR'),'4':('AmerIndian','AMR')}

ap=argparse.ArgumentParser()
ap.add_argument("--db", default=os.environ.get("FIVES_DB","5S_rDNA.db"))
ap.add_argument("--cohort-dir", default=str(Path(os.environ.get("FIVES_DATA","data"))/"results"/"wgs"/"cohort"))
ap.add_argument("--no-backup", action="store_true")
a=ap.parse_args()

def midage(s):
    try: lo,hi=s.split('-'); return (int(lo)+int(hi))//2
    except: return None
def fnum(x):
    try: return float(x)
    except: return None

subj={}
for r in csv.DictReader(open(SUBJ),delimiter='\t'):
    subj[r['SUBJID']]={'sex':{'1':'M','2':'F'}.get(r.get('SEX'),None),'age':midage(r.get('AGE','')),'race':r.get('RACE','')}
rna=set()
for r in csv.DictReader(open(SAMP),delimiter='\t'):
    if r['SMAFRZE']=='RNASEQ': rna.add('-'.join(r['SAMPID'].split('-')[:2]))

if not a.no_backup and os.path.exists(a.db):
    bk=a.db+'.bak_'+datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
    shutil.copy2(a.db,bk); print("backup:",bk)

donors=[]
for q in sorted(glob.glob(a.cohort_dir+"/*.qc.tsv")):
    rows=list(csv.DictReader(open(q),delimiter='\t'))
    if rows and rows[0].get('status')=='ok': donors.append(rows[0])
print(f"{len(donors)} ok donors to load from {a.cohort_dir}")
if not donors: sys.exit("nothing to load")

con=sqlite3.connect(a.db); con.execute("PRAGMA foreign_keys=ON"); cur=con.cursor()
cur.execute("""CREATE TABLE IF NOT EXISTS wgs_cohort_qc(
  assembly_id INTEGER PRIMARY KEY REFERENCES assembly(assembly_id),
  median_depth REAL, control_depth REAL, est_copies REAL,
  n_variants INTEGER, callable_fraction REAL, slice_reads INTEGER)""")
con.execute("BEGIN")
for r in donors:
    d=r['donor_id']; sp=subj.get(d,{}); pop,supop=RACE.get(sp.get('race'),(None,None))
    cur.execute("""INSERT INTO assembly(sample_id,cohort,population,superpopulation,sex,age,
        has_hifi,has_illumina,has_methylation,has_rnaseq,source) VALUES(?,?,?,?,?,?,0,1,0,?,?)
      ON CONFLICT(sample_id) DO UPDATE SET cohort=excluded.cohort,population=excluded.population,
        superpopulation=excluded.superpopulation,sex=excluded.sex,age=excluded.age,
        has_illumina=1,has_rnaseq=excluded.has_rnaseq,source=excluded.source""",
      (d,COHORT,pop,supop,sp.get('sex'),sp.get('age'),1 if d in rna else 0,SOURCE))
cur.execute("SELECT sample_id,assembly_id FROM assembly WHERE cohort=?", (COHORT,))
aid=dict(cur.fetchall()); ids=tuple(aid.values()); qm=",".join("?"*len(ids))
cur.execute(f"DELETE FROM read_variant WHERE modality=? AND assembly_id IN ({qm})",(MOD,*ids))
cur.execute(f"DELETE FROM wgs_cohort_qc WHERE assembly_id IN ({qm})",ids)
nv=0
for r in donors:
    d=r['donor_id']; i=aid[d]; vf=a.cohort_dir+f"/{d}.variants.tsv"
    if os.path.exists(vf):
        for v in csv.DictReader(open(vf),delimiter='\t'):
            cur.execute("INSERT INTO read_variant(assembly_id,modality,consensus_pos,ref,alt,depth,alt_depth,vaf,region) VALUES(?,?,?,?,?,?,?,?,?)",
              (i,MOD,int(v['consensus_pos']),v['ref'],v['alt'],int(v['depth']),int(v['alt_depth']),float(v['vaf']),v['region'])); nv+=1
    cur.execute("INSERT INTO wgs_cohort_qc VALUES(?,?,?,?,?,?,?)",
      (i,fnum(r['median_depth']),fnum(r['control_depth']),fnum(r['est_copies']),int(r['n_variants']),fnum(r['callable_fraction']),int(r['slice_reads'])))
con.commit()
print(f"loaded {len(donors)} donors, {nv} read_variant rows")
for label,q in [("assembly cohort",f"SELECT count(*) FROM assembly WHERE cohort='{COHORT}'"),
  ("read_variant",f"SELECT count(*) FROM read_variant WHERE modality='{MOD}' AND assembly_id IN (SELECT assembly_id FROM assembly WHERE cohort='{COHORT}')"),
  ("wgs_cohort_qc","SELECT count(*) FROM wgs_cohort_qc")]:
    print(f"  {label}: {cur.execute(q).fetchone()[0]}")
con.close()
