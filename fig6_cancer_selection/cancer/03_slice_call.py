#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# 03_slice_call.py — call 5S variants from controlled TCGA/CPTAC BAMs via the GDC slicing API.
# Author: Jonathan Bohlen
# (c) 2026 Jonathan Bohlen. Code accompanying Sengl et al. (2026),
# "Surveillance and selection of 5S ribosomal RNA genes in the human genome."
# Released under the MIT License; see LICENSE at the repository root.
# -----------------------------------------------------------------------------
"""
03_slice_call.py — call 5S variants from controlled TCGA/CPTAC BAMs via the GDC slicing API.

Pulls the 5S locus region from controlled GDC BAMs via the GDC slicing API
(`/slicing/view/{file_id}?region=...`, header X-Auth-Token), which returns a mini BAM of the 5S
locus (+ unmapped). Downstream:
  slice -> collate -> fastq -> bwa mem -B6 -O8 -L5,5 -T30 to the T2T 5S consensus
        -> view -q30 -F0x904 -> sort -> bcftools mpileup -d200000 -q30 -Q30 -> per-pos AD/DP
Per-sample variant + QC + gz pileup. Resumable (skips finished samples). Calling/pooling is a
separate post step.

Usage:
  python3 scripts/03_slice_call.py --modality wgs --projects TCGA-OV TCGA-THCA [--sample-types "Primary Tumor"] [--max 20]
  python3 scripts/03_slice_call.py --modality rna --projects TCGA-OV --bed reference/rna5s_grch38_loci.bed --deep

Token path is read from GDC_TOKEN_FILE (default ~/.gdc/token). Requires bwa, samtools, bcftools, curl.
"""
import os, sys, subprocess, statistics, argparse, time, random
from pathlib import Path

HERE = Path(__file__).resolve().parent.parent
T2T  = HERE.parent
CONS = T2T / "reference" / "5S_t2t_consensus.fa"          # 2168 bp, gene 630-748
# alignment reference + coord offset (set in main): default = full consensus; RNA can use the
# gene+30 reference (reference/5S_gene_30flank.fa, consensus 600-778) with --ref-offset 599 so reported
# positions stay in consensus coordinates and the Alu/NTS reads are dropped at the bwa step.
REF = CONS
OFFSET = 0
TOKEN= Path(os.environ.get("GDC_TOKEN_FILE", os.path.expanduser("~/.gdc/token")))
SLICE= "https://api.gdc.cancer.gov/slicing/view"
# GRCh38 windows:
CLUSTER = "chr1:228605000-228650000"     # RNA5S cluster slice
CTRL    = "chr12:50000000-50100000"      # single-copy autosomal control window for CN
GENE    = (630, 748)
MIN_AD, MIN_VAF, MIN_CALLABLE_DP = 5, 0.003, 50          # calling thresholds
GOLD_AWK = HERE / "scripts" / "gold_rna_classify.awk"    # gold per-pair RNA classifier
EXCL = set(range(1, 90)) | set(range(2034, 2169))        # WGS: edges only (790-932/974-1057 real in WGS)
THREADS = 6

def sh(c, env=None): subprocess.run(c, shell=True, executable="/bin/bash", env=env)
def o(c, env=None):  return subprocess.run(c, shell=True, executable="/bin/bash",
                                           capture_output=True, text=True, env=env).stdout
def region(p): return "gene" if GENE[0] <= p <= GENE[1] else ("nts_pre" if p < GENE[0] else "nts_post")

def load_manifest(program):
    p = HERE / "inventory" / f"{program}_seq_manifest.tsv"
    hdr = None; rows = []
    for ln in open(p):
        f = ln.rstrip("\n").split("\t")
        if hdr is None: hdr = f; continue
        rows.append(dict(zip(hdr, f)))
    return rows

def slice_bam(fid, regions, dest, token, tries=5):
    """GDC slicing API -> mini BAM of the requested regions (+ unmapped).
    Retries transient GDC errors (429/5xx/timeouts)."""
    rq  = "&".join(f"region={r}" for r in regions)
    hdr = f"{dest}.hdr"
    # drive retries in Python (no curl --retry) so we can read GDC's Retry-After header.
    # --speed-limit/--speed-time abort only on a STALLED connection (not slow-but-progressing
    # big unmapped pulls); -D dumps headers, -w prints the HTTP status code to stdout.
    cmd = (f'curl -s -S -D {hdr} -o {dest} -w "%{{http_code}}" '
           f'--connect-timeout 30 --speed-limit 1000 --speed-time 120 '
           f'-H "X-Auth-Token: {token}" "{SLICE}/{fid}?{rq}" 2>/dev/null')
    for k in range(tries):
        code = (o(cmd) or "").strip()[-3:]
        if code == "200" and os.path.exists(dest) and os.path.getsize(dest) > 0:
            try: os.remove(hdr)
            except OSError: pass
            return True
        try: os.remove(dest)                 # drop the 403/error body so it can't masquerade as a slice
        except OSError: pass
        wait = None                          # honor GDC's Retry-After throttle signal if present
        try:
            for ln in open(hdr):
                if ln.lower().startswith("retry-after:"):
                    wait = int(ln.split(":", 1)[1].strip()); break
        except (OSError, ValueError): pass
        if wait is None: wait = min(2 ** k, 30)            # no signal -> exponential backoff
        wait += random.uniform(0, min(wait, 8) + 1)        # + jitter (de-synchronize shards)
        if k < tries - 1: time.sleep(wait)
    try: os.remove(hdr)
    except OSError: pass
    return False

def call_from_bam(bam, env):
    pile = o(f'bcftools mpileup -Ou -d 200000 -q 30 -Q 30 -a FORMAT/AD,FORMAT/DP -f {REF} {bam} 2>/dev/null '
             f'| bcftools query -f %POS"\\t"%REF"\\t"%ALT"\\t["%DP"]\\t["%AD"]\\n"', env=env)
    depths, calls = {}, []
    for ln in pile.splitlines():
        p = ln.split("\t")
        if len(p) < 5: continue
        try: pos = int(p[0]) + OFFSET; dp = int(p[3])   # OFFSET maps gene-only local pos -> consensus
        except: continue
        depths[pos] = dp
        if pos in EXCL: continue
        ad = [int(x) for x in p[4].split(",") if x != "."]
        if not ad: continue
        for i, a in enumerate(p[2].split(",")):
            if a in ("<*>", ".", ""): continue
            adi = ad[i+1] if i+1 < len(ad) else 0
            vaf = adi/dp if dp else 0
            if adi >= MIN_AD and vaf >= MIN_VAF:
                calls.append((pos, p[1], a, dp, adi, round(vaf, 5), region(pos)))
    return depths, calls, pile

def process(r, args, token, env):
    fid = r["file_id"]; samp = r["case"]; st = r["sample_type"].replace(" ", "_")
    tag = f'{samp}.{st}.{args.modality}'
    base = HERE / args.variants_dir / tag
    vtsv, qtsv, pgz = f"{base}.variants.tsv", f"{base}.qc.tsv", f"{base}.pileup.tsv.gz"
    if os.path.exists(vtsv) and os.path.exists(qtsv): return "skip"
    wd = HERE / "slices" / f"_tmp_{fid[:8]}"; wd.mkdir(parents=True, exist_ok=True)
    mini = wd / "slice.bam"
    # regions: WGS -> cluster (+unmapped); RNA -> all RNA5S loci from BED (+unmapped if --deep)
    if args.modality == "rna" and args.bed:
        regions = [ln.split()[0] + ":" + ln.split()[1] + "-" + ln.split()[2]
                   for ln in open(args.bed) if ln.strip() and not ln.startswith("#")]
    else:
        regions = [CLUSTER]
    if args.deep and args.modality != "wgs":   # WGS: never slice genome-wide unmapped (huge); cluster only
        regions = regions + ["unmapped"]
    if not slice_bam(fid, regions, mini, token):
        open(qtsv, "w").write(f"sample\tstatus\n{tag}\tslice_failed\n"); sh(f"rm -rf {wd}"); return "slice_failed"
    nRNA = nDNA = nSG = ""
    if args.modality == "rna":
        # 5S RNA extraction: paired reads -> k=17 kmer
        # prefilter (bbduk) -> bwa -p realign to consensus -> name-sort -> per-pair classifier
        # (gold_rna_classify.awk): keep pairs where BOTH mates are gene-contained (start>=GS, end<=GE)
        # AND fragment span<=SPAN. Drops precursor/gDNA/Alu/spacer pairs that the start-only filter leaked.
        GS, GE, SPAN = 630 - args.mate_tol, 748 + args.mate_tol, args.span_max
        sh(f'samtools collate -u -O {mini} 2>/dev/null | samtools fastq -1 {wd}/r1.fq -2 {wd}/r2.fq '
           f'-0 /dev/null -s /dev/null -n 2>/dev/null', env)
        nreads = int(o(f"wc -l < {wd}/r1.fq") or 0) // 2
        sh(f'bbduk.sh in1={wd}/r1.fq in2={wd}/r2.fq ref={REF} k=17 rcomp=t threads={THREADS} '
           f'-Xmx1g overwrite=t outm={wd}/m.fq 2>/dev/null', env)
        sh(f'bwa mem -B 6 -O 8 -L 5,5 -T 30 -t {THREADS} -p {REF} {wd}/m.fq 2>/dev/null '
           f'| samtools view -b -F 0x900 2>/dev/null | samtools sort -n -o {wd}/pe.bam - 2>/dev/null', env)
        sh(f"samtools view {wd}/pe.bam 2>/dev/null | gawk -v gs={GS} -v ge={GE} -v sp={SPAN} "
           f"-v RNAF={wd}/rna.txt -v GPF={wd}/gp.txt -v CNT={wd}/cnt.txt -f {GOLD_AWK} 2>/dev/null", env)
        nRNA, nDNA, nSG = (o(f"cat {wd}/cnt.txt 2>/dev/null").split() + ["0", "0", "0"])[:3]
        if os.path.exists(f"{wd}/rna.txt") and os.path.getsize(f"{wd}/rna.txt") > 0:
            sh(f"samtools view -N {wd}/rna.txt -b {wd}/pe.bam 2>/dev/null | samtools sort -o {wd}/d.bam - 2>/dev/null && samtools index {wd}/d.bam 2>/dev/null", env)
        else:
            sh(f"samtools view -H {wd}/pe.bam 2>/dev/null | samtools view -b -o {wd}/d.bam - 2>/dev/null", env)
    else:
        # WGS/WXS: paired + singleton realign to consensus (no transcript containment)
        sh(f'samtools collate -u -O {mini} 2>/dev/null | samtools fastq -1 {wd}/r1.fq -2 {wd}/r2.fq '
           f'-s {wd}/s.fq -0 /dev/null -n 2>/dev/null', env)
        nreads = (int(o(f"wc -l < {wd}/r1.fq") or 0) + int(o(f"wc -l < {wd}/s.fq") or 0)) // 4
        sh(f'( bwa mem -B 6 -O 8 -L 5,5 -T 30 -t {THREADS} {REF} {wd}/r1.fq {wd}/r2.fq 2>/dev/null; '
           f'  bwa mem -B 6 -O 8 -L 5,5 -T 30 -t {THREADS} {REF} {wd}/s.fq 2>/dev/null | grep -v "^@" ) '
           f'| samtools view -b -q 30 -F 0x904 - | samtools sort -o {wd}/d.bam 2>/dev/null', env)
        sh(f"samtools index {wd}/d.bam 2>/dev/null", env)
    if args.keep_bam:
        (HERE / "bam").mkdir(exist_ok=True)
        sh(f"cp {wd}/d.bam {HERE / 'bam'}/{tag}.bam && cp {wd}/d.bam.bai {HERE / 'bam'}/{tag}.bam.bai 2>/dev/null", env)
    # guard: an empty BAM (no aligned 5S reads) must NOT masquerade as ok med=0
    nmap = int(o(f"samtools view -c {wd}/d.bam 2>/dev/null", env).strip() or 0)
    if nmap == 0:
        with open(qtsv, "w") as f:
            f.write("sample\tproject\tsample_type\tmodality\tstatus\n")
            f.write(f"{samp}\t{r['project']}\t{r['sample_type']}\t{args.modality}\tempty_bam\n")
        sh(f"rm -rf {wd}"); return f"empty_bam (0 aligned reads; slice_reads={nreads})"
    depths, calls, pile = call_from_bam(f"{wd}/d.bam", env)
    nonx = [v for k, v in depths.items() if k not in EXCL and 90 <= k <= 2033]
    med = statistics.median(nonx) if nonx else 0
    callable_frac = round(sum(1 for v in nonx if v >= MIN_CALLABLE_DP) / len(nonx), 4) if nonx else 0
    # ── CN normalization (DNA only): single-copy control window divides out WGS sequencing depth.
    #    est_copies = median 5S depth / control mean depth  (= per-haplotype copies).
    ctrl_dep, est_copies = "", ""
    if args.modality in ("wgs", "wxs"):
        cb = wd / "ctrl.bam"
        if slice_bam(fid, [CTRL], cb, token):
            sh(f"samtools index {cb} 2>/dev/null", env)
            md = o(f"samtools depth -a -r {CTRL} {cb} 2>/dev/null "
                   f"| awk '{{s+=$3;n++}}END{{if(n)printf \"%.2f\",s/n; else print 0}}'", env)
            ctrl_dep = float(md or 0)
            est_copies = round(med / ctrl_dep, 2) if ctrl_dep else ""
    with open(vtsv, "w") as f:
        f.write("sample\tproject\tsample_type\tconsensus_pos\tref\talt\tdepth\talt_depth\tvaf\tregion\n")
        for c in sorted(calls):
            f.write(f"{samp}\t{r['project']}\t{r['sample_type']}\t" + "\t".join(map(str, c)) + "\n")
    with open(qtsv, "w") as f:
        f.write("sample\tproject\tsample_type\tmodality\tmedian_depth\tcontrol_depth\test_copies\t"
                "n_variants\tcallable_fraction\tslice_reads\tnRNA\tnDNA_NTS\tnSingleton\tstatus\n")
        f.write(f"{samp}\t{r['project']}\t{r['sample_type']}\t{args.modality}\t{med}\t{ctrl_dep}\t{est_copies}\t"
                f"{len(calls)}\t{callable_frac}\t{nreads}\t{nRNA}\t{nDNA}\t{nSG}\tok\n")
    import gzip as _gz
    with _gz.open(pgz, "wt") as g:
        g.write(pile)
    sh(f"rm -rf {wd}")
    return f"ok med={med} cn={est_copies} nvar={len(calls)} reads={nreads}"

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--modality", choices=["wgs", "rna", "wxs"], required=True)
    ap.add_argument("--program", default="TCGA")
    ap.add_argument("--projects", nargs="*", help="e.g. TCGA-OV TCGA-THCA (default: all in program)")
    ap.add_argument("--sample-types", nargs="*", default=None, help='e.g. "Primary Tumor" "Blood Derived Normal"')
    ap.add_argument("--cases", nargs="*", default=None, help="restrict to these case barcodes (e.g. match WGS pilot for the RNA arm)")
    ap.add_argument("--bed", default=None, help="RNA5S loci BED (rna modality)")
    ap.add_argument("--deep", action="store_true", help="also slice region=unmapped and include")
    ap.add_argument("--mate-in-region", default=None, help="RNA transcript mode: keep only read pairs whose FULL FRAGMENT lies in LO-HI (e.g. 630-748 = 5S gene); excludes Alu + 3'/5' precursor/gDNA overhang")
    ap.add_argument("--mate-tol", type=int, default=10, help="nt a read may reach past the gene borders (gold containment window = 630-tol .. 748+tol)")
    ap.add_argument("--span-max", type=int, default=130, help="gold: max fragment span (hi-lo) for an RNA-classified pair")
    ap.add_argument("--variants-dir", default="variants", help="output subdir for variants/qc/pileup (use a fresh dir to avoid skipping existing pileups)")
    ap.add_argument("--keep-bam", action="store_true", help="retain aligned 5S BAM (pre-mate-filter) in bam/ for future re-filtering without re-slicing")
    ap.add_argument("--ref", default=None, help="alignment reference (default full consensus; for RNA use reference/5S_gene_30flank.fa)")
    ap.add_argument("--ref-offset", type=int, default=0, help="add to local positions -> consensus coords (gene+30 ref = 599)")
    ap.add_argument("--max", type=int, default=0, help="cap samples (pilot)")
    ap.add_argument("--shard", default=None, help="i/N — process only rows where index%%N==i (parallel workers)")
    args = ap.parse_args()
    global REF, OFFSET
    if args.ref: REF = args.ref
    OFFSET = args.ref_offset
    if not TOKEN.exists():
        sys.exit(f"No GDC token at {TOKEN} — see README (portal login -> Download Token).")
    token = TOKEN.read_text().strip()
    strat = {"wgs": "WGS", "rna": "RNA-Seq", "wxs": "WXS"}[args.modality]
    rows = [r for r in load_manifest(args.program) if r["strategy"] == strat]
    if strat == "RNA-Seq": rows = [r for r in rows if r["rna_bam"] == "genomic"]
    if args.projects:     rows = [r for r in rows if r["project"] in args.projects]
    if args.sample_types: rows = [r for r in rows if r["sample_type"] in args.sample_types]
    if args.cases:        rows = [r for r in rows if r["case"] in args.cases]
    if args.max:          rows = rows[:args.max]
    if args.shard:
        si, sn = (int(x) for x in args.shard.split("/"))
        rows = [r for j, r in enumerate(rows) if j % sn == si]
    (HERE / args.variants_dir).mkdir(exist_ok=True); (HERE / "slices").mkdir(exist_ok=True)
    print(f"{len(rows)} {strat} samples; modality={args.modality} deep={args.deep}", flush=True)
    env = dict(os.environ)
    for i, r in enumerate(rows, 1):
        try: res = process(r, args, token, env)
        except Exception as e: res = f"ERROR {e}"
        print(f"[{i}/{len(rows)}] {r['case']} {r['sample_type']}: {res}", flush=True)
    print("SLICE_CALL_DONE", flush=True)

if __name__ == "__main__":
    main()
