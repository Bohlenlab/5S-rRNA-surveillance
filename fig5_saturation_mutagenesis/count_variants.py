#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# count_variants.py — Per-base mismatch/deletion and dinucleotide substitution
# counts from aligned BAMs over a reference region.
# Author: Jonathan Bohlen
# (c) 2026 Jonathan Bohlen. Code accompanying Sengl et al. (2026),
# "Surveillance and selection of 5S ribosomal RNA genes in the human genome."
# Released under the MIT License; see LICENSE at the repository root.
# -----------------------------------------------------------------------------
import os, sys, glob, csv, argparse
import pysam
import numpy as np
from multiprocessing import Pool

PAIR_LIST = [a+b for a in "ACGT" for b in "ACGT"]
BASE_TO_IDX = {-1: -1, ord('A'):0, ord('C'):1, ord('G'):2, ord('T'):3,
               ord('a'):0, ord('c'):1, ord('g'):2, ord('t'):3}

def parse_args():
    p = argparse.ArgumentParser(description="Fast per-base mismatch/deletion stats + dinucleotide substitution counts.")
    p.add_argument("--fasta", default=os.path.join(os.environ.get("FIVES_REFS", "refs"), "5S_rRNA.fa"))
    p.add_argument("--chrom", default="5S_rRNA")
    p.add_argument("--start", type=int, default=50)
    p.add_argument("--end",   type=int, default=200)             # [start, end)
    p.add_argument("--bam-glob", default=os.path.join(os.environ.get("FIVES_DATA", "data"), "*.sorted.bam"))
    p.add_argument("--min-bq", type=int, default=0)              # base quality threshold
    p.add_argument("--mapq", type=int, default=0)                # discard reads with MAPQ < this
    p.add_argument("--include-duplicates", action="store_true")  # exclude duplicates by default
    p.add_argument("--snv-end-trim", type=int, default=0,        # trim for SNV counting
                   help="Ignore bases within N aligned positions of either read end for SNV coverage/mismatches.")
    p.add_argument("--threads", type=int, default=1)             # parallelize across BAMs
    p.add_argument("--verbose", action="store_true")
    return p.parse_args()

def ensure_bai(bam):
    bai = bam + ".bai"
    if not os.path.exists(bai):
        pysam.index(bam)

def md_has_mismatch_in_region(read, region_start, region_end):
    md = read.get_tag("MD") if read.has_tag("MD") else None
    if md is None:
        return True
    rpos = read.reference_start
    i, L = 0, len(md)
    while i < L:
        j = i
        while j < L and md[j].isdigit():
            j += 1
        if j > i:
            rpos += int(md[i:j]); i = j
        if i >= L: break
        c = md[i]
        if c == '^':
            i += 1; j = i
            while j < L and md[j].isalpha(): j += 1
            rpos += (j - i); i = j; continue
        if region_start <= rpos < region_end:
            return True
        rpos += 1; i += 1
    return False

def process_one_bam(args_tuple):
    bam_path, fasta_path, chrom, start, end, min_bq, mapq, include_duplicates, snv_end_trim, verbose = args_tuple
    sample = os.path.splitext(bam_path)[0]
    ensure_bai(bam_path)
    with pysam.AlignmentFile(bam_path, "rb", require_index=True) as bam:
        if chrom not in bam.references:
            if verbose: print(f"[{sample}] SKIP: '{chrom}' not in BAM references.", flush=True)
            return
        with pysam.FastaFile(fasta_path) as fa:
            ref_len = fa.get_reference_length(chrom)
            if end > ref_len:
                if verbose: print(f"[{sample}] SKIP: region end {end} exceeds contig length {ref_len}.", flush=True)
                return
            ref_slice = fa.fetch(chrom, start, end).upper()

        region_len = end - start

        # ---------- SNV coverage counts ----------
        if snv_end_trim <= 0:
            # Fast path (no end-trim): use bam.count_coverage
            A, C, G, T = bam.count_coverage(chrom, start, end, quality_threshold=min_bq, read_callback="all")
            A = np.asarray(A, dtype=np.int64)
            C = np.asarray(C, dtype=np.int64)
            G = np.asarray(G, dtype=np.int64)
            T = np.asarray(T, dtype=np.int64)
        else:
            # Trimmed path: build A/C/G/T by walking aligned bases, skipping ends
            A = np.zeros(region_len, dtype=np.int64)
            C = np.zeros(region_len, dtype=np.int64)
            G = np.zeros(region_len, dtype=np.int64)
            T = np.zeros(region_len, dtype=np.int64)
            for read in bam.fetch(chrom, start, end):
                if read.is_unmapped: continue
                if read.mapping_quality < mapq: continue
                if (not include_duplicates) and read.is_duplicate: continue
                pairs = read.get_aligned_pairs(matches_only=True, with_seq=False)
                if not pairs: continue
                qseq = read.query_sequence
                qquals = read.query_qualities
                if qseq is None: continue
                last_idx = len(pairs) - 1
                for k, (qpos, rpos) in enumerate(pairs):
                    # skip ends within N aligned positions
                    if k < snv_end_trim or (last_idx - k) < snv_end_trim:
                        continue
                    if rpos is None or qpos is None:
                        continue
                    if rpos < start or rpos >= end:
                        continue
                    if qquals is not None and qquals[qpos] < min_bq:
                        continue
                    b = BASE_TO_IDX.get(ord(qseq[qpos]), -1)
                    if b == 0: A[rpos - start] += 1
                    elif b == 1: C[rpos - start] += 1
                    elif b == 2: G[rpos - start] += 1
                    elif b == 3: T[rpos - start] += 1
                    else:
                        continue

        # ---------- deletions + perfect reads ----------
        deletions = np.zeros(region_len, dtype=np.int64)
        fully_covering = set()
        reads_with_mismatch = set()

        for read in bam.fetch(chrom, start, end):
            if read.is_unmapped: continue
            if read.mapping_quality < mapq: continue
            if (not include_duplicates) and read.is_duplicate: continue

            rpos = read.reference_start
            cig = read.cigartuples or []
            for op, length in cig:
                if op == 2:  # D
                    s = max(rpos, start); e = min(rpos + length, end)
                    if e > s:
                        deletions[s-start:e-start] += 1
                    rpos += length
                elif op in (0,7,8,3):  # M,=,X,N
                    rpos += length
                else:
                    continue

            if read.reference_start <= start and read.reference_end >= end:
                fully_covering.add(read.query_name)
                if md_has_mismatch_in_region(read, start, end):
                    reads_with_mismatch.add(read.query_name)

        perfect_reads = len(fully_covering - reads_with_mismatch)

        # ---------- dinucleotide substitution counts ----------
        n_windows = max(0, region_len - 1)
        if n_windows == 0:
            dinuc_counts = np.zeros((0,16), dtype=np.int64)
            dinuc_cov = np.zeros((0,), dtype=np.int64)
        else:
            dinuc_counts = np.zeros((n_windows,16), dtype=np.int64)
            dinuc_cov = np.zeros((n_windows,), dtype=np.int64)

            for read in bam.fetch(chrom, start, end):
                if read.is_unmapped: continue
                if read.mapping_quality < mapq: continue
                if (not include_duplicates) and read.is_duplicate: continue

                pairs = read.get_aligned_pairs(matches_only=True, with_seq=False)
                if not pairs: continue
                qseq = read.query_sequence
                qquals = read.query_qualities
                if qseq is None: continue

                for k in range(len(pairs) - 1):
                    q1, r1 = pairs[k]
                    q2, r2 = pairs[k+1]
                    if r1 is None or r2 is None:
                        continue
                    if r2 != r1 + 1:
                        continue
                    if r1 < start or r2 >= end:
                        continue
                    if qquals is not None:
                        if qquals[q1] < min_bq or qquals[q2] < min_bq:
                            continue
                    b1 = BASE_TO_IDX.get(ord(qseq[q1]), -1)
                    if b1 < 0: continue
                    b2 = BASE_TO_IDX.get(ord(qseq[q2]), -1)
                    if b2 < 0: continue
                    idx = (r1 - start)
                    dinuc_counts[idx, b1*4 + b2] += 1
                    dinuc_cov[idx] += 1

        # ---------- write outputs ----------
        out1 = f"{sample}_mismatches.csv"
        with open(out1, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=[
                "position","ref_base","total_coverage","mismatch_total","A","C","G","T","deletions"
            ])
            w.writeheader()
            for i in range(region_len):
                ref_base = ref_slice[i]
                total = int(A[i] + C[i] + G[i] + T[i])
                ref_cov = {"A":int(A[i]), "C":int(C[i]), "G":int(G[i]), "T":int(T[i])}.get(ref_base, 0)
                mism = max(0, total - ref_cov)
                w.writerow({
                    "position": start + i,
                    "ref_base": ref_base,
                    "total_coverage": total,
                    "mismatch_total": mism,
                    "A": int(A[i]), "C": int(C[i]), "G": int(G[i]), "T": int(T[i]),
                    "deletions": int(deletions[i]),
                })

        out2 = f"{sample}_region_perfect_reads.csv"
        with open(out2, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=["region_start","region_end","reads_fully_covering","perfect_reads"])
            w.writeheader()
            w.writerow({
                "region_start": start,
                "region_end": end,
                "reads_fully_covering": len(fully_covering),
                "perfect_reads": perfect_reads,
            })

        out3 = f"{sample}_dinuc_substitutions.csv"
        dinuc_fields = ["position","ref_dinuc","total_pairs"] + PAIR_LIST + ["mismatch_both","first_only","second_only","match_both"]
        with open(out3, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=dinuc_fields)
            w.writeheader()
            for i in range(max(0, region_len-1)):
                ref1 = ref_slice[i]; ref2 = ref_slice[i+1]
                row = {
                    "position": start + i,
                    "ref_dinuc": ref1 + ref2,
                    "total_pairs": int(dinuc_cov[i]),
                }
                counts_row = dinuc_counts[i] if n_windows > 0 else np.zeros(16, dtype=np.int64)
                for j, pair in enumerate(PAIR_LIST):
                    row[pair] = int(counts_row[j])

                mismatch_both = 0; first_only = 0; second_only = 0; match_both = 0
                for j, pair in enumerate(PAIR_LIST):
                    c = int(counts_row[j])
                    p0, p1 = pair[0], pair[1]
                    if p0 == ref1 and p1 == ref2: match_both += c
                    elif p0 != ref1 and p1 == ref2: first_only += c
                    elif p0 == ref1 and p1 != ref2: second_only += c
                    else: mismatch_both += c
                row["mismatch_both"] = mismatch_both
                row["first_only"] = first_only
                row["second_only"] = second_only
                row["match_both"]  = match_both
                w.writerow(row)

        if verbose:
            print(f"[{sample}] wrote {out1}, {out2}, {out3} | windows={max(0, region_len-1)}", flush=True)

def main():
    a = parse_args()
    bam_files = sorted(glob.glob(a.bam_glob))
    if not bam_files:
        print(f"No BAMs matched: {a.bam_glob}")
        return
    work = [(b, a.fasta, a.chrom, a.start, a.end, a.min_bq, a.mapq, a.include_duplicates, a.snv_end_trim, a.verbose)
            for b in bam_files]
    if a.threads > 1 and len(bam_files) > 1:
        with Pool(processes=a.threads) as pool:
            pool.map(process_one_bam, work)
    else:
        for item in work:
            process_one_bam(item)

if __name__ == "__main__":
    main()
