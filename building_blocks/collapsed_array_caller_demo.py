#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# collapsed_array_caller_demo.py — Minimal, self-contained demonstration of
# low-VAF variant calling in a collapsed tandem array.
# Author: Jonathan Bohlen
# (c) 2026 Jonathan Bohlen. Code accompanying Sengl et al. (2026),
# "Surveillance and selection of 5S ribosomal RNA genes in the human genome."
# Released under the MIT License; see LICENSE at the repository root.
# -----------------------------------------------------------------------------
"""
Concept demo for the collapsed-array short-read caller.

A tandem array of near-identical repeat copies cannot be resolved copy-by-copy
with short reads. The approach used in this study maps every read from the array
to a SINGLE consensus repeat unit, so that reads from all copies pool into one
per-position pileup, and then calls variants at low allele fraction: a variant
present in k of an individual's n copies appears at variant-allele fraction ~k/n,
which for a large array is a fraction of a percent, so a low calling threshold is
required.

This script reproduces that logic end-to-end on synthetic data in a fraction of a
second, with only the Python standard library and no external aligner or caller.
It simulates an array in which a few copies carry one substitution, pools reads
against the consensus unit, builds the pileup, calls at a 0.3% threshold, and
verifies that the planted low-fraction variant is recovered at the expected VAF
while sequencing error stays below threshold.

The production pipeline follows the same principle but maps real short reads with
bwa and calls with bcftools mpileup (see fig2_shortread_ukbb/). Run: `python collapsed_array_caller_demo.py`
"""
import random

BASES = "ACGT"
SEED = 7
N_COPIES = 100          # repeat copies in the array (one individual)
UNIT_LEN = 240          # length of the consensus repeat unit (bp)
READ_LEN = 100          # short-read length (bp)
READS_PER_COPY = 40     # reads simulated per copy
ERROR_RATE = 0.001      # per-base sequencing substitution rate
VAF_THRESHOLD = 0.003   # 0.3% calling threshold
MIN_ALT_DEPTH = 3       # minimum supporting reads


def make_consensus(rng):
    return [rng.choice(BASES) for _ in range(UNIT_LEN)]


def simulate(rng, consensus, variant_pos, variant_alt, n_variant_copies):
    """Build N_COPIES copies (n_variant_copies carry the variant) and emit reads
    against consensus coordinates. Returns a list of (offset, read_bases)."""
    copies = []
    for i in range(N_COPIES):
        seq = list(consensus)
        if i < n_variant_copies:
            seq[variant_pos] = variant_alt
        copies.append(seq)

    reads = []
    for seq in copies:
        for _ in range(READS_PER_COPY):
            offset = rng.randint(0, UNIT_LEN - READ_LEN)
            read = []
            for j in range(READ_LEN):
                b = seq[offset + j]
                if rng.random() < ERROR_RATE:                 # sequencing error
                    b = rng.choice([x for x in BASES if x != b])
                read.append(b)
            reads.append((offset, read))
    return reads


def pileup(reads):
    """Pool reads into a per-position base-count pileup on the consensus unit."""
    counts = [{b: 0 for b in BASES} for _ in range(UNIT_LEN)]
    for offset, read in reads:
        for j, b in enumerate(read):
            counts[offset + j][b] += 1
    return counts


def call_variants(counts, consensus):
    """Call alt alleles whose fraction exceeds VAF_THRESHOLD with >= MIN_ALT_DEPTH support."""
    calls = []
    for pos, c in enumerate(counts):
        depth = sum(c.values())
        if depth == 0:
            continue
        ref = consensus[pos]
        for alt in BASES:
            if alt == ref:
                continue
            ad = c[alt]
            vaf = ad / depth
            if vaf >= VAF_THRESHOLD and ad >= MIN_ALT_DEPTH:
                calls.append((pos, ref, alt, vaf, ad, depth))
    return calls


def main():
    rng = random.Random(SEED)
    consensus = make_consensus(rng)
    variant_pos = 123
    variant_alt = next(b for b in BASES if b != consensus[variant_pos])
    n_variant_copies = 3                      # 3 of 100 copies -> expected VAF 3%
    expected_vaf = n_variant_copies / N_COPIES

    reads = simulate(rng, consensus, variant_pos, variant_alt, n_variant_copies)
    counts = pileup(reads)
    calls = call_variants(counts, consensus)

    print(f"array copies={N_COPIES}  reads={len(reads)}  unit={UNIT_LEN} bp  "
          f"threshold={VAF_THRESHOLD:.1%}")
    print(f"planted: pos {variant_pos} {consensus[variant_pos]}>{variant_alt} "
          f"in {n_variant_copies}/{N_COPIES} copies (expected VAF {expected_vaf:.1%})")
    print("calls:")
    for pos, ref, alt, vaf, ad, depth in calls:
        print(f"  pos {pos}  {ref}>{alt}  VAF {vaf:.3%}  (alt {ad}/{depth})")

    planted = [c for c in calls if c[0] == variant_pos and c[2] == variant_alt]
    assert planted, "planted low-fraction variant was not recovered"
    assert abs(planted[0][3] - expected_vaf) < 0.01, "recovered VAF off expectation"
    assert all(c[0] == variant_pos for c in calls), "false-positive call above threshold"
    print("\nOK: planted 3% variant recovered; no false positives above threshold.")


if __name__ == "__main__":
    main()
