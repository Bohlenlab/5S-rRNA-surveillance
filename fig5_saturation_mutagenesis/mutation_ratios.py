#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# mutation_ratios.py — Per-variant replicate mutation-frequency ratios versus a
# normalization sample, with pairwise and all-vs-all correlations.
# Author: Jonathan Bohlen
# (c) 2026 Jonathan Bohlen. Code accompanying Sengl et al. (2026),
# "Surveillance and selection of 5S ribosomal RNA genes in the human genome."
# Released under the MIT License; see LICENSE at the repository root.
# -----------------------------------------------------------------------------
"""
Process per-position mismatch tables for multiple biological replicates,
normalize to a reference (normalization) file, and output per-pair tables
of variant ratios for replicate 1 and replicate 2, plus correlation summaries.

Also writes all-vs-all correlation matrices across all samples/replicates:
  - correlation_matrix_pearson.csv
  - correlation_matrix_spearman.csv

Input file format (CSV, inferred from example):
- position (int)
- ref_base (str of 'A','C','G','T')
- total_coverage (int)
- optional: mismatch_total (int)
- A, C, G, T (int counts supporting each base at this position; includes ref)
- optional: deletions (int)

What the script does:
1) Reads a normalization file (same format) and computes per-position mutation
   frequencies for A/C/G/T **excluding wildtype** (i.e., alt != ref_base).
   Frequency = count_alt / denominator, where denominator is configurable
   (default: total_coverage). A small pseudocount avoids division-by-zero.
2) For every CSV in an input directory (matching --file-glob, default: *.csv),
   computes the same frequencies and then the ratio (sample_freq / norm_freq)
   for each (position, ref, alt) within an indicated range. Positions are also
   renumbered so that --new-start-pos becomes position 1.
3) Files are paired by biological replicate using a regex (default: r'_(\\d+)_').
   A pair is formed if we find replicate 1 and replicate 2 for the same "base name"
   (i.e., the filename with the matched replicate token removed). For each pair,
   we output a table with columns:
      variant_id (POS_ref_alt), position, ref, alt, ratio_rep1, ratio_rep2, ratio_mean
   where POS is the **renumbered** position by default (configurable).
4) We also write a summary CSV with Pearson and Spearman correlations for each 1↔2 pair.
5) Finally, we compute **all-vs-all** Pearson and Spearman correlation matrices across
   every sample/replicate, using pairwise NaN handling for overlapping variants only.

Usage:
    python process_rep_mutation_ratios.py \
        --input-dir /path/to/folder \
        --norm-file /path/to/normalization.csv \
        --range 50:140 \
        --new-start-pos 50 \
        --outdir outputs

Common tweaks:
- If your replicate indicator isn't like \"..._1_.../..._2_...\", change --rep-regex.
- If you'd like variant_id to use original genome positions, set --id-pos original.
- If you prefer a different frequency denominator, set --denominator accordingly.
"""
import argparse
import glob
import os
import re
from typing import Dict, Tuple

import numpy as np
import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(
        description="Compute per-variant replicate ratios vs. normalization and output per-pair tables + summary + all-vs-all matrices."
    )
    p.add_argument("--input-dir", required=True, help="Folder containing sample CSVs.")
    p.add_argument("--norm-file", required=True, help="Normalization CSV (same format).")
    p.add_argument("--range", required=True, help="Position range as START:END (inclusive). Uses ORIGINAL positions.")
    p.add_argument("--new-start-pos", type=int, required=True,
                   help="Original position that will be renumbered to 1.")
    p.add_argument("--outdir", default="outputs", help="Directory for outputs (created if missing).")
    p.add_argument("--file-glob", default="*.csv", help="Glob for input files, relative to --input-dir.")
    p.add_argument("--rep-regex", default=r"_(\d+)_",
                   help=r"Regex (with a capturing group) that extracts the replicate number, default: r'_(\d+)_'.")
    p.add_argument("--denominator", choices=["coverage", "mismatch", "nonref_total"], default="coverage",
                   help="Frequency denominator. 'coverage' = total_coverage; 'mismatch' = mismatch_total; "
                        "'nonref_total' = sum of alt counts (A/C/G/T except ref).")
    p.add_argument("--pseudocount", type=float, default=1e-9,
                   help="Small value added to normalization frequency to avoid division by zero.")
    p.add_argument("--id-pos", choices=["renumbered", "original"], default="renumbered",
                   help="Which position to use inside variant_id (POS_ref_alt).")
    p.add_argument("--exclude-norm-zero", action="store_true",
                   help="If set, drop variants where normalization frequency is 0 (after pseudocount).")
    return p.parse_args()


def load_table(path: str) -> pd.DataFrame:
    df = pd.read_csv(path)
    # Accept a flexible schema and synthesize missing columns when possible
    required_min = {"position", "ref_base", "total_coverage", "A", "C", "G", "T"}
    missing_min = required_min - set(df.columns)
    if missing_min:
        raise ValueError(f"{path} is missing required columns: {missing_min} (must at least have {sorted(required_min)})")

    # Ensure integer dtypes where sensible
    for col in ["position", "total_coverage"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce").astype("Int64")

    # Synthesize mismatch_total if absent: sum of non-ref counts
    if "mismatch_total" not in df.columns:
        base_ix = df["ref_base"].map({"A": 0, "C": 1, "G": 2, "T": 3}).astype("Int64")
        counts = df[["A", "C", "G", "T"]].to_numpy(dtype=float)
        idx = np.arange(len(df))
        ref_counts = counts[idx, base_ix.fillna(0).astype(int)]
        nonref_sum = counts.sum(axis=1) - ref_counts
        df["mismatch_total"] = nonref_sum.astype(int)

    # Provide deletions column if missing
    if "deletions" not in df.columns:
        df["deletions"] = 0

    return df


def compute_freqs(df: pd.DataFrame, denom_mode: str) -> pd.DataFrame:
    """
    Return long-form DataFrame with columns:
      position, ref, alt, count_alt, denom, freq
    where alt != ref.
    """
    bases = ["A", "C", "G", "T"]

    # Pre-compute ref counts efficiently
    base_ix = df["ref_base"].map({"A": 0, "C": 1, "G": 2, "T": 3}).astype("Int64")
    counts = df[["A", "C", "G", "T"]].to_numpy(dtype=float)
    idx = np.arange(len(df))
    ref_counts = counts[idx, base_ix.fillna(0).astype(int)]
    nonref_sum = counts.sum(axis=1) - ref_counts

    # Compute denominator by mode
    if denom_mode == "coverage":
        denom = df["total_coverage"].astype(float).to_numpy()
    elif denom_mode == "mismatch":
        denom = df["mismatch_total"].replace(0, np.nan).astype(float).to_numpy()
    elif denom_mode == "nonref_total":
        denom = pd.Series(nonref_sum).replace(0, np.nan).astype(float).to_numpy()
    else:
        raise ValueError("Unknown denominator mode")

    rows = []
    # use numpy for faster looping
    positions = df["position"].astype(int).to_numpy()
    refs = df["ref_base"].astype(str).to_numpy()
    for i in range(len(df)):
        pos = int(positions[i])
        ref = refs[i]
        d = denom[i]
        if not (d == d and d not in (0, None)):  # NaN or zero
            # still emit rows with NaN freq to preserve keys
            for alt_j, alt in enumerate(["A", "C", "G", "T"]):
                if alt == ref:
                    continue
                rows.append({"position": pos, "ref": ref, "alt": alt,
                             "count_alt": counts[i, alt_j], "denom": d, "freq": np.nan})
            continue
        for alt_j, alt in enumerate(["A", "C", "G", "T"]):
            if alt == ref:
                continue
            count_alt = counts[i, alt_j]
            freq = count_alt / d
            rows.append({"position": pos, "ref": ref, "alt": alt,
                         "count_alt": count_alt, "denom": d, "freq": freq})
    return pd.DataFrame(rows)


def apply_range_and_renumber(df: pd.DataFrame, start: int, end: int, new_start_pos: int) -> pd.DataFrame:
    """Filter to [start, end] (inclusive) on ORIGINAL positions and add renumbered_position."""
    df = df[(df["position"] >= start) & (df["position"] <= end)].copy()
    df["renumbered_position"] = df["position"] - new_start_pos + 1
    return df


def variant_id_series(df: pd.DataFrame, use_renum: bool) -> pd.Series:
    if use_renum:
        pos_series = df["renumbered_position"]
    else:
        pos_series = df["position"]
    return pos_series.astype(int).astype(str) + "_" + df["ref"] + "_" + df["alt"]


def base_key_from_filename(fname: str, rep_regex: str) -> Tuple[str, str]:
    """
    Return (base_key, replicate) where base_key is filename with the replicate token removed.
    Example: name = "Input_5S_1_trimFinal_aligned.sorted_mismatches.csv", rep_regex = r'_(\\d+)_'
             -> replicate='1', base_key='Input_5S__trimFinal_aligned.sorted_mismatches.csv'
    """
    m = re.search(rep_regex, fname)
    if not m:
        return fname, None
    rep = m.group(1)
    start, end = m.span(0)
    base_key = fname[:start] + fname[end:]
    return base_key, rep


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    # Parse range
    try:
        r_start, r_end = args.range.split(":")
        r_start, r_end = int(r_start), int(r_end)
    except Exception:
        raise SystemExit("--range must look like START:END (e.g., 50:140)")

    # Load normalization and compute freqs within range
    norm_df = load_table(args.norm_file)
    norm_freqs = compute_freqs(norm_df, args.denominator)
    norm_freqs = apply_range_and_renumber(norm_freqs, r_start, r_end, args.new_start_pos)

    # For easy lookup: dict[(pos, ref, alt)] -> (freq, renum_pos)
    norm_key = {}
    for _, row in norm_freqs.iterrows():
        key = (int(row["position"]), row["ref"], row["alt"])
        norm_key[key] = (row["freq"], int(row["renumbered_position"]))

    # Collect per-sample variant ratio tables
    all_files = sorted(glob.glob(os.path.join(args.input_dir, args.file_glob)))
    # Exclude the normalization file from samples if it sits in the same dir/glob
    all_files = [f for f in all_files if os.path.abspath(f) != os.path.abspath(args.norm_file)]
    if not all_files:
        raise SystemExit("No sample files found with the provided --input-dir and --file-glob.")

    # Map base_key -> {rep: filepath}
    rep_map: Dict[str, Dict[str, str]] = {}
    for f in all_files:
        fname = os.path.basename(f)
        base_key, rep = base_key_from_filename(fname, args.rep_regex)
        rep_map.setdefault(base_key, {})
        if rep is not None:
            rep_map[base_key][rep] = f
        else:
            rep_map[base_key][os.path.splitext(fname)[0]] = f

    # Process each file to compute ratios vs normalization
    per_file_ratio: Dict[str, pd.DataFrame] = {}
    long_for_matrix = []  # collect rows for pivot (variant_id, sample_label, ratio)
    for f in all_files:
        df = load_table(f)
        freqs = compute_freqs(df, args.denominator)
        freqs = apply_range_and_renumber(freqs, r_start, r_end, args.new_start_pos)

        # Join to normalization (by original pos/ref/alt)
        # Vectorized merge for speed
        freqs["_key"] = freqs["position"].astype(int).astype(str) + "_" + freqs["ref"] + "_" + freqs["alt"]
        norm_freqs["_key"] = norm_freqs["position"].astype(int).astype(str) + "_" + norm_freqs["ref"] + "_" + norm_freqs["alt"]
        merged = pd.merge(freqs, norm_freqs[["position", "ref", "alt", "freq", "renumbered_position", "_key"]],
                          on=["_key"], suffixes=("", "_norm"), how="left")

        # Use norm freq from merged, fall back to 0 when missing (with pseudocount later)
        norm_freq = merged["freq_norm"]
        sample_freq = merged["freq"]
        norm_adj = norm_freq.fillna(0).astype(float) + args.pseudocount
        ratio = sample_freq.astype(float) / norm_adj

        out_df = pd.DataFrame({
            "position": merged["position_norm"].fillna(merged["position"]).astype("Int64"),
            "renumbered_position": merged["renumbered_position"],
            "ref": merged["ref"],
            "alt": merged["alt"],
            "sample_freq": sample_freq,
            "norm_freq": norm_freq,
            "ratio": ratio
        })

        per_file_ratio[f] = out_df

        # For all-vs-all: build variant_id consistently (use renumbered unless --id-pos original)
        use_renum = (args.id_pos == "renumbered")
        out_df = out_df.copy()
        if use_renum:
            out_df["variant_id"] = out_df["renumbered_position"].astype(int).astype(str) + "_" + out_df["ref"] + "_" + out_df["alt"]
        else:
            out_df["variant_id"] = out_df["position"].astype(int).astype(str) + "_" + out_df["ref"] + "_" + out_df["alt"]

        sample_label = os.path.basename(f)
        long_for_matrix.append(out_df[["variant_id", "ratio"]].assign(sample=sample_label))

    # ==== Per-pair output and pairwise correlation summary (as before) ====
    summary_rows = []
    for base_key, reps in rep_map.items():
        f1 = reps.get("1")
        f2 = reps.get("2")
        if not (f1 and f2):
            continue

        r1 = per_file_ratio.get(f1, pd.DataFrame())
        r2 = per_file_ratio.get(f2, pd.DataFrame())

        # Prepare identifiers
        for df_part in (r1, r2):
            df_part["variant_id_renum"] = (
                df_part["renumbered_position"].astype(int).astype(str) + "_" +
                df_part["ref"] + "_" + df_part["alt"]
            )
            df_part["variant_id_orig"] = (
                df_part["position"].astype(int).astype(str) + "_" +
                df_part["ref"] + "_" + df_part["alt"]
            )

        id_col = "variant_id_renum" if args.id_pos == "renumbered" else "variant_id_orig"
        keep_cols = ["position", "renumbered_position", "ref", "alt"]

        r1_small = r1[[id_col, *keep_cols, "ratio"]].rename(columns={"ratio": "ratio_rep1"})
        r2_small = r2[[id_col, *keep_cols, "ratio"]].rename(columns={"ratio": "ratio_rep2"})

        merged = pd.merge(r1_small, r2_small, on=[id_col, *keep_cols], how="outer")
        out = merged.copy()
        out["variant_id"] = out[id_col]
        out["position"] = out["renumbered_position"].astype("Int64")

        # Mean across replicates
        out["ratio_mean"] = out[["ratio_rep1", "ratio_rep2"]].astype(float).mean(axis=1, skipna=True)

        out_final = out[["variant_id", "position", "ref", "alt", "ratio_rep1", "ratio_rep2", "ratio_mean"]].sort_values(
            by=["position", "ref", "alt"], kind="stable"
        )

        # Correlations on rows where both replicates have values
        r1_vals = out_final["ratio_rep1"].astype(float)
        r2_vals = out_final["ratio_rep2"].astype(float)
        mask = r1_vals.notna() & r2_vals.notna()
        n_common = int(mask.sum())
        if n_common >= 2:
            pearson_r = r1_vals[mask].corr(r2_vals[mask], method="pearson")
            spearman_rho = r1_vals[mask].corr(r2_vals[mask], method="spearman")
        else:
            pearson_r = float("nan")
            spearman_rho = float("nan")

        pair_name = re.sub(r"[^\w.-]+", "_", base_key.rstrip("_")).rstrip("_")
        out_path = os.path.join(args.outdir, f"{pair_name}_replicate_ratios.csv")
        out_final.to_csv(out_path, index=False)
        print(f"Wrote {out_path}  (rows={len(out_final)})")

        summary_rows.append({
            "pair_name": pair_name,
            "n_common_variants": n_common,
            "pearson_r": pearson_r,
            "spearman_rho": spearman_rho
        })

    # Write correlation summary if we computed any pairs
    if summary_rows:
        summary_df = pd.DataFrame(summary_rows)
        summary_path = os.path.join(args.outdir, "replicate_correlation_summary.csv")
        summary_df.to_csv(summary_path, index=False)
        print(f"Wrote {summary_path}  (pairs={len(summary_df)})")

    # ==== All-vs-all correlation matrices (efficient) ====
    if long_for_matrix:
        long_df = pd.concat(long_for_matrix, ignore_index=True)
        # Pivot to wide: rows = variant_id, cols = sample label, values = ratio
        wide = long_df.pivot_table(index="variant_id", columns="sample", values="ratio", aggfunc="mean")
        # Pearson (pairwise complete observations)
        pearson_mat = wide.corr(method="pearson", min_periods=2)
        spearman_mat = wide.corr(method="spearman", min_periods=2)
        pearson_path = os.path.join(args.outdir, "correlation_matrix_pearson.csv")
        spearman_path = os.path.join(args.outdir, "correlation_matrix_spearman.csv")
        pearson_mat.to_csv(pearson_path)
        spearman_mat.to_csv(spearman_path)
        print(f"Wrote {pearson_path} and {spearman_path}  (n_samples={wide.shape[1]})")

    print("Done.")


if __name__ == "__main__":
    main()
