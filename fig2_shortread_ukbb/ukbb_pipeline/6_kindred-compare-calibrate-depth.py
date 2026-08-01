#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# 6_kindred-compare-calibrate-depth.py — for each trio, compute the fraction of a
# child's 5S variants also detected in each parent across multiple AD cutoffs.
# Author: Jonathan Bohlen
# (c) 2026 Jonathan Bohlen. Code accompanying Sengl et al. (2026),
# "Surveillance and selection of 5S ribosomal RNA genes in the human genome."
# Released under the MIT License; see LICENSE at the repository root.
# -----------------------------------------------------------------------------
"""
6_kindred-compare-calibrate-depth.py

For each child in a trio table, compute the fraction of the child's variants
(defined by POS,REF,ALT with AD >= cutoff) that are also detected in:
  - parent1,
  - parent2,
  - either parent,
  - both parents.

Multi-cutoff design:
- Multiple AD thresholds via -m/--cutoffs.
- Optional rarity filter applied globally before comparisons:
    --max-sample-count N and --rare-min-ad R:
      keep only variants that appear (AD >= R) in < N samples total.
- Optional per-child detailed tables at one chosen cutoff:
    --shared-cutoff C and --shared-depth D
    (only variants with child AD >= D and parent AD >= D are kept in the table).

Output:
- A wide summary TSV with columns like:
    child_id, parent1_id, parent2_id,
    pct_in_parent1_AD>10, pct_in_parent2_AD>10, pct_in_either_AD>10, pct_in_both_AD>10,
    n_child_AD>10, n_in_parent1_AD>10, ... (repeated for each cutoff)

Usage:
  python 6_kindred-compare-calibrate-depth.py trios.tsv /path/to/variants_dir /path/to/output_dir \
      -m 10 20 50 --shared-cutoff 20 --shared-depth 10 --max-sample-count 24 --rare-min-ad 10 -t 8
"""

import argparse
import os
import glob
from concurrent.futures import ProcessPoolExecutor, as_completed
from collections import defaultdict
import itertools
import pandas as pd
from tqdm import tqdm

# ---------------- File discovery ----------------
def find_sample_path(variants_dir: str, sample_id: str) -> str | None:
    exact = os.path.join(variants_dir, f"{sample_id}.tsv")
    if os.path.exists(exact):
        return exact
    matches = glob.glob(os.path.join(variants_dir, f"{sample_id}*.tsv"))
    return matches[0] if matches else None

def build_sample_path_map(variants_dir: str, ids) -> dict[str, str]:
    mp = {}
    for sid in ids:
        p = find_sample_path(variants_dir, sid)
        if p:
            mp[sid] = p
    return mp

# ---------------- Parsing & caching ----------------
def _parse_sample_file(path: str) -> pd.DataFrame:
    """
    Return long-form per-sample rows with explicit AD per ALT:
      columns: POS (int64), REF (str), ALT (str), DP (int64), AD (int64)
    """
    df = pd.read_csv(
        path, sep="\t", header=None,
        names=["POS","REF","ALT_raw","DP","AD_raw"],
        dtype=str
    )
    if df.empty:
        return pd.DataFrame(columns=["POS","REF","ALT","DP","AD"])
    df["POS"] = pd.to_numeric(df["POS"], errors="coerce").fillna(0).astype("int64")
    df["DP"]  = pd.to_numeric(df["DP"],  errors="coerce").fillna(0).astype("int64")

    df["ALT_list"] = df["ALT_raw"].str.split(",")
    df["AD_list"]  = df["AD_raw"].str.replace('"','', regex=False).str.split(",")

    df = df[["POS","REF","DP","ALT_list","AD_list"]].reset_index(drop=False).rename(columns={"index":"row_id"})
    df = df.explode("ALT_list", ignore_index=True)
    df["alt_idx"] = df.groupby("row_id").cumcount()

    def pick_ad(row):
        ads = row["AD_list"]
        i = row["alt_idx"] + 1  # index 0 is ref depth
        try:
            return int(ads[i]) if i < len(ads) else 0
        except Exception:
            return 0

    df["AD"] = df.apply(pick_ad, axis=1)
    df = df.rename(columns={"ALT_list":"ALT"})[["POS","REF","ALT","DP","AD"]]
    return df

# process-local cache: sample_id -> parsed DF
_CACHE: dict[str, pd.DataFrame] = {}

def load_sample_df(sample_id: str, sample_paths: dict[str,str]) -> pd.DataFrame:
    df = _CACHE.get(sample_id)
    if df is None:
        df = _parse_sample_file(sample_paths[sample_id])
        _CACHE[sample_id] = df
    return df

# ---------------- Rarity pre-pass (parallel) ----------------
def _present_variant_keys_for_count(sample_id: str, sample_paths: dict[str,str], min_ad: int):
    df = load_sample_df(sample_id, sample_paths)
    present = df.loc[df["AD"] >= min_ad, ["POS","REF","ALT"]].drop_duplicates()
    return [tuple(x) for x in present.itertuples(index=False, name=None)]

def build_variant_counts_parallel(sample_ids, sample_paths, min_ad: int, workers: int):
    counts = defaultdict(int)
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futs = {ex.submit(_present_variant_keys_for_count, sid, sample_paths, min_ad): sid
                for sid in sample_ids if sid in sample_paths}
        with tqdm(total=len(futs), desc="rarity pass", unit="sample", dynamic_ncols=True) as bar:
            for fut in as_completed(futs):
                for key in fut.result():
                    counts[key] += 1
                bar.update(1)
    return counts

# ---------------- Core per-(child, cutoff) task ----------------
def process_child_cutoff(child_id: str, p1_id: str, p2_id: str, cutoff: int,
                         shared_cutoff: int | None, shared_depth: int,
                         output_dir: str | None, sample_paths: dict[str,str],
                         rare_key_set: set[str] | None):
    df_child = load_sample_df(child_id, sample_paths).copy()
    df_p1    = load_sample_df(p1_id,   sample_paths).copy() if p1_id in sample_paths else pd.DataFrame(columns=df_child.columns)
    df_p2    = load_sample_df(p2_id,   sample_paths).copy() if p2_id in sample_paths else pd.DataFrame(columns=df_child.columns)

    # Filter by cutoff
    df_child = df_child[df_child["AD"] >= cutoff]
    df_p1    = df_p1[df_p1["AD"] >= cutoff]
    df_p2    = df_p2[df_p2["AD"] >= cutoff]

    # Optional rarity filter (applied to the child's keys only)
    if rare_key_set is not None and not df_child.empty:
        child_keys = (df_child["POS"].astype(str) + "|" +
                      df_child["REF"].astype(str) + "|" +
                      df_child["ALT"].astype(str))
        df_child = df_child[child_keys.isin(rare_key_set)]

    # Unique child variant set
    child_keys_df = df_child[["POS","REF","ALT"]].drop_duplicates()
    n_child = len(child_keys_df)

    if n_child == 0:
        return {
            "child_id": child_id, "parent1_id": p1_id, "parent2_id": p2_id, "cutoff": cutoff,
            "n_child": 0, "n_in_parent1": 0, "n_in_parent2": 0, "n_in_either": 0, "n_in_both": 0,
            "pct_in_parent1": float("nan"), "pct_in_parent2": float("nan"),
            "pct_in_either": float("nan"), "pct_in_both": float("nan"),
        }

    # Build parent presence sets at this cutoff
    p1_keys = set(map(tuple, df_p1[["POS","REF","ALT"]].drop_duplicates().itertuples(index=False, name=None)))
    p2_keys = set(map(tuple, df_p2[["POS","REF","ALT"]].drop_duplicates().itertuples(index=False, name=None)))
    child_keys = set(map(tuple, child_keys_df.itertuples(index=False, name=None)))

    in_p1 = child_keys & p1_keys
    in_p2 = child_keys & p2_keys
    in_both = in_p1 & in_p2
    in_either = in_p1 | in_p2

    n_p1 = len(in_p1)
    n_p2 = len(in_p2)
    n_both = len(in_both)
    n_either = len(in_either)

    # Optional per-child table at shared_cutoff: keep variants meeting shared_depth in child and parent(s)
    if shared_cutoff is not None and cutoff == shared_cutoff and output_dir is not None:
        # gather detailed rows from child with AD>=shared_depth and mark parent depths at >=shared_depth
        child_sd = load_sample_df(child_id, sample_paths)
        p1_sd = load_sample_df(p1_id, sample_paths) if p1_id in sample_paths else pd.DataFrame(columns=child_sd.columns)
        p2_sd = load_sample_df(p2_id, sample_paths) if p2_id in sample_paths else pd.DataFrame(columns=child_sd.columns)

        child_sd = child_sd[child_sd["AD"] >= shared_depth][["POS","REF","ALT","AD"]].rename(columns={"AD":"AD_child"})
        p1_sd    = p1_sd[p1_sd["AD"] >= shared_depth][["POS","REF","ALT","AD"]].rename(columns={"AD":"AD_parent1"})
        p2_sd    = p2_sd[p2_sd["AD"] >= shared_depth][["POS","REF","ALT","AD"]].rename(columns={"AD":"AD_parent2"})

        # Only consider variants that are in the child set at the main cutoff (after rarity filter)
        det = child_keys_df.merge(child_sd, on=["POS","REF","ALT"], how="left")
        det = det.merge(p1_sd, on=["POS","REF","ALT"], how="left")
        det = det.merge(p2_sd, on=["POS","REF","ALT"], how="left")

        out_dir = os.path.join(output_dir, f"shared_AD{cutoff}")
        os.makedirs(out_dir, exist_ok=True)
        out_path = os.path.join(out_dir, f"{child_id}_{p1_id}_{p2_id}_AD{cutoff}.tsv")
        det.to_csv(out_path, sep="\t", index=False)

    return {
        "child_id": child_id,
        "parent1_id": p1_id,
        "parent2_id": p2_id,
        "cutoff": cutoff,
        "n_child": n_child,
        "n_in_parent1": n_p1,
        "n_in_parent2": n_p2,
        "n_in_either": n_either,
        "n_in_both": n_both,
        "pct_in_parent1": n_p1 / n_child if n_child else float("nan"),
        "pct_in_parent2": n_p2 / n_child if n_child else float("nan"),
        "pct_in_either": n_either / n_child if n_child else float("nan"),
        "pct_in_both": n_both / n_child if n_child else float("nan"),
    }

# ---------------- Main ----------------
def main():
    ap = argparse.ArgumentParser(description="Trio parent presence across multiple AD cutoffs (parallel).")
    ap.add_argument("trio_file", help="TSV with columns including IID, parent1_ID, parent2_ID")
    ap.add_argument("variants_dir", help="Directory of per-sample variant TSVs")
    ap.add_argument("output_dir", help="Directory to write outputs")
    ap.add_argument("-m","--cutoffs", nargs="+", type=int, required=True,
                    help="List of AD thresholds (e.g. -m 10 20 50)")
    ap.add_argument("--shared-cutoff", type=int, default=None,
                    help="Cutoff at which to write per-child detailed tables")
    ap.add_argument("--shared-depth", type=int, default=10,
                    help="Minimum AD in child/parent to include in shared table")
    ap.add_argument("--max-sample-count", type=int, default=None,
                    help="Keep only variants that appear in < N samples (rarity filter)")
    ap.add_argument("--rare-min-ad", type=int, default=10,
                    help="AD threshold for counting presence in rarity pass")
    ap.add_argument("-t","--threads", type=int, default=1,
                    help="Parallel workers")
    ap.add_argument("--require-both-in-sample", action="store_true",
                    help="Only evaluate rows where child and both parents have TSVs")
    args = ap.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    trios = pd.read_csv(args.trio_file, sep="\t", dtype=str)
    needed = ["IID","parent1_ID","parent2_ID"]
    miss = [c for c in needed if c not in trios.columns]
    if miss:
        raise ValueError(f"Trio file must contain columns {needed}; missing: {miss}")

    rows = trios[needed].dropna()
    all_ids = pd.unique(rows.values.ravel())
    sample_paths = build_sample_path_map(args.variants_dir, all_ids)

    if args.require_both_in_sample:
        rows = rows[rows.apply(lambda r: r["IID"] in sample_paths
                                         and r["parent1_ID"] in sample_paths
                                         and r["parent2_ID"] in sample_paths, axis=1)]
    if rows.empty:
        raise SystemExit("No trios to process after filtering / file lookup.")

    # Rarity pass (global, once), if requested
    if args.max_sample_count is not None:
        workers = args.threads if args.threads and args.threads > 0 else max(1, (os.cpu_count() or 1)-1)
        counts = build_variant_counts_parallel(all_ids, sample_paths, args.rare_min_ad, workers)
        rare_list = [k for k, c in counts.items() if c < args.max_sample_count]
        rare_key_set = set(f"{pos}|{ref}|{alt}" for (pos, ref, alt) in rare_list)
        print(f"Rarity filter: {len(rare_key_set)} variants remain (<{args.max_sample_count} samples at AD ≥ {args.rare_min_ad}).")
    else:
        rare_key_set = None

    # Build tasks: (child, parents, cutoff)
    tasks = [(r.IID, r.parent1_ID, r.parent2_ID, cutoff) for r, cutoff in
             itertools.product(rows.itertuples(index=False), args.cutoffs)]

    workers = args.threads if args.threads and args.threads > 0 else max(1, (os.cpu_count() or 1)-1)
    results = []
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futs = {
            ex.submit(process_child_cutoff, child, p1, p2, cutoff,
                      args.shared_cutoff, args.shared_depth, args.output_dir,
                      sample_paths, rare_key_set): (child, cutoff)
            for (child, p1, p2, cutoff) in tasks
        }
        with tqdm(total=len(futs), desc="trio analysis", unit="task", dynamic_ncols=True) as bar:
            for fut in as_completed(futs):
                try:
                    results.append(fut.result())
                except Exception as e:
                    child, cutoff = futs[fut]
                    tqdm.write(f"[ERROR] {child} AD>{cutoff}: {e}")
                bar.update(1)

    if not results:
        raise SystemExit("No results produced.")

    # Long to wide: one row per child, columns per cutoff
    df = pd.DataFrame(results)
    # keep unique parent IDs per child (in case duplicates)
    parents = rows.drop_duplicates(subset=["IID"]).set_index("IID")[["parent1_ID","parent2_ID"]]

    # Build wide blocks for each metric
    def wide(metric):
        w = df.pivot(index="child_id", columns="cutoff", values=metric)
        w.columns = [f"{metric}_AD>{c}" for c in w.columns]
        return w

    metrics = ["pct_in_parent1","pct_in_parent2","pct_in_either","pct_in_both",
               "n_child","n_in_parent1","n_in_parent2","n_in_either","n_in_both"]
    blocks = [wide(m) for m in metrics]
    final = pd.concat([parents, *blocks], axis=1)
    final = final.reset_index().rename(columns={"index":"child_id"})

    out_summary = os.path.join(args.output_dir, "trio_parent_presence_summary.tsv")
    final.to_csv(out_summary, sep="\t", index=False, float_format="%.6f")
    print(f"✅ Summary written to {out_summary}")

if __name__ == "__main__":
    main()
