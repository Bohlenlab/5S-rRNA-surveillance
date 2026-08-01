#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# 4_variant_summary.py — summarize per-sample 5S variant TSVs into a per-position
# table of mean/median alt depth and carrier count (parallel I/O).
# Author: Jonathan Bohlen
# (c) 2026 Jonathan Bohlen. Code accompanying Sengl et al. (2026),
# "Surveillance and selection of 5S ribosomal RNA genes in the human genome."
# Released under the MIT License; see LICENSE at the repository root.
# -----------------------------------------------------------------------------
"""
4_variant_summary.py: Summarize variants across samples by position/outcome,
average depth, median depth, and carrier count — with parallel I/O.

Usage:
    python 4_variant_summary.py /path/to/variants_dir MIN_AD /path/to/output.tsv [--workers N]
"""
import argparse
import os
import glob
import pandas as pd
from concurrent.futures import ProcessPoolExecutor, as_completed
from functools import partial
from tqdm import tqdm

def parse_one_file(path: str, min_ad: int):
    """
    Returns a per-sample dataframe:
        POS, REF, ALT, AD
    One row per ALT in this sample with AD > min_ad.
    """
    # Read raw TSV: POS, REF, ALT_raw, DP, AD_raw (no header)
    df = pd.read_csv(
        path, sep='\t', header=None,
        names=['POS','REF','ALT_raw','DP','AD_raw'],
        dtype=str
    )

    if df.empty:
        return pd.DataFrame(columns=['POS','REF','ALT','AD'])

    # Convert numeric
    df['POS'] = pd.to_numeric(df['POS'], errors='coerce').fillna(0).astype('int64')

    # Preprocess lists
    alt_lists = df['ALT_raw'].str.split(',').tolist()
    ad_lists  = df['AD_raw'].str.replace('"','', regex=False).str.split(',').tolist()

    recs = []
    for pos, ref, alts, ads in zip(df['POS'].tolist(), df['REF'].tolist(), alt_lists, ad_lists):
        if not alts or not ads or len(ads) < 2:
            continue
        for i, alt in enumerate(alts):
            try:
                ad = int(ads[i+1])  # AD_list[0] is ref depth
            except (IndexError, ValueError):
                ad = 0
            if ad > min_ad:
                recs.append((pos, ref, alt, ad))

    if not recs:
        return pd.DataFrame(columns=['POS','REF','ALT','AD'])

    return pd.DataFrame(recs, columns=['POS','REF','ALT','AD'])


def main(variants_dir: str, min_ad: int, output_file: str, workers: int | None):
    paths = sorted(glob.glob(os.path.join(variants_dir, '*.tsv')))
    if not paths:
        print(f"No .tsv files found in {variants_dir}")
        return

    if workers is None or workers <= 0:
        cpu = os.cpu_count() or 1
        workers = max(1, cpu - 1)

    mapper = partial(parse_one_file, min_ad=min_ad)

    per_file_results = []
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futures = {ex.submit(mapper, p): p for p in paths}
        with tqdm(total=len(futures), unit="file", dynamic_ncols=True) as bar:
            for fut in as_completed(futures):
                p = futures[fut]
                try:
                    df_small = fut.result()
                    if not df_small.empty:
                        per_file_results.append(df_small)
                except Exception as e:
                    tqdm.write(f"[ERROR] {os.path.basename(p)}: {e}")
                bar.update(1)

    if not per_file_results:
        print("No carrier variants found above threshold.")
        pd.DataFrame(columns=['POS','REF','ALT','avg_AD','median_AD','n_carriers']).to_csv(
            output_file, sep='\t', index=False
        )
        return

    combined = pd.concat(per_file_results, ignore_index=True)

    # Group by variant and calculate mean, median, carrier count
    summary = combined.groupby(['POS','REF','ALT'], as_index=False).agg(
        avg_AD=('AD','mean'),
        median_AD=('AD','median'),
        n_carriers=('AD','size')  # each AD row = one carrier
    )

    summary = summary.sort_values(['POS','REF','ALT'])

    summary.to_csv(output_file, sep='\t', index=False, float_format='%.6f')
    print(f"[OK] {len(paths)} files → {len(summary)} variant rows written to {output_file}")


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description="Summarize variant calls across samples (parallel).")
    ap.add_argument('variants_dir', help='Directory with per-sample variant TSVs')
    ap.add_argument('min_ad', type=int, help='Minimum alt depth threshold (AD > MIN_AD)')
    ap.add_argument('output_file', help='TSV path for summary output')
    ap.add_argument('--workers', type=int, default=None,
                    help='Number of parallel workers (default: CPU-1)')
    args = ap.parse_args()
    main(args.variants_dir, args.min_ad, args.output_file, args.workers)
