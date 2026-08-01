#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# 0_separation_script.py — extract 5S rDNA (and optionally other) reads from CRAM
# files into per-sample FASTQ using region/BED-targeted samtools view, in parallel.
# Author: Jonathan Bohlen
# (c) 2026 Jonathan Bohlen. Code accompanying Sengl et al. (2026),
# "Surveillance and selection of 5S ribosomal RNA genes in the human genome."
# Released under the MIT License; see LICENSE at the repository root.
# -----------------------------------------------------------------------------
import os
import sys
import subprocess
import argparse
import traceback
import gzip
import shutil
from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm

# ----------------------------------------------------------------------
# CONFIG
# ----------------------------------------------------------------------
_SCRIPT_DIR   = os.path.dirname(os.path.abspath(__file__))
_PIPELINE_DIR = os.path.dirname(_SCRIPT_DIR)
_INPUT_DIR    = os.path.join(_PIPELINE_DIR, "input")
_RESULTS_DIR  = os.path.join(_PIPELINE_DIR, "results")

CRAM_DIR   = os.environ.get("CRAM_DIR",  os.path.join(_INPUT_DIR,   "cram"))
REF        = os.environ.get("REF", "")
OUT_5S     = os.environ.get("OUT_5S",   os.path.join(_RESULTS_DIR, "fastq_5S"))
BED_5S     = os.environ.get("BED_5S",   os.path.join(_INPUT_DIR,   "5S_regions.bed"))
OUT_OTHER   = os.environ.get("OUT_OTHER", os.path.join(_RESULTS_DIR, "fastq_other"))
BED_OTHER   = os.environ.get("BED_OTHER", os.path.join(_INPUT_DIR,   "other_regions.bed"))

_RESUME    = os.environ.get("RESUME",    "0") == "1"
_SAM_T     = int(os.environ.get("STEP0_SAM_T", "0"))
_OTHER_ONLY = os.environ.get("OTHER_ONLY", "0") == "1"   # skip 5S, extract other only
_TQDM_POS  = int(os.environ.get("TQDM_POSITION", "0"))  # bar row (0=top, 1=below)

if not REF:
    sys.exit("[ERROR] REF is not set.")
if not _OTHER_ONLY and not os.path.exists(BED_5S):
    sys.exit(f"[ERROR] 5S BED file not found: {BED_5S}")
if _OTHER_ONLY and not os.path.exists(BED_OTHER):
    sys.exit(f"[ERROR] other BED file not found: {BED_OTHER}")

# Parse 5S BED → positional region string (only needed in 5S mode).
if not _OTHER_ONLY:
    with open(BED_5S) as _f:
        _rows = [l.split() for l in _f if l.strip() and not l.startswith("#")]
    _REGION = f"{_rows[0][0]}:{_rows[0][1]}-{_rows[0][2]}" if len(_rows) == 1 else None
    os.makedirs(OUT_5S, exist_ok=True)
else:
    _REGION = None
    os.makedirs(OUT_OTHER, exist_ok=True)

# ----------------------------------------------------------------------
# SAM COMMAND TEMPLATES
#
# Both paths use the same optimised flags:
#   -F 0x900           drop supplementary + secondary before decode
#   decode_md=0        skip MD/NM reconstruction (expensive for CRAM)
#   required_fields    decode only QNAME + SEQ + QUAL (0x601)
#
# 5S path  — positional region string: single index seek, maximum speed.
# other path — -M -L BED: multi-region index access across many regions.
#
# The two paths are never run together in one process invocation; they are
# run concurrently as separate background jobs so the slower multi-region scan
# (~125 seeks/CRAM) does not delay the 5S pipeline.
#
# Note on brace escaping: these strings are used with .format(), so literal
# {{ and }} are needed to produce { and } in the final shell command (the
# awk action block delimiters).  \" and \\n produce " and \n in the shell.
# ----------------------------------------------------------------------
_SAM_5S = (
    "samtools view -@ {sam_t} -F 0x900"
    " --input-fmt-option decode_md=0"
    " --input-fmt-option required_fields=0x601"
    " --reference '{ref}' '{cram}' {region} 2>/dev/null"
    " | awk '!/^@/{{print \"@\"$1\"\\n\"$10\"\\n+\\n\"$11}}'"
)

_SAM_OTHER = (
    "samtools view -@ {sam_t} -F 0x900"
    " --input-fmt-option decode_md=0"
    " --input-fmt-option required_fields=0x601"
    " -M -L '{bed}' --reference '{ref}' '{cram}' 2>/dev/null"
    " | awk '!/^@/{{print \"@\"$1\"\\n\"$10\"\\n+\\n\"$11}}'"
)


# ----------------------------------------------------------------------
# PER-SAMPLE WORKER
# ----------------------------------------------------------------------

def process_cram(cram_path):
    sample = os.path.basename(cram_path).replace(".cram", "")

    if _OTHER_ONLY:
        out = os.path.join(OUT_OTHER, f"{sample}.other.fq.gz")
        if _RESUME and os.path.exists(out):
            return sample
        cmd = _SAM_OTHER.format(sam_t=_SAM_T, bed=BED_OTHER, ref=REF, cram=cram_path)
        label = "other"
    else:
        out = os.path.join(OUT_5S, f"{sample}.5S.fq.gz")
        if _RESUME and os.path.exists(out):
            return sample
        region_arg = f"'{_REGION}'" if _REGION else f"-M -L '{BED_5S}'"
        cmd = _SAM_5S.format(sam_t=_SAM_T, ref=REF, cram=cram_path, region=region_arg)
        label = "5S"

    proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE,
                            executable="/bin/bash")
    with gzip.open(out, "wb", compresslevel=1) as f:
        shutil.copyfileobj(proc.stdout, f, length=1024 * 1024)
    proc.wait()
    if proc.returncode != 0:
        raise RuntimeError(
            f"{label} extraction failed for {sample} (exit {proc.returncode})"
        )
    return sample


# ----------------------------------------------------------------------
# PARALLEL EXECUTION
# ----------------------------------------------------------------------

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=4)
    args = parser.parse_args()

    cram_files = sorted([
        os.path.join(CRAM_DIR, f)
        for f in os.listdir(CRAM_DIR)
        if f.endswith(".cram")
    ])

    mode = "other-only" if _OTHER_ONLY else "5S"
    print(f"[INFO] Mode     : {mode}")
    print(f"[INFO] Found {len(cram_files)} CRAM files")
    print(f"[INFO] Workers  : {args.workers}  SAM_T: +{_SAM_T}")
    if _OTHER_ONLY:
        print(f"[INFO] other BED : {BED_OTHER}  →  {OUT_OTHER}")
    else:
        print(f"[INFO] 5S  BED  : {BED_5S}  →  {OUT_5S}")
    if _RESUME:
        print("[INFO] Resume mode: skipping samples whose output already exists")

    bar_desc = "other" if _OTHER_ONLY else "5S  "
    failed = []
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        futures = {ex.submit(process_cram, path): path for path in cram_files}
        with tqdm(total=len(futures), unit="file", dynamic_ncols=True,
                  desc=bar_desc, position=_TQDM_POS, leave=True) as bar:
            for fut in as_completed(futures):
                try:
                    fut.result()
                except Exception:
                    sample = os.path.basename(futures[fut])
                    tqdm.write(f"[ERROR] {sample}\n{traceback.format_exc()}")
                    failed.append(sample)
                bar.update(1)

    if failed:
        print(f"\n[WARN] {len(failed)} file(s) failed:")
        for s in failed:
            print(f"  {s}")
    print("[ALL DONE] FASTQ extraction complete.")
