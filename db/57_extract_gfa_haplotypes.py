#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# 57_extract_gfa_haplotypes.py — reconstruct per-haplotype assembly FASTAs from a Minigraph-Cactus full GFA.
# Author: Jonathan Bohlen
# (c) 2026 Jonathan Bohlen. Code accompanying Sengl et al. (2026),
# "Surveillance and selection of 5S ribosomal RNA genes in the human genome."
# Released under the MIT License; see LICENSE at the repository root.
# -----------------------------------------------------------------------------
"""
57 — Extract per-haplotype assembly FASTA from a Minigraph-Cactus *full* GFA.

The full graph is lossless: every input haplotype is stored as a W-line (walk)
over the S-line segments. Reconstructing a walk reproduces the original contig
sequence exactly (forward segments as-is, '<' segments reverse-complemented).

Single pass: S-lines precede W-lines in MC GFA output, so all segments are loaded
into memory and walks are reconstructed on the fly.

Usage:
  57_extract_gfa_haplotypes.py <graph.gfa[.gz]> <out_dir> [--limit N]
                               [--samples S1,S2,...] [--exclude CHM13,GRCh38,_MINIGRAPH_]
                               [--list-samples]

One FASTA per sample-haplotype is written: <out_dir>/<sample>.<hap>.fa
Headers: ><sample>#<hap>#<seqid>:<start>-<end>
"""
import sys, os, gzip, argparse, time

_RC = bytes.maketrans(b"ACGTNacgtnRYKMSWBDHVrykmswbdhv",
                      b"TGCANtgcanYRMKSWVHDByrmkswvhdb")
def revcomp(s):
    return s.translate(_RC)[::-1]

def opener(path):
    return gzip.open(path, "rb") if path.endswith(".gz") else open(path, "rb")

def parse_walk(walk):
    """Yield (segment_id_bytes, is_reverse) from a GFA walk string like b'>12>3<7'."""
    i, n = 0, len(walk)
    while i < n:
        orient = walk[i:i+1]          # b'>' or b'<'
        j = i + 1
        while j < n and walk[j:j+1] not in (b">", b"<"):
            j += 1
        yield walk[i+1:j], (orient == b"<")
        i = j

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("gfa")
    ap.add_argument("outdir")
    ap.add_argument("--limit", type=int, default=None,
                    help="reconstruct only the first N distinct samples (after excludes)")
    ap.add_argument("--samples", default=None, help="comma list of samples to keep")
    ap.add_argument("--exclude", default="CHM13,GRCh38,GRCh38_hg38,_MINIGRAPH_,CHM13v2",
                    help="comma list of sample names to skip (references)")
    ap.add_argument("--list-samples", action="store_true",
                    help="just scan W-lines and print the sample list, then exit")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    keep = set(a.samples.split(",")) if a.samples else None
    excl = set(x for x in a.exclude.split(",") if x)

    # --- optional: quick sample listing (still must stream, but skips loading segs) ---
    if a.list_samples:
        seen = {}
        t0 = time.time()
        with opener(a.gfa) as fh:
            for line in fh:
                if line[:2] == b"W\t":
                    smp = line.split(b"\t", 2)[1].decode()
                    seen[smp] = seen.get(smp, 0) + 1
        for s in sorted(seen):
            print(f"{s}\t{seen[s]}")
        print(f"# {len(seen)} samples, {time.time()-t0:.0f}s", file=sys.stderr)
        return

    seg = {}
    t0 = time.time()
    n_seg = 0
    chosen = []           # order of accepted samples (for --limit)
    chosen_set = set()
    open_files = {}       # (sample,hap) -> file handle
    bp_written = 0
    walks_written = 0

    def accept(smp):
        if smp in excl:
            return False
        if keep is not None:
            return smp in keep
        if a.limit is not None:
            if smp in chosen_set:
                return True
            if len(chosen_set) >= a.limit:
                return False
            chosen_set.add(smp); chosen.append(smp)
            return True
        return True

    with opener(a.gfa) as fh:
        for line in fh:
            tag = line[:2]
            if tag == b"S\t":
                f = line.rstrip(b"\n").split(b"\t")
                seg[f[1]] = f[2]
                n_seg += 1
                if n_seg % 10_000_000 == 0:
                    print(f"[seg] {n_seg:,} segments  {time.time()-t0:.0f}s", flush=True)
            elif tag == b"W\t":
                f = line.rstrip(b"\n").split(b"\t")
                # W <sample> <hapidx> <seqid> <start> <end> <walk>
                smp, hap, seqid, start, end, walk = f[1], f[2], f[3], f[4], f[5], f[6]
                smp_s = smp.decode()
                if not accept(smp_s):
                    continue
                parts = []
                for sid, rev in parse_walk(walk):
                    s = seg.get(sid)
                    if s is None:
                        sys.exit(f"ERROR: segment {sid!r} referenced in walk not found "
                                 f"(graph not fully loaded / not a full graph?)")
                    parts.append(revcomp(s) if rev else s)
                seq = b"".join(parts)
                key = (smp_s, hap.decode())
                fpath = os.path.join(a.outdir, f"{key[0]}.{key[1]}.fa")
                fhndl = open_files.get(key)
                if fhndl is None:
                    fhndl = open_files[key] = open(fpath, "ab")
                hdr = f">{smp_s}#{hap.decode()}#{seqid.decode()}:{start.decode()}-{end.decode()}\n"
                fhndl.write(hdr.encode())
                w = memoryview(seq)
                for i in range(0, len(seq), 60):
                    fhndl.write(w[i:i+60]); fhndl.write(b"\n")
                bp_written += len(seq); walks_written += 1
    for fhndl in open_files.values():
        fhndl.close()
    print(f"[done] {n_seg:,} segments loaded; wrote {walks_written} contigs "
          f"({bp_written:,} bp) across {len(open_files)} haplotypes in "
          f"{time.time()-t0:.0f}s", flush=True)
    print("[done] samples:", ",".join(sorted({k[0] for k in open_files})))

if __name__ == "__main__":
    main()
