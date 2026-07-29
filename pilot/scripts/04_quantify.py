#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Step B pilot — 04: quantify H3K27ac signal across the pilot bigWigs -> region x sample matrix.

Wraps deepTools `multiBigwigSummary` (env: atac_hdac).

Two modes:
  bins  (default) : genome-wide fixed bins (no reference bed needed) — the standard, dependency-free
                    batch-assessment substrate. Good enough to answer "batch vs lineage".
  bed             : quantify over an H3K27ac peak reference (SE-relevant). Needs --ref-bed, e.g. built
                    locally by 00_make_reference_regions.sh from the AllCell bed.

Runs on HPC compute (SLURM) or locally. Output results/signal.tab feeds 05_normalize_pca.py.
"""
import argparse, glob, os, subprocess, sys


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--bw-dir", default="../data/bigwigs")
    ap.add_argument("--selection", default="../data/selection.tsv")
    ap.add_argument("--out-prefix", default="../results/signal")
    ap.add_argument("--mode", choices=["bins", "bed"], default="bins")
    ap.add_argument("--bin-size", type=int, default=10000)
    ap.add_argument("--ref-bed", default="../data/reference_regions.bed", help="bed mode only")
    ap.add_argument("--threads", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", "4")))
    a = ap.parse_args()

    bws = sorted(glob.glob(os.path.join(a.bw_dir, "SRX*.bw")))
    if not bws:
        sys.exit(f"[04] no bigWigs in {a.bw_dir} — run 03_download_bigwigs.sh --go first")
    labels = [os.path.splitext(os.path.basename(b))[0] for b in bws]
    os.makedirs(os.path.dirname(os.path.abspath(a.out_prefix)), exist_ok=True)
    npz, tab = a.out_prefix + ".npz", a.out_prefix + ".tab"

    subcmd = "bins" if a.mode == "bins" else "BED-file"   # deepTools subcommand names
    cmd = ["multiBigwigSummary", subcmd, "-b", *bws, "--labels", *labels,
           "-o", npz, "--outRawCounts", tab, "-p", str(a.threads)]
    if a.mode == "bins":
        cmd += ["--binSize", str(a.bin_size)]
    else:
        if not os.path.exists(a.ref_bed):
            sys.exit(f"[04] --ref-bed {a.ref_bed} missing (run 00_make_reference_regions.sh locally)")
        cmd += ["--BED", a.ref_bed]

    print("[04] " + " ".join(cmd[:6]) + f" ... ({len(bws)} bigWigs, mode={a.mode})", file=sys.stderr)
    subprocess.run(cmd, check=True)
    print(f"[04] wrote {tab} and {npz}", file=sys.stderr)


if __name__ == "__main__":
    main()
