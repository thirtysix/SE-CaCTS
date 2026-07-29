#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Coarsen every retained 100bp genome archive to `--factor`x coarser bins (default 10 -> 1kb), IN A JOB.

Derives the coarse archive from the RETAINED 100bp one with no re-pull. 1kb is ~1000x finer than the
megabase-scale CN inference the archive exists for (PULL_DESIGN §8.5.6), at 1/10 the size (~72 GB -> ~10 GB).
Writes <SRX>.bin<newbp>.f16.gz + .json alongside the originals; does NOT delete the 100bp (verify first).

  python 52_coarsen_archives.py --out-root out --factor 10 [--jobs N]
"""
from __future__ import annotations
import argparse, glob, os, sys
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.environ.get("CNROSE") or os.path.join(HERE, "..", "..", "cnrose"))
from cnrose.io import read_binned_archive, coarsen_archive, write_binned_archive   # noqa: E402

FACTOR = 10

def one(path):
    v, meta = read_binned_archive(path)
    newbin = meta["binsize"] * FACTOR
    outp = path.replace(f".bin{meta['binsize']}.", f".bin{newbin}.")
    if os.path.exists(outp):
        return (path, 0, 0)
    cv, cl = coarsen_archive(v, meta["layout"], FACTOR)
    write_binned_archive(cv, cl, outp, binsize=newbin)
    return (path, os.path.getsize(path), os.path.getsize(outp))

def main():
    global FACTOR
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out-root", required=True)
    ap.add_argument("--factor", type=int, default=10)
    ap.add_argument("--jobs", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", 4)))
    a = ap.parse_args()
    FACTOR = a.factor
    arcs = sorted(glob.glob(os.path.join(a.out_root, "*", "*.bin100.f16.gz")))
    print(f"[52] {len(arcs)} archives -> {100*a.factor}bp, {a.jobs} workers", flush=True)
    ob = nb = done = 0
    with Pool(a.jobs) as pool:
        for i, (p, o, n) in enumerate(pool.imap_unordered(one, arcs, chunksize=8), 1):
            ob += o; nb += n; done += (n > 0)
            if i % 400 == 0:
                print(f"[52]   {i}/{len(arcs)}", flush=True)
    print(f"[52] done: {done} written ({len(arcs)-done} already existed)")
    if ob:
        print(f"[52] size: {ob/1e9:.1f} GB (100bp) -> {nb/1e9:.1f} GB ({100*a.factor}bp)  = {ob/nb:.1f}x smaller")

if __name__ == "__main__":
    main()
