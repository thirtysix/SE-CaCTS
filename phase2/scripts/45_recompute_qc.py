#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Recompute qc.json for already-pulled samples FROM THE RETAINED ARTIFACTS — no bigWig, no re-download.

This is the retention set paying off (PULL_DESIGN.md §8.5): the genome-wide archive + the grid f32 column +
the grid BED are exactly the inputs `qc_from_retained` needs, so a QC-code fix (e.g. the union-grid
dynamic-range / FRiP fixes, 2026-07-21) can be applied to samples pulled under the old code without touching
ChIP-Atlas again. bigwig_header / bigwig_bytes and the caller-added n_peaks/n_super/tangent_cutoff are
preserved from the existing qc.json.

  python 45_recompute_qc.py --out-root <pull out/> --grid <grid.20.bed> --grid-label grid.20
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "cnrose"))
sys.path.insert(0, os.path.join(os.environ.get("CNROSE", ""), "") or os.path.join(HERE, "..", "..", "cnrose"))
from cnrose.io import qc_from_retained, read_bed3, read_binned_archive       # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out-root", required=True, help="the pull output tree (out/, sharded out/<last2>/)")
    ap.add_argument("--grid", required=True, help="the grid BED matching --grid-label")
    ap.add_argument("--grid-label", default="grid.20")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--force", action="store_true",
                    help="recompute even samples already written by the fixed code (default: skip them, so "
                         "this is safe to run while a pull is in progress — only pre-fix qc.json is touched)")
    a = ap.parse_args()

    grid = read_bed3(a.grid, canonical_only=True)
    n_ok = n_skip = 0
    for qpath in sorted(glob.glob(os.path.join(a.out_root, "*", "*.qc.json"))):
        srx = os.path.basename(qpath)[:-len(".qc.json")]
        d = os.path.dirname(qpath)
        # skip samples already produced by the fixed code (they have the _raw field), unless --force. Makes
        # this safe to run mid-pull: only pre-fix qc.json is rewritten.
        if not a.force:
            try:
                if "dynamic_range_p99_over_median_raw" in json.load(open(qpath)):
                    n_skip += 1
                    continue
            except (OSError, ValueError):
                pass
        f32 = os.path.join(d, f"{srx}.{a.grid_label}.f32")
        arch = os.path.join(d, f"{srx}.bin100.f16.gz")
        if not (os.path.exists(f32) and os.path.exists(arch)):
            print(f"  skip {srx}: missing {a.grid_label}.f32 or archive")
            n_skip += 1
            continue
        gsig = np.fromfile(f32, dtype="<f4")
        binned, meta = read_binned_archive(arch)
        stats = qc_from_retained(binned, meta["binsize"], grid_regions=grid, grid_signal=gsig,
                                 layout=meta["layout"])
        old = json.load(open(qpath))
        # preserve bigWig-only + caller-added fields; overwrite the recomputed statistics
        keep = {k: old[k] for k in ("bigwig_bytes", "bigwig_header", "n_peaks", "n_regions", "n_super",
                                    "tangent_cutoff", "n_super_cn", "cn_cutoff") if k in old}
        merged = {**keep, **stats, "qc_recomputed_from_retained": True}
        if a.dry_run:
            print(f"  {srx}: dyn_range {old.get('dynamic_range_p99_over_median')} -> "
                  f"{stats.get('dynamic_range_p99_over_median')}, "
                  f"FRiP {old.get('frip_proxy')} -> {stats.get('frip_proxy')}")
        else:
            tmp = qpath + ".tmp"
            json.dump(merged, open(tmp, "w"), indent=1)
            os.replace(tmp, qpath)
        n_ok += 1
    print(f"[45] recomputed {n_ok} qc.json{' (dry-run)' if a.dry_run else ''}; {n_skip} skipped")


if __name__ == "__main__":
    main()
