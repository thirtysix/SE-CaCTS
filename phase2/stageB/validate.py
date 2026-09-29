#!/usr/bin/env python3
"""Stage B validation gate (ROADMAP B2): our ChIP-Atlas-v1 emulation against ChIP-Atlas's own v1 files for the
same experiments. No production run without a pass.

Per experiment: QC tuple (reads, mapped %, duplicate %, peaks) vs the experiment list; RPM step (smallest
bigWig value = 1e6 / reads after dedup); bigWig at 1 kb bins (Pearson on log1p) and the atlas grid.20 column
(Spearman, max relative difference) vs ChIP-Atlas's bigWig; bed20 and SE calls (bp Jaccard, recovery) vs the
pull's v1 files. Thresholds: bins r >= 0.99, grid rho >= 0.99, bed20 Jaccard >= 0.95, SE Jaccard >= 0.9.

  python validate.py --validate-dir $S/validate --manifest $S/manifest.validate.tsv --pull-out $W/out \
      --bed20 $W/data/bed20 --grid $W/data/grid.20.bed --qc-table qc_v1.tsv --cnrose $W/cnrose \
      --compare compare_v1_v2.py --out validate.tsv
"""
import argparse
import importlib.util
import json
import os
import subprocess
import sys

import numpy as np
import pandas as pd


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    for k in ("--validate-dir", "--manifest", "--pull-out", "--bed20", "--grid", "--qc-table", "--cnrose",
              "--compare", "--out"):
        ap.add_argument(k, required=True)
    a = ap.parse_args()
    sys.path.insert(0, a.cnrose)
    from cnrose.io import read_bed3, read_f32, scan_bigwig
    from scipy.stats import spearmanr
    spec = importlib.util.spec_from_file_location("cmp", a.compare)
    cmp = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(cmp)

    man = pd.read_csv(a.manifest, sep="\t")
    qc = pd.read_csv(a.qc_table, sep="\t").set_index("srx")
    grid = read_bed3(a.grid)
    rows = []
    for srx, layout in zip(man.srx, man.layout):
        d = os.path.join(a.validate_dir, srx)
        if not os.path.exists(os.path.join(d, f"{srx}.se.bed")):
            rows.append({"srx": srx, "status": "not run"}); continue
        v1bw = os.path.join(a.validate_dir, "v1", f"{srx}.bw")
        os.makedirs(os.path.dirname(v1bw), exist_ok=True)
        if not os.path.exists(v1bw):
            subprocess.check_call(["curl", "-sfL", "-o", v1bw, f"https://chip-atlas.dbcls.jp/data/hg38/eachData/bw/{srx}.bw"])
        r = {"srx": srx, "layout": layout}
        mine = json.load(open(os.path.join(d, f"{srx}.stageB_qc.json")))
        for k in ("reads", "pct_mapped", "pct_dup", "n_peaks"):
            r[f"{k}_ours"], r[f"{k}_v1"] = mine[k], qc.loc[srx, k] if srx in qc.index else np.nan
        (g1,), b1, _ = scan_bigwig(v1bw, [grid], binsize=1000)
        (g2,), b2, _ = scan_bigwig(os.path.join(d, f"{srx}.bw"), [grid], binsize=1000)
        x, y = np.log1p(b1.astype(float)), np.log1p(b2.astype(float))
        nz = (b1 > 0) | (b2 > 0)
        r["bins_r"] = float(np.corrcoef(x[nz], y[nz])[0, 1])
        r["scale_ratio"] = float(b2.sum() / b1.sum())
        anyg = (g1 > 0) | (g2 > 0)
        r["grid_rho"] = float(spearmanr(g1[anyg], g2[anyg]).correlation)
        r["grid_max_rel"] = float(np.max(np.abs(g2 - g1)[anyg] / np.maximum(g1[anyg], 1e-6)))
        r["grid_median_rel"] = float(np.median(np.abs(g2 - g1)[anyg] / np.maximum(g1[anyg], 1e-6)))
        p1 = read_bed3(os.path.join(a.bed20, f"{srx}.20.bed"), canonical_only=True)
        p2 = read_bed3(os.path.join(d, f"{srx}.20.bed"), canonical_only=True)
        for k, v in cmp.interval_metrics(p1, p2).items():
            r[f"peaks_{k}"] = v
        s1 = read_bed3(os.path.join(a.pull_out, srx[-2:], f"{srx}.se.bed"), canonical_only=True)
        s2 = read_bed3(os.path.join(d, f"{srx}.se.bed"), canonical_only=True)
        for k, v in cmp.interval_metrics(s1, s2).items():
            r[f"se_{k}"] = v
        for k in ("timing.align", "timing.post"):
            fn = os.path.join(d, k)
            if os.path.exists(fn):
                for kv in open(fn).read().split():
                    kk, vv = kv.split("=")
                    r[kk] = int(vv)
        r["passes"] = bool(r["bins_r"] >= 0.99 and r["grid_rho"] >= 0.99 and r["peaks_bp_jaccard"] >= 0.95
                           and r["se_bp_jaccard"] >= 0.9)
        rows.append(r)
    t = pd.DataFrame(rows)
    t.to_csv(a.out, sep="\t", index=False)
    pd.set_option("display.width", 250)
    cols = ["srx", "layout", "reads_ours", "reads_v1", "pct_mapped_ours", "pct_mapped_v1", "pct_dup_ours", "pct_dup_v1",
            "n_peaks_ours", "n_peaks_v1", "bins_r", "scale_ratio", "grid_rho", "grid_median_rel", "peaks_bp_jaccard",
            "se_bp_jaccard", "align_s", "post_s", "passes"]
    print(t[[c for c in cols if c in t]].round(4).to_string(index=False))
    print(f"\n[validate] {int(t.get('passes', pd.Series(dtype=bool)).fillna(False).sum())}/{len(t)} pass")


if __name__ == "__main__":
    main()
