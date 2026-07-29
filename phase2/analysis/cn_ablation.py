#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""CN ablation: which specific-SE calls are driven by copy-number amplification rather than biology?

Compares two `score_pilot.py` runs that differ in EXACTLY one step — scoring-time CN correction applied
vs `--no-cn`. Both arms evaluate CN and report it as `cn_mean`, so amplification is visible in either.

Two readouts:

  1. IDENTITY-GENE STABILITY. `hierarchy_summary.genes` carries each curated master-TF's best rank over ALL
     SEs (not just the reported top-N), so the rank shift between arms is free. A real lineage identity SE
     should be roughly rank-stable under correction; one that collapses when CN is divided out was riding
     an amplicon.
  2. AMPLICON-DRIVEN CALLS. SEs in the UNCORRECTED top-N that the corrected arm drops, ranked by cn_mean.
     This is the pilot's BCAT1 / 19q13 / BCAS1 pattern (ROADMAP Phase 4) generalized to the full atlas.

  ~/miniconda3/envs/atac_hdac/bin/python phase2/analysis/cn_ablation.py \
      --corrected phase2/scores/atlas.s3 --uncorrected phase2/scores/atlas.s3.nocn
"""
from __future__ import annotations

import argparse
import os

import numpy as np
import pandas as pd

LEVELS = ["OncotreeLineage", "OncotreePrimaryDisease", "OncotreeSubtype", "line"]


def parse_genes(cell):
    """'SOX17:6;PAX8:90' -> {'SOX17': 6, 'PAX8': 90}"""
    out = {}
    if isinstance(cell, str) and cell:
        for part in cell.split(";"):
            if ":" in part:
                g, r = part.rsplit(":", 1)
                try:
                    out[g] = int(r)
                except ValueError:
                    pass
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--corrected", required=True, help="output prefix of the CN-CORRECTED run")
    ap.add_argument("--uncorrected", required=True, help="output prefix of the --no-cn run")
    ap.add_argument("--out", default=None, help="write the amplicon table here (default <corrected>.cn_ablation.tsv)")
    ap.add_argument("--cn-amp", type=float, default=1.3, help="cn_mean above which a locus counts as amplified")
    a = ap.parse_args()
    out_path = a.out or f"{a.corrected}.cn_ablation.tsv"

    cs = pd.read_csv(f"{a.corrected}.hierarchy_summary.tsv", sep="\t")
    us = pd.read_csv(f"{a.uncorrected}.hierarchy_summary.tsv", sep="\t")
    key = ["level", "group"]
    m = cs[key + ["n_lines", "genes"]].merge(us[key + ["genes"]], on=key, suffixes=("_cor", "_unc"))

    print("=" * 100)
    print("1. IDENTITY-GENE RANK STABILITY UNDER CN CORRECTION  (uncorrected -> corrected)")
    print("=" * 100)
    rows = []
    for _, r in m.iterrows():
        gc, gu = parse_genes(r["genes_cor"]), parse_genes(r["genes_unc"])
        for g in set(gc) | set(gu):
            rc, ru = gc.get(g), gu.get(g)
            if rc is None or ru is None:
                continue
            rows.append(dict(level=r["level"], group=r["group"], n_lines=r["n_lines"], gene=g,
                             rank_unc=ru, rank_cor=rc, log2_shift=float(np.log2(rc / ru))))
    ident = pd.DataFrame(rows)
    if not ident.empty:
        # only the ones that MATTER: a gene that is near the top in at least one arm
        top = ident[(ident.rank_cor <= 50) | (ident.rank_unc <= 50)].copy()
        top = top.sort_values("log2_shift")
        print(f"  {len(ident)} identity-gene x group pairs; {len(top)} reach rank<=50 in some arm\n")
        print(f"  {'level':<22} {'group':<32} {'gene':<8} {'unc':>7} {'cor':>7}  shift")
        for _, r in top.iterrows():
            arrow = "IMPROVED by correction" if r.log2_shift < -0.5 else (
                    "DEMOTED by correction " if r.log2_shift > 0.5 else "stable                ")
            print(f"  {r.level:<22} {str(r.group)[:32]:<32} {r.gene:<8} {r.rank_unc:>7} {r.rank_cor:>7}  "
                  f"{r.log2_shift:+5.2f}  {arrow}")

    print()
    print("=" * 100)
    print(f"2. AMPLICON-DRIVEN CALLS — in the UNCORRECTED top-N, dropped after correction, cn_mean > {a.cn_amp}")
    print("=" * 100)
    amp = []
    for lev in LEVELS:
        fc = f"{a.corrected}.{lev}.top_specific.tsv"
        fu = f"{a.uncorrected}.{lev}.top_specific.tsv"
        if not (os.path.exists(fc) and os.path.exists(fu)):
            continue
        C, U = pd.read_csv(fc, sep="\t"), pd.read_csv(fu, sep="\t")
        if "cn_mean" not in U.columns:
            print(f"  [skip] {lev}: no cn_mean column (re-run the scorer)")
            continue
        kept = set(zip(C["group"], C["se"]))
        U = U.assign(dropped=[(g, s) not in kept for g, s in zip(U["group"], U["se"])])
        hit = U[U["dropped"] & (U["cn_mean"] > a.cn_amp)]
        amp.append(hit)
        n_drop = int(U["dropped"].sum())
        print(f"  {lev:<24} top-N rows={len(U):>5}  dropped by correction={n_drop:>5}  "
              f"of which amplified={len(hit):>4}")
    if amp:
        A = pd.concat(amp, ignore_index=True).sort_values("cn_mean", ascending=False)
        A.to_csv(out_path, sep="\t", index=False)
        print(f"\n  top amplicon-driven calls (wrote {len(A)} rows -> {out_path}):\n")
        print(f"  {'level':<22} {'group':<28} {'se':<12} {'gene':<10} {'cn':>6} {'rank_unc':>9}")
        for _, r in A.head(20).iterrows():
            print(f"  {r.level:<22} {str(r.group)[:28]:<28} {r.se:<12} {str(r.nearest_gene):<10} "
                  f"{r.cn_mean:>6.2f} {int(r['rank']):>9}")

    print()
    print("=" * 100)
    print("3. CN CONTEXT OF THE CORRECTED CALLS  (are the surviving top hits in neutral CN?)")
    print("=" * 100)
    for lev in LEVELS:
        f = f"{a.corrected}.{lev}.top_specific.tsv"
        if not os.path.exists(f):
            continue
        C = pd.read_csv(f, sep="\t")
        if "cn_mean" not in C.columns:
            continue
        print(f"  {lev:<24} cn_mean median={C.cn_mean.median():.3f}  "
              f"amplified(>{a.cn_amp})={100.0 * (C.cn_mean > a.cn_amp).mean():5.1f}%  "
              f"deleted(<0.7)={100.0 * (C.cn_mean < 0.7).mean():5.1f}%")


if __name__ == "__main__":
    main()
