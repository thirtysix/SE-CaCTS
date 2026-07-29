#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""CALL-BASED CN ablation, on the HONEST (permutation) null.

The original `cn_ablation.py` counts amplicon-driven SEs by RANK — "in the uncorrected top-15, dropped
after correction". That is rank-based and null-invariant, but it is not a statement about CALLS. With the
permutation FDR available for BOTH arms we can now ask the cleaner question:

    which SE x group tests PASS the permutation FDR without CN correction but FAIL it with correction?

Those are **amplicon-driven calls** — significant only because amplification inflated the signal. The
converse set (fail uncorrected, pass corrected) are **rescued calls** — real specificity that an
amplicon-inflated context was masking. Both arms report `cn_mean` (CN is always evaluated; `--no-cn` only
skips APPLYING it), so the copy number at each disputed locus is read directly.

Restricted to the levels the panel supports — OncotreeLineage and OncotreePrimaryDisease. Subtype has ~1
call and cell-line permutation is degenerate (gotcha 72), so a call-based ablation is undefined there; the
per-line amplicon story (MYCN in SK-N-BE(2), etc.) remains the RANK-based `cn_ablation.py`.

  ~/miniconda3/envs/atac_hdac/bin/python phase2/analysis/cn_ablation_calls.py \
      --corrected phase2/scores/atlas.s3.perm --uncorrected phase2/scores/atlas.s3.perm.nocn
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np
import pandas as pd

_ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
sys.path.insert(0, _ROOT)
from secacts_env import DATAROOT, cache_path                      # noqa: E402  (local paths live in .env)
HERE = os.path.dirname(os.path.abspath(__file__))
SECACTS = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(SECACTS, "cnrose"))

from cnrose.cn.depmap import load_gene_coords, DepMapGeneCN     # noqa: E402

LEVELS = ["OncotreeLineage", "OncotreePrimaryDisease"]
GENE_CACHE = cache_path("gene_coords.GRCh38.106.tsv")


def nearest_gene_fn():
    """Same universe as score_pilot's nearest_gene column: protein-coding genes present in DepMap."""
    gc = load_gene_coords(f"{DATAROOT}/0.human_genome/Homo_sapiens.GRCh38.106.chr.gtf.gz", cache_path=GENE_CACHE)
    prov = DepMapGeneCN(f"{DATAROOT}/DepMap/2026q1/OmicsCNGeneWGS.csv", gc)
    gidx = {}
    for g in prov.usable:
        gidx.setdefault(gc[g][0], []).append((g, (gc[g][1] + gc[g][2]) // 2))
    for c in gidx:
        gidx[c].sort(key=lambda t: t[1])
    mids = {c: np.array([m for _, m in v]) for c, v in gidx.items()}

    def nearest(chrom, mid):
        v = gidx.get(chrom)
        if not v:
            return "", 0
        i = int(np.searchsorted(mids[chrom], mid))
        best, bd = "", None
        for j in (i - 1, i, i + 1):
            if 0 <= j < len(v):
                d = abs(v[j][1] - mid)
                if bd is None or d < bd:
                    best, bd = v[j][0], d
        return best, (bd or 0) // 1000
    return nearest


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--corrected", default=f"{SECACTS}/phase2/scores/atlas.s3.perm")
    ap.add_argument("--uncorrected", default=f"{SECACTS}/phase2/scores/atlas.s3.perm.nocn")
    ap.add_argument("--catalog", default=f"{SECACTS}/phase2/results/atlas.s3.union_catalog.bed.gz")
    ap.add_argument("--cn-amp", type=float, default=1.3)
    ap.add_argument("--out", default=f"{SECACTS}/phase2/scores/atlas.s3.perm.cn_ablation_calls.tsv")
    a = ap.parse_args()

    import gzip
    coords = {}
    with gzip.open(a.catalog, "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            coords[f[3]] = (f[0], (int(f[1]) + int(f[2])) // 2)
    nearest = nearest_gene_fn()

    rows = []
    print(f"{'level':<24} {'uncorr':>7} {'corr':>7} {'ampl-driven':>12} {'(cn>' + str(a.cn_amp) + ')':>9} "
          f"{'rescued':>8} {'stable':>7}")
    for lev in LEVELS:
        fu = f"{a.uncorrected}.{lev}.specific.tsv.gz"
        fc = f"{a.corrected}.{lev}.specific.tsv.gz"
        if not (os.path.exists(fu) and os.path.exists(fc)):
            print(f"  [skip] {lev}: missing dump(s)")
            continue
        U = pd.read_csv(fu, sep="\t")
        C = pd.read_csv(fc, sep="\t")
        uset = set(zip(U["group"], U["se"]))
        cset = set(zip(C["group"], C["se"]))
        ampl = uset - cset                       # pass uncorrected, fail corrected
        resc = cset - uset                       # fail uncorrected, pass corrected
        stable = uset & cset
        Ui = U.set_index(["group", "se"])
        ampl_cn = Ui.loc[list(ampl), "cn_mean"] if ampl else pd.Series(dtype=float)
        n_amp_hi = int((ampl_cn > a.cn_amp).sum())
        print(f"{lev:<24} {len(uset):>7,} {len(cset):>7,} {len(ampl):>12,} {n_amp_hi:>9,} "
              f"{len(resc):>8,} {len(stable):>7,}")
        for (grp, se) in ampl:
            cn = float(Ui.loc[(grp, se), "cn_mean"])
            g, d = nearest(*coords[se]) if se in coords else ("", 0)
            rows.append(dict(level=lev, kind="amplicon_driven", group=grp, se=se, cn_mean=cn,
                             nearest_gene=g, dist_kb=d))
        for (grp, se) in resc:
            cn = float(C.set_index(["group", "se"]).loc[(grp, se), "cn_mean"])
            g, d = nearest(*coords[se]) if se in coords else ("", 0)
            rows.append(dict(level=lev, kind="rescued", group=grp, se=se, cn_mean=cn, nearest_gene=g, dist_kb=d))

    R = pd.DataFrame(rows)
    if not R.empty:
        R.to_csv(a.out, sep="\t", index=False)
        print("\n=== cn_mean by outcome (amplicon-driven should be HIGH; rescued/stable near-neutral) ===")
        for k in ["amplicon_driven", "rescued"]:
            s = R[R.kind == k]["cn_mean"]
            if len(s):
                print(f"  {k:<16} n={len(s):>4}  median cn={s.median():.3f}  "
                      f"share cn>{a.cn_amp}={100 * (s > a.cn_amp).mean():.0f}%")
        amp = R[(R.kind == "amplicon_driven")].sort_values("cn_mean", ascending=False)
        if len(amp):
            print(f"\n=== top amplicon-driven CALLS (pass uncorrected FDR, fail corrected) ===")
            print(f"  {'level':<22} {'group':<28} {'gene':<10} {'cn':>6}")
            for _, r in amp.head(15).iterrows():
                print(f"  {r.level:<22} {str(r.group)[:28]:<28} {str(r.nearest_gene):<10} {r.cn_mean:>6.2f}")
        print(f"\n[cn_ablation_calls] wrote {a.out} ({len(R)} rows)")
    else:
        print("\n[cn_ablation_calls] no disputed calls at the supported levels")


if __name__ == "__main__":
    main()
