#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Is correction strength a property of the CN SOURCE? A paired test on lines that have both.

v2 corrects 282 lines with DepMap WGS copy number and 104 with WES (CMP pureCN or DepMap MC_WES). Per line,
the Spearman of SE signal against SE copy number drops after correction; if a source is noisier, dividing
by it over-corrects and pushes that rho NEGATIVE. Comparing sources across DIFFERENT lines confounds the
source with the lines (and their lineages). Here every line is scored against BOTH sources, so the
difference is the source alone:

  per line  rho_raw / rho_corrected with WGS CN, and with CMP WES CN; plus r(WGS CN, CMP CN) over SEs

  ~/miniconda3/envs/atac_hdac/bin/python phase2/analysis/cn_source_paired.py
"""
from __future__ import annotations

import gzip
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, pearsonr, wilcoxon

HERE = os.path.dirname(os.path.abspath(__file__))
SECACTS = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, SECACTS)
sys.path.insert(0, os.path.join(SECACTS, "cnrose"))
from secacts_env import DATAROOT, cache_path                                  # noqa: E402
from cnrose.cn.depmap import load_gene_coords, DepMapGeneCN                    # noqa: E402
from cnrose.cn.cmp import CellModelPassportsWesCN                              # noqa: E402
from cnrose.cn.base import correct                                             # noqa: E402

RES = os.environ.get("SECACTS_RES", os.path.join(SECACTS, "phase2", "results_v2"))
OUT = os.environ.get("PAIRED_OUT", os.path.join(SECACTS, "phase2", "scores_v2", "cn_source_paired.tsv"))


def main():
    ps = pd.read_csv(os.path.join(SECACTS, "phase2/data/pull_set.v2.tsv"), sep="\t")
    wgs_lines = ps[ps.cn_provider == "depmap_wgs"].drop_duplicates("key").set_index("key")
    S = pd.read_csv(os.path.join(RES, "atlas.s3.se_signal.tsv.gz"), sep="\t", index_col=0)
    srx_key = dict(zip(ps.srx, ps.key))
    keys = pd.Series([srx_key.get(c) for c in S.columns], index=S.columns)
    keys = keys[keys.isin(wgs_lines.index)]
    L = S[keys.index].T.groupby(keys).mean().T                                     # SE x line (WGS lines)
    del S
    coords = {}
    with gzip.open(os.path.join(RES, "atlas.s3.union_catalog.bed.gz"), "rt") as fh:
        for line in fh:
            f = line.split("\t")
            coords[f[3]] = (f[0], int(f[1]), int(f[2]))
    se = list(L.index)

    gc = load_gene_coords(os.path.join(DATAROOT, "0.human_genome/Homo_sapiens.GRCh38.106.chr.gtf.gz"),
                          cache_path=cache_path("gene_coords.GRCh38.106.tsv"))
    wgs = DepMapGeneCN(os.path.join(DATAROOT, "DepMap/2026q1/OmicsCNGeneWGS.csv"), gc)
    wgs.preload(list(L.columns))
    cmp = CellModelPassportsWesCN(os.path.join(DATAROOT, "CellModelPassports/WES_pureCN_CNV_genes_latest.csv.gz"),
                                  os.path.join(DATAROOT, "CellModelPassports/model_list_20240110.csv"),
                                  cache_dir=cache_path("cmp_wes"))
    cvcl = wgs_lines["cvcl"].to_dict()
    # the per-model cache holds every line CMP WES covers (a full-file scan wrote it); a line with no cache
    # file has no CMP WES copy number, so skip it rather than re-scanning the ~1 GB source for nothing
    both = [k for k in L.columns if cmp.resolve(cvcl[k]) and os.path.exists(cmp._cache_file(cmp.resolve(cvcl[k])))]
    L = L[both]
    cmp.preload([cvcl[k] for k in both])

    rows = []
    for k in L.columns:
        tw, tc = wgs.track(k), cmp.track(cvcl[k])
        if tw is None or tc is None:
            continue
        cw = np.fromiter((tw.region_cn(*coords[s]) for s in se), float, len(se))
        cc = np.fromiter((tc.region_cn(*coords[s]) for s in se), float, len(se))
        x = L[k].values.astype(float)
        rec = dict(key=k, line=wgs_lines.loc[k, "cell"], lineage=wgs_lines.loc[k, "lineage"],
                   r_cn=round(float(pearsonr(np.log2(cw), np.log2(cc))[0]), 4))
        for tag, c in (("wgs", cw), ("cmp", cc)):
            y = np.maximum(correct(x, c, model="log2offset", floor=0.1), 0.0)
            rec[f"rho_raw_{tag}"] = round(float(spearmanr(x, c)[0]), 4)
            rec[f"rho_cor_{tag}"] = round(float(spearmanr(y, c)[0]), 4)
            rec[f"sd_log2cn_{tag}"] = round(float(np.std(np.log2(c))), 4)
        rows.append(rec)
        if len(rows) % 50 == 0:
            print(f"[paired] {len(rows)} lines", file=sys.stderr, flush=True)
    D = pd.DataFrame(rows)
    D.to_csv(OUT, sep="\t", index=False)

    print(f"\n[paired] {len(D)} lines carry both DepMap WGS and CMP WES copy number")
    print(f"  r(log2 CN WGS, log2 CN CMP) over SEs: median {D.r_cn.median():.3f}")
    for m in ("rho_raw", "rho_cor", "sd_log2cn"):
        a, b = D[f"{m}_wgs"], D[f"{m}_cmp"]
        p = wilcoxon(a, b).pvalue
        print(f"  {m:<10} WGS median {a.median():+.3f} | CMP median {b.median():+.3f} | "
              f"paired diff (CMP-WGS) median {np.median(b - a):+.3f}  Wilcoxon p={p:.2g}")
    print(f"[paired] wrote {OUT}")


if __name__ == "__main__":
    main()
