#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Phase 1 — 21: validate CCLE 2019 SNP6 CN (cnrose.cn.ccle) against DepMap WGS before using it (FINDINGS §21).

Same test the CMP and MC_WES backends passed (CN_COVERAGE.md §5): every line with both sources, log2 CN at
5 Mb bins over chr1-22, Pearson r per line; then which held-back lines it would reach (own arrays, or a
Cellosaurus relative's, from phase2/analysis/out/cn_routes.tsv).

  python3 21_ccle_snp6_validate.py --seg data_cna_hg19.seg --clinical data_clinical_sample.txt \
      --depmap-wgs OmicsCNGeneWGS.csv --gene-coords gene_coords.GRCh38.106.tsv --cnrose <cnrose root> \
      --out ccle_snp6_vs_wgs.tsv
"""
import argparse
import sys

import numpy as np

HG38 = {"chr1": 248956422, "chr2": 242193529, "chr3": 198295559, "chr4": 190214555, "chr5": 181538259,
        "chr6": 170805979, "chr7": 159345973, "chr8": 145138636, "chr9": 138394717, "chr10": 133797422,
        "chr11": 135086622, "chr12": 133275309, "chr13": 114364328, "chr14": 107043718, "chr15": 101991189,
        "chr16": 90338345, "chr17": 83257441, "chr18": 80373285, "chr19": 58617616, "chr20": 64444167,
        "chr21": 46709983, "chr22": 50818468}
COARSE = 5_000_000


def coarse_log2(track):
    out = []
    for c, L in HG38.items():
        for s in range(0, L - COARSE, COARSE):
            v = track.region_cn(c, s, s + COARSE)
            out.append(np.log2(v) if v and v > 0 else np.nan)
    return np.array(out)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    for k in ("--seg", "--clinical", "--depmap-wgs", "--gene-coords", "--cnrose", "--out"):
        ap.add_argument(k, required=True)
    ap.add_argument("--chain", help="hg19ToHg38.over.chain.gz (UCSC)")
    ap.add_argument("--export-keys", help="comma-separated ModelIDs: write their hg38 segments to --export and stop")
    ap.add_argument("--export", help="TSV(.gz) of key, chrom, start, end, ratio for cnrose.cn.segfile")
    a = ap.parse_args()
    sys.path.insert(0, a.cnrose)
    import pandas as pd
    from cnrose.cn.ccle import Ccle2019SnpCN
    from cnrose.cn.depmap import DepMapGeneCN, load_gene_coords

    ccle = Ccle2019SnpCN(a.seg, a.clinical, chain=a.chain)
    if a.export_keys:
        rows = []
        for k in a.export_keys.split(","):
            t = ccle.track(k)
            if t is None:
                print(f"[21] {k}: no CCLE 2019 segments", file=sys.stderr)
                continue
            for c, (st, en, ra) in t._chrom.items():
                rows += [(k, c, int(s), int(e), float(r)) for s, e, r in zip(st, en, ra)]
        pd.DataFrame(rows, columns=["key", "chrom", "start", "end", "ratio"]).to_csv(a.export, sep="\t", index=False)
        print(f"[21] exported {len({r[0] for r in rows})} lines, {len(rows)} segments -> {a.export}", file=sys.stderr)
        return
    wgs = DepMapGeneCN(a.depmap_wgs, load_gene_coords(None, cache_path=a.gene_coords))
    both = sorted(set(ccle.sample_for) & set(wgs.index) if hasattr(wgs, "index") else set(ccle.sample_for))
    wgs.preload(both)
    rows = []
    for mid in both:
        t2 = wgs.track(mid)
        if t2 is None:
            continue
        t1 = ccle.track(mid)
        if t1 is None:
            continue
        x, y = coarse_log2(t1), coarse_log2(t2)
        ok = np.isfinite(x) & np.isfinite(y)
        if ok.sum() > 50:
            rows.append((mid, float(np.corrcoef(x[ok], y[ok])[0, 1]), int(ok.sum()), ccle.dropped.get(mid, 0)))
    d = pd.DataFrame(rows, columns=["model_id", "r", "n_bins", "segments_dropped_in_liftover"])
    d.to_csv(a.out, sep="\t", index=False)
    print(f"[21] CCLE 2019 SNP6 vs DepMap WGS: {len(d)} shared lines; median r {d.r.median():.3f}, "
          f"mean {d.r.mean():.3f}, r < 0.5: {int((d.r < 0.5).sum())}; median segments dropped in liftover "
          f"{d.segments_dropped_in_liftover.median():.0f}", file=sys.stderr)


if __name__ == "__main__":
    main()
