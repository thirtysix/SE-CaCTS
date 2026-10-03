#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Phase 1 — 26: measured copy number from the Childhood Cancer Model Atlas for lines scored on input-inferred CN
(v3.2; phase1/data/cn_gap_sources_2026-10.md, "cheapest route" 2).

Source: `CCMA_cnv_042024.csv` (Sun et al. 2023, Cancer Cell; Mendeley Data doi:10.17632/rnfs539pfw.3, CC BY 4.0):
WGS CNV calls as gene x sample rows that carry their segment (hg19 start/end, sv_type, norm_rd = read-depth ratio
to normal). Only ALTERED segments are listed, so any region not listed is neutral.

Output: a segment file for cnrose `SegmentFileCN` (key chrom start end ratio, hg38), one key per line:
  - altered segments lifted to hg38 by their endpoints (pyliftover, UCSC chain), as the CCLE SNP6 backend does;
    segments whose ends land on different chromosomes, the minus strand or out of order are dropped and counted
  - every gap filled with ratio 1.0 EXPLICITLY: CNTrack.region_cn falls back to the nearest segment within 3 Mb,
    which would spread an altered segment's ratio into the neutral DNA around it
  - all pieces <= 1 Mb, so region_cn's 3 Mb scan bound still finds them; overlapping lifted segments keep the
    first (CCMA segments do not overlap on hg19)
  - chrY left neutral: its norm_rd is a ratio to a near-zero normal (3.6e11 and 1.2e12 for the two DIPG lines), and
    chrY is not scored (v3.1 item 5)

    python3 phase1/scripts/26_ccma_cn.py --ccma $DATAROOT/CCMA/CCMA_cnv_042024.csv \
        --out phase2/data/cn_ccma_wgs.hg38.tsv.gz
"""
import argparse
import os
import sys

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..", "..")))
from secacts_env import DATAROOT                                          # noqa: E402

HG38 = {"chr1": 248956422, "chr2": 242193529, "chr3": 198295559, "chr4": 190214555, "chr5": 181538259,
        "chr6": 170805979, "chr7": 159345973, "chr8": 145138636, "chr9": 138394717, "chr10": 133797422,
        "chr11": 135086622, "chr12": 133275309, "chr13": 114364328, "chr14": 107043718, "chr15": 101991189,
        "chr16": 90338345, "chr17": 83257441, "chr18": 80373285, "chr19": 58617616, "chr20": 64444167,
        "chr21": 46709983, "chr22": 50818468, "chrX": 156040895, "chrY": 57227415}
# CCMA sample -> pull-set key (Cellosaurus: SU-DIPG-IV CVCL_IT39, SU-DIPG-XIII CVCL_IT41)
SAMPLES = {"SU_DIPG_4": "CVCL_IT39", "SU_DIPG_13": "CVCL_IT41"}


def pieces(chrom, s, e, r, piece):
    while s < e:
        t = min(e, s + piece)
        yield chrom, s, t, r
        s = t


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ccma", default=os.path.join(DATAROOT, "CCMA/CCMA_cnv_042024.csv"))
    ap.add_argument("--chain", default=os.path.expanduser("~/.pyliftover/hg19ToHg38.over.chain.gz"))
    ap.add_argument("--piece", type=int, default=1_000_000)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    from pyliftover import LiftOver
    lo = LiftOver(a.chain)

    d = pd.read_csv(a.ccma, usecols=["sample", "chr", "start", "end", "sv_type", "norm_rd"], dtype={"chr": str})
    d = d[d["sample"].isin(SAMPLES)].drop_duplicates(["sample", "chr", "start", "end"])
    rows = []
    for smp, g in d.groupby("sample"):
        key, segs, drop = SAMPLES[smp], [], 0
        for c, s, e, r in zip(g["chr"], g["start"], g["end"], g["norm_rd"]):
            ch = "chr" + c.replace("chr", "")
            if ch == "chrY":
                continue
            p, q = lo.convert_coordinate(ch, int(s) - 1), lo.convert_coordinate(ch, int(e) - 1)
            if not p or not q or p[0][0] != q[0][0] or p[0][2] != "+" or q[0][2] != "+" or q[0][1] <= p[0][1]:
                drop += 1
                continue
            segs.append((p[0][0], int(p[0][1]), int(q[0][1]) + 1, float(r)))
        segs.sort()
        out, n_alt = [], 0
        for ch, L in HG38.items():
            pos = 0
            for c2, s2, e2, r in (x for x in segs if x[0] == ch):
                s2 = max(s2, pos)
                if e2 <= s2:
                    continue                                              # overlap after liftover: keep the first
                out += pieces(ch, pos, s2, 1.0, a.piece)
                out += pieces(ch, s2, e2, r, a.piece)
                n_alt += 1
                pos = e2
            out += pieces(ch, pos, L, 1.0, a.piece)
        rows += [(key, *x) for x in out]
        alt_bp = sum(e - s for _, s, e, r in out if r != 1.0)
        print(f"[ccma] {smp} -> {key}: {len(g)} altered segments, {drop} dropped in liftover, {n_alt} kept "
              f"({alt_bp / 1e6:.0f} Mb altered); {len(out)} pieces", file=sys.stderr)
    pd.DataFrame(rows, columns=["key", "chrom", "start", "end", "ratio"]).to_csv(a.out, sep="\t", index=False)
    print(f"[ccma] wrote {a.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
