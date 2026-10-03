#!/usr/bin/env python3
"""Count aligned reads per 50 kb bin from a SAM stream, for measured copy number from low-pass WGS (v3.2;
phase1/scripts/27_wgs_cn.py, phase2/roihu/wgs_cn.slurm).

Reads with MAPQ >= --min-mapq that are not unmapped, secondary or supplementary are counted at their leftmost
position. Output: one float32 array per chromosome (chr1-22, chrX), length chrom_len // bin as in the ChIP-input
bins (phase1/scripts/20_cn_input_calibration.py), plus `_reads` (reads counted) and `_seen` (alignment records).

    bowtie2 ... | python3 phase2/wgs_bin.py --out <key>.npz
"""
import argparse
import sys

import numpy as np

HG38 = {"chr1": 248956422, "chr2": 242193529, "chr3": 198295559, "chr4": 190214555, "chr5": 181538259,
        "chr6": 170805979, "chr7": 159345973, "chr8": 145138636, "chr9": 138394717, "chr10": 133797422,
        "chr11": 135086622, "chr12": 133275309, "chr13": 114364328, "chr14": 107043718, "chr15": 101991189,
        "chr16": 90338345, "chr17": 83257441, "chr18": 80373285, "chr19": 58617616, "chr20": 64444167,
        "chr21": 46709983, "chr22": 50818468, "chrX": 156040895}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--bin", type=int, default=50_000)
    ap.add_argument("--min-mapq", type=int, default=10)
    a = ap.parse_args()
    cnt = {c: np.zeros(L // a.bin, dtype=np.float64) for c, L in HG38.items()}
    seen = used = 0
    for line in sys.stdin:
        if line[0] == "@":
            continue
        f = line.split("\t", 5)
        seen += 1
        flag = int(f[1])
        if flag & 0x904 or int(f[4]) < a.min_mapq:          # unmapped, secondary, supplementary; low MAPQ
            continue
        c = cnt.get(f[2])
        if c is None:
            continue
        i = (int(f[3]) - 1) // a.bin
        if i < len(c):
            c[i] += 1
            used += 1
    np.savez_compressed(a.out, **{c: v.astype(np.float32) for c, v in cnt.items()},
                        _reads=np.array([used]), _seen=np.array([seen]))
    print(f"[wgs_bin] {seen} records, {used} reads counted -> {a.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
