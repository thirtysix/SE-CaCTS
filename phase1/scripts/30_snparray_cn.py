#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Phase 1 — 30: measured copy number from Affymetrix SNP-array CEL files on GEO (SNP6, CytoScan HD), for lines
scored on input-inferred CN (v3.2; FINDINGS §61).

rawcopy (rawcopy.org; conda env `rawcopy`, see phase1/data/README or FINDINGS §61 for the install) normalises each
CEL against its built-in reference of normals and segments it (hg19). Segments -> hg38 by their endpoints
(pyliftover) -> one value per 50 kb bin (the segment over the bin centre) -> ratio 2**log2, written next to the
WGS ratios with an index.tsv row (source `array_snp6` / `array_cytoscan`), so 27_wgs_cn.py evaluate calibrates
each platform on lines of the SAME GEO series that have DepMap WGS, and the per-source gate decides.

    python3 phase1/scripts/30_snparray_cn.py fetch
    python3 phase1/scripts/30_snparray_cn.py run          # rawcopy, one CEL at a time (laptop: one worker)
    python3 phase1/scripts/30_snparray_cn.py bins --out phase2/analysis/out/wgs_cn/ratio
"""
import argparse
import os
import subprocess
import sys
import urllib.request

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
CACHE = os.path.join(ROOT, ".cache", "geo_cel")
RAW = os.path.join(ROOT, "phase2", "analysis", "out", "rawcopy")
RENV = os.path.expanduser("~/miniconda3/envs/rawcopy")
GEO = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={}&targ=self&form=text&view=brief"
HG38_LEN = {"chr1": 248956422, "chr2": 242193529, "chr3": 198295559, "chr4": 190214555, "chr5": 181538259,
            "chr6": 170805979, "chr7": 159345973, "chr8": 145138636, "chr9": 138394717, "chr10": 133797422,
            "chr11": 135086622, "chr12": 133275309, "chr13": 114364328, "chr14": 107043718, "chr15": 101991189,
            "chr16": 90338345, "chr17": 83257441, "chr18": 80373285, "chr19": 58617616, "chr20": 64444167,
            "chr21": 46709983, "chr22": 50818468, "chrX": 156040895}
# (name, key, role, GSM, platform); calibration lines come from the targets' own series and have DepMap WGS
SAMPLES = [
    # SNP6, GSE22208 (lymphoma lines)
    ("OCI-LY4", "CVCL_8801", "target", "GSM552452", "snp6"), ("OCI-LY8", "CVCL_8803", "target", "GSM552454", "snp6"),
    ("SU-DHL-7", "CVCL_4380", "target", "GSM552466", "snp6"),
    ("DB", "ACH-000334", "calib", "GSM552401", "snp6"), ("HT_snp6", "ACH-000914", "calib", "GSM552431", "snp6"),
    ("KARPAS422_snp6", "ACH-000315", "calib", "GSM552444", "snp6"), ("SUDHL4", "ACH-000365", "calib", "GSM552463", "snp6"),
    ("WSUNHL", "ACH-001709", "calib", "GSM552470", "snp6"),
    # SNP6, GSE22305 (melanoma), GSE34211, GSE50253 (ALCL)
    ("SK-MEL-147", "CVCL_3876", "target", "GSM555181", "snp6"),
    ("A375", "ACH-000219", "calib", "GSM555173", "snp6"), ("SKMEL19", "ACH-002005", "calib", "GSM555175", "snp6"),
    ("NCI-H2107", "CVCL_1527", "target", "GSM845544", "snp6"),
    ("JB6", "CVCL_H633", "target", "GSM1217150", "snp6"), ("Mac-1", "CVCL_H631", "target", "GSM1217151", "snp6"),
    # CytoScan HD, GSE209728 / GSE101533 (neuroblastoma)
    ("CLB-Pe", "CVCL_9534", "target", "GSM6387482", "cytoscan"), ("NBL-S_cyto", "ACH-001361", "target", "GSM2705773", "cytoscan"),
    ("SKNAS_cyto", "ACH-000260", "calib", "GSM6387470", "cytoscan"), ("SKNBE2", "ACH-000312", "calib", "GSM6387471", "cytoscan"),
    ("KELLY", "ACH-000259", "calib", "GSM6387476", "cytoscan"), ("IMR32_cyto", "ACH-000310", "calib", "GSM6387477", "cytoscan"),
    ("NGP_cyto", "ACH-001366", "calib", "GSM6387479", "cytoscan"), ("SHSY5Y_cyto", "ACH-001188", "calib", "GSM6387475", "cytoscan"),
]


def cmd_fetch(a):
    os.makedirs(CACHE, exist_ok=True)
    for name, key, role, gsm, plat in SAMPLES:
        fn = os.path.join(CACHE, f"{gsm}.CEL.gz")
        if os.path.exists(fn):
            continue
        with urllib.request.urlopen(GEO.format(gsm), timeout=120) as r:
            urls = [l.split(" = ", 1)[1] for l in r.read().decode().splitlines()
                    if l.startswith("!Sample_supplementary_file") and l.upper().endswith(".CEL.GZ")]
        if not urls:
            print(f"[fetch] {name} {gsm}: no CEL", file=sys.stderr)
            continue
        urllib.request.urlretrieve(urls[0].replace("ftp://", "https://"), fn + ".tmp")
        os.replace(fn + ".tmp", fn)
        print(f"[fetch] {name} {gsm}: {os.path.getsize(fn) / 1e6:.0f} MB", file=sys.stderr)


def cmd_run(a):
    """rawcopy per CEL (its own folder), single core; skips samples already done."""
    os.makedirs(RAW, exist_ok=True)
    env = dict(os.environ, PATH=f"{RENV}/bin:{os.environ['PATH']}", OMP_NUM_THREADS="1")
    for name, key, role, gsm, plat in SAMPLES:
        if a.only and name not in a.only.split(","):
            continue
        out = os.path.join(RAW, gsm)
        if os.path.exists(os.path.join(out, "segments.txt")) or any(
                f == "segments.txt" for _, _, fs in os.walk(out) for f in fs):
            continue
        src = os.path.join(CACHE, f"{gsm}.CEL.gz")
        if not os.path.exists(src):
            continue
        d = os.path.join(RAW, "_in", gsm)
        os.makedirs(d, exist_ok=True)
        cel = os.path.join(d, f"{gsm}.CEL")
        if not os.path.exists(cel):
            subprocess.run(f"gunzip -c {src} > {cel}", shell=True, check=True)
        r = (f'suppressMessages(library(rawcopy)); rawcopy(CELfiles.or.directory="{cel}", outdir="{RAW}", cores=1)')
        p = subprocess.run([f"{RENV}/bin/Rscript", "-e", r], env=env, capture_output=True, text=True)
        print(f"[run] {name} {gsm}: exit {p.returncode}", file=sys.stderr)
        if p.returncode:
            print(p.stderr[-1500:], file=sys.stderr)
        os.remove(cel)


def find_segments(gsm):
    for dp, _, fs in os.walk(RAW):
        if "segments.txt" in fs and gsm in dp:
            return os.path.join(dp, "segments.txt")
    return None


def cmd_bins(a):
    from pyliftover import LiftOver
    lo = LiftOver(os.path.expanduser("~/.pyliftover/hg19ToHg38.over.chain.gz"))
    rows = []
    for name, key, role, gsm, plat in SAMPLES:
        fn = find_segments(gsm)
        if fn is None:
            print(f"[bins] {name} {gsm}: no rawcopy segments", file=sys.stderr)
            continue
        s = pd.read_csv(fn, sep="\t")
        ccol = [c for c in s.columns if c.lower().startswith("chrom")][0]
        vcol = [c for c in s.columns if c.lower() in ("value", "log2ratio", "logratio")][0]
        out = {c: np.full(L // 50_000, np.nan) for c, L in HG38_LEN.items()}
        drop = 0
        for c, st, en, v in zip(s[ccol], s["Start"], s["End"], s[vcol]):
            ch = "chr" + str(c).replace("chr", "").replace("23", "X")
            if ch not in out or not np.isfinite(v):
                continue
            p, q = lo.convert_coordinate(ch, int(st) - 1), lo.convert_coordinate(ch, int(en) - 1)
            if not p or not q or p[0][0] != ch or q[0][0] != ch or q[0][1] <= p[0][1]:
                drop += 1
                continue
            i0, i1 = int(p[0][1] + 25_000) // 50_000, int(q[0][1] - 25_000) // 50_000   # bins whose centre is inside
            out[ch][max(i0, 0): min(i1, len(out[ch]) - 1) + 1] = 2.0 ** v
        np.savez_compressed(os.path.join(a.out, f"{name}.npz"), **{c: v.astype(np.float32) for c, v in out.items()},
                            _reads=np.array([len(s)]))
        rows.append({"name": name, "key": key, "role": role, "reads": len(s), "source": f"array_{plat}"})
        print(f"[bins] {name}: {len(s)} segments, {drop} dropped in liftover", file=sys.stderr)
    idx_f = os.path.join(a.out, "index.tsv")
    idx = pd.read_csv(idx_f, sep="\t")
    idx = idx[~idx.name.isin([r["name"] for r in rows])]
    pd.concat([idx, pd.DataFrame(rows)], ignore_index=True).to_csv(idx_f, sep="\t", index=False)


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("fetch")
    r = sub.add_parser("run")
    r.add_argument("--only")
    b = sub.add_parser("bins")
    b.add_argument("--out", required=True)
    a = ap.parse_args()
    {"fetch": cmd_fetch, "run": cmd_run, "bins": cmd_bins}[a.cmd](a)


if __name__ == "__main__":
    main()
