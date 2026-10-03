#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Phase 1 — 29: measured copy number from array CGH on GEO, for lines scored on input-inferred CN (v3.2;
phase1/data/cn_gap_sources_2026-10.md "GEO arrays"; FINDINGS §59).

Array CGH sample tables on GEO already carry the processed log2 ratio of the line against normal reference DNA
(VALUE), so no raw-file processing is needed. Per sample: probe log2 ratios -> hg38 (pyliftover, the platform's
build) -> mean per 50 kb bin -> ratio 2**mean, written as `<name>.npz` next to the WGS ratios with an index.tsv row
(source `array_cgh`), so `BinnedInputCN(segment=True, slope=1, zero_is_deletion=False)` reads them like the WGS
tracks and `27_wgs_cn.py evaluate` scores them. Calibration rows are lines on the same array series that have
DepMap WGS (independent data).

Platforms:
  nimblegen_hg18  NimbleGen HG18 whole-genome tiling (e.g. GPL15436, GSE43272): probe ids CHR01FS000032108 encode
                  chromosome and hg18 position
  agilent_<build> Agilent CGH Feature Extraction file (GEO supplementary; the GSM table is empty, e.g. GSE22694):
                  SystematicName chr:start-end, LogRatio = log10(red / green); `flip` when the line is green (Cy3)

    python3 phase1/scripts/29_array_cn.py --out phase2/analysis/out/wgs_cn/ratio
"""
import argparse
import io
import os
import re
import sys
import urllib.request

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, ROOT)
from secacts_env import DATAROOT                                          # noqa: E402

HG38_LEN = {"chr1": 248956422, "chr2": 242193529, "chr3": 198295559, "chr4": 190214555, "chr5": 181538259,
            "chr6": 170805979, "chr7": 159345973, "chr8": 145138636, "chr9": 138394717, "chr10": 133797422,
            "chr11": 135086622, "chr12": 133275309, "chr13": 114364328, "chr14": 107043718, "chr15": 101991189,
            "chr16": 90338345, "chr17": 83257441, "chr18": 80373285, "chr19": 58617616, "chr20": 64444167,
            "chr21": 46709983, "chr22": 50818468, "chrX": 156040895}
CHAINS = {"hg17": os.path.join(DATAROOT, "0.human_genome/hg17ToHg38.over.chain.gz"),
          "hg18": os.path.join(DATAROOT, "0.human_genome/hg18ToHg38.over.chain"),
          "hg19": os.path.expanduser("~/.pyliftover/hg19ToHg38.over.chain.gz")}
# (name as scored, key, role, GSM, platform). GSE43272: DLBCL lines vs male control DNA, NimbleGen 3x720K HG18.
SAMPLES = [("TMD8", "CVCL_A442", "target", "GSM1059801", "nimblegen_hg18"),
           ("U-2932", "CVCL_1896", "target", "GSM1059802", "nimblegen_hg18"),
           ("HBL-1", "CVCL_M572", "target", "GSM1059804", "nimblegen_hg18"),
           ("OCI-LY1_acgh", "CVCL_1879", "check", "GSM1059803", "nimblegen_hg18"),   # cross-check of the cSRA track
           ("BJAB", "ACH-001447", "calib", "GSM1059805", "nimblegen_hg18"),
           ("HT", "ACH-000914", "calib", "GSM1059806", "nimblegen_hg18"),
           ("KARPAS422", "ACH-000315", "calib", "GSM1059807", "nimblegen_hg18"),
           # GSE22694: gastric lines (Cy3) vs spleen DNA (Cy5), Agilent 244K alpha (GPL4544); build chosen on the
           # calibration lines (--agilent-build)
           ("ISt-1", "CVCL_F645", "target", "GSM562416", "agilent_flip"),
           ("AGS", "ACH-000880", "calib", "GSM562392", "agilent_flip"),
           ("KATOIII", "ACH-000793", "calib", "GSM562404", "agilent_flip"),
           ("NCIN87", "ACH-000427", "calib", "GSM562405", "agilent_flip"),
           ("SNU16", "ACH-000581", "calib", "GSM562407", "agilent_flip")]
GEO = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={}&targ=self&form=text&view=data"


def fetch(gsm, cache):
    fn = os.path.join(cache, f"{gsm}.tsv.gz")
    if not os.path.exists(fn):
        with urllib.request.urlopen(GEO.format(gsm), timeout=300) as r:
            txt = r.read().decode()
        rows = "\n".join(l for l in txt.splitlines() if l and l[0] not in "!#^")
        pd.read_csv(io.StringIO(rows), sep="\t").to_csv(fn, sep="\t", index=False)
    return pd.read_csv(fn, sep="\t")


def fetch_fe(gsm, cache):
    fn = os.path.join(cache, f"{gsm}.fe.txt.gz")
    if not os.path.exists(fn):
        with urllib.request.urlopen(GEO.format(gsm).replace("view=data", "view=brief"), timeout=120) as r:
            url = [l.split(" = ", 1)[1] for l in r.read().decode().splitlines()
                   if l.startswith("!Sample_supplementary_file") and "_CGH" in l][0]
        urllib.request.urlretrieve(url.replace("ftp://", "https://"), fn)
    rows = []
    import gzip
    with gzip.open(fn, "rt", errors="replace") as fh:
        hdr = None
        for l in fh:
            f = l.rstrip("\n").split("\t")
            if f[0] == "FEATURES":
                hdr = {k: i for i, k in enumerate(f)}
            elif hdr and f[0] == "DATA" and f[hdr["ControlType"]] == "0":
                rows.append((f[hdr["SystematicName"]], f[hdr["LogRatio"]]))
    return pd.DataFrame(rows, columns=["loc", "lr10"])


def probes(d, platform, build=None):
    if platform.startswith("agilent"):
        m = d["loc"].str.extract(r"^(chr[0-9XY]+):0*(\d+)-0*(\d+)$")
        lr = pd.to_numeric(d["lr10"], errors="coerce") * np.log2(10)
        if platform.endswith("flip"):
            lr = -lr
        out = pd.DataFrame({"chrom": m[0], "pos": (pd.to_numeric(m[1]) + pd.to_numeric(m[2])) // 2, "lr": lr})
        return out.dropna(), build
    if platform == "nimblegen_hg18":
        m = d["ID_REF"].astype(str).str.extract(r"^CHR(\d+|X|Y)FS0*(\d+)$")
        chrom = m[0].str.lstrip("0").replace({"23": "X", "24": "Y"})
        return pd.DataFrame({"chrom": "chr" + chrom, "pos": pd.to_numeric(m[1]), "lr": pd.to_numeric(d["VALUE"],
                             errors="coerce")}).dropna(), "hg18"
    raise SystemExit(f"unknown platform {platform}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True, help="the WGS ratio dir (index.tsv is appended to)")
    ap.add_argument("--cache", default=os.path.join(ROOT, ".cache/geo_acgh"))
    ap.add_argument("--min-probes", type=int, default=3, help="probes needed in a 50 kb bin")
    ap.add_argument("--agilent-build", default="hg18", choices=["hg17", "hg18"])
    ap.add_argument("--only", help="comma list of names")
    ap.add_argument("--suffix", default="", help="appended to output names (build tests)")
    a = ap.parse_args()
    os.makedirs(a.cache, exist_ok=True)
    from pyliftover import LiftOver
    lifts = {}
    rows = []
    for name, key, role, gsm, platform in SAMPLES:
        if a.only and name not in a.only.split(","):
            continue
        d = fetch_fe(gsm, a.cache) if platform.startswith("agilent") else fetch(gsm, a.cache)
        p, asm = probes(d, platform, a.agilent_build)
        name = name + a.suffix
        if asm != "hg38":
            lo = lifts.setdefault(asm, LiftOver(CHAINS[asm]))
            hit = [lo.convert_coordinate(c, int(x) - 1) for c, x in zip(p.chrom, p.pos)]
            ok = np.array([bool(h) and h[0][0] == c for h, c in zip(hit, p.chrom)])
            p = p[ok].assign(pos=[int(h[0][1]) for h, k in zip(hit, ok) if k])
        out = {}
        for c, L in HG38_LEN.items():
            n = L // 50_000
            q = p[(p.chrom == c) & (p.pos // 50_000 < n)]
            i = (q.pos.values // 50_000).astype(int)
            s = np.bincount(i, weights=q.lr.values, minlength=n)
            k = np.bincount(i, minlength=n)
            v = np.full(n, np.nan)
            m = k >= a.min_probes
            v[m] = 2.0 ** (s[m] / k[m])
            out[c] = v.astype(np.float32)
        np.savez_compressed(os.path.join(a.out, f"{name}.npz"), **out, _reads=np.array([len(p)]))
        rows.append({"name": name, "key": key, "role": role, "reads": len(p),
                     "source": "array_" + platform.split("_")[0]})          # array_nimblegen, array_agilent
        print(f"[array] {name} {gsm}: {len(p):,} probes on hg38, bins with data "
              f"{sum(np.isfinite(v).sum() for v in out.values()):,}", file=sys.stderr)
    idx_f = os.path.join(a.out, "index.tsv")
    idx = pd.read_csv(idx_f, sep="\t")
    idx = idx[~idx.name.isin([r["name"] for r in rows])]
    pd.concat([idx, pd.DataFrame(rows)], ignore_index=True).to_csv(idx_f, sep="\t", index=False)


if __name__ == "__main__":
    main()
