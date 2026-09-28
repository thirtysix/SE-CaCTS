#!/usr/bin/env python3
"""ChIP-Atlas pipeline v1 vs v2 on the same experiments (FINDINGS 2026-09-28 §13, §15).

v1 = ChIP-Atlas's published hg38 files, exactly as our atlas used them. v2 = their v2 container
(fastp -> bwa-mem2 -> MACS3, BAMPE for pairs) run by us on the same raw reads. Every metric is split by
library layout, because v1 aligned paired-end data as single reads while v2 aligns pairs.

  sample     one experiment -> <cmp>/<SRX>.metrics.json          (Roihu: bigWigs live on scratch)
  summarize  all metrics.json -> summary.tsv + medians by layout   (anywhere)

Both bigWigs go through the same code path, cnrose.io.scan_bigwig, so the grid values and 1 kb bins are
in identical units (integral of coverage; uncovered bases = 0 for bins). cnrose must be on PYTHONPATH.
"""
from __future__ import annotations

import argparse
import glob
import gzip
import hashlib
import json
import os
import sys

import numpy as np

HLA = ("chr6", 28_510_120, 33_480_577)          # GRCh38 MHC, the densest alt-haplotype region
BIN = 1000


# ------------------------------------------------------------------ intervals
def by_chrom(ivs):
    d = {}
    for c, s, e in ivs:
        d.setdefault(c, []).append((s, e))
    return {c: np.array(sorted(v), dtype=np.int64).reshape(-1, 2) for c, v in d.items()}


def merge(a):
    if len(a) == 0:
        return a
    out = [list(a[0])]
    for s, e in a[1:]:
        if s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return np.array(out, dtype=np.int64)


def frac_hit(A, B):
    """Fraction of A's intervals overlapping any interval of B."""
    hit = n = 0
    for c, a in A.items():
        n += len(a)
        b = merge(B.get(c, np.empty((0, 2), np.int64)))
        if not len(b):
            continue
        j = np.searchsorted(b[:, 1], a[:, 0], side="right")      # first b ending after a starts
        ok = j < len(b)
        hit += int((ok & (b[np.minimum(j, len(b) - 1), 0] < a[:, 1])).sum())
    return hit / n if n else float("nan")


def bp_jaccard(A, B):
    inter = union = 0
    for c in set(A) | set(B):
        a = merge(A.get(c, np.empty((0, 2), np.int64)))
        b = merge(B.get(c, np.empty((0, 2), np.int64)))
        la, lb = int((a[:, 1] - a[:, 0]).sum()) if len(a) else 0, int((b[:, 1] - b[:, 0]).sum()) if len(b) else 0
        i = j = ov = 0
        while i < len(a) and j < len(b):
            ov += max(0, min(a[i, 1], b[j, 1]) - max(a[i, 0], b[j, 0]))
            if a[i, 1] < b[j, 1]:
                i += 1
            else:
                j += 1
        inter += ov
        union += la + lb - ov
    return inter / union if union else float("nan")


def interval_metrics(a_ivs, b_ivs):
    A, B = by_chrom(a_ivs), by_chrom(b_ivs)
    return {"n_v1": len(a_ivs), "n_v2": len(b_ivs), "v1_recovered_by_v2": frac_hit(A, B),
            "v2_recovered_by_v1": frac_hit(B, A), "bp_jaccard": bp_jaccard(A, B)}


# ------------------------------------------------------------------ stats
def spearman(x, y):
    from scipy.stats import spearmanr
    return float(spearmanr(x, y).correlation)


def top_overlap(x, y, k=500):
    return len(set(np.argsort(-x)[:k]) & set(np.argsort(-y)[:k])) / k


def catalog_rows(catalog, grid):
    """Grid rows under each atlas SE (the same lookup aggregate.py uses; grid rows are sorted, disjoint)."""
    g = {}
    for i, (c, s, e) in enumerate(grid):
        g.setdefault(c, ([], [], []))
        g[c][0].append(s); g[c][1].append(e); g[c][2].append(i)
    g = {c: tuple(np.asarray(v) for v in t) for c, t in g.items()}
    rows = []
    for c, s, e in catalog:
        if c not in g:
            rows.append(np.empty(0, int)); continue
        gs, ge, gi = g[c]
        rows.append(gi[np.searchsorted(ge, s, side="right"):np.searchsorted(gs, e, side="left")])
    return rows


def se_sum(grid_vals, rows):
    return np.array([grid_vals[r].sum() if r.size else 0.0 for r in rows])


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 22), b""):
            h.update(chunk)
    return h.hexdigest()


def read_bed_gz(path, canonical):
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as fh:
        return [(f[0], int(f[1]), int(f[2])) for f in (l.rstrip("\n").split("\t") for l in fh)
                if f[0] in canonical]


def alt_primary_regions(path, canonical):
    """Primary-assembly spans that have an alt haplotype in UCSC hg38 (altLocations.txt.gz)."""
    out = []
    with gzip.open(path, "rt") as fh:
        for l in fh:
            f = l.rstrip("\n").split("\t")
            if f[1] in canonical:
                out.append((f[1], int(f[2]), int(f[3])))
    return out


def bin_mask(layout, regions):
    off = {L["chrom"]: (L["offset"], L["nbins"]) for L in layout}
    m = np.zeros(sum(L["nbins"] for L in layout), bool)
    for c, s, e in regions:
        if c in off:
            base, nb = off[c]
            m[base + max(0, s // BIN): base + min(nb, (e - 1) // BIN + 1)] = True
    return m


# ------------------------------------------------------------------ one sample
def cmd_sample(a):
    from cnrose.io import CANONICAL, read_bed3, read_f32, scan_bigwig
    import pyBigWig

    srx, v1d, v2d, cd = a.srx, a.v1_dir, a.v2_dir, a.cmp_dir
    grid = read_bed3(a.grid)
    out = {"srx": srx, "layout": a.layout}

    # provenance: is the v1 bigWig we fetched today the one the atlas pulled in July?
    prov = json.load(open(f"{v1d}/{srx}.provenance.json"))
    out["v1_bw_same_as_atlas"] = sha256(a.v1_bw) == prov["sha256"]

    # QC
    v2s = open(f"{v2d}/{srx}.stats.tsv").read().rstrip("\n").split("\t")
    out["v2_stats_tsv"] = v2s
    fp = glob.glob(f"{v2d}/*fastp*.json")
    if fp:
        j = json.load(open(fp[0]))
        out["fastp"] = {"reads_before": j["summary"]["before_filtering"]["total_reads"],
                        "reads_after": j["summary"]["after_filtering"]["total_reads"],
                        "read_len_after": j["summary"]["after_filtering"].get("read1_mean_length"),
                        "dup_rate": j.get("duplication", {}).get("rate"),
                        "insert_size_peak": j.get("insert_size", {}).get("peak")}
    for tag, p in (("v1", f"{v1d}/{srx}.qc.json"), ("v2", f"{cd}/{srx}.v2.qc.json")):
        if os.path.exists(p):
            q = json.load(open(p))
            out[f"{tag}_qc"] = {k: q.get(k) for k in ("dynamic_range_p99_over_median", "frip_proxy",
                                                       "zero_fraction", "n_peaks", "n_super")}

    # peaks (canonical chromosomes, as the atlas uses)
    p1 = read_bed3(a.v1_peaks, canonical_only=True)
    p2 = read_bed3(f"{v2d}/{srx}.20_peaks.narrowPeak", canonical_only=True)
    out["peaks"] = interval_metrics(p1, p2)
    out["peaks"]["v2_noncanonical"] = len(read_bed3(f"{v2d}/{srx}.20_peaks.narrowPeak")) - len(p2)

    # coverage: one scan per bigWig, same code path -> grid (exact) + 1 kb bins
    (g1,), b1, lay1 = scan_bigwig(a.v1_bw, [grid], binsize=BIN)
    (g2,), b2, lay2 = scan_bigwig(f"{v2d}/{srx}.bw", [grid], binsize=BIN)
    assert [L["chrom"] for L in lay1] == [L["chrom"] for L in lay2], "chrom layouts differ"
    stored = read_f32(f"{v1d}/{srx}.grid.20.f32", len(grid))
    out["v1_grid_rescan_max_rel_err"] = float(np.max(np.abs(g1 - stored) / np.maximum(np.abs(stored), 1e-6)))

    x, y = b1.astype(np.float64), b2.astype(np.float64)
    tot = {}
    for tag, p, b in (("v1", a.v1_bw, x), ("v2", f"{v2d}/{srx}.bw", y)):
        bw = pyBigWig.open(p); hdr = bw.header(); bw.close()
        tot[tag] = {"canonical": float(b.sum() * BIN), "all": float(hdr["sumData"])}
    nz = (x > 0) | (y > 0)
    sig = (x > np.median(x[x > 0])) & (y > 0)
    lr = np.log2(y[sig] / x[sig]); lr -= np.median(lr)
    alt = bin_mask(lay1, alt_primary_regions(a.alt, CANONICAL))
    hla = bin_mask(lay1, [HLA])
    v1sig = x > np.median(x[x > 0])
    out["coverage_1kb"] = {
        "spearman": spearman(x[nz], y[nz]),
        "pearson_log1p": float(np.corrcoef(np.log1p(x[nz]), np.log1p(y[nz]))[0, 1]),
        "v2_over_v1_total": tot["v2"]["canonical"] / tot["v1"]["canonical"],
        "v2_noncanonical_frac": 1 - tot["v2"]["canonical"] / tot["v2"]["all"],
        "v1_noncanonical_frac": 1 - tot["v1"]["canonical"] / tot["v1"]["all"],
        "alt_bins_frac_genome": float(alt.mean()),
        "alt_log2_v2_v1_centered": float(np.median(lr[alt[sig]])),
        "nonalt_log2_v2_v1_centered": float(np.median(lr[~alt[sig]])),
        "hla_log2_v2_v1_centered": float(np.median(lr[hla[sig]])) if hla[sig].any() else None,
        "alt_dropout": float((y[v1sig & alt] == 0).mean()),        # v1 signal, v2 nothing
        "nonalt_dropout": float((y[v1sig & ~alt] == 0).mean()),
    }

    # grid + atlas SE level (raw grid sums over the v2 atlas union catalog, as aggregate.py builds it)
    anyg = (g1 > 0) | (g2 > 0)
    out["grid"] = {"spearman": spearman(g1[anyg], g2[anyg]), "top500": top_overlap(g1, g2)}
    rows = catalog_rows(read_bed_gz(a.catalog, CANONICAL), grid)
    s1, s2 = se_sum(g1, rows), se_sum(g2, rows)
    out["atlas_se"] = {"spearman": spearman(s1, s2), "top500": top_overlap(s1, s2),
                       "top100": top_overlap(s1, s2, 100)}

    # natural replicates: other experiments of the same line, scored the same way (v1 files from the pull)
    if a.replicates and os.path.exists(a.replicates):
        reps = [l.rstrip("\n").split("\t") for l in open(a.replicates)][1:]
        reps = [r for r in reps if r[0] == srx]
        rr = []
        for _, rep, same_study, rep_layout in reps:
            f = f"{a.pull_out}/{rep[-2:]}/{rep}.grid.20.f32"
            if not os.path.exists(f):
                continue
            s3 = se_sum(read_f32(f, len(grid)).astype(np.float64), rows)
            rr.append({"rep": rep, "same_study": same_study == "1", "rep_layout": rep_layout,
                       "rho_v1": spearman(s1, s3), "rho_v2": spearman(s2, s3)})
        out["replicates"] = rr

    # SE calls: v1 (the pull), v2 (v2 bigWig + v2 peaks), and v2 bigWig on v1 peaks (coverage effect alone)
    se1 = read_bed3(f"{v1d}/{srx}.se.bed", canonical_only=True)
    out["se"] = interval_metrics(se1, read_bed3(f"{cd}/{srx}.v2.se.bed", canonical_only=True))
    out["se_v2bw_v1peaks"] = interval_metrics(se1, read_bed3(f"{cd}/{srx}.v2bw_v1pk.se.bed",
                                                              canonical_only=True))

    with open(a.out, "w") as fh:
        json.dump(out, fh, indent=1)
    print(f"[cmp] {srx} {a.layout}: SE rho {out['atlas_se']['spearman']:.3f}, "
          f"SE recovery {out['se']['v1_recovered_by_v2']:.2f}, peaks recovery "
          f"{out['peaks']['v1_recovered_by_v2']:.2f}, alt {out['coverage_1kb']['alt_log2_v2_v1_centered']:+.2f}")


# ------------------------------------------------------------------ summary
def cmd_summarize(a):
    import pandas as pd
    from scipy.stats import mannwhitneyu

    rows = []
    for p in sorted(glob.glob(os.path.join(a.cmp_dir, "*.metrics.json"))):
        m = json.load(open(p))
        r = {"srx": m["srx"], "layout": m["layout"], "v1_bw_ok": m["v1_bw_same_as_atlas"],
             "v1_rescan_err": m["v1_grid_rescan_max_rel_err"]}
        for blk in ("peaks", "se", "se_v2bw_v1peaks", "coverage_1kb", "grid", "atlas_se"):
            for k, v in m.get(blk, {}).items():
                r[f"{blk}.{k}"] = v
        for tag in ("v1_qc", "v2_qc"):
            for k, v in (m.get(tag) or {}).items():
                r[f"{tag}.{k}"] = v
        if m.get("fastp"):
            r["fastp.insert_size_peak"] = m["fastp"]["insert_size_peak"]
            r["fastp.dup_rate"] = m["fastp"]["dup_rate"]
        reps = [x for x in m.get("replicates", []) if not x["same_study"]]
        if reps:
            r["rep.n_cross_study"] = len(reps)
            r["rep.rho_v1_median"] = float(np.median([x["rho_v1"] for x in reps]))
            r["rep.rho_v2_median"] = float(np.median([x["rho_v2"] for x in reps]))
            for lay in ("SINGLE", "PAIRED"):
                sub = [x for x in reps if x["rep_layout"] == lay]
                if sub:
                    r[f"rep.{lay}.rho_v1"] = float(np.median([x["rho_v1"] for x in sub]))
                    r[f"rep.{lay}.rho_v2"] = float(np.median([x["rho_v2"] for x in sub]))
        rows.append(r)
    df = pd.DataFrame(rows).set_index("srx")
    df.to_csv(os.path.join(a.cmp_dir, "summary.tsv"), sep="\t")
    num = df.select_dtypes("number")
    med = num.groupby(df["layout"]).median().T
    pv = {c: mannwhitneyu(num.loc[df.layout == "SINGLE", c].dropna(),
                          num.loc[df.layout == "PAIRED", c].dropna()).pvalue
          if num[c].notna().groupby(df.layout).sum().min() >= 3 else np.nan for c in num.columns}
    med["MWU_p"] = pd.Series(pv)
    pd.set_option("display.width", 200)
    print(f"{len(df)} samples ({(df.layout == 'SINGLE').sum()} SE, {(df.layout == 'PAIRED').sum()} PE); "
          f"v1 bigWig identical to the atlas's: {int(df.v1_bw_ok.sum())}/{len(df)}; "
          f"max v1 re-scan error {df.v1_rescan_err.max():.2g}")
    print(med.round(3).to_string())


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("sample")
    s.add_argument("--srx", required=True)
    s.add_argument("--layout", required=True)
    s.add_argument("--v1-bw", required=True, help="ChIP-Atlas v1 bigWig, fetched again")
    s.add_argument("--v1-dir", required=True, help="the pull's shard dir holding <SRX>.se.bed/.grid.20.f32")
    s.add_argument("--v1-peaks", required=True, help="ChIP-Atlas v1 bed20")
    s.add_argument("--v2-dir", required=True, help="pipeline-v2.sh --outdir for this sample")
    s.add_argument("--cmp-dir", required=True, help="where the v2 cnrose calls were written")
    s.add_argument("--grid", required=True)
    s.add_argument("--catalog", required=True, help="atlas union catalog (atlas.s3.union_catalog.bed.gz)")
    s.add_argument("--alt", required=True, help="UCSC hg38 altLocations.txt.gz")
    s.add_argument("--replicates", help="replicates.tsv (srx, rep, same_study, rep_layout)")
    s.add_argument("--pull-out", help="the pull's out/ root (for replicate grid columns)")
    s.add_argument("--out", required=True)
    m = sub.add_parser("summarize")
    m.add_argument("--cmp-dir", required=True)
    a = ap.parse_args(argv)
    {"sample": cmd_sample, "summarize": cmd_summarize}[a.cmd](a)


if __name__ == "__main__":
    main()
