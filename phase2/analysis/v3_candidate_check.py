#!/usr/bin/env python3
"""Signal verification of the pulled v3 candidates (ROADMAP Phase 3c, step A3b).

Runs on Roihu, where the per-sample grid columns live. For every candidate experiment:
  QC        MACS2 peaks (the >= 2,000 S3norm gate), dynamic range, SE count (cnrose qc.json)
  H3K27ac?  fraction of its bed20 peaks on the fixed grid, and the correlation of its SE profile with the
            atlas-wide mean profile (input, IgG or another mark correlate poorly)
  identity  correlation with every scored line's centroid; rank of the labelled line (scored lines) or
            of the labelled lineage among lineage centroids built without that line (new lines)
Profiles are log1p SE signal over the v2 union catalog (grid rows summed, as aggregate.py does); Pearson on
log1p is close to invariant to S3norm's per-sample power transform, so raw candidate sums compare with the
S3-normalised atlas. Baselines are the same statistics for the atlas's own experiments, leave-one-out.

  python v3_candidate_check.py --work /scratch/$PROJ/se-cacts/phase2 --meta-dir <dir> --out <tsv>
"""
from __future__ import annotations

import argparse
import json
import os

import numpy as np
import pandas as pd


def z(M, axis=0):
    M = M - M.mean(axis=axis, keepdims=True)
    return M / np.maximum(M.std(axis=axis, keepdims=True), 1e-12)


def read_bed(path, cols=3):
    return pd.read_csv(path, sep="\t", header=None, usecols=range(cols), names=["c", "s", "e"][:cols],
                       dtype={"c": str, "s": np.int64, "e": np.int64}, comment="#")


def catalog_rows(cat, grid):
    g = {c: (d.s.values, d.e.values, d.index.values) for c, d in grid.groupby("c")}
    out = []
    for c, s, e in cat.itertuples(index=False):
        if c not in g:
            out.append(np.empty(0, int)); continue
        gs, ge, gi = g[c]
        out.append(gi[np.searchsorted(ge, s, side="right"):np.searchsorted(gs, e, side="left")])
    return out


def in_grid_fraction(peaks, grid_by_chrom):
    hit = 0
    for c, d in peaks.groupby("c"):
        if c not in grid_by_chrom:
            continue
        gs, ge = grid_by_chrom[c]
        j = np.searchsorted(ge, d.s.values, side="right")
        ok = j < len(gs)
        hit += int((ok & (gs[np.minimum(j, len(gs) - 1)] < d.e.values)).sum())
    return hit / len(peaks) if len(peaks) else np.nan


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--work", required=True, help="phase2 work dir on Roihu (out/, data/, results_v2/)")
    ap.add_argument("--meta-dir", required=True,
                    help="dir with expansion_v3_candidates.tsv, pull_set.v2.tsv, lineage_resolved.tsv")
    ap.add_argument("--new-lines", help="cvcl, cell_line, lineage, OncotreePrimaryDisease for lines not yet in "
                    "lineage_resolved.tsv (from DepMap Model.csv by RRID)")
    ap.add_argument("--candidates", help="TSV with srx, cvcl, line_status (scored_v2|new), route, and optional "
                    "key (expected scored line), cell_line, lineage; default: <meta-dir>/expansion_v3_candidates.tsv "
                    "(qc_pass rows)")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    W, M = a.work, a.meta_dir

    cand = pd.read_csv(a.candidates or f"{M}/expansion_v3_candidates.tsv", sep="\t")
    if "qc_pass" in cand:
        cand = cand[cand.qc_pass == True]
    cand = cand.set_index("srx")
    ps = pd.read_csv(f"{M}/pull_set.v2.tsv", sep="\t").set_index("srx")
    lin = pd.read_csv(f"{M}/lineage_resolved.tsv", sep="\t").drop_duplicates("cvcl").set_index("cvcl")
    disease = pd.Series(dtype=str)
    if a.new_lines:
        nl = pd.read_csv(a.new_lines, sep="\t").drop_duplicates("cvcl").set_index("cvcl")
        lin = pd.concat([lin[~lin.index.isin(nl.index)], nl[["cell_line", "lineage"]]])
        disease = nl.OncotreePrimaryDisease

    grid = read_bed(f"{W}/data/grid.20.bed")
    cat = read_bed(f"{W}/results_v2/atlas.s3.union_catalog.bed")
    rows = catalog_rows(cat, grid)
    gbc = {c: (d.s.values, d.e.values) for c, d in grid.groupby("c")}

    sig = f"{W}/results_v2/atlas.s3.se_signal.tsv"
    hdr = open(sig).readline().rstrip("\n").split("\t")
    S = pd.read_csv(sig, sep="\t", index_col=0, dtype={c: np.float32 for c in hdr[1:]})
    assert len(S) == len(cat), f"catalog {len(cat)} rows != signal {len(S)}"
    keep = (S.values > 0).any(axis=1)
    cols = [c for c in S.columns if c in ps.index]
    X = np.log1p(S.loc[keep, cols].values.astype(np.float64))            # SEs x atlas experiments
    rows = [r for r, k in zip(rows, keep) if k]
    n = X.shape[0]
    print(f"[check] atlas {len(cols)} experiments x {n} SEs; {len(cand)} candidates")

    keys = ps.loc[cols, "key"].values
    lines = pd.unique(keys)
    line_lineage = ps.loc[cols].groupby("key").lineage.first()
    line_name = ps.loc[cols].groupby("key").cell.first()
    C = np.stack([X[:, keys == k].mean(1) for k in lines], axis=1)       # raw line centroids
    nper = np.array([(keys == k).sum() for k in lines])
    Zc = z(C)
    gz = z(X.mean(1, keepdims=True))[:, 0]
    lineages = sorted(line_lineage.unique())
    lmask = {g: np.array([line_lineage[k] == g for k in lines]) for g in lineages}
    Lsum = {g: C[:, lmask[g]].sum(1) for g in lineages}                   # sum of line centroids
    Lz = {g: z((Lsum[g] / lmask[g].sum())[:, None])[:, 0] for g in lineages}
    line_idx = {k: i for i, k in enumerate(lines)}

    def lineage_rank(x_z, own_lineage, own_line=None):
        """Rank (1 = best) of own_lineage among lineage centroids; own line left out of its centroid."""
        rs = {g: float(Lz[g] @ x_z / n) for g in lineages}
        if own_line in line_idx and line_lineage.get(own_line) == own_lineage:
            k = lmask[own_lineage].sum()
            if k < 2:
                rs.pop(own_lineage)
            else:
                loo = (Lsum[own_lineage] - C[:, line_idx[own_line]]) / (k - 1)
                rs[own_lineage] = float(z(loo[:, None])[:, 0] @ x_z / n)
        if own_lineage not in rs:
            return np.nan, max(rs, key=rs.get)
        order = sorted(rs, key=rs.get, reverse=True)
        return order.index(own_lineage) + 1, order[0]

    # ---- baselines on the atlas's own experiments
    Zx = z(X)
    r_global_atlas = gz @ Zx / n
    R = Zx.T @ Zc / n                                                    # experiments x lines
    base_line_rank, base_lin_rank = [], []
    rng = np.random.default_rng(0)
    for j, k in enumerate(keys):
        li = int(np.where(lines == k)[0][0])
        if nper[li] >= 2:
            loo = z(((C[:, li] * nper[li] - X[:, j]) / (nper[li] - 1))[:, None])[:, 0]
            r = R[j].copy(); r[li] = loo @ Zx[:, j] / n
            base_line_rank.append(int((r > r[li]).sum()) + 1)
    for j in rng.choice(len(cols), size=min(600, len(cols)), replace=False):
        base_lin_rank.append(lineage_rank(Zx[:, j], line_lineage[keys[j]], own_line=keys[j])[0])
    base_line_rank = np.array(base_line_rank)
    base_lin_rank = np.array([x for x in base_lin_rank if x == x])
    thr_global = float(np.percentile(r_global_atlas, 1))
    thr_line = int(np.percentile(base_line_rank, 95))
    thr_lin = int(np.percentile(base_lin_rank, 95))
    print(f"[check] baselines: r(global) atlas p1 {thr_global:.3f} (median {np.median(r_global_atlas):.3f}); "
          f"own-line rank LOO median {np.median(base_line_rank):.0f}, top-1 {np.mean(base_line_rank == 1):.0%}, "
          f"p95 {thr_line}; own-lineage rank (line left out) top-1 {np.mean(base_lin_rank == 1):.0%}, p95 {thr_lin}")

    # ---- candidates
    cvcl2key = ps.drop_duplicates("cvcl").set_index("cvcl").key
    out = []
    for srx, r in cand.iterrows():
        d = f"{W}/out/{srx[-2:]}"
        row = {"srx": srx, "cvcl": r.cvcl, "line_status": r.line_status, "route": r.route,
               "cell_line": r.get("cell_line") if isinstance(r.get("cell_line"), str) else lin.cell_line.get(r.cvcl, ""),
               "lineage": r.get("lineage") if isinstance(r.get("lineage"), str) else lin.lineage.get(r.cvcl, "")}
        if not os.path.exists(f"{d}/{srx}.done"):
            row["status"] = "not pulled"; out.append(row); continue
        q = json.load(open(f"{d}/{srx}.qc.json"))
        row.update(n_peaks=q.get("n_peaks"), n_super=q.get("n_super"),
                   dyn_range=q.get("dynamic_range_p99_over_median"), frip=q.get("frip_proxy"))
        pk = f"{W}/data/bed20/{srx}.20.bed"
        row["in_grid"] = in_grid_fraction(read_bed(pk), gbc) if os.path.exists(pk) else np.nan
        g = np.fromfile(f"{d}/{srx}.grid.20.f32", dtype="<f4").astype(np.float64)
        x = np.log1p(np.array([g[i].sum() if i.size else 0.0 for i in rows]))
        xz = z(x[:, None])[:, 0]
        row["r_global"] = float(gz @ xz / n)
        rl = Zc.T @ xz / n
        best = int(np.argmax(rl))
        row["best_line"], row["best_line_r"] = line_name[lines[best]], float(rl[best])
        row["best_line_lineage"] = line_lineage[lines[best]]
        key = r.get("key") if isinstance(r.get("key"), str) else cvcl2key.get(r.cvcl)
        if key is not None and key in set(lines):
            li = int(np.where(lines == key)[0][0])
            row["own_line_rank"] = int((rl > rl[li]).sum()) + 1
            row["own_line_r"] = float(rl[li])
        own_lin = row["lineage"]
        if own_lin in lineages:
            row["own_lineage_rank"], row["best_lineage"] = lineage_rank(xz, own_lin, own_line=key)
        flags = ["non-cancer"] if disease.get(r.cvcl) == "Non-Cancerous" else []
        if (row["n_peaks"] or 0) < 2000:
            flags.append("peaks<2000")
        if row["r_global"] < thr_global:
            flags.append("not-H3K27ac-like")
        if row.get("own_line_rank", 0) > thr_line:
            flags.append("identity?")
        if r.line_status == "new" and row.get("own_lineage_rank", 0) > thr_lin:
            flags.append("lineage?")
        row["flags"] = ",".join(flags)
        row["status"] = "pass" if not flags else "flag"
        out.append(row)

    df = pd.DataFrame(out)
    df.to_csv(a.out, sep="\t", index=False, float_format="%.4g")
    ok = df[df.status != "not pulled"]
    print(f"[check] {len(ok)} checked: {int((ok.status == 'pass').sum())} pass, {int((ok.status == 'flag').sum())} flagged")
    for f in ["non-cancer", "peaks<2000", "not-H3K27ac-like", "identity?", "lineage?"]:
        print(f"  {f:17s} {int(ok['flags'].fillna('').str.contains(f, regex=False).sum())}")
    print(f"  in-grid peak fraction median {ok.in_grid.median():.3f} (p10 {ok.in_grid.quantile(.1):.3f}); "
          f"r(global) median {ok.r_global.median():.3f} (atlas {np.median(r_global_atlas):.3f})")
    sc = ok[ok.line_status == "scored_v2"]
    if len(sc):
        print(f"  scored-line replicates: own line top-1 {np.mean(sc.own_line_rank == 1):.0%}, "
              f"top-5 {np.mean(sc.own_line_rank <= 5):.0%} (atlas LOO top-1 {np.mean(base_line_rank == 1):.0%})")
    nw = ok[ok.line_status == "new"]
    if len(nw):
        print(f"  new lines: {nw.cvcl.nunique()} lines, own lineage top-1 {np.mean(nw.own_lineage_rank == 1):.0%} "
              f"(atlas baseline {np.mean(base_lin_rank == 1):.0%})")
        pas = nw[nw.status == "pass"]
        print(f"  new lines with >= 1 passing experiment: {pas.cvcl.nunique()}; with >= 2 studies unknown here")


if __name__ == "__main__":
    main()
