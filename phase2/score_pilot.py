#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Phase-4 HIERARCHICAL SE specificity scoring (pilot): score each super-enhancer's lineage specificity at
multiple resolutions of the DepMap Oncotree hierarchy — cancer type (OncotreeLineage), primary disease,
subtype, and individual cell line — reusing ../pyCaCTS (build_rep_matrix + JSD + empirical FDR).

Pipeline (ROADMAP Phase 4, all local): SE x sample signal (aggregate.py) -> per-sample SCORING-TIME CN
correction (symmetric; 2nd call site of cnrose.cn.correct) -> `--norm` across samples -> collapse
replicate SRX to cell line -> for each hierarchy level, per-group MEAN (pyCaCTS build_rep_matrix) -> CaCTS
JSD (lower = more specific) -> empirical-null FDR -> recovery of known master-TF SEs per group.

Normalization: the recommended path is `aggregate.py --norm s3norm` (fitted on the fine grid, the correct
resolution for a nonlinear transform) followed by `--norm none` here. `--norm quantile` is the pilot
baseline, retained for comparison — see `normalize()` and s3norm.py.

Specificity threshold: FDR <= 0.10 under a per-group empirical null with **global BH** (`--fdr-scope global`,
the default; see specificity.py). Sharing one testing budget across all SE x group tests is what makes
per-group counts comparable — `--fdr-scope pergroup` reproduces pyCaCTS's behaviour, under which a group can
return zero calls. Note this changes only WHICH SEs pass, never the JSD rankings.

Specificity is inherently multi-resolution: an SE can be pan-lineage-specific (e.g. ESR1 in Breast) or
subtype-specific (GATA1 in erythroid CML vs SPI1 in AML). 13 samples is a proof-of-mechanism. atac_hdac env.
"""
from __future__ import annotations

import argparse
import gzip
import os
import sys
import time

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
SECACTS = os.path.dirname(HERE)
_ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
sys.path.insert(0, _ROOT)
from secacts_env import DATAROOT, cache_path                      # noqa: E402  (local paths live in .env)
sys.path.insert(0, os.path.join(DATAROOT, "002.AI_projects", "pyCaCTS"))
sys.path.insert(0, os.path.join(SECACTS, "cnrose"))

from pycacts.score import cacts_score_matrix, rank_specific        # noqa: E402
from pycacts.grouping import build_rep_matrix                      # noqa: E402
from specificity import fdr_matrix                                 # noqa: E402
from permutation import permutation_fdr                            # noqa: E402
from cnrose.cn.depmap import load_gene_coords, DepMapGeneCN        # noqa: E402
from cnrose.cn.base import correct                                 # noqa: E402

# master-TF / identity genes keyed by keywords in the Oncotree group label (specific rules first).
IDENTITY_RULES = [
    (("serous ovarian", "ovarian epithelial", "ovary", "fallopian"), ["PAX8", "SOX17", "WT1", "MECOM", "GATA6"]),
    (("ductal carcinoma", "breast"),               ["ESR1", "FOXA1", "GATA3", "TFAP2C", "SPDEF", "GRHL2", "XBP1"]),
    (("colon", "colorectal", "bowel"),             ["CDX2", "HNF4A", "ASCL2", "VDR", "KLF5", "CDX1", "TCF7L2"]),
    (("t-cell", "t-lymphoblastic", "lymphoid"),    ["TAL1", "LMO2", "RUNX1", "TLX1", "TLX3", "LEF1", "TCF7", "MYB"]),
    (("acute myeloid",),                           ["SPI1", "CEBPA", "IRF8", "MEF2C", "RUNX1", "MYB"]),
    (("chronic myeloid", "myeloproliferative"),    ["GATA1", "GATA2", "TAL1", "KLF1", "MYB"]),   # K562 = erythroid
    (("myeloid",),                                 ["SPI1", "CEBPA", "GATA1", "GATA2", "RUNX1", "MYB"]),
]
LEVELS = ["OncotreeLineage", "OncotreePrimaryDisease", "OncotreeSubtype", "line"]
# pilot SRX -> DepMap ModelID for lines absent from pull_set (no WGS CN, but Oncotree-annotated)
EXTRA_MODEL = {"SRX10809652": "ACH-000768"}   # MDA-MB-231


def identity_for(label):
    ll = str(label).lower()
    for keys, genes in IDENTITY_RULES:
        if any(k in ll for k in keys):
            return genes
    return []


def quantile_normalize(M):
    ranks = np.argsort(np.argsort(M, axis=0), axis=0)
    ref = np.sort(M, axis=0).mean(axis=1)
    return ref[ranks]


def normalize(M, how, samples):
    """Cross-sample normalization of the SE x sample matrix.

    'none'     — the matrix is already normalized upstream. This is the RECOMMENDED Phase-2 path:
                 `aggregate.py --norm s3norm` fits S3norm on the FINE GRID and re-sums, which is the
                 correct resolution for a nonlinear transform (Σf(xᵢ) ≠ f(Σxᵢ)).
    'quantile' — the pilot baseline. Kept for comparison only: forcing identical distributions
                 manufactures false positives (DESIGN.md §normalization).
    's3norm'   — S3norm fitted directly on SE sums. An APPROXIMATION of the grid-level fit above; use
                 when only the SE matrix is available.
    """
    if how == "none":
        return M
    if how == "quantile":
        return quantile_normalize(M)
    if how == "s3norm":
        from s3norm import s3norm_matrix
        return s3norm_matrix(M, ref="medoid", srx=samples)[0].astype(float)
    raise ValueError(f"unknown normalization: {how!r}")


def load_coords(bed):
    """Union-catalog BED -> {se_id: (chrom, start, end)}. Transparently reads .gz (the at-scale
    atlas catalogs are gzipped by scripts/51_compress_results.sh)."""
    opener = gzip.open if bed.endswith(".gz") else open
    coords = {}
    with opener(bed, "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            coords[f[3]] = (f[0], int(f[1]), int(f[2]))
    return coords


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--signal", default=os.path.join(SECACTS, "phase2/rehearse/pilot.se_signal.tsv"))
    ap.add_argument("--catalog", default=os.path.join(SECACTS, "phase2/rehearse/pilot.union_catalog.bed"))
    ap.add_argument("--pull-set", default=os.path.join(SECACTS, "phase2/data/pull_set.tsv"))
    ap.add_argument("--model", default=os.path.join(DATAROOT, "DepMap/2026q1/Model.csv"))
    ap.add_argument("--cn-gene-csv", default=os.path.join(DATAROOT, "DepMap/2026q1/OmicsCNGeneWGS.csv"))
    ap.add_argument("--gtf", default=os.path.join(DATAROOT, "0.human_genome/Homo_sapiens.GRCh38.106.chr.gtf.gz"))
    ap.add_argument("--gene-cache", default=cache_path("gene_coords.GRCh38.106.tsv"))
    ap.add_argument("--out", default=os.path.join(SECACTS, "phase2/rehearse/pilot_scores"))
    ap.add_argument("--topn", type=int, default=15)
    ap.add_argument("--fdr", type=float, default=0.10, help="empirical-null FDR threshold for 'specific'")
    ap.add_argument("--norm", default="quantile", choices=["quantile", "s3norm", "none"],
                    help="cross-sample normalization (default quantile = the pilot baseline). Use 'none' "
                         "when the signal matrix came from `aggregate.py --norm s3norm` (grid-level fit).")
    ap.add_argument("--fdr-method", default="analytic", choices=["analytic", "permutation"],
                    help="'analytic' (default) fits a normal to each group's JSD column — FAST but "
                         "MIS-CALIBRATED (gotcha 70: it over-calls ~10x because JSD is bounded and "
                         "left-skewed). 'permutation' measures the null by shuffling which line carries "
                         "which group label; degenerate at 'line' level, which it will skip.")
    ap.add_argument("--n-perm", type=int, default=50, help="permutations for --fdr-method permutation")
    ap.add_argument("--dump-specific", type=float, default=None, metavar="FDR",
                    help="also write EVERY SE x group test at or below this FDR, long-format and gzipped, to "
                         "<out>.<level>.specific.tsv.gz. The top_specific tables only carry --topn rows per "
                         "group, which is an optimistic sample for any downstream analysis (Phase 6).")
    ap.add_argument("--no-cn", action="store_true",
                    help="skip APPLYING scoring-time CN correction — the ablation arm. CN is still evaluated "
                         "and reported as cn_mean, so the two arms differ in exactly one step. Compare a "
                         "corrected run against a --no-cn run to label amplification-driven SEs.")
    ap.add_argument("--fdr-null", default="pergroup", choices=["pergroup", "global"],
                    help="empirical-null calibration (see specificity.py). 'global' pools all groups and is "
                         "over-conservative for tight groups — it drops SOX17 below significance.")
    ap.add_argument("--fdr-scope", default="global", choices=["global", "pergroup"],
                    help="Benjamini-Hochberg scope. 'global' (default) shares ONE testing budget across all "
                         "SE x group tests, making per-group counts comparable; 'pergroup' reproduces "
                         "pyCaCTS.empirical_fdr, under which a group can return zero calls.")
    a = ap.parse_args()

    M = pd.read_csv(a.signal, sep="\t", index_col=0)               # SE x SRX (uncorrected)
    coords = load_coords(a.catalog)
    samples = list(M.columns)
    se_ids = list(M.index)

    ps = pd.read_csv(a.pull_set, sep="\t")
    srx_model = dict(zip(ps["srx"], ps["model_id"]))
    srx_model.update(EXTRA_MODEL)
    model = pd.read_csv(a.model, index_col="ModelID")
    name_of = model["StrippedCellLineName"].to_dict()

    # scoring-time CN over each SE region, per sample (symmetric per-copy, floor 0.1)
    gene_coords = load_gene_coords(a.gtf, cache_path=a.gene_cache)
    prov = DepMapGeneCN(a.cn_gene_csv, gene_coords)
    prov.preload([srx_model.get(s) for s in samples if isinstance(srx_model.get(s), str)])
    raw = M.values.astype(float)
    cn = np.ones_like(raw)
    # CN is a property of the CELL LINE, not the experiment, so evaluate region_cn once per ModelID and
    # reuse it across that line's replicate SRX. Exactly equivalent (region_cn is pure, tracks are cached),
    # but at atlas scale it is the difference between ~14M and ~92M calls (2,136 samples / 324 models).
    cn_by_model, t0 = {}, time.time()
    for j, s in enumerate(samples):
        mid = srx_model.get(s)
        if not isinstance(mid, str):
            continue
        if mid not in cn_by_model:
            tr = prov.track(mid)
            cn_by_model[mid] = None if tr is None else np.fromiter(
                (tr.region_cn(*coords[i]) for i in se_ids), dtype=float, count=len(se_ids))
            if len(cn_by_model) % 25 == 0:
                print(f"[score]   CN tracks: {len(cn_by_model)} models, {time.time() - t0:.0f}s",
                      file=sys.stderr, flush=True)
        if cn_by_model[mid] is not None:
            cn[:, j] = cn_by_model[mid]
    n_cn = sum(v is not None for v in cn_by_model.values())
    print(f"[score] scoring-time CN: {n_cn}/{len(cn_by_model)} models with a DepMap track "
          f"({time.time() - t0:.0f}s)", file=sys.stderr, flush=True)
    # CN is always EVALUATED (it is reported per call as cn_mean, so any hit can be checked against the
    # copy number at its locus); --no-cn only skips APPLYING it. The two arms differ in exactly one step,
    # which is what makes the ablation interpretable.
    if a.no_cn:
        print("[score] scoring-time CN correction DISABLED (--no-cn: the ablation arm)",
              file=sys.stderr, flush=True)
        corrected = raw
    else:
        corrected = correct(raw, cn, model="log2offset", floor=0.1)
        # log2offset is (sig+eps)/cn - eps, which goes NEGATIVE when cn > sig+eps — i.e. a near-zero-signal SE
        # sitting in a high-CN region. Signal is non-negative by construction, and JSD is a divergence between
        # distributions, so a negative cell makes pycacts emit NaN (score.py propagates it by design). Clip here
        # rather than in cnrose.cn.correct, which is validated bit-for-bit against ROSE2 at calling time.
        n_neg = int((corrected < 0).sum())
        if n_neg:
            print(f"[score] CN correction produced {n_neg:,} negative cells "
                  f"({100.0 * n_neg / corrected.size:.4f}%, min {corrected.min():.3f}); clipping to 0",
                  file=sys.stderr, flush=True)
            corrected = np.maximum(corrected, 0.0)

    # normalise (batch) then collapse replicate SRX -> cell line (ModelID)
    col_model = pd.Series([srx_model.get(s) for s in samples], index=samples)
    def to_lines(mat):
        nm = pd.DataFrame(normalize(mat, a.norm, samples), index=se_ids, columns=samples)
        keep = col_model.dropna()
        return nm[keep.index].T.groupby(keep).mean().T             # SE x ModelID
    lines_cor = to_lines(corrected)          # the uncorrected collapse was computed but never read
    # per-line mean CN at each SE, collapsed the same way — reported as cn_mean so an amplicon-driven call
    # is visible in the output rather than needing a separate investigation. NOT normalized (it is a ratio).
    keep0 = col_model.dropna()
    cn_lines = pd.DataFrame(cn, index=se_ids, columns=samples)[keep0.index].T.groupby(keep0).mean().T

    # protein-coding annotation universe + nearest / identity-window helpers
    gidx = {}
    for g in prov.usable:
        gidx.setdefault(gene_coords[g][0], []).append(g)
    def nearest(chrom, mid):
        best, bd = None, None
        for g in gidx.get(chrom, ()):
            _, s, e = gene_coords[g]
            d = abs((s + e) // 2 - mid)
            if bd is None or d < bd:
                best, bd = g, d
        return best, (bd or 0)
    def identity_near(genes, chrom, s, e, w=100_000):
        for g in genes:
            p = gene_coords.get(g)
            if p and p[0] == chrom and not (p[2] < s - w or p[1] > e + w):
                return g
        return None

    print(f"[score] {len(se_ids)} SEs; {len(samples)} samples -> {lines_cor.shape[1]} cell lines; "
          f"norm={a.norm}\n")
    summary = []
    for level in LEVELS:
        rep, gsize = build_rep_matrix(lines_cor, model, level, min_group_n=1)
        rep.columns = [str(c) for c in rep.columns]
        jsd = cacts_score_matrix(rep)
        cn_rep, _ = build_rep_matrix(cn_lines, model, level, min_group_n=1)   # mean CN per group, same grouping
        cn_rep.columns = [str(c) for c in cn_rep.columns]
        # one FDR matrix per level; global BH shares the testing budget across groups (specificity.py)
        if a.fdr_method == "permutation" and level != "line":
            FDR = np.power(10.0, permutation_fdr(jsd, lines_cor, model, level, n_perm=a.n_perm,
                                                 scope=a.fdr_scope))
        else:
            if a.fdr_method == "permutation":
                print("  [perm] SKIPPING permutation at 'line' level (degenerate — permuting labels only "
                      "renames single-line groups); falling back to the analytic null.",
                      file=sys.stderr, flush=True)
            FDR = np.power(10.0, fdr_matrix(jsd, null=a.fdr_null, scope=a.fdr_scope))
        # label each group with cell name (line level) and its identity gene set
        rows = []
        print(f"================  {level}  ({rep.shape[1]} groups)  ================")
        for grp in rep.columns:
            disp = name_of.get(grp, grp) if level == "line" else grp
            genes = identity_for(model.loc[grp, "OncotreeSubtype"] if level == "line" else grp)
            s = jsd[grp].sort_values()
            fdr = FDR[grp]
            hits, best = [], {}
            for i, se in enumerate(s.index):
                c, ss, ee = coords[se]
                idg = identity_near(genes, c, ss, ee)
                if idg:
                    hits.append((i + 1, idg, se))
                    if idg not in best or i + 1 < best[idg]:
                        best[idg] = i + 1
            top15 = sum(1 for r, _, _ in hits if r <= 15)
            gtxt = ", ".join(f"{g}(#{best[g]})" for g in sorted(best, key=lambda g: best[g])[:5])
            n = int(gsize.get(grp, 1))
            n05 = int((fdr <= 0.05).sum()); n10 = int((fdr <= 0.10).sum()); n25 = int((fdr <= 0.25).sum())
            print(f"  {disp:<34} n={n}  specific(FDR≤.05/.10/.25)={n05:>4}/{n10:>4}/{n25:>4}  "
                  f"id@15={top15}  {gtxt or '—'}")
            rows.append(dict(level=level, group=disp, n_lines=n, n_spec_fdr05=n05, n_spec_fdr10=n10,
                             n_spec_fdr25=n25, top15_identity=top15,
                             genes=";".join(f"{g}:{best[g]}" for g in best)))
        summary.extend(rows)
        # write per-level top-specific SEs
        recs = []
        for grp in rep.columns:
            disp = name_of.get(grp, grp) if level == "line" else grp
            genes = identity_for(model.loc[grp, "OncotreeSubtype"] if level == "line" else grp)
            tbl = rank_specific(jsd, grp).head(a.topn)
            fdr = FDR[grp]
            for _, r in tbl.iterrows():
                se = r["tf"]; c, ss, ee = coords[se]
                gene, dist = nearest(c, (ss + ee) // 2)
                idg = identity_near(genes, c, ss, ee)
                recs.append(dict(level=level, group=disp, se=se, rank=int(r["rank"]),
                                 jsd=r["cacts_score"], fdr=float(fdr[se]), nearest_gene=gene,
                                 dist_kb=dist // 1000, identity=idg or "",
                                 cn_mean=round(float(cn_rep[grp].get(se, float("nan"))), 3)))
        pd.DataFrame(recs).to_csv(f"{a.out}.{level}.top_specific.tsv", sep="\t", index=False)

        if a.dump_specific is not None:
            # every SE x group at or below the cutoff, long-format. Rank is within-group by ascending JSD,
            # matching rank_specific, so this is a strict superset of the top_specific table.
            order = np.argsort(jsd.values, axis=0, kind="stable")
            rank_of = np.empty_like(order)
            np.put_along_axis(rank_of, order, np.arange(1, jsd.shape[0] + 1)[:, None].repeat(jsd.shape[1], 1), axis=0)
            keep = np.argwhere(FDR.values <= a.dump_specific)
            se_arr, grp_arr = np.array(se_ids), np.array(rep.columns, dtype=object)
            out = pd.DataFrame({
                "level": level,
                "group": [name_of.get(g, g) if level == "line" else g for g in grp_arr[keep[:, 1]]],
                "se": se_arr[keep[:, 0]],
                "rank": rank_of[keep[:, 0], keep[:, 1]],
                "jsd": jsd.values[keep[:, 0], keep[:, 1]],
                "fdr": FDR.values[keep[:, 0], keep[:, 1]],
                "cn_mean": np.round(cn_rep.values[keep[:, 0], keep[:, 1]], 3),
            }).sort_values(["group", "rank"])
            path = f"{a.out}.{level}.specific.tsv.gz"
            out.to_csv(path, sep="\t", index=False, compression="gzip")
            print(f"  [dump] {len(out):,} SE x group tests at FDR<={a.dump_specific} -> {path}",
                  file=sys.stderr, flush=True)
        print()

    pd.DataFrame(summary).to_csv(a.out + ".hierarchy_summary.tsv", sep="\t", index=False)
    print(f"[score] wrote {a.out}.<level>.top_specific.tsv + .hierarchy_summary.tsv", file=sys.stderr)


if __name__ == "__main__":
    main()
