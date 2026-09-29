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
from permutation import permutation_fdr, permutation_fdr_units, rep_agg   # noqa: E402
from cnrose.cn.depmap import load_gene_coords, DepMapGeneCN, DepMapMcWesCN   # noqa: E402
from cnrose.cn.cmp import CellModelPassportsWesCN                  # noqa: E402
from cnrose.cn.inferred import BinnedInputCN, load_blacklist       # noqa: E402
from cnrose.cn.segfile import SegmentFileCN                        # noqa: E402
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
    ap.add_argument("--cmp-wes", default=os.path.join(DATAROOT, "CellModelPassports/WES_pureCN_CNV_genes_latest.csv.gz"),
                    help="CMP WES pureCN 2025 — CN for pull-set rows with cn_provider=cmp_wes")
    ap.add_argument("--cmp-model-list", default=os.path.join(DATAROOT, "CellModelPassports/model_list_20240110.csv"))
    ap.add_argument("--mc-wes", default=os.path.join(DATAROOT, "DepMap/2026q1/OmicsCNGeneMC_WES.csv"),
                    help="DepMap MC_WES — CN for pull-set rows with cn_provider=depmap_mc_wes")
    ap.add_argument("--model-condition", default=os.path.join(DATAROOT, "DepMap/2026q1/ModelCondition.csv"))
    ap.add_argument("--ccle-segments", default=os.path.join(SECACTS, "phase2/data/cn_ccle_snp6.hg38.tsv.gz"),
                    help="CCLE 2019 SNP6 lifted to hg38 — CN for pull-set rows with cn_provider=ccle_snp6 (FINDINGS §23)")
    ap.add_argument("--inferred-bins", help="dir of <SRX>.npz input bins (phase1 script 20) — CN for rows with "
                    "cn_provider=input_inferred (FINDINGS §24)")
    ap.add_argument("--inferred-map", default=os.path.join(SECACTS, "phase2/analysis/out/cn_admissions.tsv"),
                    help="TSV with key and cn_input (the admitted input SRX per line)")
    ap.add_argument("--inferred-slope", type=float, default=0.81)
    ap.add_argument("--blacklist", help="ENCODE hg38 blacklist v2 BED(.gz), for input_inferred")
    ap.add_argument("--lines-meta", default=os.path.join(SECACTS, "phase1/data/lineage_resolved.tsv"),
                    help="Oncotree labels (by CVCL) for lines with no DepMap ModelID; their pull-set key is the CVCL")
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
    ap.add_argument("--levels", default=",".join(LEVELS),
                    help=f"comma-separated hierarchy levels to score (default all: {','.join(LEVELS)})")
    ap.add_argument("--shuffle-labels", type=int, default=None, metavar="SEED",
                    help="CALIBRATION TEST: permute the Oncotree labels across the scored lines (one shared "
                         "permutation, so the lineage/disease/subtype nesting and every group size are kept) "
                         "before scoring. Nothing real then exists, so a calibrated FDR must call ~nothing.")
    ap.add_argument("--exclude-keys", default="",
                    help="comma-separated line keys (ModelID / CVCL) to drop before scoring — sensitivity runs")
    ap.add_argument("--cn-diagnostic", action="store_true",
                    help="write <out>.cn_by_line.tsv: per line, its CN source and the Spearman of SE signal vs "
                         "SE copy number before and after correction. Under-correction (a noisier CN source) "
                         "shows as a residual positive rho — the lineage-varying-strength confound.")
    ap.add_argument("--save-matrices", action="store_true",
                    help="also write the full SE x group JSD and FDR matrices per group level (not 'line') to "
                         "<out>.<level>.{jsd,fdr}.tsv.gz — for figures that need values beyond the call dumps")
    ap.add_argument("--agg", default="mean", choices=["mean", "q25", "min"],
                    help="how a group combines its members before JSD. 'mean' (default) is the original score; "
                         "'q25' scores the group CONSENSUS — the value >=75%% of members reach (the minimum for "
                         "n <= 4) — so one member's private SE cannot make a small group look specific. Use it "
                         "at the small-group levels (subtype, line); see permutation.py.")
    ap.add_argument("--line-members", default="lines", choices=["lines", "experiments", "studies"],
                    help="members of a cell-line group: 'lines' (one profile, permutation degenerate) or "
                         "'experiments' (the line's replicate experiments; the null permutes experiments "
                         "between lines, so a call means specific AND reproducible across its experiments), "
                         "or 'studies' (experiments averaged within each study first, so a call must replicate "
                         "across independent studies and a study's shared batch effects cannot pose as one)")
    ap.add_argument("--srx-study", default=os.path.join(SECACTS, "phase2/data/srx_study.tsv"),
                    help="experiment -> study accession (for --line-members studies)")
    ap.add_argument("--null-strata", default="none", choices=["none", "lineage", "disease", "adaptive"],
                    help="consensus levels: shuffle labels only among relatives, so a call means specific "
                         "relative to the group's lineage (or disease) rather than to the whole panel. "
                         "'adaptive' uses the primary disease when it holds >= --strata-min lines, else the "
                         "lineage. Writes <out>.<level>.strata.tsv (the comparison set per group).")
    ap.add_argument("--strata-min", type=int, default=4)
    ap.add_argument("--min-members", type=int, default=2,
                    help="with a consensus --agg, groups with fewer members are not called (FDR set to 1): one "
                         "member has no consensus, and its 'specificity' measures how unlike a random line it is "
                         "(highest for tumour types with no relatives in the panel), not a group property")
    ap.add_argument("--keep-frac", type=float, default=0.05,
                    help="share of the permutation null kept per group (only the left tail is ever read)")
    a = ap.parse_args()
    levels = [l for l in a.levels.split(",") if l]
    bad = [l for l in levels if l not in LEVELS]
    if bad:
        ap.error(f"unknown level(s) {bad}; choose from {LEVELS}")

    M = pd.read_csv(a.signal, sep="\t", index_col=0)               # SE x SRX (uncorrected)
    coords = load_coords(a.catalog)
    samples = list(M.columns)
    se_ids = list(M.index)

    ps = pd.read_csv(a.pull_set, sep="\t")
    # v2 pull sets carry `key` (ModelID, or CVCL for a line DepMap does not hold) and `cn_provider`;
    # a v1 pull set is all DepMap WGS, keyed by ModelID.
    if "key" not in ps.columns:
        ps["key"], ps["cn_provider"] = ps["model_id"], "depmap_wgs"
    srx_model = dict(zip(ps["srx"], ps["key"]))
    srx_model.update(EXTRA_MODEL)
    excl = {k for k in a.exclude_keys.split(",") if k}
    if excl:
        drop = [c for c in M.columns if srx_model.get(c) in excl]
        M = M.drop(columns=drop)
        print(f"[score] --exclude-keys: dropped {len(drop)} samples from "
              f"{len({srx_model[c] for c in drop})}/{len(excl)} requested lines", file=sys.stderr, flush=True)
        samples = list(M.columns)
    # An SE with zero signal in EVERY sample is not a test: it has no quantified territory (v2: union loci
    # that fall wholly outside the fixed grid). Left in, pycacts scores the all-zero row as a uniform
    # profile (JSD 0.347) that clears the permutation null in every group — 2,108 spurious lineage calls in
    # the uncorrected v2 arm — and scoring-time correction turns its zeros into CN-derived pseudo-signal.
    zero = (M.values == 0).all(axis=1)
    if zero.any():
        M = M.loc[~zero]
        print(f"[score] dropped {int(zero.sum())} SE loci with zero signal in every sample (not tests)",
              file=sys.stderr, flush=True)
    se_ids = list(M.index)
    model = pd.read_csv(a.model, index_col="ModelID")
    # lines outside DepMap get their Oncotree labels from the Cellosaurus-NCIt crosswalk (phase1 script 14)
    extra = sorted(set(k for k in ps["key"].dropna() if k not in model.index))
    if extra:
        lm = pd.read_csv(a.lines_meta, sep="\t").set_index("cvcl").loc[extra]
        add = pd.DataFrame({"StrippedCellLineName": lm["cell_line"], "OncotreeLineage": lm["lineage"],
                            "OncotreePrimaryDisease": lm["primary_disease"], "OncotreeSubtype": lm["subtype"]},
                           index=pd.Index(extra, name="ModelID"))
        model = pd.concat([model, add])
        print(f"[score] {len(extra)} line(s) outside DepMap, labelled from {os.path.basename(a.lines_meta)}",
              file=sys.stderr, flush=True)
    name_of = model["StrippedCellLineName"].to_dict()

    # scoring-time CN over each SE region, per sample (symmetric per-copy, floor 0.1). Each line reads the
    # source it was admitted on (phase1 script 15): DepMap WGS -> CMP WES pureCN 2025 -> DepMap MC_WES.
    gene_coords = load_gene_coords(a.gtf, cache_path=a.gene_cache)
    prov = DepMapGeneCN(a.cn_gene_csv, gene_coords)
    lines_ps = ps.dropna(subset=["key"]).drop_duplicates("key").set_index("key")
    src_of = lines_ps["cn_provider"].to_dict()
    # a line scored on a relative's CN names that relative in `cn_cvcl` (phase1 script 23); CMP keys by CVCL
    cvcl_of = (lines_ps["cn_cvcl"].fillna(lines_ps["cvcl"]) if "cn_cvcl" in lines_ps else lines_ps["cvcl"]).to_dict()
    used = {srx_model.get(s) for s in samples if isinstance(srx_model.get(s), str)}
    by_src = {p: [k for k in used if src_of.get(k, "depmap_wgs") == p]
              for p in ("depmap_wgs", "cmp_wes", "depmap_mc_wes", "ccle_snp6", "input_inferred")}
    prov.preload(by_src["depmap_wgs"])
    cmp_prov = mcw_prov = None
    if by_src["cmp_wes"]:
        cmp_prov = CellModelPassportsWesCN(a.cmp_wes, a.cmp_model_list, cache_dir=cache_path("cmp_wes"))
        cmp_prov.preload([cvcl_of[k] for k in by_src["cmp_wes"]])
    if by_src["depmap_mc_wes"]:
        mcw_prov = DepMapMcWesCN(a.mc_wes, a.model_condition, gene_coords)
        mcw_prov.preload(by_src["depmap_mc_wes"])
    ccle_prov = SegmentFileCN(a.ccle_segments, "ccle_snp6") if by_src["ccle_snp6"] else None
    inf_prov = None
    if by_src["input_inferred"]:
        if not (a.inferred_bins and a.blacklist):
            sys.exit("[score] input_inferred lines need --inferred-bins and --blacklist")
        im = pd.read_csv(a.inferred_map, sep="\t").dropna(subset=["cn_input"])
        inf_prov = BinnedInputCN(a.inferred_bins, dict(zip(im["key"], im["cn_input"])),
                                 blacklist=load_blacklist(a.blacklist), slope=a.inferred_slope)

    def track_for(key):
        src = src_of.get(key, "depmap_wgs")
        if src == "ccle_snp6":
            return ccle_prov.track(key)
        if src == "input_inferred":
            return inf_prov.track(key)
        if src == "cmp_wes":
            return cmp_prov.track(cvcl_of[key])
        if src == "depmap_mc_wes":
            return mcw_prov.track(key)
        return prov.track(key)
    print("[score] CN source per line: " + ", ".join(f"{p}={len(v)}" for p, v in by_src.items()),
          file=sys.stderr, flush=True)
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
            tr = track_for(mid)
            cn_by_model[mid] = None if tr is None else np.fromiter(
                (tr.region_cn(*coords[i]) for i in se_ids), dtype=float, count=len(se_ids))
            if len(cn_by_model) % 25 == 0:
                print(f"[score]   CN tracks: {len(cn_by_model)} models, {time.time() - t0:.0f}s",
                      file=sys.stderr, flush=True)
        if cn_by_model[mid] is not None:
            cn[:, j] = cn_by_model[mid]
    n_cn = sum(v is not None for v in cn_by_model.values())
    print(f"[score] scoring-time CN: {n_cn}/{len(cn_by_model)} models with a CN track "
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
    nm_samples, unit_line = None, None
    if a.line_members in ("experiments", "studies"):
        # the per-experiment profiles, normalised exactly as to_lines does before it collapses them
        keep_ = col_model.dropna()
        nm_samples = pd.DataFrame(normalize(corrected, a.norm, samples), index=se_ids, columns=samples)[keep_.index]
        if a.line_members == "studies":
            # one member per (line, study): a study's experiments are averaged, so a line needs >= 2 studies
            # to be testable and within-study batch effects count once
            st = pd.read_csv(a.srx_study, sep="\t").set_index("srx")["study"]
            unit = pd.Series([f"{col_model[c]}|{st.get(c, c)}" for c in nm_samples.columns], index=nm_samples.columns)
            nm_samples = nm_samples.T.groupby(unit).mean().T
            unit_line = pd.Series({u: u.split("|")[0] for u in nm_samples.columns})
            print(f"[score] line members = studies: {nm_samples.shape[1]} (line, study) units over "
                  f"{len(set(u.split('|')[0] for u in nm_samples.columns))} lines", file=sys.stderr, flush=True)
    def to_lines(mat):
        nm = pd.DataFrame(normalize(mat, a.norm, samples), index=se_ids, columns=samples)
        keep = col_model.dropna()
        return nm[keep.index].T.groupby(keep).mean().T             # SE x ModelID
    lines_cor = to_lines(corrected)          # the uncorrected collapse was computed but never read
    # per-line mean CN at each SE, collapsed the same way — reported as cn_mean so an amplicon-driven call
    # is visible in the output rather than needing a separate investigation. NOT normalized (it is a ratio).
    keep0 = col_model.dropna()
    cn_lines = pd.DataFrame(cn, index=se_ids, columns=samples)[keep0.index].T.groupby(keep0).mean().T
    if a.cn_diagnostic:
        from scipy.stats import spearmanr
        raw_lines = pd.DataFrame(raw, index=se_ids, columns=samples)[keep0.index].T.groupby(keep0).mean().T
        diag = []
        for k in lines_cor.columns:
            c = cn_lines[k].values
            if np.allclose(c, 1.0):
                continue
            diag.append(dict(key=k, line=name_of.get(k, k), cn_source=src_of.get(k, "depmap_wgs"),
                             n_samples=int((keep0 == k).sum()), frac_amp=round(float((c > 1.3).mean()), 4),
                             rho_raw=round(float(spearmanr(raw_lines[k].values, c)[0]), 4),
                             rho_corrected=round(float(spearmanr(lines_cor[k].values, c)[0]), 4)))
        pd.DataFrame(diag).to_csv(f"{a.out}.cn_by_line.tsv", sep="\t", index=False)
        del raw_lines
        print(f"[score] wrote {a.out}.cn_by_line.tsv ({len(diag)} lines)", file=sys.stderr, flush=True)

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

    if a.shuffle_labels is not None:
        # one permutation of whole label ROWS among the scored lines: group sizes and the hierarchy's
        # nesting survive exactly, only the line -> group assignment is destroyed
        cols = ["OncotreeLineage", "OncotreePrimaryDisease", "OncotreeSubtype"]
        scored = [k for k in lines_cor.columns if k in model.index]
        perm = np.random.default_rng(a.shuffle_labels).permutation(len(scored))
        model = model.copy()
        model.loc[scored, cols] = model.loc[scored, cols].values[perm]
        print(f"[score] CALIBRATION: Oncotree labels SHUFFLED across {len(scored)} lines "
              f"(seed {a.shuffle_labels})", file=sys.stderr, flush=True)
    print(f"[score] {len(se_ids)} SEs; {len(samples)} samples -> {lines_cor.shape[1]} cell lines; "
          f"norm={a.norm}\n")
    summary = []
    for level in levels:
        exp_units = level == "line" and a.line_members in ("experiments", "studies")
        use_units = exp_units or a.agg != "mean"
        if use_units:
            # generic unit -> group aggregation (permutation.rep_agg): units are cell lines, or, for a
            # cell-line group with --line-members experiments, that line's experiments
            if exp_units:
                Xdf = nm_samples
                lab = (unit_line if unit_line is not None else col_model).loc[nm_samples.columns].values
                if a.shuffle_labels is not None:        # calibration: experiments reassigned to random lines
                    lab = np.random.default_rng(a.shuffle_labels).permutation(lab)
                    print(f"[score] CALIBRATION: experiment -> line labels SHUFFLED (seed {a.shuffle_labels})",
                          file=sys.stderr, flush=True)
            else:
                ok = [k for k in lines_cor.columns if k in model.index]
                labs = (pd.Series(ok, index=ok) if level == "line" else model.loc[ok, level]).dropna()
                Xdf, lab = lines_cor[labs.index], labs.values.astype(str)
            groups = sorted(set(lab))
            codes = np.array([groups.index(x) for x in lab]) if len(groups) < 50 else \
                pd.Categorical(lab, categories=groups).codes
            Xv = Xdf.values.astype(np.float32)
            rep = pd.DataFrame(rep_agg(Xv, np.asarray(codes), len(groups), a.agg), index=se_ids, columns=groups)
            gsize = pd.Series(lab).value_counts().reindex(groups)
            print(f"[score] {level}: {len(groups)} groups aggregated by {a.agg} over "
                  f"{a.line_members if exp_units else 'lines'}", file=sys.stderr, flush=True)
        else:
            rep, gsize = build_rep_matrix(lines_cor, model, level, min_group_n=1)
        rep.columns = [str(c) for c in rep.columns]
        jsd = cacts_score_matrix(rep)
        cn_rep, _ = build_rep_matrix(cn_lines, model, level, min_group_n=1)   # mean CN per group, same grouping
        cn_rep.columns = [str(c) for c in cn_rep.columns]
        # one FDR matrix per level; global BH shares the testing budget across groups (specificity.py)
        if a.fdr_method == "permutation" and use_units and (level != "line" or exp_units):
            strata = None
            if a.null_strata != "none":
                # each unit's line, then that line's parent group; stratum sizes counted in lines
                unit_key = (unit_line if (exp_units and unit_line is not None) else
                            (col_model if exp_units else pd.Series(Xdf.columns, index=Xdf.columns))).loc[Xdf.columns]
                lin = model["OncotreeLineage"].reindex(unit_key.values).fillna("NA").values
                dis = model["OncotreePrimaryDisease"].reindex(unit_key.values).fillna("NA").values
                if a.null_strata == "lineage":
                    strata = lin
                elif a.null_strata == "disease":
                    strata = dis
                else:
                    per_line = pd.DataFrame({"k": unit_key.values, "d": dis}).drop_duplicates("k")
                    n_d = per_line.d.value_counts()
                    strata = np.where(pd.Series(dis).map(n_d).fillna(0).values >= a.strata_min, dis, lin)
                sdf = pd.DataFrame({"group": lab, "stratum": strata, "line": unit_key.values})
                n_lines = sdf.drop_duplicates("line").groupby("stratum").size()
                out_s = sdf.drop_duplicates("group")[["group", "stratum"]].assign(
                    lines_in_stratum=lambda d: d.stratum.map(n_lines).values)
                out_s["group"] = [name_of.get(g, g) if level == "line" else g for g in out_s.group]
                out_s.to_csv(f"{a.out}.{level}.strata.tsv", sep="\t", index=False)
                print(f"[score] {level}: null shuffles within {len(n_lines)} {a.null_strata} strata",
                      file=sys.stderr, flush=True)
            FDR = np.power(10.0, permutation_fdr_units(jsd, Xv, lab, a.agg, n_perm=a.n_perm,
                                                       keep_frac=a.keep_frac, scope=a.fdr_scope,
                                                       strata=strata))
            if a.agg != "mean":
                small = [g for g in FDR.columns if int(gsize.get(g, 0)) < a.min_members]
                FDR[small] = 1.0
                print(f"[score] {len(small)} {level} groups with < {a.min_members} members left uncalled "
                      f"(rankings only)", file=sys.stderr, flush=True)
        elif a.fdr_method == "permutation" and level != "line":
            FDR = np.power(10.0, permutation_fdr(jsd, lines_cor, model, level, n_perm=a.n_perm,
                                                 scope=a.fdr_scope))
        else:
            if a.fdr_method == "permutation":
                print("  [perm] SKIPPING permutation at 'line' level (degenerate — permuting labels only "
                      "renames single-line groups); falling back to the analytic null.",
                      file=sys.stderr, flush=True)
            FDR = np.power(10.0, fdr_matrix(jsd, null=a.fdr_null, scope=a.fdr_scope))
        if a.save_matrices and level != "line":
            jsd.round(5).to_csv(f"{a.out}.{level}.jsd.tsv.gz", sep="\t", compression="gzip")
            FDR.to_csv(f"{a.out}.{level}.fdr.tsv.gz", sep="\t", compression="gzip", float_format="%.4g")
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
