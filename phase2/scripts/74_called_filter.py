#!/usr/bin/env python3
"""A derived copy of a scores directory in which a specific super-enhancer must be a super-enhancer of its group.

Scoring measures H3K27ac signal at every union-catalogue locus (an SE in ANY experiment of the panel) in every
line, so a locus can be "specific to VCaP" because VCaP has the most signal there even though no VCaP experiment
calls an SE at it (FINDINGS §28, §39: 30% of v3 lineage calls). This keeps a call only when at least one
experiment of the group (lineage, primary disease, subtype, or the line itself) calls an SE that overlaps the
locus: the group's member loci in atlas.s3.se_presence, widened to every catalogue locus overlapping one of them
(the catalogue nests small loci inside large ones).

FDRs are kept as computed over every catalogue locus. Ranks are renumbered 1..n within each group, which is
exact: FDR is monotone in JSD within a group (ranks are contiguous in every group of v3), so the kept calls are
the top of the group's called loci. The JSD/FDR matrices are masked (NaN) where the locus is not an SE of the
group. Shuffled-label arms are copied unchanged (calibration of the test, not the definition; 0 calls).

    python3 phase2/scripts/74_called_filter.py --scores phase2/scores_v3 --results phase2/results_v3 \\
        --pull-set phase2/data/pull_set.v3.tsv --out phase2/scores_v3c
"""
import argparse, os, shutil, sys
import numpy as np
import pandas as pd
import scipy.sparse as sp

SECACTS = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
sys.path.insert(0, SECACTS)
from secacts_env import DATAROOT                                         # noqa: E402

LEVELS = ["OncotreeLineage", "OncotreePrimaryDisease", "OncotreeSubtype", "line"]
EXTRA_MODEL = {"SRX10809652": "ACH-000768"}                             # as score_pilot.py
MAIN_ARMS = ["atlas.s3.perm", "atlas.s3.perm.nocn", "atlas.s3.perm.noinf"]
LINE_ARMS = ["atlas.s3.lines.all", "atlas.s3.lines.all.nocn", "atlas.s3.lines.sub", "atlas.s3.lines.dis",
             "atlas.s3.lines.lin"]


def labels(ps, model_csv, lines_meta):
    """Per line key: the four group names score_pilot.py uses (Model.csv, else the Cellosaurus crosswalk)."""
    model = pd.read_csv(model_csv, index_col="ModelID")
    extra = sorted(set(k for k in ps["key"].dropna() if k not in model.index))
    if extra:
        lm = pd.read_csv(lines_meta, sep="\t").set_index("cvcl").loc[extra]
        model = pd.concat([model, pd.DataFrame({
            "StrippedCellLineName": lm["cell_line"], "OncotreeLineage": lm["lineage"],
            "OncotreePrimaryDisease": lm["primary_disease"], "OncotreeSubtype": lm["subtype"]},
            index=pd.Index(extra, name="ModelID"))])
    keys = sorted(ps["key"].dropna().unique())
    L = model.loc[keys, ["StrippedCellLineName", "OncotreeLineage", "OncotreePrimaryDisease", "OncotreeSubtype"]]
    return L.rename(columns={"StrippedCellLineName": "line"})


def overlap_graph(cat):
    """Sparse adjacency of catalogue loci that overlap each other."""
    idx = {s: i for i, s in enumerate(cat.se)}
    rows, cols = [], []
    for _, d in cat.groupby("chrom"):
        d = d.sort_values("start")
        st, en, ii = d.start.values, d.end.values, d.se.map(idx).values
        for k in range(len(d)):
            j = k + 1
            while j < len(d) and st[j] < en[k]:
                if min(en[k], en[j]) > max(st[k], st[j]):
                    rows += [ii[k], ii[j]]; cols += [ii[j], ii[k]]
                j += 1
    n = len(cat)
    return sp.csr_matrix((np.ones(len(rows), dtype=np.float32), (rows, cols)), shape=(n, n))


def coverage(P, A, srx_group, groups):
    """bool[locus, group]: an SE called by an experiment of the group overlaps the locus."""
    col = {g: j for j, g in enumerate(groups)}
    r = [i for i, g in enumerate(srx_group) if g in col]
    E = sp.csr_matrix((np.ones(len(r), dtype=np.float32), (r, [col[srx_group[i]] for i in r])),
                      shape=(P.shape[1], len(groups)))
    M = (P @ E) > 0
    return ((M + (A @ M)) > 0).toarray()


def filter_calls(df, cov, se_idx, gcol):
    if df.empty:
        return df, 0
    i = df["se"].map(se_idx)
    j = df["group"].map(gcol)
    ok = i.notna() & j.notna()
    keep = np.zeros(len(df), bool)
    keep[ok.values] = cov[i[ok].astype(int).values, j[ok].astype(int).values]
    out = df[keep].sort_values(["group", "rank"]).copy()
    out["rank"] = out.groupby("group").cumcount() + 1
    return out, int((~ok).sum())


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--scores", required=True)
    ap.add_argument("--results", required=True)
    ap.add_argument("--pull-set", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--model", default=os.path.join(DATAROOT, "DepMap/2026q1/Model.csv"))
    ap.add_argument("--lines-meta", default=os.path.join(SECACTS, "phase1/data/lineage_resolved.tsv"))
    a = ap.parse_args()
    os.makedirs(os.path.join(a.out, "out_lines"), exist_ok=True)

    ps = pd.read_csv(a.pull_set, sep="\t")
    lab = labels(ps, a.model, a.lines_meta)
    cat = pd.read_csv(os.path.join(a.results, "atlas.s3.union_catalog.bed.gz"), sep="\t", header=None,
                      usecols=[0, 1, 2, 3], names=["chrom", "start", "end", "se"])
    se_idx = {s: i for i, s in enumerate(cat.se)}
    pres = pd.read_csv(os.path.join(a.results, "atlas.s3.se_presence.tsv.gz"), sep="\t", index_col=0)
    pres = pres.reindex(cat.se).fillna(0)
    P = sp.csr_matrix((pres.values > 0).astype(np.float32))
    A = overlap_graph(cat)
    srx_key = dict(zip(ps["srx"], ps["key"])); srx_key.update(EXTRA_MODEL)
    keys = [srx_key.get(s) for s in pres.columns]
    print(f"[74] {len(cat):,} loci, {P.shape[1]} experiments ({sum(k is None for k in keys)} without a line), "
          f"{A.nnz // 2:,} overlapping locus pairs", file=sys.stderr)
    inferred = set(ps.loc[ps.cn_provider == "input_inferred", "key"])

    def cov_for(level, arm):
        excl = inferred if ".noinf" in arm else set()
        srx_group = [lab.at[k, "line" if level == "line" else level] if (k in lab.index and k not in excl) else None
                     for k in keys]
        groups = sorted({g for g in srx_group if isinstance(g, str)})
        return coverage(P, A, srx_group, groups), {g: j for j, g in enumerate(groups)}

    cache = {}

    def cov_cached(level, arm):
        k = (level, ".noinf" in arm)
        if k not in cache:
            cache[k] = cov_for(level, arm)
        return cache[k]

    def do_arm(src_dir, dst_dir, arm, levels):
        H = pd.read_csv(os.path.join(src_dir, f"{arm}.hierarchy_summary.tsv"), sep="\t")
        H = H.drop(columns=[c for c in ("top15_identity", "genes") if c in H.columns])
        for lev in levels:
            f = os.path.join(src_dir, f"{arm}.{lev}.specific.tsv.gz")
            if not os.path.exists(f):
                continue
            cov, gcol = cov_cached(lev, arm)
            S = pd.read_csv(f, sep="\t")
            F, unk = filter_calls(S, cov, se_idx, gcol)
            F.to_csv(os.path.join(dst_dir, f"{arm}.{lev}.specific.tsv.gz"), sep="\t", index=False)
            T = pd.read_csv(os.path.join(src_dir, f"{arm}.{lev}.top_specific.tsv"), sep="\t")
            TF, _ = filter_calls(T, cov, se_idx, gcol)
            TF.to_csv(os.path.join(dst_dir, f"{arm}.{lev}.top_specific.tsv"), sep="\t", index=False)
            n05 = F[F.fdr <= 0.05].groupby("group").size(); n10 = F.groupby("group").size()
            m = H.level == lev
            H.loc[m, "n_spec_fdr05"] = H.loc[m, "group"].map(n05).fillna(0).astype(int).values
            H.loc[m, "n_spec_fdr10"] = H.loc[m, "group"].map(n10).fillna(0).astype(int).values
            H.loc[m, "n_spec_fdr25"] = np.nan                                  # not recoverable from <= 0.10 rows
            print(f"[74] {arm} {lev}: {len(S):,} -> {len(F):,} calls ({len(F) / max(len(S), 1):.0%})"
                  + (f"; {unk} rows with an unknown group or locus dropped" if unk else ""), file=sys.stderr)
            sp_ = os.path.join(src_dir, f"{arm}.{lev}.strata.tsv")
            if os.path.exists(sp_):
                shutil.copy2(sp_, os.path.join(dst_dir, os.path.basename(sp_)))
        H.to_csv(os.path.join(dst_dir, f"{arm}.hierarchy_summary.tsv"), sep="\t", index=False)

    for arm in MAIN_ARMS:
        if os.path.exists(os.path.join(a.scores, f"{arm}.hierarchy_summary.tsv")):
            do_arm(a.scores, a.out, arm, LEVELS)
    has_lines = os.path.isdir(os.path.join(a.scores, "out_lines"))                    # v2 has no per-line arms
    for arm in LINE_ARMS if has_lines else []:
        do_arm(os.path.join(a.scores, "out_lines"), os.path.join(a.out, "out_lines"), arm, ["line"])
    # the main arm's full matrices: NaN where the locus is not an SE of the group
    for lev in LEVELS[:3]:
        cov, gcol = cov_cached(lev, "atlas.s3.perm")
        for kind in ("jsd", "fdr"):
            f = os.path.join(a.scores, f"atlas.s3.perm.{lev}.{kind}.tsv.gz")
            if not os.path.exists(f):
                continue
            X = pd.read_csv(f, sep="\t", index_col=0)
            ri = X.index.map(se_idx).to_numpy(dtype=float)
            for c in X.columns:
                if c in gcol:
                    ok = ~np.isnan(ri)
                    mask = np.zeros(len(X), bool)
                    mask[ok] = cov[ri[ok].astype(int), gcol[c]]
                    X.loc[~mask, c] = np.nan
                else:
                    X[c] = np.nan
            X.to_csv(os.path.join(a.out, os.path.basename(f)), sep="\t")
    # unchanged: shuffled-label arms (0 calls; calibration of the test) and per-line CN
    for f in os.listdir(a.scores):
        if (".shuffle." in f and f.startswith(("atlas.s3.perm", "atlas.s3.analytic"))) or f == "atlas.s3.perm.cn_by_line.tsv":
            shutil.copy2(os.path.join(a.scores, f), os.path.join(a.out, f))
    for f in os.listdir(os.path.join(a.scores, "out_lines")) if has_lines else []:
        if ".shuffle." in f:
            shutil.copy2(os.path.join(a.scores, "out_lines", f), os.path.join(a.out, "out_lines", f))
    print(f"[74] -> {a.out}; now run the derive step (cn_ablation_calls, concordance_bridge2) on it", file=sys.stderr)


if __name__ == "__main__":
    main()
