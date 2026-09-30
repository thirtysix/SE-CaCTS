#!/usr/bin/env python3
"""Candidate lineage master TFs supported by independent layers (a quick grab for novel biology).

For every lineage and every TF (pyCaCTS's curated 1,671):
  SE          a lineage-specific super-enhancer (v3.0.1: FDR <= 0.10 and an SE of the lineage) within 100 kb
  expression  pyCaCTS calls the TF lineage-specific in DepMap expression (category "specific", FDR <= 0.10)
  dependency  DepMap CRISPR (Chronos): lines of the lineage depend on it (mean <= -0.5) more than other lines
              (difference <= -0.3), with at least 3 screened lines in the lineage
A TF with SE + expression is a two-layer candidate; + dependency, three. The 22 pre-registered master regulators
(poster Fig. 4) are marked so the rest can be read as candidates; novelty still needs a literature check.

    ~/miniconda3/envs/atac_hdac/bin/python phase2/analysis/novel_tf_candidates.py --scores phase2/scores_v3c
"""
import argparse, os, re, ast, sys
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
SECACTS = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, SECACTS); sys.path.insert(0, os.path.join(SECACTS, "cnrose"))
from secacts_env import DATAROOT, cache_path                                       # noqa: E402
from cnrose.cn.depmap import load_gene_coords                                       # noqa: E402

PYCACTS = os.path.join(os.path.dirname(SECACTS), "pyCaCTS")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--scores", default=os.path.join(SECACTS, "phase2/scores_v3c"))
    ap.add_argument("--catalog", default=os.path.join(SECACTS, "phase2/results_v3/atlas.s3.union_catalog.bed.gz"))
    ap.add_argument("--window", type=int, default=100_000)
    ap.add_argument("--out", default=os.path.join(SECACTS, "phase2/analysis/out/novel_tf_candidates.tsv"))
    a = ap.parse_args()

    tfs = [t.strip() for t in open(os.path.join(PYCACTS, "data/CaCTS_merged_1671_TFs.txt")) if t.strip() and t.strip() != "NameTF"]
    src = open(os.path.join(SECACTS, "phase2/figures/poster_figures.py")).read()
    prereg = {g for _, gs in ast.literal_eval(re.search(r"IDENTITY = (\[.*?\])\n", src, re.S).group(1)) for g in gs}

    # SE layer: best-ranked lineage-specific SE within the window of each TF
    gc = load_gene_coords(None, cache_path=cache_path("gene_coords.GRCh38.106.tsv"))
    cat = pd.read_csv(a.catalog, sep="\t", header=None, usecols=[0, 1, 2, 3], names=["chrom", "start", "end", "se"])
    S = pd.read_csv(f"{a.scores}/atlas.s3.perm.OncotreeLineage.specific.tsv.gz", sep="\t").merge(cat, on="se")
    se_hits = []
    for tf in tfs:
        if tf not in gc:
            continue
        c, s, e = gc[tf]
        d = S[(S.chrom == c) & (S.end >= s - a.window) & (S.start <= e + a.window)]
        for g, x in d.groupby("group"):
            b = x.sort_values("rank").iloc[0]
            dist = max(0, max(b.start, s) - min(b.end, e))
            se_hits.append((g, tf, int(b["rank"]), float(b.fdr), int(dist // 1000), len(x)))
    SE = pd.DataFrame(se_hits, columns=["lineage", "tf", "se_rank", "se_fdr", "dist_kb", "n_se"])

    # expression layer: pyCaCTS lineage-specific TFs
    X = pd.read_csv(os.path.join(PYCACTS, "results/divisions/mtfs_by_OncotreeLineage.tsv"), sep="\t")
    X = X[(X.category == "specific") & (X.fdr_log10 <= -1)]
    expr = {(r.group, r.tf): (int(r.jsd_rank), float(r.fdr_log10)) for r in X.itertuples()}

    # dependency layer: DepMap CRISPR, lineage vs the rest
    D = os.path.join(DATAROOT, "DepMap/2026q1")
    head = pd.read_csv(os.path.join(D, "CRISPRGeneEffect.csv"), nrows=0).columns
    sym = {c.split(" (")[0]: c for c in head[1:]}
    want = [sym[t] for t in set(SE.tf) if t in sym]
    G = pd.read_csv(os.path.join(D, "CRISPRGeneEffect.csv"), usecols=[head[0]] + want, index_col=0)
    G.columns = [c.split(" (")[0] for c in G.columns]
    lin = pd.read_csv(os.path.join(D, "Model.csv"), index_col="ModelID")["OncotreeLineage"].reindex(G.index)

    rows = []
    for r in SE.itertuples():
        ex = expr.get((r.lineage, r.tf))
        dep_in = dep_d = n_in = None
        if r.tf in G.columns:
            m = lin == r.lineage
            n_in = int(m.sum())
            if n_in >= 3:
                dep_in = float(G.loc[m, r.tf].mean()); dep_d = dep_in - float(G.loc[~m, r.tf].mean())
        dep = dep_in is not None and dep_in <= -0.5 and dep_d <= -0.3
        rows.append(dict(lineage=r.lineage, tf=r.tf, se_rank=r.se_rank, se_fdr=round(r.se_fdr, 4), dist_kb=r.dist_kb,
                         expr_rank=ex[0] if ex else None, expr_log10fdr=ex[1] if ex else None,
                         crispr_mean=None if dep_in is None else round(dep_in, 2),
                         crispr_diff=None if dep_d is None else round(dep_d, 2), crispr_lines=n_in,
                         layers=1 + (ex is not None) + dep, preregistered=r.tf in prereg))
    R = pd.DataFrame(rows).sort_values(["layers", "se_rank"], ascending=[False, True])
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    R.to_csv(a.out, sep="\t", index=False)
    two = R[R.layers >= 2]
    print(f"TF x lineage pairs with a specific SE: {len(R):,}; with expression too: {len(two):,} "
          f"({int(two.preregistered.sum())} pre-registered); all three layers: {int((R.layers == 3).sum())}")
    print("\nThree layers (SE + expression + dependency):")
    print(R[R.layers == 3].to_string(index=False))
    print("\nTwo layers, not pre-registered, best SE rank per lineage (top 3):")
    print(two[~two.preregistered].sort_values("se_rank").groupby("lineage").head(3)
          .sort_values(["lineage", "se_rank"]).to_string(index=False))


if __name__ == "__main__":
    main()
