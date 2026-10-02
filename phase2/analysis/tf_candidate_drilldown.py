#!/usr/bin/env python3
"""Drill-down for TF candidates (FINDINGS §41): which lines drive each lineage-specific SE, and is it the TF's SE.

For each (lineage, TF): every v3.0.1 lineage call within 100 kb of the TF, with the genes listed for that locus
(docs/data/se_genes.json: edge-to-body distance, 0 = overlap), the lineage's lines ranked by mean S3norm signal at
the locus (uncorrected; CN shown), how many of each line's experiments call an SE overlapping it, the line's
atlas-wide signal rank, and the TF's DepMap expression in that line.

    python3 phase2/analysis/tf_candidate_drilldown.py Bone:VAX1 Uterus:POU3F3 Lung:ST18 Skin:ZNF536
"""
import argparse, gzip, importlib.util, json, os, sys
import numpy as np
import pandas as pd

SECACTS = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
sys.path.insert(0, SECACTS)
spec = importlib.util.spec_from_file_location("cf", os.path.join(SECACTS, "phase2/scripts/74_called_filter.py"))
cf = importlib.util.module_from_spec(spec); spec.loader.exec_module(cf)
sys.path.insert(0, os.path.join(SECACTS, "cnrose"))
from secacts_env import DATAROOT, cache_path                             # noqa: E402
from cnrose.cn.depmap import load_gene_coords                            # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("pairs", nargs="+", help="Lineage:TF")
    ap.add_argument("--scores", default=os.path.join(SECACTS, "phase2/scores_v3c"))
    ap.add_argument("--results", default=os.path.join(SECACTS, "phase2/results_v3"))
    ap.add_argument("--pull-set", default=os.path.join(SECACTS, "phase2/data/pull_set.v3.tsv"))
    ap.add_argument("--window", type=int, default=100_000)
    ap.add_argument("--top", type=int, default=8)
    ap.add_argument("--out", default=os.path.join(SECACTS, "phase2/analysis/out/tf_candidate_drilldown.tsv"))
    a = ap.parse_args()

    cat = pd.read_csv(f"{a.results}/atlas.s3.union_catalog.bed.gz", sep="\t", header=None, usecols=[0, 1, 2, 3],
                      names=["chrom", "start", "end", "se"]).set_index("se")
    S = pd.read_csv(f"{a.scores}/atlas.s3.perm.OncotreeLineage.specific.tsv.gz", sep="\t").join(cat, on="se")
    gc = load_gene_coords(None, cache_path=cache_path("gene_coords.GRCh38.106.tsv"))
    genes = json.load(open(os.path.join(SECACTS, "docs/data/se_genes.json")))

    targets = []
    for p in a.pairs:
        lin, tf = p.split(":")
        c, s, e = gc[tf]
        d = S[(S.group == lin) & (S.chrom == c) & (S.end >= s - a.window) & (S.start <= e + a.window)]
        targets += [(lin, tf, r) for r in d.itertuples()]
    loci = sorted({r.se for _, _, r in targets})
    print(f"{len(targets)} calls on {len(loci)} loci", file=sys.stderr)

    def rows_of(path):
        with gzip.open(path, "rt") as fh:
            head = fh.readline().rstrip("\n").split("\t")
            keep = []
            for line in fh:
                k = line[:line.index("\t")]
                if k in loci:
                    keep.append(line.rstrip("\n").split("\t"))
        return pd.DataFrame([r[1:] for r in keep], index=[r[0] for r in keep], columns=head[1:]).astype(float)
    sig = rows_of(f"{a.results}/atlas.s3.se_signal.tsv.gz")
    pres_all = pd.read_csv(f"{a.results}/atlas.s3.se_presence.tsv.gz", sep="\t", index_col=0)

    ps = pd.read_csv(a.pull_set, sep="\t")
    ps = ps[ps.srx.isin(sig.columns)].copy()
    ps.loc[ps.srx.eq("SRX10809652"), "key"] = "ACH-000768"
    L = cf.labels(ps, os.path.join(DATAROOT, "DepMap/2026q1/Model.csv"),
                  os.path.join(SECACTS, "phase1/data/lineage_resolved.tsv"))
    key_of = dict(zip(ps.srx, ps.key))
    G = cf.overlap_graph(cat.reset_index())
    se_idx = {s: i for i, s in enumerate(cat.index)}
    pres_all = pres_all.reindex(columns=sig.columns, fill_value=0)

    # expression (DepMap log2 TPM+1) of each TF per line
    D = os.path.join(DATAROOT, "DepMap/2026q1")
    eh = pd.read_csv(os.path.join(D, "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv"), nrows=0).columns
    tfs = sorted({tf for _, tf, _ in targets})
    ecols = [c for c in eh if c.split(" (")[0] in tfs]
    E = pd.read_csv(os.path.join(D, "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv"),
                    usecols=["ModelID", "IsDefaultEntryForModel"] + ecols)
    E = E[E.IsDefaultEntryForModel == "Yes"].drop(columns="IsDefaultEntryForModel").set_index("ModelID")
    E.columns = [c.split(" (")[0] for c in E.columns]

    out = []
    for lin, tf, r in targets:
        i = se_idx[r.se]
        nb = set(G[i].indices) | {i}
        nb_ids = [cat.index[j] for j in nb]
        called = pres_all.reindex(nb_ids).fillna(0).max(axis=0)          # experiment calls an overlapping SE
        line_sig = sig.loc[r.se].groupby(sig.columns.map(key_of)).mean()
        line_called = called.groupby(called.index.map(key_of)).agg(["sum", "size"])
        atlas_rank = line_sig.rank(ascending=False)
        members = L.index[L.OncotreeLineage == lin]
        m = line_sig.reindex(members).dropna().sort_values(ascending=False)
        gl = genes.get(r.se, [])
        gtxt = "; ".join(f"{g}:{d}kb" for g, d, *_ in gl[:6])
        for k in m.index[:a.top]:
            out.append(dict(lineage=lin, tf=tf, se=r.se, locus=f"{r.chrom}:{r.start}-{r.end}",
                            size_kb=round((r.end - r.start) / 1000, 1), rank=r.rank, fdr=round(r.fdr, 4),
                            cn_mean=r.cn_mean, genes=gtxt, line=L.loc[k, "line"],
                            subtype=L.loc[k, "OncotreeSubtype"], signal=round(m[k], 2),
                            atlas_rank=int(atlas_rank[k]), of_lines=len(line_sig),
                            exps_calling=f"{int(line_called.loc[k, 'sum'])}/{int(line_called.loc[k, 'size'])}",
                            tf_expr=round(float(E.loc[k, tf]), 2) if (k in E.index and tf in E.columns) else None))
    T = pd.DataFrame(out)
    T.to_csv(a.out, sep="\t", index=False)
    with pd.option_context("display.width", 250, "display.max_colwidth", 60, "display.max_rows", 500):
        for (lin, tf, se), d in T.groupby(["lineage", "tf", "se"], sort=False):
            h = d.iloc[0]
            print(f"\n== {tf} in {lin}: {se} {h.locus} ({h.size_kb} kb) rank {h['rank']} FDR {h.fdr} CN {h.cn_mean}"
                  f"\n   genes: {h.genes}")
            print(d[["line", "subtype", "signal", "atlas_rank", "of_lines", "exps_calling", "tf_expr"]].to_string(index=False))


if __name__ == "__main__":
    main()
