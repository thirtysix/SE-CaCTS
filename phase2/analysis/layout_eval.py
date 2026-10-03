#!/usr/bin/env python3
"""Retention of a build's calls when the within-line library-layout effect is removed (arm `layout`, FINDINGS §16).

Per group: calls kept (same locus, FDR <= 0.10 in both), new calls, and the group's share of paired-end experiments,
with the Spearman correlation between retention and PE share (v2: +0.82, the signature of a level batch effect).

    python3 phase2/analysis/layout_eval.py --base phase2/scores_v31f --arm phase2/scores_v31f_layout \\
        --pe phase2/analysis/out/layout_by_lineage_calls.v31f.tsv --out phase2/analysis/out/v31f_robust/layout_retention
"""
import argparse, os
import pandas as pd
from scipy import stats


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", required=True, help="dir with atlas.s3.perm.<level>.specific.tsv.gz")
    ap.add_argument("--arm", required=True, help="dir with atlas.s3.perm.layoutcorr.<level>.specific.tsv.gz")
    ap.add_argument("--pe", required=True, help="layout_batch.py per-lineage table (pe_frac)")
    ap.add_argument("--fdr", type=float, default=0.10)
    ap.add_argument("--out", required=True, help="output prefix")
    a = ap.parse_args()
    pe = pd.read_csv(a.pe, sep="\t", index_col=0)["pe_frac"]
    for lev in ("OncotreeLineage", "OncotreePrimaryDisease"):
        b = pd.read_csv(os.path.join(a.base, f"atlas.s3.perm.{lev}.specific.tsv.gz"), sep="\t")
        r = pd.read_csv(os.path.join(a.arm, f"atlas.s3.perm.layoutcorr.{lev}.specific.tsv.gz"), sep="\t")
        b, r = b[b.fdr <= a.fdr], r[r.fdr <= a.fdr]
        rows = []
        for g in sorted(set(b.group) | set(r.group)):
            bs, rs = set(b.se[b.group == g]), set(r.se[r.group == g])
            rows.append({"group": g, "base": len(bs), "layoutcorr": len(rs), "kept": len(bs & rs), "new": len(rs - bs),
                         "retention": round(len(bs & rs) / len(bs), 3) if bs else float("nan"),
                         "pe_frac": pe.get(g, float("nan"))})
        T = pd.DataFrame(rows)
        T.to_csv(f"{a.out}.{lev}.tsv", sep="\t", index=False)
        tot_b, kept, new = T.base.sum(), T.kept.sum(), T.new.sum()
        print(f"[{lev}] {tot_b} -> {T.layoutcorr.sum()} calls; kept {kept} ({kept / tot_b:.1%}), new {new}")
        if lev == "OncotreeLineage":
            t = T.dropna(subset=["pe_frac", "retention"])
            t = t[t.base >= 20]
            rho = stats.spearmanr(t.pe_frac, t.retention)
            lo = t[t.pe_frac < 0.25]
            print(f"  retention vs PE share: Spearman {rho.correlation:+.2f} (p={rho.pvalue:.2g}, {len(t)} lineages "
                  f"with >= 20 calls); PE share < 25%: {lo.kept.sum()}/{lo.base.sum()} kept "
                  f"({lo.kept.sum() / max(lo.base.sum(), 1):.1%}), the rest {t.kept.sum() - lo.kept.sum()}/"
                  f"{t.base.sum() - lo.base.sum()} ({(t.kept.sum() - lo.kept.sum()) / max(t.base.sum() - lo.base.sum(), 1):.1%})")
            print(T.sort_values("retention").to_string(index=False))


if __name__ == "__main__":
    main()
