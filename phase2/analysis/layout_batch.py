#!/usr/bin/env python3
"""Is library layout (single- vs paired-end) a hidden batch effect in the atlas?

ChIP-Atlas v1 aligned paired-end runs as single reads (bowtie2 on the reads, MACS2 without BAMPE), so the
two layouts went through different effective pipelines. Layout is also unevenly spread across lineages,
which is what would let a batch effect turn into an apparent lineage call. Four checks:

  1. profile concordance: same line, different study, same vs different layout (Spearman of SE signal)
  2. S3norm exponent B by layout, with read length, depth, instrument and line held fixed
  3. layout-sensitive SEs: within each line that has both layouts, mean PE - mean SE (log signal),
     t-test across lines, BH; calibrated by random sign flips per line
  4. do lineage calls lean on layout-sensitive SEs in PE-rich lineages?

Run with atac_hdac:  python phase2/analysis/layout_batch.py
"""
from __future__ import annotations

import itertools
import os

import numpy as np
import pandas as pd
from scipy import stats

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
P = lambda *a: os.path.join(ROOT, *a)
INST = r"(NovaSeq|NextSeq|HiSeq X|HiSeq 4000|HiSeq 2500|HiSeq 2000|HiSeq 3000|Genome Analyzer|BGISEQ|DNBSEQ)"


def bh(p):
    p = np.asarray(p); o = np.argsort(p); n = len(p)
    q = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    out = np.empty(n); out[o] = np.minimum(q, 1); return out


def ols(d, y, cols, fe=()):
    X = [np.ones(len(d))] + [d[c].values for c in cols]
    for f in fe:
        X += [(d[f] == l).astype(float).values for l in sorted(d[f].unique())[1:]]
    X = np.column_stack(X)
    b, *_ = np.linalg.lstsq(X, d[y].values, rcond=None)
    r = d[y].values - X @ b
    dof = len(d) - np.linalg.matrix_rank(X)
    se = np.sqrt(np.diag(np.linalg.pinv(X.T @ X)) * (r @ r / dof))
    pv = 2 * stats.t.sf(np.abs(b / se), dof)
    return {c: (b[i + 1], pv[i + 1]) for i, c in enumerate(cols)}


def main():
    lay = pd.read_csv(P("phase2/data/srx_layout.tsv"), sep="\t").set_index("srx")
    ps = pd.read_csv(P("phase2/data/pull_set.v2.tsv"), sep="\t").set_index("srx")
    S = pd.read_csv(P("phase2/results_v2/atlas.s3.se_signal.tsv.gz"), sep="\t", index_col=0)
    S = S.loc[(S.values > 0).any(axis=1)]
    cols = [c for c in S.columns if c in lay.index and c in ps.index]
    meta = pd.DataFrame({"line": ps.loc[cols, "key"], "lineage": ps.loc[cols, "lineage"],
                         "layout": lay.loc[cols, "layout"], "study": lay.loc[cols, "study"]}, index=cols)
    X = np.log1p(S[cols].values)
    print(f"[layout] {len(cols)} experiments ({(meta.layout == 'PAIRED').mean():.0%} paired), {X.shape[0]} SEs")

    # 1. concordance, cross-study pairs only (a study usually has one layout, so same-study pairs would
    #    measure the study, not the layout)
    both = meta.groupby("line").layout.nunique().pipe(lambda s: s[s == 2].index)
    sub = meta[meta.line.isin(both)]
    R = stats.rankdata(X[:, [cols.index(c) for c in sub.index]], axis=0)
    R = (R - R.mean(0)) / R.std(0)
    rho = R.T @ R / R.shape[0]
    ix = {s: i for i, s in enumerate(sub.index)}
    pairs = [(ln, meta.layout[a] == meta.layout[b], rho[ix[a], ix[b]])
             for ln, g in sub.groupby("line") for a, b in itertools.combinations(g.index, 2)
             if meta.study[a] != meta.study[b]]
    pr = pd.DataFrame(pairs, columns=["line", "same", "rho"])
    per = pr.groupby(["line", "same"]).rho.mean().unstack().dropna()
    d = per[True] - per[False]
    print(f"[1] cross-study pairs on {len(per)} lines: same-layout minus cross-layout rho median {d.median():+.4f} "
          f"({(d > 0).sum()}/{len(d)} lines higher; Wilcoxon p={stats.wilcoxon(per[True], per[False]).pvalue:.2g}); "
          f"median cross-study rho {pr.rho.median():.3f}")

    # 2. S3norm exponent
    p3 = pd.read_csv(P("phase2/results_v2/atlas.s3.s3norm_params.tsv.gz"), sep="\t").set_index("sample")
    dd = p3.join(lay).join(ps[["key"]]).dropna(subset=["layout", "read_len"])
    dd["inst"] = dd.instrument.str.extract(INST, expand=False).fillna("other")
    dd["PE"] = (dd.layout == "PAIRED").astype(float)
    dd["log2_readlen"] = np.log2(dd.read_len.clip(20, 300))
    dd["log2_reads"] = np.log2(dd.reads.clip(1e6))
    ref = p3.index[p3.is_ref][0]
    print(f"[2] S3norm reference {ref} is {lay.loc[ref, 'layout']}; B median SE {dd.B[dd.PE == 0].median():.3f}, "
          f"PE {dd.B[dd.PE == 1].median():.3f}")
    for fe in ((), ("inst",), ("inst", "key"), ("inst", "key", "study")):
        co = ols(dd, "B", ["PE", "log2_readlen", "log2_reads"], fe)
        print(f"    B ~ PE + log2 read length + log2 reads | fixed effects {','.join(fe) or '-':17s}: "
              + "  ".join(f"{c} {v[0]:+.3f} (p={v[1]:.1g})" for c, v in co.items()))

    # 3. layout-sensitive SEs within lines
    lines = list(both)
    D = np.stack([X[:, [cols.index(c) for c in sub.index[(sub.line == ln) & (sub.layout == "PAIRED")]]].mean(1)
                  - X[:, [cols.index(c) for c in sub.index[(sub.line == ln) & (sub.layout == "SINGLE")]]].mean(1)
                  for ln in lines], axis=1)                            # SEs x lines
    t, p = stats.ttest_1samp(D, 0, axis=1)
    q = bh(np.nan_to_num(p, nan=1.0))
    hits = q < 0.05
    rng = np.random.default_rng(0)
    null = []
    for _ in range(20):
        tn, pn = stats.ttest_1samp(D * rng.choice([-1, 1], size=D.shape[1]), 0, axis=1)
        null.append(int((bh(np.nan_to_num(pn, nan=1.0)) < 0.05).sum()))
    eff = D.mean(1)
    print(f"[3] {len(lines)} lines with both layouts: {hits.sum()} SEs layout-sensitive at FDR<0.05 "
          f"({(hits & (t > 0)).sum()} higher in PE, {(hits & (t < 0)).sum()} higher in SE); "
          f"sign-flip null: median {int(np.median(null))}, max {max(null)} of 20")
    if hits.any():
        print(f"    effect of those SEs: median |PE-SE| {np.median(np.abs(eff[hits])):.2f} log units "
              f"(x{np.exp(np.median(np.abs(eff[hits]))):.2f}); all SEs median |PE-SE| {np.median(np.abs(eff)):.2f}")
    sens = pd.DataFrame({"t": t, "q": q, "pe_minus_se": eff}, index=S.index)

    # 3b. the same test on SHAPE: remove each experiment's global level (median log signal over SEs) first.
    #     A pure level shift moves every SE of a sample together; only a shape difference can make one SE
    #     look specific to a PE-rich or SE-rich group.
    Dc = D - np.median(D, axis=0, keepdims=True)
    tc, pc = stats.ttest_1samp(Dc, 0, axis=1)
    qc = bh(np.nan_to_num(pc, nan=1.0))
    hc = qc < 0.05
    effc = Dc.mean(1)
    print(f"[3b] level removed: {hc.sum()} SEs layout-sensitive in shape ({(hc & (tc > 0)).sum()} PE-up, "
          f"{(hc & (tc < 0)).sum()} SE-up); median |effect| of those {np.median(np.abs(effc[hc])) if hc.any() else 0:.2f} "
          f"log units (x{np.exp(np.median(np.abs(effc[hc]))) if hc.any() else 1:.2f}); median |level shift| per line "
          f"{np.median(np.abs(np.median(D, axis=0))):.2f}, direction SE-up in {(np.median(D, axis=0) < 0).sum()}/{len(lines)} lines")
    sens["t_shape"], sens["q_shape"], sens["pe_minus_se_shape"] = tc, qc, effc

    # 4. lineage calls vs layout-sensitive SEs
    calls = pd.read_csv(P("phase2/scores_v2/atlas.s3.perm.OncotreeLineage.specific.tsv.gz"), sep="\t")
    calls = calls[calls.fdr <= 0.1]
    pe_frac = meta.groupby("lineage").layout.apply(lambda s: (s == "PAIRED").mean())
    rows = []
    for g, c in calls.groupby("group"):
        s = sens.reindex(c.se)
        rows.append((g, pe_frac.get(g, np.nan), len(c), ((s.q_shape < 0.05) & (s.t_shape > 0)).mean(),
                     ((s.q_shape < 0.05) & (s.t_shape < 0)).mean(), s.pe_minus_se_shape.mean()))
    L = pd.DataFrame(rows, columns=["lineage", "pe_frac", "n_calls", "frac_PE_up", "frac_SE_up",
                                    "mean_shape_effect"]).set_index("lineage")
    L["net_PE_up"] = L.frac_PE_up - L.frac_SE_up
    r = stats.spearmanr(L.pe_frac, L.net_PE_up)
    r2 = stats.spearmanr(L.pe_frac, L.mean_shape_effect)
    print(f"[4] {len(calls)} lineage calls; shape-sensitive among them {(sens.reindex(calls.se).q_shape < 0.05).mean():.1%} "
          f"(all SEs {(sens.q_shape < 0.05).mean():.1%}). If layout drove calls, PE-rich lineages would call PE-up SEs: "
          f"Spearman(PE fraction, net PE-up share) {r.correlation:+.2f} (p={r.pvalue:.2g}); "
          f"Spearman(PE fraction, mean shape effect) {r2.correlation:+.2f} (p={r2.pvalue:.2g}); {len(L)} lineages")
    print(L.sort_values("pe_frac").round(3).to_string())
    out = P("phase2/analysis/out"); os.makedirs(out, exist_ok=True)
    sens.to_csv(os.path.join(out, "layout_sensitive_ses.tsv.gz"), sep="\t", float_format="%.4g")
    L.to_csv(os.path.join(out, "layout_by_lineage_calls.tsv"), sep="\t", float_format="%.4g")


if __name__ == "__main__":
    main()
