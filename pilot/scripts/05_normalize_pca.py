#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Step B pilot — 05: normalize the region x sample matrix and test batch (study) vs lineage.

For each of {raw-log, quantile, reference(median-of-ratios)} normalizations:
  - PCA of samples (numpy SVD; no sklearn),
  - PC1/PC2 scatter coloured by LINEAGE, annotated by cell line (batch-probe replicates flagged),
  - silhouette of samples in PC-space by LINEAGE vs by STUDY/cell-line.

Decision read (see ../PILOT.md): normalization should RAISE lineage silhouette and LOWER study silhouette,
and the cross-study replicate of the same cell line should co-cluster. Env: atac_hdac (numpy/pandas/
scipy/matplotlib). Runs locally or on HPC.
"""
import argparse, os, sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# ---------- normalizations ----------
def log_raw(X):
    return np.log1p(X)


def quantile_norm(X):
    """Quantile-normalize columns (samples) to a common (mean-of-sorted) reference distribution."""
    Xl = np.log1p(X)
    order = np.argsort(Xl, axis=0)
    ranks_sorted = np.take_along_axis(Xl, order, axis=0)
    ref = ranks_sorted.mean(axis=1)
    out = np.empty_like(Xl)
    for j in range(Xl.shape[1]):
        out[order[:, j], j] = ref
    return out


def reference_norm(X):
    """DESeq-style median-of-ratios size factors, then log."""
    Xf = X.astype(float)
    with np.errstate(divide="ignore"):
        loggeo = np.mean(np.log(np.where(Xf > 0, Xf, np.nan)), axis=1)
    keep = np.isfinite(loggeo)
    geo = np.exp(loggeo[keep])
    sf = np.array([np.median((Xf[keep, j] / geo)[np.isfinite(Xf[keep, j] / geo)])
                   for j in range(Xf.shape[1])])
    sf[~np.isfinite(sf) | (sf <= 0)] = 1.0
    return np.log1p(Xf / sf[None, :])


NORMS = {"raw_log": log_raw, "quantile": quantile_norm, "reference": reference_norm}


# ---------- PCA + silhouette ----------
def pca_scores(M, k=None):
    """M: samples x features. Return PC scores (samples x k) and explained-variance ratio."""
    Mc = M - M.mean(axis=0, keepdims=True)
    U, S, _ = np.linalg.svd(Mc, full_matrices=False)
    scores = U * S
    evr = (S ** 2) / (S ** 2).sum()
    if k:
        scores = scores[:, :k]
    return scores, evr


def silhouette(scores, labels):
    """Mean silhouette over samples in Euclidean PC-space. Singleton-only labels -> nan."""
    labels = np.asarray(labels)
    n = len(labels)
    from scipy.spatial.distance import squareform, pdist
    D = squareform(pdist(scores))
    sil = []
    for i in range(n):
        same = (labels == labels[i]) & (np.arange(n) != i)
        if same.sum() == 0:
            continue  # i is a singleton in its cluster
        a = D[i, same].mean()
        b = min(D[i, labels == L].mean() for L in set(labels) if L != labels[i])
        sil.append((b - a) / max(a, b) if max(a, b) > 0 else 0.0)
    return float(np.mean(sil)) if sil else float("nan")


def load_matrix(tab):
    df = pd.read_csv(tab, sep="\t")
    df.columns = [c.strip().strip("'\"").lstrip("#").strip() for c in df.columns]
    meta = [c for c in df.columns[:3]]
    samples = [c for c in df.columns if c not in meta]
    X = df[samples].to_numpy(float)
    # drop non-informative regions (any-NaN or all-zero)
    good = ~np.isnan(X).any(axis=1) & (X.sum(axis=1) > 0)
    return X[good], samples


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--signal", default="../results/signal.tab")
    ap.add_argument("--selection", default="../data/selection.tsv")
    ap.add_argument("--outdir", default="../results")
    ap.add_argument("--top-var", type=int, default=20000, help="keep N most-variable regions for PCA")
    a = ap.parse_args()

    X, samples = load_matrix(a.signal)
    sel = pd.read_csv(a.selection, sep="\t").set_index("srx")
    lineage = [sel.loc[s, "lineage"] if s in sel.index else "?" for s in samples]
    cline = [sel.loc[s, "cell_line"] if s in sel.index else s for s in samples]

    # Batch is measured by the cross-study replicate probe: the same cell line profiled in different
    # studies (roles 'batch-probe' + its 'lineage' rep). If those co-cluster, cross-study batch is mild.
    # (A GSE-level silhouette is undefined here because each experiment is its own GSE by design.)
    probe_line = None
    if "role" in sel.columns and (sel["role"] == "batch-probe").any():
        probe_line = sel.loc[sel.index[sel["role"] == "batch-probe"][0], "cell_line"]
    probe_idx = [i for i, s in enumerate(samples)
                 if s in sel.index and probe_line and sel.loc[s, "cell_line"] == probe_line]

    from scipy.spatial.distance import pdist, squareform

    def probe_cohesion(scores):
        """mean within-probe PC-distance / median over all pairs. <1 = replicates tighter than average."""
        if len(probe_idx) < 2:
            return float("nan")
        D = squareform(pdist(scores))
        within = [D[probe_idx[i], probe_idx[j]]
                  for i in range(len(probe_idx)) for j in range(i + 1, len(probe_idx))]
        allp = D[np.triu_indices(len(scores), 1)]
        return float(np.mean(within) / np.median(allp))

    os.makedirs(a.outdir, exist_ok=True)
    rows = []
    fig, axes = plt.subplots(1, len(NORMS), figsize=(6 * len(NORMS), 5.4), squeeze=False)
    palette = plt.cm.tab10(np.linspace(0, 1, len(set(lineage))))
    lin2col = {L: palette[i] for i, L in enumerate(sorted(set(lineage)))}

    for ax, (name, fn) in zip(axes[0], NORMS.items()):
        M = fn(X).T  # samples x features
        # restrict to most-variable regions (computed on this normalization)
        v = M.var(axis=0)
        idx = np.argsort(v)[-min(a.top_var, M.shape[1]):]
        scores, evr = pca_scores(M[:, idx], k=min(len(samples) - 1, 10))
        sil_lin = silhouette(scores, lineage)
        coh = probe_cohesion(scores)
        rows.append(dict(normalization=name, sil_lineage=round(sil_lin, 3),
                         probe_cohesion=round(coh, 3) if coh == coh else float("nan"),
                         probe_clusters=("yes" if coh == coh and coh < 1 else "no"),
                         pc1_var=round(float(evr[0]), 3), pc2_var=round(float(evr[1]), 3)))
        for i, s in enumerate(samples):
            ax.scatter(scores[i, 0], scores[i, 1], color=lin2col[lineage[i]], s=90,
                       edgecolor="k", linewidth=0.5, zorder=3)
            tag = cline[i] + ("*" if sel.loc[s, "role"] == "batch-probe" and s in sel.index else "")
            ax.annotate(tag, (scores[i, 0], scores[i, 1]), fontsize=7,
                        xytext=(4, 3), textcoords="offset points")
        coh_s = f"{coh:.2f}" if coh == coh else "n/a"
        ax.set_title(f"{name}\nsil(lineage)={sil_lin:.2f}  {probe_line or 'probe'} cohesion={coh_s}", fontsize=11)
        ax.set_xlabel(f"PC1 ({evr[0]*100:.0f}%)")
        ax.set_ylabel(f"PC2 ({evr[1]*100:.0f}%)")

    handles = [plt.Line2D([0], [0], marker="o", ls="", mfc=lin2col[L], mec="k", label=L)
               for L in sorted(set(lineage))]
    fig.legend(handles=handles, loc="lower center", ncol=len(handles), fontsize=8, frameon=False)
    fig.suptitle("Step B pilot — H3K27ac batch(study) vs lineage under 3 normalizations "
                 "(* = cross-study replicate)", fontsize=12)
    fig.tight_layout(rect=(0, 0.06, 1, 0.96))
    png = os.path.join(a.outdir, "pca_batch_vs_lineage.png")
    fig.savefig(png, dpi=140)

    summ = pd.DataFrame(rows)
    summ.to_csv(os.path.join(a.outdir, "pca_summary.tsv"), sep="\t", index=False)
    print(summ.to_string(index=False), file=sys.stderr)
    print(f"[05] wrote {png} and pca_summary.tsv", file=sys.stderr)
    print("[05] read: pick the normalization with the highest sil(lineage) AND probe_cohesion<1 "
          "(same-line, different-study replicates cluster => batch is normalizable). Carry it into Phase 2.",
          file=sys.stderr)


if __name__ == "__main__":
    main()
