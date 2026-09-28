#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Poster figure panels for SE-CaCTS (v2 atlas), print-ready: SVG + PDF + 300-dpi PNG per panel.

Sized for an A0 portrait poster with ~260 mm columns: every panel is drawn at its final physical size, so
the type sizes below are the printed sizes (16-18 pt text reads at ~1.5 m). Titles are left to the poster
layout; panels carry axes, legends and direct labels only.

Palette (validated all-pairs with the dataviz validator, light surface): teal = corrected / specific /
DepMap WGS, red = amplicon-driven / removed, orange = CMP WES, violet = DepMap WES.

  ~/miniconda3/envs/atac_hdac/bin/python phase2/figures/poster_figures.py [panel ...]
"""
from __future__ import annotations

import gzip
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                                      # noqa: E402
import numpy as np                                                   # noqa: E402
import pandas as pd                                                  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
SECACTS = os.path.abspath(os.path.join(HERE, "..", ".."))
RES = os.environ.get("POSTER_RES", os.path.join(SECACTS, "phase2", "results_v2"))
SC = os.environ.get("POSTER_SC", os.path.join(SECACTS, "phase2", "scores_v2"))
SC1 = os.path.join(SECACTS, "phase2", "scores")
OUT = os.environ.get("POSTER_OUT", os.path.join(SECACTS, "poster", "figures"))

TEAL, RED, ORANGE, VIOLET = "#008a7e", "#b4443a", "#d4731c", "#4a3aa7"
INK, MUTED, FAINT, GRID, SURF = "#12222a", "#5a6b73", "#9aa5a9", "#e3e8ea", "#ffffff"
MM = 1 / 25.4

plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Liberation Sans", "Arial", "DejaVu Sans"],
    "font.size": 16, "axes.labelsize": 17, "xtick.labelsize": 15, "ytick.labelsize": 15,
    "legend.fontsize": 15, "text.color": INK, "axes.labelcolor": INK, "xtick.color": MUTED,
    "ytick.color": MUTED, "axes.edgecolor": FAINT, "axes.linewidth": 0.8, "axes.spines.top": False,
    "axes.spines.right": False, "axes.grid": False, "grid.color": GRID, "grid.linewidth": 0.8,
    "xtick.major.size": 0, "ytick.major.size": 0, "xtick.major.pad": 6, "ytick.major.pad": 6,
    "legend.frameon": False, "svg.fonttype": "none", "pdf.fonttype": 42, "figure.dpi": 100,
    "savefig.facecolor": SURF,
})


def save(fig, name):
    os.makedirs(OUT, exist_ok=True)
    for ext in ("svg", "pdf", "png"):
        fig.savefig(os.path.join(OUT, f"{name}.{ext}"), dpi=300 if ext == "png" else None,
                    bbox_inches="tight", pad_inches=0.08)
    plt.close(fig)
    print(f"[fig] {name}")


def hgrid(ax):
    ax.xaxis.grid(True); ax.set_axisbelow(True)


def vgrid(ax):
    ax.yaxis.grid(True); ax.set_axisbelow(True)


# ---------------------------------------------------------------------------------------------- data
def v2_panel():
    """Lines in the scored v2 atlas (post-QC-gate), with CN source and lineage."""
    ps = pd.read_csv(os.path.join(SECACTS, "phase2/data/pull_set.v2.tsv"), sep="\t")
    kept = set(pd.read_csv(os.path.join(RES, "atlas.s3.s3norm_params.tsv.gz"), sep="\t")["sample"])
    ps = ps[ps.srx.isin(kept)]
    lin = pd.read_csv(os.path.join(SECACTS, "phase1/data/lineage_resolved.tsv"), sep="\t").set_index("cvcl")
    L = ps.groupby("key").agg(cvcl=("cvcl", "first"), cn=("cn_provider", "first"), n=("srx", "size"))
    L["lineage"] = L["cvcl"].map(lin["lineage"])
    return L


# ------------------------------------------------------------------------------------------ panels
def fig_panel_expansion():
    """Lines per lineage in the scored atlas, split by the copy-number source that admitted them."""
    L = v2_panel()
    order = ["depmap_wgs", "cmp_wes", "depmap_mc_wes"]
    lab = {"depmap_wgs": "DepMap WGS (v1 panel)", "cmp_wes": "+ CMP WES (new)", "depmap_mc_wes": "+ DepMap WES (new)"}
    col = {"depmap_wgs": TEAL, "cmp_wes": ORANGE, "depmap_mc_wes": VIOLET}
    T = L.groupby(["lineage", "cn"]).size().unstack(fill_value=0).reindex(columns=order, fill_value=0)
    T = T.loc[T.sum(axis=1).sort_values().index]
    fig, ax = plt.subplots(figsize=(250 * MM, 190 * MM))
    y = np.arange(len(T))
    left = np.zeros(len(T))
    for c in order:
        v = T[c].values
        ax.barh(y, v, left=left, height=0.68, color=col[c], label=lab[c], edgecolor=SURF, linewidth=1.2)
        left += v
    for yi, tot, new in zip(y, left, T[["cmp_wes", "depmap_mc_wes"]].sum(axis=1).values):
        ax.text(tot + 0.6, yi, f"{int(tot)}" + (f"  (+{int(new)})" if new else ""), va="center",
                fontsize=13, color=MUTED)
    ax.set_yticks(y, T.index, fontsize=14)
    ax.set_xlabel("cell lines in the scored atlas")
    ax.set_xlim(0, left.max() * 1.18)
    hgrid(ax); ax.spines["left"].set_visible(False)
    n1, n2 = int(T["depmap_wgs"].sum()), int(T.values.sum())
    ax.legend(loc="lower right", handlelength=1.1, borderaxespad=0.2,
              title=f"{n1} → {n2} lines (+{100 * (n2 - n1) / n1:.0f}%)", title_fontsize=15)
    save(fig, "fig1_panel_expansion")


LEVEL_LABEL = {"OncotreeLineage": "Lineage", "OncotreePrimaryDisease": "Primary disease",
               "OncotreeSubtype": "Subtype"}


def fig_resolution():
    """Calls per group against lines per group, v1 vs v2, one small multiple per hierarchy level."""
    H1 = pd.read_csv(os.path.join(SC1, "atlas.s3.perm.hierarchy_summary.tsv"), sep="\t")
    H2 = pd.read_csv(os.path.join(SC, "atlas.s3.perm.hierarchy_summary.tsv"), sep="\t")
    levs = list(LEVEL_LABEL)
    fig, axes = plt.subplots(1, 3, figsize=(250 * MM, 105 * MM), sharey=True)
    rng = np.random.default_rng(0)
    for ax, lev in zip(axes, levs):
        for H, lab, face, edge, z in ((H1, "v1 · 282 lines", "none", FAINT, 2),
                                      (H2, "v2 · 386 lines", TEAL, SURF, 3)):
            h = H[H.level == lev]
            x = h["n_lines"].values * np.exp(rng.uniform(-0.06, 0.06, len(h)))      # de-overlap
            ax.scatter(x, h["n_spec_fdr10"].values, s=46, facecolor=face, edgecolor=edge if face != "none" else FAINT,
                       linewidth=1.3 if face == "none" else 1.0, zorder=z, label=lab)
        h1, h2 = H1[H1.level == lev], H2[H2.level == lev]
        ax.set_xscale("log"); ax.set_yscale("symlog", linthresh=1, linscale=0.6)
        ax.set_xticks([1, 2, 5, 10, 20, 50], ["1", "2", "5", "10", "20", "50"])
        ax.set_title(f"{LEVEL_LABEL[lev]}", fontsize=16, color=INK, loc="left", pad=30)
        ax.text(0.0, 1.02, f"groups with calls: {int((h1.n_spec_fdr10 > 0).sum())}/{len(h1)}"
                           f" → {int((h2.n_spec_fdr10 > 0).sum())}/{len(h2)}",
                transform=ax.transAxes, va="bottom", fontsize=12.5, color=MUTED)
        vgrid(ax); ax.xaxis.grid(True)
    axes[0].set_ylabel("specific SEs per group\n(FDR ≤ 0.10)")
    axes[1].set_xlabel("cell lines in the group")
    axes[0].set_yticks([0, 1, 10, 100, 1000], ["0", "1", "10", "100", "1,000"])
    axes[0].set_ylim(-0.3, 5000)
    h, l = axes[0].get_legend_handles_labels()
    fig.legend(h, l, loc="upper right", ncol=2, bbox_to_anchor=(1.0, 1.13), handletextpad=0.3, columnspacing=1.4)
    fig.tight_layout(w_pad=1.2)
    save(fig, "fig2_resolution")


def _count_calls(prefix, lev):
    f = f"{prefix}.{lev}.specific.tsv.gz"
    if not os.path.exists(f):
        return None
    return int((pd.read_csv(f, sep="\t")["fdr"] <= 0.10).sum())


def fig_calibration():
    """The calibration test: the same procedure on real vs shuffled labels, analytic vs permutation null."""
    lev = "OncotreeLineage"
    H = pd.read_csv(os.path.join(SC, "atlas.s3.perm.hierarchy_summary.tsv"), sep="\t")
    n_se = len(pd.read_csv(os.path.join(SC, f"atlas.s3.perm.{lev}.fdr.tsv.gz"), sep="\t", usecols=[0]))
    n_tests = n_se * int((H.level == lev).sum())
    rows = [("Normal-approx. null", "real labels", _count_calls(os.path.join(SC, "atlas.s3.analytic"), lev)),
            ("Normal-approx. null", "shuffled labels", _count_calls(os.path.join(SC, "atlas.s3.analytic.shuffle"), lev)),
            ("Permutation null", "real labels", int(H[H.level == lev]["n_spec_fdr10"].sum())),
            ("Permutation null", "shuffled labels", _count_calls(os.path.join(SC, "atlas.s3.perm.shuffle"), lev))]
    fig, ax = plt.subplots(figsize=(250 * MM, 85 * MM))
    ylab, y = [], []
    for i, (null, lab, n) in enumerate(rows):
        yi = (1 - i // 2) * 2.6 + (1 - i % 2)
        col = (TEAL if lab == "real labels" else (RED if null.startswith("Normal") else FAINT))
        ax.barh(yi, n, height=0.72, color=col)
        pct = 100 * n / n_tests
        ax.text(max(n, 0) + n_tests * 0.0015, yi, f"{n:,}  ({pct:.2f}% of tests)" if n else "0 calls",
                va="center", fontsize=14, color=INK if (lab == "shuffled labels") else MUTED,
                fontweight="bold" if lab == "shuffled labels" else "normal")
        ylab.append(f"{lab}"); y.append(yi)
    ax.set_yticks(y, ylab)
    for yy, t in ((4.1, "Normal-approximation null"), (1.5, "Label-permutation null (B = 1,000)")):
        ax.text(0, yy, t, fontsize=15, color=INK, fontweight="bold", va="bottom", transform=ax.get_yaxis_transform())
    ax.set_xlabel(f"SE × lineage tests called specific at FDR ≤ 0.10  (of {n_tests:,})")
    ax.set_xlim(0, max(r[2] for r in rows) * 1.45)
    ax.set_ylim(-0.6, 4.6)
    hgrid(ax); ax.spines["left"].set_visible(False)
    ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{int(v):,}"))
    save(fig, "fig3_calibration")


def _ablation_sets(lev):
    U = pd.read_csv(os.path.join(SC, f"atlas.s3.perm.nocn.{lev}.specific.tsv.gz"), sep="\t")
    C = pd.read_csv(os.path.join(SC, f"atlas.s3.perm.{lev}.specific.tsv.gz"), sep="\t")
    U, C = U[U.fdr <= 0.10], C[C.fdr <= 0.10]
    u, c = set(zip(U.group, U.se)), set(zip(C.group, C.se))
    return U, C, u, c


def fig_cn_ablation():
    """Left: calls without vs with CN correction. Right: CN at the disputed loci — removed vs rescued."""
    A = pd.read_csv(os.path.join(SC, "atlas.s3.perm.cn_ablation_calls.tsv"), sep="\t")
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(250 * MM, 105 * MM), gridspec_kw={"width_ratios": [1, 1.55]})
    levs = ["OncotreeLineage", "OncotreePrimaryDisease"]
    for i, lev in enumerate(levs):
        U, C, u, c = _ablation_sets(lev)
        yb = (1 - i) * 2.4
        ax.barh(yb + 0.45, len(u), height=0.7, color=FAINT)
        ax.barh(yb - 0.45, len(c), height=0.7, color=TEAL)
        ax.text(len(u) + 80, yb + 0.45, f"{len(u):,}", va="center", fontsize=14, color=MUTED)
        ax.text(len(c) + 80, yb - 0.45, f"{len(c):,}", va="center", fontsize=14, color=INK, fontweight="bold")
        ax.text(0, yb + 1.05, LEVEL_LABEL[lev], fontsize=15, fontweight="bold", va="bottom")
    ax.set_yticks([2.85, 1.95, 0.45, -0.45], ["uncorrected", "CN-corrected"] * 2)
    ax.set_xlabel("specific SEs (FDR ≤ 0.10)")
    ax.set_ylim(-1.1, 3.9); hgrid(ax); ax.spines["left"].set_visible(False)
    ax.set_xlim(0, ax.get_xlim()[1] * 1.25)
    ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{int(v):,}"))

    rng = np.random.default_rng(1)
    kinds = [("rescued", "rescued by correction", TEAL, 0.0), ("amplicon_driven", "removed by correction", RED, 1.0)]
    for k, lab, col, yb in kinds:
        v = A[A.kind == k]["cn_mean"].clip(lower=0.3)
        yj = yb + rng.uniform(-0.22, 0.22, len(v))
        bx.scatter(v, yj, s=9, color=col, alpha=0.35, edgecolor="none", zorder=3)
        med, amp = v.median(), 100 * (v > 1.3).mean()
        bx.text(0.31, yb - 0.26, f"{lab}\nn = {len(v):,} · median CN {med:.2f} · {amp:.0f}% at CN > 1.3",
                fontsize=13, color=INK, va="top", linespacing=1.15)
    amp = (A[A.kind == "amplicon_driven"].sort_values("cn_mean", ascending=False)
           .drop_duplicates("nearest_gene").head(4).sort_values("cn_mean"))
    for i, r in enumerate(amp.itertuples()):
        bx.annotate(r.nearest_gene, (r.cn_mean, 1.2), xytext=(0, 14 + 22 * (i % 3)), textcoords="offset points",
                    ha="center", fontsize=13, color=INK, fontstyle="italic",
                    arrowprops=dict(arrowstyle="-", color=FAINT, lw=0.8, shrinkA=0, shrinkB=2))
    bx.axvline(1.0, color=FAINT, lw=1, zorder=1)
    bx.set_xscale("log")
    bx.set_xticks([0.5, 1, 2, 5, 10, 20, 50, 100], ["0.5", "1", "2", "5", "10", "20", "50", "100"])
    bx.set_yticks([]); bx.spines["left"].set_visible(False)
    bx.set_ylim(-0.72, 1.75)
    bx.set_xlabel("mean copy number at the locus (group)")
    hgrid(bx)
    fig.tight_layout(w_pad=4.0)
    save(fig, "fig4_cn_ablation")


# pre-specified identity TFs per lineage, from the literature — shown whether or not they pass
IDENTITY = [("Ovary/Fallopian Tube", ["PAX8", "SOX17", "MECOM"]), ("Bowel", ["CDX2", "HNF4A"]),
            ("Breast", ["ESR1", "GATA3", "FOXA1"]), ("Myeloid", ["SPI1", "CEBPA"]),
            ("Lymphoid", ["IKZF1", "PAX5"]), ("Skin", ["SOX10", "MITF"]),
            ("Peripheral Nervous System", ["PHOX2B", "HAND2"]), ("Lung", ["NKX2-1", "ASCL1"]),
            ("Liver", ["HNF1A"]), ("Prostate", ["AR"]), ("Kidney", ["PAX2"]), ("Head and Neck", ["TP63"])]


def fig_identity():
    """Known lineage master TFs: the lineage-best SE near each gene, and its FDR in every lineage."""
    sys.path.insert(0, os.path.join(SECACTS, "cnrose"))
    sys.path.insert(0, SECACTS)
    from secacts_env import cache_path
    from cnrose.cn.depmap import load_gene_coords
    gc = load_gene_coords(None, cache_path=cache_path("gene_coords.GRCh38.106.tsv"))
    lev = "OncotreeLineage"
    J = pd.read_csv(os.path.join(SC, f"atlas.s3.perm.{lev}.jsd.tsv.gz"), sep="\t", index_col=0)
    F = pd.read_csv(os.path.join(SC, f"atlas.s3.perm.{lev}.fdr.tsv.gz"), sep="\t", index_col=0)
    cat = pd.read_csv(os.path.join(RES, "atlas.s3.union_catalog.bed.gz"), sep="\t", header=None,
                      usecols=[0, 1, 2, 3], names=["chrom", "start", "end", "se"]).set_index("se")
    cols = [g for g, _ in IDENTITY if g in J.columns]
    rows, lab = [], []
    for grp, genes in IDENTITY:
        if grp not in J.columns:
            continue
        for g in genes:
            if g not in gc:
                continue
            c, s, e = gc[g]
            near = cat[(cat.chrom == c) & (cat.end >= s - 100_000) & (cat.start <= e + 100_000)].index   # score_pilot identity_near
            near = [x for x in near if x in J.index and np.isfinite(J.loc[x, grp])]
            if not near:
                continue
            best = min(near, key=lambda x: J.loc[x, grp])
            rows.append(-np.log10(F.loc[best, cols].astype(float).clip(lower=1e-12)).values)
            rank = int((J[grp] < J.loc[best, grp]).sum()) + 1
            lab.append((g, grp, rank, float(F.loc[best, grp])))
    M = np.array(rows)
    fig, ax = plt.subplots(figsize=(250 * MM, 16 * MM + 9.2 * MM * len(rows)))
    from matplotlib.colors import LinearSegmentedColormap
    cmap = LinearSegmentedColormap.from_list("teal", ["#ffffff", "#bfe3de", "#5fb8ad", TEAL, "#0a4f49"])
    vmax = 3.0
    im = ax.imshow(np.minimum(M, vmax), cmap=cmap, vmin=0, vmax=vmax, aspect="auto", interpolation="nearest")
    for i, (g, grp, rank, fdr) in enumerate(lab):
        j = cols.index(grp)
        ax.add_patch(plt.Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False, edgecolor=INK, lw=1.6))
    ax.set_xticks(range(len(cols)), cols, rotation=35, ha="right", fontsize=13)
    ax.set_yticks(range(len(lab)), [f"{g}  #{r:,}" + ("" if f <= 0.10 else "  (n.s.)") for g, _, r, f in lab],
                  fontsize=13.5)
    for t, (g, *_ ) in zip(ax.get_yticklabels(), lab):
        t.set_fontstyle("italic")
    ax.tick_params(length=0)
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.set_xticks(np.arange(-0.5, len(cols)), minor=True); ax.set_yticks(np.arange(-0.5, len(lab)), minor=True)
    ax.grid(which="minor", color=SURF, linewidth=1.5); ax.tick_params(which="minor", length=0)
    cb = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.02, ticks=[0, 1, 2, 3])
    cb.ax.set_yticklabels(["1", "0.1", "0.01", "≤0.001"], fontsize=12)
    cb.set_label("FDR (lineage)", fontsize=14); cb.outline.set_visible(False)
    save(fig, "fig5_identity")
    return pd.DataFrame(lab, columns=["gene", "lineage", "rank", "fdr"])


def fig_concordance():
    """Genes beside group-specific SEs are group-specific in expression; the link decays with distance."""
    C = pd.read_csv(os.path.join(SC, "atlas.s3.perm.concordance2.summary.tsv"), sep="\t")
    P = pd.read_csv(os.path.join(SC, "atlas.s3.perm.concordance2.pairs.tsv.gz"), sep="\t")
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(250 * MM, 95 * MM), gridspec_kw={"width_ratios": [1, 1.3]})
    x = np.arange(len(C))
    w = 0.26
    for k, (col, lab, c) in enumerate((("per_pair", "SE-adjacent gene", TEAL), ("shuffled", "same gene, shuffled group", FAINT),
                                       ("background", "all genes (background)", "#cfd6d8"))):
        v = 100 * C[col].values
        ax.bar(x + (k - 1) * w, v, width=w * 0.92, color=c, label=lab)
    for xi, r in zip(x, C.itertuples()):
        ax.text(xi - w, 100 * r.per_pair + 0.5, f"{100 * r.per_pair:.0f}%  ({r.enrichment:.1f}×)", ha="left",
                fontsize=13, color=INK)
    ax.set_xticks(x, [LEVEL_LABEL[l] for l in C.level])
    ax.set_ylabel("% of pairs with the gene\nspecific to the SE's group")
    ax.set_ylim(0, max(100 * C.per_pair.max() * 1.75, 10))
    vgrid(ax); ax.legend(loc="upper right", fontsize=12.5, handlelength=1.0, ncol=1, borderaxespad=0)
    ax.spines["bottom"].set_visible(False)

    P["kb"] = P["dist_bp"] / 1000
    bins = [0, 10, 25, 50, 100, 250]
    labs = ["<10", "10–25", "25–50", "50–100", ">100"]
    P["bin"] = pd.cut(P["kb"], bins=bins + [1e9], labels=labs + ["x"], right=False)
    P = P[P["bin"] != "x"]
    for lev, mk in (("OncotreeLineage", "o"), ("OncotreePrimaryDisease", "s")):
        d = P[P.level == lev].groupby("bin", observed=True).agg(c=("concordant", "mean"), s=("shuffled", "mean"))
        bx.plot(range(len(d)), 100 * d["c"], marker=mk, ms=8, lw=2, color=TEAL,
                label=f"{LEVEL_LABEL[lev]}", mec=SURF, mew=1.5)
        bx.plot(range(len(d)), 100 * d["s"], marker=mk, ms=7, lw=1.6, color=FAINT, mec=SURF, mew=1.2)
    bx.text(len(labs) - 1, 100 * P[P.level == "OncotreeLineage"]["shuffled"].mean() + 1.2, "shuffled control",
            ha="right", fontsize=13, color=MUTED)
    bx.set_xticks(range(len(labs)), labs)
    bx.set_xlabel("SE → gene distance (kb)")
    bx.set_ylabel("% concordant")
    bx.set_ylim(0, None); vgrid(bx)
    bx.legend(loc="upper right", fontsize=13)
    fig.tight_layout(w_pad=2.2)
    save(fig, "fig6_concordance")


def fig_cn_sources():
    """Correction strength by CN source: signal-vs-CN correlation per line, before and after correction."""
    D = pd.read_csv(os.path.join(SC, "atlas.s3.perm.cn_by_line.tsv"), sep="\t")
    order = [("depmap_wgs", "DepMap WGS", TEAL), ("cmp_wes", "CMP WES", ORANGE), ("depmap_mc_wes", "DepMap WES", VIOLET)]
    fig, ax = plt.subplots(figsize=(250 * MM, 90 * MM))
    rng = np.random.default_rng(2)
    for i, (k, lab, col) in enumerate(order):
        d = D[D.cn_source == k]
        for j, (c, alpha) in enumerate((("rho_raw", 0.35), ("rho_corrected", 0.9))):
            xs = i * 3 + j + rng.uniform(-0.22, 0.22, len(d))
            ax.scatter(xs, d[c], s=22, color=col if j else FAINT, alpha=alpha, edgecolor="none", zorder=3)
            m = d[c].median()
            ax.plot([i * 3 + j - 0.32, i * 3 + j + 0.32], [m, m], color=INK, lw=2.2, zorder=4)
        ax.text(i * 3 + 0.5, ax.get_ylim()[1] if False else 0.47, f"{lab}\n(n = {len(d)})", ha="center",
                va="top", fontsize=14)
    ax.set_xticks([i * 3 + j for i in range(3) for j in range(2)], ["raw", "corrected"] * 3, fontsize=13.5)
    ax.axhline(0, color=FAINT, lw=1)
    ax.set_ylabel("per-line Spearman ρ,\nSE signal vs copy number")
    ax.set_ylim(-0.35, 0.5); vgrid(ax); ax.spines["bottom"].set_visible(False)
    save(fig, "fig7_cn_sources")


PANELS = {"expansion": fig_panel_expansion, "resolution": fig_resolution, "calibration": fig_calibration,
          "ablation": fig_cn_ablation, "identity": fig_identity, "concordance": fig_concordance,
          "cnsources": fig_cn_sources}

if __name__ == "__main__":
    want = sys.argv[1:] or list(PANELS)
    for k in want:
        PANELS[k]()
