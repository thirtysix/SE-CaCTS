#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Poster figure panels for SE-CaCTS (v2 by default; v3 via POSTER_* below), print-ready: SVG + PDF + 300-dpi PNG per panel.

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
# the atlas drawn and the one it is compared with (fig2): v2 vs v1 by default, v3 vs v2 via 72_finish_v3.sh
PS = os.environ.get("POSTER_PS", os.path.join(SECACTS, "phase2", "data", "pull_set.v2.tsv"))
PREV_SC = os.environ.get("POSTER_PREV_SC", SC1)
LABEL, PREV_LABEL = os.environ.get("POSTER_LABEL", "v2"), os.environ.get("POSTER_PREV_LABEL", "v1")

TEAL, RED, ORANGE, VIOLET = "#008a7e", "#b4443a", "#d4731c", "#4a3aa7"
OCHRE, BLUE = "#b8930f", "#3f6fb0"
INK, MUTED, FAINT, GRID, SURF = "#12222a", "#5a6b73", "#9aa5a9", "#e3e8ea", "#ffffff"
RAMP = ["#ffffff", "#bfe3de", "#5fb8ad", TEAL, "#0a4f49"]
# poster palette variants: a roles JSON from ClaudeSkills/conference-poster/scripts/palette_roles.py replaces the
# data colours (primary, secondary, the five CN-source categories, the FDR ramp) and the greys; unset = the above
PALETTE = os.environ.get("POSTER_PALETTE")
if PALETTE:
    import json
    _R = json.load(open(PALETTE))
    TEAL, RED = _R["fig_primary"], _R["fig_secondary"]
    _, ORANGE, VIOLET, OCHRE, BLUE = _R["fig_categories"][:5]
    INK, MUTED, FAINT, GRID = _R["fig_ink"], _R["fig_muted"], _R["fig_faint"], _R["fig_grid"]
    RAMP = _R["fig_ramp"]
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


def _cell_text(rgba):
    """SURF or INK, whichever has the higher WCAG contrast on a heatmap cell (palette variants only)."""
    lum = lambda c: sum(w * (v / 12.92 if v <= 0.04045 else ((v + 0.055) / 1.055) ** 2.4)      # noqa: E731
                        for w, v in zip((0.2126, 0.7152, 0.0722), matplotlib.colors.to_rgb(c)))
    ratio = lambda a, b: (max(lum(a), lum(b)) + 0.05) / (min(lum(a), lum(b)) + 0.05)          # noqa: E731
    return SURF if ratio(SURF, rgba) >= ratio(INK, rgba) else INK


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
# CN sources in the order they were added to the atlas, with the colour each keeps in every panel
PROVIDERS = [("depmap_wgs", "DepMap WGS", TEAL), ("cmp_wes", "CMP WES", ORANGE), ("depmap_mc_wes", "DepMap WES", VIOLET),
             ("ccle_snp6", "CCLE SNP6", OCHRE), ("input_inferred", "inferred from ChIP input", BLUE)]


def atlas_panel():
    """Lines in the scored atlas (post-QC-gate), with CN source and lineage."""
    ps = pd.read_csv(PS, sep="\t")
    kept = set(pd.read_csv(os.path.join(RES, "atlas.s3.s3norm_params.tsv.gz"), sep="\t")["sample"])
    ps = ps[ps.srx.isin(kept)]
    lin = pd.read_csv(os.path.join(SECACTS, "phase1/data/lineage_resolved.tsv"), sep="\t").set_index("cvcl")
    L = ps.groupby("key").agg(cvcl=("cvcl", "first"), cn=("cn_provider", "first"), n=("srx", "size"),
                              ps_lineage=("lineage", "first"))
    L["lineage"] = L["cvcl"].map(lin["lineage"]).fillna(L["ps_lineage"])
    return L


# ------------------------------------------------------------------------------------------ panels
def fig_panel_expansion():
    """Lines per lineage in the scored atlas, split by the copy-number source that admitted them."""
    L = atlas_panel()
    order = [k for k, _, _ in PROVIDERS if k in set(L["cn"])]
    lab = {k: n for k, n, _ in PROVIDERS}
    col = {k: c for k, _, c in PROVIDERS}
    T = L.groupby(["lineage", "cn"]).size().unstack(fill_value=0).reindex(columns=order, fill_value=0)
    T = T.loc[T.sum(axis=1).sort_values().index]
    fig, ax = plt.subplots(figsize=(205 * MM, 140 * MM))
    y = np.arange(len(T))
    left = np.zeros(len(T))
    for c in order:
        v = T[c].values
        ax.barh(y, v, left=left, height=0.68, color=col[c], label=lab[c], edgecolor=SURF, linewidth=1.2)
        left += v
    for yi, tot in zip(y, left):
        ax.text(tot + 0.6, yi, f"{int(tot)}", va="center", fontsize=12.5, color=MUTED)
    ax.set_yticks(y, T.index, fontsize=13)
    ax.set_ylim(-0.6, len(T) - 0.4)
    ax.set_xlabel("cell lines")
    ax.set_xlim(0, left.max() * 1.08)
    hgrid(ax); ax.spines["left"].set_visible(False)
    ax.legend(loc="lower right", handlelength=1.1, borderaxespad=0.2, fontsize=13.5,
              title="copy-number source", title_fontsize=14, alignment="left")
    save(fig, "fig1_panel_expansion")


LEVEL_LABEL = {"OncotreeLineage": "Lineage", "OncotreePrimaryDisease": "Primary disease",
               "OncotreeSubtype": "Subtype"}


def fig_resolution():
    """Calls per group against lines per group, previous vs current atlas, one small multiple per level."""
    H1 = pd.read_csv(os.path.join(PREV_SC, "atlas.s3.perm.hierarchy_summary.tsv"), sep="\t")
    H2 = pd.read_csv(os.path.join(SC, "atlas.s3.perm.hierarchy_summary.tsv"), sep="\t")
    nl = lambda H: int(H.loc[H.level == "OncotreeLineage", "n_lines"].sum())         # noqa: E731
    levs = list(LEVEL_LABEL)
    fig, axes = plt.subplots(1, 3, figsize=(250 * MM, 105 * MM), sharey=True)
    rng = np.random.default_rng(0)
    for ax, lev in zip(axes, levs):
        for H, lab, face, edge, z in ((H1, f"{PREV_LABEL} · {nl(H1)} lines", "none", FAINT, 2),
                                      (H2, f"{LABEL} · {nl(H2)} lines", TEAL, SURF, 3)):
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


DEL_CN = 0.3   # below this group-mean CN, dividing by CN (floored at 0.1) inflates noise: deletions, absent chrY


def _ablation_pair(lev):
    """(uncorrected, corrected) distinct (group, SE) calls at FDR <= 0.10. Cell line comes from the consensus-of-
    studies arm, which has a permutation null; the mean arm falls back to the analytic null at line level."""
    # cell line: the per-line "vs all" comparison (untestable lines left out of BH, FINDINGS §37); older atlases: conss
    lines_all = os.path.join(SC, "out_lines", "atlas.s3.lines.all.line.specific.tsv.gz")
    stem = ((("out_lines/atlas.s3.lines.all" if os.path.exists(lines_all) else "atlas.s3.conss"), "line")
            if lev == "line" else ("atlas.s3.perm", lev))
    u = os.path.join(SC, f"{stem[0]}.nocn.{stem[1]}.specific.tsv.gz")
    c = os.path.join(SC, f"{stem[0]}.{stem[1]}.specific.tsv.gz")
    if not (os.path.exists(u) and os.path.exists(c)):
        return None
    n = lambda f: len({(g, e) for g, e, q in pd.read_csv(f, sep="\t", usecols=["group", "se", "fdr"])   # noqa: E731
                       .itertuples(index=False) if q <= 0.10})
    return n(u), n(c)


ABL_LABEL = dict(LEVEL_LABEL, line="Cell line")


def fig_cn_ablation(levs=("OncotreeLineage", "OncotreePrimaryDisease"), name="fig4_cn_ablation", width=420,
                    base_mm=38, per_level_mm=48):
    """Left: calls without vs with CN correction per level. Right: CN at the disputed loci — removed vs rescued,
    with the rescues inside deep deletions (CN < 0.3) shown apart as the correction artifacts they are."""
    A = pd.read_csv(os.path.join(SC, "atlas.s3.perm.cn_ablation_calls.tsv"), sep="\t")
    levs = [l for l in levs if _ablation_pair(l) is not None]
    k = len(levs)
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(width * MM, (base_mm + per_level_mm * k) * MM),
                                 gridspec_kw={"width_ratios": [1, 1.7]})
    ticks, labs = [], []
    per_level = k > 2          # counts differ by orders of magnitude across levels: scale bars within each level
    for i, lev in enumerate(levs):
        nu, nc = _ablation_pair(lev)
        sc = 1 / max(nu, nc, 1) if per_level else 1
        yb = (k - 1 - i) * 2.4
        ax.barh(yb + 0.45, nu * sc, height=0.7, color=FAINT)
        ax.barh(yb - 0.45, nc * sc, height=0.7, color=TEAL)
        ax.text(nu * sc, yb + 0.45, f"  {nu:,}", va="center", fontsize=14, color=MUTED)
        ax.text(nc * sc, yb - 0.45, f"  {nc:,}", va="center", fontsize=14, color=INK, fontweight="bold")
        ax.text(0, yb + 1.05, ABL_LABEL[lev], fontsize=15, fontweight="bold", va="bottom")
        ticks += [yb + 0.45, yb - 0.45]; labs += ["uncorrected", "CN-corrected"]
    ax.set_yticks(ticks, labs)
    ax.set_ylim(-1.1, (k - 1) * 2.4 + 1.5); ax.spines["left"].set_visible(False)
    if per_level:
        ax.set_xlim(0, 1.45); ax.set_xticks([]); ax.spines["bottom"].set_visible(False)
        ax.set_xlabel("specific super-enhancers (FDR ≤ 0.10)")
    else:
        ax.set_xlabel("specific SEs (FDR ≤ 0.10)"); hgrid(ax)
        ax.set_xlim(0, ax.get_xlim()[1] * 1.3)
        ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{int(v):,}"))

    rng = np.random.default_rng(1)
    dele = (A.kind == "rescued") & (A.cn_mean < DEL_CN)
    kinds = [(A[(A.kind == "rescued") & ~dele], TEAL, 0.0), (A[A.kind == "amplicon_driven"], RED, 1.0)]
    for d, col, yb in kinds:
        v = d["cn_mean"]
        bx.scatter(v, yb + rng.uniform(-0.3, 0.3, len(v)), s=10, color=col, alpha=0.35, edgecolor="none", zorder=3)
    v = A[dele]["cn_mean"].clip(lower=0.05)          # deletion "rescues": open circles on grey, not counted
    bx.scatter(v, 0.0 + rng.uniform(-0.3, 0.3, len(v)), s=18, facecolor="none", edgecolor=MUTED, lw=0.9, zorder=4)
    bx.axvspan(0.045, DEL_CN, color=GRID, alpha=0.55, zorder=0, lw=0)
    bx.axvline(1.0, color=FAINT, lw=1, zorder=1)
    bx.set_xscale("log"); bx.set_xlim(0.045, 150)
    bx.set_xticks([0.1, 0.3, 1, 3, 10, 30, 100], ["0.1", "0.3", "1", "3", "10", "30", "100"])
    bx.set_yticks([0, 1], ["added", "removed"], fontsize=16)
    for t, c in zip(bx.get_yticklabels(), (TEAL, RED)):
        t.set_color(c); t.set_fontweight("bold")
    bx.spines["left"].set_visible(False)
    # the far right of "removed": the highest-CN amplicons, named (as in the first version of this panel)
    amp = (A[A.kind == "amplicon_driven"].sort_values("cn_mean", ascending=False)
           .drop_duplicates("nearest_gene").head(6).sort_values("cn_mean"))
    for i, r in enumerate(amp.itertuples()):
        bx.annotate(r.nearest_gene, (r.cn_mean, 1.32), xytext=(0, 12 + 20 * (i % 3)), textcoords="offset points",
                    ha="center", fontsize=13, color=INK, fontstyle="italic",
                    arrowprops=dict(arrowstyle="-", color=FAINT, lw=0.8, shrinkA=0, shrinkB=2))
    bx.set_ylim(-0.55, 1.95)
    bx.set_xlabel("mean copy number of the group at the locus")
    hgrid(bx)
    fig.tight_layout(w_pad=4.0)
    save(fig, name)


def fig_cn_ablation_all():
    """fig4 at lineage, primary disease and subtype (the user's choice, 2026-10-01; subtype counts depend on how
    BH is pooled, FINDINGS §37, which the presenter states verbally)."""
    fig_cn_ablation(("OncotreeLineage", "OncotreePrimaryDisease", "OncotreeSubtype"), "fig4b_cn_ablation_levels",
                    base_mm=30, per_level_mm=34)                    # compact: the poster's panel 03 also holds Fig. 7


# pre-specified identity TFs per lineage, from the literature — shown whether or not they pass
IDENTITY = [("Ovary/Fallopian Tube", ["PAX8", "SOX17", "MECOM"]), ("Bowel", ["CDX2", "HNF4A"]),
            ("Breast", ["ESR1", "GATA3", "FOXA1"]), ("Myeloid", ["SPI1", "CEBPA"]),
            ("Lymphoid", ["IKZF1", "PAX5"]), ("Skin", ["SOX10", "MITF"]),
            ("Peripheral Nervous System", ["PHOX2B", "HAND2"]), ("Lung", ["NKX2-1", "ASCL1"]),
            ("Liver", ["HNF1A"]), ("Prostate", ["AR"]), ("Kidney", ["PAX2"]), ("Head and Neck", ["TP63"])]


def identity_table(sc=None, res=None, select="jsd", labels=None, links=None, level="OncotreeLineage", groups=None):
    """Known lineage master TFs: one SE per gene (within 100 kb, the best in the gene's own lineage) and its JSD rank,
    FDR and (with `labels`, a fused_labels.groups table) CN status in every lineage.

    select="jsd" takes the lowest JSD, as the poster did. From v3.1 a lineage is tested only where its own experiments
    call an SE (FDR 1 elsewhere), so the lowest JSD can be an SE the lineage never calls (HNF4A in Bowel: rank 2,
    FDR 1); select="fdr" takes the lowest FDR, ties by JSD, which is the SE a reader finds in the atlas.
    links = {gene: {se: info}}: SEs linked to the gene another way (Hi-C contact, 77_stage_master_tfs.py), candidates
    beside the 100 kb ones; a record whose SE is outside 100 kb carries that info as `link`.
    level/groups: the same at a deeper Oncotree level. groups = [(group, its lineage), ...] are the columns (default:
    the IDENTITY lineages themselves); a gene's own groups (`own`) are its lineage's, and its SE is the best over them."""
    sys.path.insert(0, os.path.join(SECACTS, "cnrose"))
    sys.path.insert(0, SECACTS)
    from secacts_env import cache_path
    from cnrose.cn.depmap import load_gene_coords
    gc = load_gene_coords(None, cache_path=cache_path("gene_coords.GRCh38.106.tsv"))
    sc, res, lev = sc or SC, res or RES, level
    J = pd.read_csv(os.path.join(sc, f"atlas.s3.perm.{lev}.jsd.tsv.gz"), sep="\t", index_col=0)
    F = pd.read_csv(os.path.join(sc, f"atlas.s3.perm.{lev}.fdr.tsv.gz"), sep="\t", index_col=0)
    cat = pd.read_csv(os.path.join(res, "atlas.s3.union_catalog.bed.gz"), sep="\t", header=None,
                      usecols=[0, 1, 2, 3], names=["chrom", "start", "end", "se"]).set_index("se")
    cn = {}
    if labels:
        L = pd.read_csv(labels, sep="\t", usecols=["level", "group", "se", "cn_status"]).query("level == @lev")
        cn = dict(zip(zip(L.group, L.se), L.cn_status))
    groups = [(g, ln) for g, ln in (groups or [(g, g) for g, _ in IDENTITY]) if g in J.columns]
    cols = [g for g, _ in groups]
    key = {"jsd": lambda own: lambda x: J.loc[x, own].min(),
           "fdr": lambda own: lambda x: (F.loc[x, own].min(), J.loc[x, own].min())}[select]
    recs = []
    for grp, genes in IDENTITY:
        own = [c for c, ln in groups if ln == grp]
        if not own:
            continue
        ok = lambda x: x in J.index and np.isfinite(J.loc[x, own].astype(float)).any()       # noqa: E731
        for g in genes:
            if g not in gc:
                continue
            c, s, e = gc[g]
            near = cat[(cat.chrom == c) & (cat.end >= s - 100_000) & (cat.start <= e + 100_000)].index   # score_pilot identity_near
            near = [x for x in near if ok(x)]
            win, extra = set(near), (links or {}).get(g, {})
            near += [x for x in extra if x not in win and ok(x)]
            if not near:
                continue
            best = min(near, key=key(own))
            recs.append(dict(gene=g, lineage=grp, own=own, se=best, chrom=cat.loc[best, "chrom"],
                             start=int(cat.loc[best, "start"]), end=int(cat.loc[best, "end"]),
                             rank=[int((J[c] < J.loc[best, c]).sum()) + 1 for c in cols],
                             fdr=F.loc[best, cols].astype(float).values, cn=[cn.get((c, best)) for c in cols],
                             link=None if best in win else extra[best]))
    return cols, recs


def fig_identity(width=215, table=None, name="fig5_identity", cn_mark=False, contrast_text=False, groups=None,
                 fdr_label="FDR (lineage)", cell_fs=10.5, rot=40):
    """Draw identity_table(): -log10 FDR per lineage, the rank where FDR < 1 (* = FDR <= 0.10), the gene's own lineage
    outlined. cn_mark appends † to CN-unmasked calls (pass only with copy-number correction); contrast_text picks each
    label's colour by WCAG contrast (always on for palette variants; the poster default is a fixed threshold). A row
    whose SE is linked by Hi-C rather than within 100 kb is labelled with its distance. groups (deeper levels, as in
    identity_table): columns are bracketed by lineage and the outline spans the gene's own lineage's groups."""
    cols, recs = table or identity_table()
    rows = [-np.log10(pd.Series(r["fdr"]).clip(lower=1e-12)).values for r in recs]
    ranks, fdrs = [r["rank"] for r in recs], [r["fdr"] for r in recs]
    own_j = [[cols.index(o) for o in r.get("own", [r["lineage"]])] for r in recs]
    best_j = [min(oj, key=lambda j: r["fdr"][j]) for oj, r in zip(own_j, recs)]
    lab = [(r["gene"], r["lineage"], r["rank"][j], float(r["fdr"][j])) for r, j in zip(recs, best_j)]
    M = np.array(rows)
    fig, ax = plt.subplots(figsize=(width * MM, 30 * MM + 8.0 * MM * len(rows)))
    from matplotlib.colors import LinearSegmentedColormap
    cmap = LinearSegmentedColormap.from_list("teal", RAMP)
    vmax = 3.0
    im = ax.imshow(np.minimum(M, vmax), cmap=cmap, vmin=0, vmax=vmax, aspect="auto", interpolation="nearest")
    for i, oj in enumerate(own_j):
        j0, j1 = min(oj), max(oj)
        ax.add_patch(plt.Rectangle((j0 - 0.46, i - 0.44), j1 - j0 + 0.92, 0.88, fill=False, edgecolor=INK, lw=1.8, zorder=5))
    for i in range(len(lab)):
        for j in range(len(cols)):
            if fdrs[i][j] < 1.0:
                dag = "†" if cn_mark and fdrs[i][j] <= 0.10 and recs[i]["cn"][j] == "unmasked" else ""
                ax.text(j, i, f"{ranks[i][j]:,}" + ("*" if fdrs[i][j] <= 0.10 else "") + dag, ha="center", va="center", fontsize=cell_fs,
                        fontweight="bold" if fdrs[i][j] <= 0.10 else "normal",
                        color=_cell_text(cmap(min(M[i, j], vmax) / vmax)) if PALETTE or contrast_text else SURF if M[i, j] > 1.4 else INK)
    ax.set_xticks(range(len(cols)), cols, rotation=rot, ha="right", rotation_mode="anchor", fontsize=12.5)
    ax.set_yticks(range(len(lab)), [g + (f" (Hi-C, {r['link']['kb']:.0f} kb)" if r.get("link") else "")
                                     for (g, *_), r in zip(lab, recs)], fontsize=13.5)
    for t, (g, *_ ) in zip(ax.get_yticklabels(), lab):
        t.set_fontstyle("italic")
    ax.tick_params(length=0)
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.set_xticks(np.arange(-0.5, len(cols)), minor=True); ax.set_yticks(np.arange(-0.5, len(lab)), minor=True)
    ax.grid(which="minor", color=SURF, linewidth=1.5); ax.tick_params(which="minor", length=0)
    if groups:                                                      # lineage brackets above the columns
        lin = dict(groups)
        blocks = []
        for j, c in enumerate(cols):
            if blocks and blocks[-1][0] == lin[c]:
                blocks[-1][2] = j
            else:
                blocks.append([lin[c], j, j])
        for ln, j0, j1 in blocks:
            ax.plot([j0 - 0.42, j1 + 0.42], [-0.75, -0.75], color=MUTED, lw=1.2, clip_on=False)
            ax.text((j0 + j1) / 2, -0.95, ln, rotation=90, ha="left", va="center", rotation_mode="anchor",
                    fontsize=12.5, color=INK, clip_on=False)
            if j0 > 0:
                ax.axvline(j0 - 0.5, color=FAINT, lw=1.0, zorder=4)
    cb = fig.colorbar(im, ax=ax, fraction=0.035 * 215 / width, pad=0.02 * 215 / width, ticks=[0, 1, 2, 3])
    cb.ax.set_yticklabels(["1", "0.1", "0.01", "≤0.001"], fontsize=12)
    cb.set_label(fdr_label, fontsize=14); cb.outline.set_visible(False)
    save(fig, name)
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
    order = [p for p in PROVIDERS if p[0] in set(D.cn_source)]
    fig, ax = plt.subplots(figsize=(250 * MM, (90 if len(order) <= 3 else 105) * MM))
    rng = np.random.default_rng(2)
    for i, (k, lab, col) in enumerate(order):
        d = D[D.cn_source == k]
        for j, (c, alpha) in enumerate((("rho_raw", 0.35), ("rho_corrected", 0.9))):
            xs = i * 3 + j + rng.uniform(-0.22, 0.22, len(d))
            ax.scatter(xs, d[c], s=22, color=col if j else FAINT, alpha=alpha, edgecolor="none", zorder=3)
            m = d[c].median()
            ax.plot([i * 3 + j - 0.32, i * 3 + j + 0.32], [m, m], color=INK, lw=2.2, zorder=4)
        ax.text(i * 3 + 0.5, 0.66, f"{lab}\n(n = {len(d)})".replace("inferred from ", "inferred from\n"), ha="center",
                va="top", fontsize=13.5, linespacing=1.1)
    ax.set_xticks([i * 3 + j for i in range(len(order)) for j in range(2)], ["raw", "corr."] * len(order), fontsize=13.5)
    ax.axhline(0, color=FAINT, lw=1)
    ax.set_ylabel("per-line Spearman ρ,\nSE signal vs copy number")
    ax.set_ylim(-0.35, 0.68 if len(order) > 3 else 0.5); vgrid(ax); ax.spines["bottom"].set_visible(False)
    save(fig, "fig7_cn_sources")


PANELS = {"expansion": fig_panel_expansion, "resolution": fig_resolution, "calibration": fig_calibration,
          "ablation": fig_cn_ablation, "ablation_levels": fig_cn_ablation_all, "identity": fig_identity, "concordance": fig_concordance,
          "cnsources": fig_cn_sources}

if __name__ == "__main__":
    want = sys.argv[1:] or list(PANELS)
    for k in want:
        PANELS[k]()
