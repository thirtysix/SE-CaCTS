"""Per-sample pipeline: bigWig + constituent peaks -> SE calls (+ grid quantification, archive, QC).

This is the unit the SLURM array / local loop maps over (DESIGN.md §9). Because the pull DELETES each
bigWig after this function returns (PULL_DESIGN.md §8.5), this is the only chance to emit anything that
needs it: multiple exact grid columns, the genome-wide binned archive, and the QC numbers.
"""
from __future__ import annotations

import json
import os

import numpy as np

from .io import scan_bigwig, read_bed3, write_f32, write_binned_archive, qc_stats
from .stitch import stitch
from .callsuper import call_super
from .cn.base import correct


def _rank(signal):
    """ROSE enhancerRank: 1 = highest signal."""
    order = np.argsort(-signal, kind="stable")
    rank = np.empty(len(signal), dtype=int)
    rank[order] = np.arange(1, len(signal) + 1)
    return rank


def _write_catalog(out_prefix, suffix, regions, signal, is_super, extra_cols=None):
    """Write {out_prefix}{suffix}.se.bed + .enhancers.tsv. extra_cols: list of (header, values)."""
    rank = _rank(signal)
    ex = extra_cols or []
    with open(f"{out_prefix}{suffix}.enhancers.tsv", "w") as fh:
        fh.write("CHROM\tSTART\tSTOP\tNUM_LOCI\tCONSTITUENT_SIZE\tSIGNAL"
                 + "".join(f"\t{h}" for h, _ in ex) + "\tenhancerRank\tisSuper\n")
        for i, (r, sig, rk, sup) in enumerate(zip(regions, signal, rank, is_super)):
            extra = "".join(f"\t{vals[i]:.6g}" for _, vals in ex)
            fh.write(f"{r['chrom']}\t{r['start']}\t{r['end']}\t{r['num_loci']}\t"
                     f"{r['constituent_size']}\t{sig:.6g}{extra}\t{rk}\t{sup}\n")
    n_super = 0
    with open(f"{out_prefix}{suffix}.se.bed", "w") as fh:
        for i, (r, sig, rk, sup) in enumerate(zip(regions, signal, rank, is_super)):
            if sup:
                n_super += 1
                fh.write(f"{r['chrom']}\t{r['start']}\t{r['end']}\tSE_{i}\t{sig:.6g}\t{rk}\n")
    return n_super


def _grid_label(path):
    """Stable output label from a grid filename: .../grid.20.bed -> 'grid.20'."""
    b = os.path.basename(path)
    return b[:-4] if b.endswith(".bed") else b


def call_sample(bw_path, peaks_bed, out_prefix, grid_bed=None, window=12500, agg="sum",
                canonical_only=True, cn_provider=None, cn_key=None, correct_at="none",
                cn_model="log2offset", cn_agg="wlen", cn_beta=1.0, cn_floor=1.0,
                signal_format="tsv", bins=None, qc=False, chunk=5_000_000):
    """Call super-enhancers for one sample and write outputs.

    Always writes the AGNOSTIC catalog:
      {out}.enhancers.tsv / {out}.se.bed   — CN-agnostic stitched regions + SEs
      {out}.signal.tsv                     — grid quantification (uncorrected), if grid_bed given

    `grid_bed` accepts a list: every grid is quantified in the same pass while the bigWig is open, which is
    free against a download-bound pipeline and keeps a looser-threshold future EXACT (PULL_DESIGN.md §8.5.2).
    `signal_format='f32'` writes {out}.<gridlabel>.f32 (7.5x smaller, exact); 'tsv' keeps the historical
    .signal.tsv for the first grid. `bins=100` writes the genome-wide archive; `qc=True` writes {out}.qc.json.
    If a cn_provider + cn_key are given and correct_at=='calling', ALSO writes the CN-CORRECTED catalog
    (LILY-style pre-call correction → a different SE set; DESIGN.md §3.2, dual catalogs):
      {out}.cn.enhancers.tsv / {out}.cn.se.bed   (with REGION_CN + corrected SIGNAL columns)
    The grid fine matrix stays UNCORRECTED (CN is separable → applied at scoring downstream, DESIGN.md §3.2).

    cn_floor: clamp CN to >= cn_floor before correcting. Default 1.0 = AMPLIFY-ONLY (only demote amplified
    regions; leave deletions untouched). This matters at CALLING time: symmetric ÷ratio boosts deleted-region
    signal into false SEs (rehearsal finding, 2026-07-19). Symmetric per-copy correction (cn_floor≈0.1) is for
    the SCORING stage, not calling.
    """
    peaks = read_bed3(peaks_bed, canonical_only=canonical_only)
    regions = stitch(peaks, window=window)
    coords = [(r["chrom"], r["start"], r["end"]) for r in regions]

    # ONE pass over the bigWig for everything it is needed for (PULL_DESIGN.md §8.5): the calling signal,
    # every grid column, and the genome-wide archive. bw.intervals() dominates the cost, so re-reading per
    # artifact would triple it. QC bins ride along at the archive resolution, or 1 kb if no archive.
    grid_paths = ([grid_bed] if isinstance(grid_bed, str) else list(grid_bed or []))
    grids = [read_bed3(g, canonical_only=canonical_only) for g in grid_paths]
    qc_bin = bins or (1000 if qc else None)
    sigs, archive, layout = scan_bigwig(bw_path, [coords] + grids, binsize=qc_bin, agg=agg,
                                        canonical_only=canonical_only, chunk=chunk)
    signal, grid_sigs = sigs[0], sigs[1:]

    os.makedirs(os.path.dirname(os.path.abspath(out_prefix)) or ".", exist_ok=True)

    cutoff, is_super = call_super(signal)
    n_super = _write_catalog(out_prefix, "", regions, signal, is_super)
    out = {"n_regions": len(regions), "n_super": n_super, "cutoff": cutoff}

    if cn_provider is not None and correct_at == "calling":
        track = cn_provider.track(cn_key) if cn_key else None
        if track is None:
            out["cn_status"] = f"no CN for key {cn_key!r} ({cn_provider.name}); corrected catalog skipped"
        else:
            region_cn = np.array([track.region_cn(c, s, e, agg=cn_agg) for c, s, e in coords])
            csig = correct(signal, region_cn, model=cn_model, beta=cn_beta, floor=cn_floor)
            cn_cutoff, is_super_cn = call_super(csig)
            n_super_cn = _write_catalog(out_prefix, ".cn", regions, csig, is_super_cn,
                                        extra_cols=[("REGION_CN", region_cn)])
            out.update(n_super_cn=n_super_cn, cn_cutoff=cn_cutoff, cn_source=track.source,
                       cn_only_agnostic=int((is_super & ~is_super_cn).sum()),
                       cn_only_corrected=int((~is_super & is_super_cn).sum()))

    # ---- durable outputs (PULL_DESIGN.md §8.5)
    first_grid = grids[0] if grids else None
    first_gsig = grid_sigs[0] if grid_sigs else None
    for gi, (gpath, grid, gsig) in enumerate(zip(grid_paths, grids, grid_sigs)):
        label = _grid_label(gpath)
        if signal_format in ("f32", "both"):
            write_f32(f"{out_prefix}.{label}.f32", gsig)
        if signal_format in ("tsv", "both"):
            # the primary grid keeps the historical .signal.tsv name (aggregate.py reads it)
            path = out_prefix + ".signal.tsv" if gi == 0 else f"{out_prefix}.{label}.signal.tsv"
            with open(path, "w") as fh:
                fh.write("CHROM\tSTART\tSTOP\tSIGNAL\n")
                for (c, s, e), v in zip(grid, gsig):
                    fh.write(f"{c}\t{s}\t{e}\t{v:.6g}\n")
        out.setdefault("grids", {})[label] = len(grid)

    if bins and archive is not None:
        meta = write_binned_archive(archive, layout, f"{out_prefix}.bin{bins}.f16.gz", binsize=bins)
        out["archive_bins"] = meta["total_bins"]

    if qc and archive is not None:
        stats = qc_stats(bw_path, archive, qc_bin, grid_regions=first_grid, grid_signal=first_gsig,
                         layout=layout)
        stats.update(n_peaks=len(peaks), n_regions=len(regions), n_super=n_super, tangent_cutoff=cutoff)
        if "n_super_cn" in out:
            stats.update(n_super_cn=out["n_super_cn"], cn_cutoff=out["cn_cutoff"])
        with open(out_prefix + ".qc.json", "w") as fh:
            json.dump(stats, fh, indent=1)
        out["qc"] = {k: stats[k] for k in ("dynamic_range_p99_over_median", "frip_proxy") if k in stats}

    return out
