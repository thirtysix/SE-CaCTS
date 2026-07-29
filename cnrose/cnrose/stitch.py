"""Constituent stitching — a faithful port of ROSE2's LocusCollection.stitchCollection.

ROSE2 (utils.py) reads a BED into Locus objects with `Locus.start = BED start`, `Locus.end = BED end`
(bedToGFF copies the coords verbatim — no shift), then single-linkage-merges: each locus is expanded by
±window and any overlapping locus is absorbed, iterating until nothing new overlaps. For loci sorted by
start, `overlaps([cS-W, cE+W], B)` reduces to `B.start <= cE + W` (the right-side gap condition), so the
whole thing is a sorted sweep that merges the next peak whenever `next.start <= cluster_end + window`.
Merged region = [min start, max end] of the constituents. `-t 0` (no TSS handling) is the H3K27ac default
and the only mode implemented here (TSS-break is a stub for later).

Validated against rose2.utils.LocusCollection.stitchCollection(window, 'both') on the 13 pilot peak sets
(see ../tests/validate_vs_rose2.py).
"""
from __future__ import annotations


def stitch(peaks, window=12500, tss=None):
    """Stitch constituent peaks into candidate super-enhancer regions.

    Parameters
    ----------
    peaks : iterable of (chrom, start, end)
        Constituent intervals (BED coords: 0-based start, end-exclusive — used verbatim, as ROSE does).
    window : int
        Max stitch distance in bp (ROSE `-s`, default 12500 for H3K27ac).
    tss : None
        TSS-exclusion is not implemented (H3K27ac SE calling uses `-t 0`); passing a value raises.

    Returns
    -------
    list of dict: {chrom, start, end, num_loci, constituent_size}
        constituent_size = summed bp of the constituent peaks (ROSE's CONSTITUENT_SIZE, col 6).
    """
    if tss:
        raise NotImplementedError("TSS-exclusion (-t > 0) is not implemented; H3K27ac SE calling uses -t 0")

    by_chrom = {}
    for c, s, e in peaks:
        s, e = int(s), int(e)
        if e < s:
            s, e = e, s
        by_chrom.setdefault(c, []).append((s, e))

    out = []
    for chrom in sorted(by_chrom):
        ivs = sorted(by_chrom[chrom])          # by start, then end
        cs, ce = ivs[0]
        const = [(cs, ce)]
        for s, e in ivs[1:]:
            if s <= ce + window:               # ROSE merge condition (B.start <= clusterEnd + W)
                if e > ce:
                    ce = e
                const.append((s, e))
            else:
                out.append(_region(chrom, cs, ce, const))
                cs, ce, const = s, e, [(s, e)]
        out.append(_region(chrom, cs, ce, const))
    return out


def _region(chrom, start, end, const):
    return {
        "chrom": chrom,
        "start": start,
        "end": end,
        "num_loci": len(const),
        "constituent_size": sum(e - s for s, e in const),
    }
