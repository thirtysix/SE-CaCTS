"""CN abstractions — the two pieces that make copy-number source-agnostic and separable (DESIGN.md §3).

  CNTrack   : a canonical per-line CN representation — segments carrying a LINEAR ratio vs. neutral
              (1.0 = neutral), plus ploidy. `region_cn` collapses an SE spanning heterogeneous segments
              to one number (fork C: length-weighted mean over the SE's overlap with each segment).
  CNProvider: the source-agnostic interface — `track(key) -> CNTrack | None`. Backends (DepMap, Progenetix,
              CMP, input-inference, …) each just emit a CNTrack; new sources are additive (DESIGN.md §4).
  correct   : the separable transform — one function called at calling-time AND scoring-time (DESIGN.md §3.2).
              log2-offset is primary (β free, matches the Phase-4 regression); ÷ratio is the cross-check.
"""
from __future__ import annotations

from abc import ABC, abstractmethod

import numpy as np


class CNTrack:
    """Segment-level copy number for one cell line: intervals -> linear ratio vs neutral (1.0 = neutral)."""

    def __init__(self, segments, ploidy=2.0, source="unknown"):
        by = {}
        for c, s, e, r in segments:
            by.setdefault(c, []).append((int(s), int(e), float(r)))
        self._chrom = {}
        for c, ivs in by.items():
            ivs.sort()
            self._chrom[c] = (
                np.array([i[0] for i in ivs], dtype=np.int64),
                np.array([i[1] for i in ivs], dtype=np.int64),
                np.array([i[2] for i in ivs], dtype=np.float64),
            )
        self.ploidy = float(ploidy)
        self.source = source

    def __len__(self):
        return sum(len(v[0]) for v in self._chrom.values())

    def region_cn(self, chrom, start, end, agg="wlen", max_nearest=3_000_000, max_seg_len=3_000_000):
        """One CN ratio for [start, end). agg='wlen' (length-weighted mean over segment overlaps),
        'max' (conservative amplicon-killer), or 'median'. Falls back to the nearest segment within
        max_nearest bp (CN is locally smooth), else 1.0 (neutral) when a chrom/region has no CN.

        Overlap-robust: intervals are sorted by start but may overlap (gene-level tracks have nested
        genes, so `end` is NOT monotonic). We bound the candidate window by start position
        (start-max_seg_len .. end) and filter by actual overlap, rather than searchsorting on `end`.
        """
        v = self._chrom.get(chrom)
        if v is None:
            return 1.0
        st, en, ra = v
        hi = np.searchsorted(st, end, side="left")                    # segs with start < end
        lo = np.searchsorted(st, start - max_seg_len, side="left")    # window floor (bounds the scan)
        if hi > lo:
            idx = np.arange(lo, hi)
            ov = np.minimum(en[idx], end) - np.maximum(st[idx], start)
            m = ov > 0
            if m.any():
                ov, rr = ov[m], ra[idx][m]
                if agg == "max":
                    return float(rr.max())
                if agg == "median":
                    return float(np.median(rr))
                return float((rr * ov).sum() / ov.sum())               # wlen (default)
        # no overlap -> nearest segment by start position
        gaps = []
        if hi - 1 >= 0:
            gaps.append((abs(start - en[hi - 1]), ra[hi - 1]))
        if hi < len(st):
            gaps.append((abs(st[hi] - end), ra[hi]))
        gaps = [(g, r) for g, r in gaps if g <= max_nearest]
        return float(min(gaps)[1]) if gaps else 1.0


class CNProvider(ABC):
    """Source-agnostic CN lookup. Subclasses resolve the join and emit a canonical CNTrack."""

    @abstractmethod
    def track(self, key):
        """Return a CNTrack for the cell line identified by `key`, or None if this source has no CN for it."""
        raise NotImplementedError

    @property
    def name(self):
        return type(self).__name__


class ChainProvider(CNProvider):
    """Resolve per line in priority order (DESIGN.md §4): first backend with a track wins."""

    def __init__(self, providers):
        self.providers = list(providers)

    def track(self, key):
        for p in self.providers:
            t = p.track(key)
            if t is not None and len(t):
                return t
        return None


def correct(signal, cn, model="log2offset", beta=1.0, floor=0.1, eps=1.0):
    """Separable CN correction of a signal vector (DESIGN.md §6.1).

    model='log2offset' (primary): 2^(log2(signal+eps) - beta*log2(cn)) - eps  — subtract a CN offset in
        log space; beta is a free covariate weight (fixed 1.0 at calling-time; fit by regression at scoring).
    model='divide'    (cross-check): signal / cn**beta  — the naive Su/Chen ÷CN-ratio form.
    `cn` is a linear ratio vs neutral (1.0 = neutral), floored at `floor` to avoid blow-ups on deep deletions.
    At beta=1 the two models coincide (up to the eps regularisation) — expected; the offset's advantage is at
    scoring, where beta is estimated. Returns corrected signal in the same (linear) space as the input.
    """
    sig = np.asarray(signal, dtype=float)
    cnr = np.clip(np.asarray(cn, dtype=float), floor, None)
    if model == "divide":
        return sig / np.power(cnr, beta)
    if model == "log2offset":
        return np.exp2(np.log2(sig + eps) - beta * np.log2(cnr)) - eps
    raise ValueError(f"unknown correction model: {model!r}")
