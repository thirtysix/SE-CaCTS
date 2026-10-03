"""CN inferred from matched ChIP input DNA (DESIGN.md §4, tier 2 — "the actual modeling task").

For lines no measured-CN resource covers, the input-control track is a copy-number readout: input DNA is
nominally unenriched genomic DNA, so binned coverage tracks local copy number.

**Everything is read REMOTELY via bigWig zoom levels.** bigWigs store precomputed multi-resolution
summaries; pyBigWig serves `stats(..., nBins=N)` out of them over HTTP range requests. Measured
2026-08-15 against ChIP-Atlas: a whole genome at 50 kb costs ~43 s and a few MB, versus 271-1615 MB to
download the file. Same trick supplies GC from UCSC's gc5Base. For 157 lines that is ~2 h of querying
instead of a ~120 GB pull — no CSC job needed. (GOTCHA 13's "~90-570 MB each" both understates current
sizes and implies a download that binned work does not require.)

A naive first pilot (bin + median-normalise, nothing else) scored median Pearson r=0.09 against DepMap
WGS while Spearman was 0.54-0.67 — the CN signal was present but swamped by artifact bins: A-375 showed
an inferred 197x bin where its true maximum is 1.4x. The three corrections every real tool applies
(GC, blacklist, outlier handling) are what close that gap, and each is separately switchable here so its
contribution can be measured rather than assumed.

Note on winsorisation: it is in direct tension with focal-amplicon detection, which is the whole point
of the correction. A MYCN amplicon spans ~6 of 55,000 bins (0.01%), so clipping at the 99th or even
99.9th percentile would erase it. `winsor_hi` therefore defaults to 99.99 and the blacklist is expected
to do most of the artifact removal.
"""
from __future__ import annotations

import gzip

import numpy as np

from .base import CNTrack, CNProvider

CANON = [f"chr{i}" for i in range(1, 23)]
GC_BIGWIG = "https://hgdownload.soe.ucsc.edu/gbdb/hg38/bbi/gc5BaseBw/gc5Base.bw"
CHIP_ATLAS_BW = "https://chip-atlas.dbcls.jp/data/hg38/eachData/bw/{srx}.bw"


def load_blacklist(path):
    """ENCODE hg38 blacklist v2 BED (plain or .gz) -> {chrom: [(start, end), ...]}."""
    op = gzip.open if str(path).endswith(".gz") else open
    out = {}
    with op(path, "rt") as fh:
        for line in fh:
            f = line.split("\t")
            if len(f) >= 3:
                out.setdefault(f[0], []).append((int(f[1]), int(f[2])))
    return out


def bin_grid(chrom_lengths, chroms=CANON, bin_size=50_000):
    """{chrom: (starts, ends)} on a fixed grid — shared by the input track, GC and the blacklist."""
    grid = {}
    for c in chroms:
        n = int(chrom_lengths.get(c, 0) // bin_size)
        if n < 4:
            continue
        s = np.arange(n, dtype=np.int64) * bin_size
        grid[c] = (s, s + bin_size)
    return grid


def blacklist_mask(grid, blacklist):
    """{chrom: bool array} — True where a bin overlaps any blacklisted interval."""
    mask = {}
    for c, (st, en) in grid.items():
        m = np.zeros(len(st), dtype=bool)
        for a, b in blacklist.get(c, []):
            m |= (st < b) & (en > a)
        mask[c] = m
    return mask


def query_binned(url, grid, retries=2, mode="covered"):
    """{chrom: float array} mean-over-SPAN per bin, from zoom levels.

    `value_mode="covered"` (DEFAULT, and MEASURED BEST) uses `type="mean"`; "span" uses `sum`/bin_size. pyBigWig's mean averages over **covered bases
    only** (see project memory `bigwig-mean-is-covered-bases-only`), and ChIP-Atlas input bigWigs are
    sparse: measured median covered fraction 0.445, with 199 of 200 bins under 90% covered. Under the
    mean convention a half-covered bin reads at full depth, overstating it ~2.2x on average and up to
    6.1x — and the distortion scales with local coverage sparsity, so it injects noise rather than a
    constant factor. For copy number the fixed-lattice convention (uncovered = 0, i.e. mean over span)
    is the correct one, exactly as the genome-wide archive does it.
    """
    import pyBigWig
    bw = None
    for _ in range(retries + 1):
        try:
            bw = pyBigWig.open(url)
            break
        except Exception:
            bw = None
    if bw is None:
        return None
    out = {}
    try:
        have = bw.chroms()
        for c, (st, en) in grid.items():
            if c not in have:
                continue
            if mode == "span":
                span = int(en[0] - st[0])
                raw = bw.stats(c, int(st[0]), int(en[-1]), nBins=len(st), type="sum")
                v = np.array([x if x is not None else np.nan for x in raw], dtype=float) / span
            else:
                raw = bw.stats(c, int(st[0]), int(en[-1]), nBins=len(st), type="mean")
                v = np.array([x if x is not None else np.nan for x in raw], dtype=float)
            out[c] = v
    finally:
        bw.close()
    return out or None


def gc_correct(values, gc, n_strata=20, min_per_stratum=20):
    """Divide each bin by the median of bins sharing its GC stratum.

    Median-based rather than a LOWESS fit: monotone-safe, no fitting failure mode on sparse data, and
    robust to the artifact bins this runs alongside.

    NaN-GC bins are left UNCORRECTED, not lumped together. Casting NaN to int yields INT_MIN, which
    np.clip then folds into stratum 0 — so ~4.7% of bins genome-wide (11% on chr21, where the GC track
    has gaps) were being divided by one arbitrary stratum's median. That bug made GC correction *harm*
    the result: median Pearson fell 0.682 -> 0.472 in the 2026-08-15 pilot.
    """
    v = np.asarray(values, float)
    g = np.asarray(gc, float)
    ok = np.isfinite(v) & np.isfinite(g)
    if ok.sum() < 100:
        return v
    lo, hi = np.nanmin(g[ok]), np.nanmax(g[ok])
    strata = np.full(g.shape, -1, dtype=np.int64)          # -1 == "no GC, do not correct"
    strata[ok] = np.clip(((g[ok] - lo) / (hi - lo + 1e-9) * n_strata).astype(np.int64), 0, n_strata - 1)
    out = v.copy()
    for s in range(n_strata):
        m = ok & (strata == s)
        if m.sum() >= min_per_stratum:
            med = np.nanmedian(v[m])
            if med > 0:
                out[m] = v[m] / med
    return out


def _best_block(x, widths, min_size):
    """Best contiguous BLOCK vs the rest, by drop in squared error (circular binary segmentation).

    Plain binary segmentation cannot find short segments: isolating a 6-bin amplicon inside 2,000 bins
    scores a first-split gain of ~0.3 against a penalty of ~1.4, so nothing splits and the amplicon is
    erased. That is exactly why Olshen & Venkatraman's CBS tests a block (two changepoints) rather than
    one. Scanning a geometric ladder of widths keeps it O(n x |widths|) instead of CBS's O(n^2).
    """
    n = len(x)
    if n < 2 * min_size:
        return None, 0.0
    cs = np.concatenate([[0.0], np.cumsum(x)])
    total = cs[-1]
    best = (None, 0.0)
    for w in widths:
        if w < min_size or w > n - min_size:
            continue
        sb = cs[w:] - cs[:-w]                       # sum of every window of width w
        gain = sb ** 2 / w + (total - sb) ** 2 / (n - w) - total ** 2 / n
        i = int(np.argmax(gain))
        if gain[i] > best[1]:
            best = ((i, i + w), float(gain[i]))
    return best


def segment(values, min_size=2, beta=3.0, max_depth=12, max_widths=24):
    """Piecewise-constant segmentation of a log-ratio vector; returns per-bin segment MEDIANS.

    Circular binary segmentation over a geometric width ladder, penalty BIC-flavoured with the noise
    scale estimated robustly from first differences (MAD). This is the step HMCan supplies with an HMM
    and DNAcopy with CBS; without it the track is raw per-bin coverage, and that bin-level noise is what
    makes per-line agreement swing from 0.15 to 0.96.

    `min_size=2` is deliberate: a MYCN amplicon occupies ~6 of 55,000 bins, so a segmenter that cannot
    return short segments would erase precisely the signal the correction exists to catch.
    """
    x = np.asarray(values, float)
    ok = np.isfinite(x)
    if ok.sum() < 4 * min_size:
        return x
    xf = x[ok]
    d = np.diff(xf)
    sigma = 1.4826 * np.median(np.abs(d - np.median(d))) / np.sqrt(2) if len(d) else 0.0
    if not np.isfinite(sigma) or sigma <= 0:
        sigma = float(np.std(xf)) or 1.0
    penalty = beta * sigma ** 2 * np.log(max(len(xf), 2))

    widths = sorted({int(round(min_size * 1.6 ** i)) for i in range(max_widths)} | {min_size})
    widths = [w for w in widths if w < len(xf)]

    edges = {0, len(xf)}
    todo = [(0, len(xf), 0)]
    while todo:
        a, b, depth = todo.pop()
        if depth >= max_depth or b - a < 2 * min_size:
            continue
        blk, gain = _best_block(xf[a:b], [w for w in widths if w < b - a], min_size)
        if blk is None or gain <= penalty:
            continue
        i, j = a + blk[0], a + blk[1]
        edges.update((i, j))
        todo += [(a, i, depth + 1), (i, j, depth + 1), (j, b, depth + 1)]

    e = sorted(edges)
    out_f = np.empty_like(xf)
    for a, b in zip(e[:-1], e[1:]):
        if b > a:
            out_f[a:b] = np.median(xf[a:b])
    out = x.copy()
    out[ok] = out_f
    return out


class ChipInputInferredCN(CNProvider):
    """CN from a matched ChIP input-control bigWig, with GC / blacklist / winsorisation switchable.

    `srx_for` maps a line key -> an input SRX (prefer the deepest; shallow inputs give noisy bins).
    GC and the bin grid are fetched once per instance and reused across every line.
    """

    def __init__(self, srx_for, blacklist=None, use_gc=False, segment_cn=False, seg_beta=3.0,
                 value_mode="covered",
                 winsor_hi=99.99, winsor_lo=0.0,
                 zero_is_deletion=True, bin_size=50_000, chroms=CANON,
                 url_template=CHIP_ATLAS_BW, gc_url=GC_BIGWIG, min_finite_frac=0.3):
        self.srx_for = dict(srx_for)
        self.blacklist = blacklist or {}
        self.use_gc = use_gc
        self.value_mode = value_mode
        self.segment_cn = segment_cn
        self.seg_beta = seg_beta
        self.winsor_hi = winsor_hi
        self.winsor_lo = winsor_lo
        self.zero_is_deletion = zero_is_deletion
        self.bin_size = bin_size
        self.chroms = list(chroms)
        self.url_template = url_template
        self.gc_url = gc_url
        self.min_finite_frac = min_finite_frac
        self._grid = self._gc = self._bmask = None
        self._cache = {}

    def _ensure_grid(self):
        if self._grid is not None:
            return
        import pyBigWig
        bw = pyBigWig.open(self.gc_url)
        try:
            lengths = bw.chroms()
        finally:
            bw.close()
        self._grid = bin_grid(lengths, self.chroms, self.bin_size)
        self._bmask = blacklist_mask(self._grid, self.blacklist)
        self._gc = query_binned(self.gc_url, self._grid, mode=self.value_mode) if self.use_gc else None

    def raw_bins(self, key):
        """Per-bin mean input coverage for one line, before any correction (for variant comparison)."""
        self._ensure_grid()
        srx = self.srx_for.get(key)
        return query_binned(self.url_template.format(srx=srx), self._grid,
                            mode=self.value_mode) if srx else None

    def build_from_bins(self, cov, srx="?"):
        """Apply blacklist -> GC -> winsorise -> median-normalise to pre-fetched bins."""
        self._ensure_grid()
        segs, pooled = [], []
        per_chrom = {}
        for c, v in cov.items():
            if c not in self._grid:
                continue
            v = v.copy()
            keep = ~self._bmask[c]
            if self.zero_is_deletion:
                v[np.isnan(v) & keep] = 0.0          # after blacklisting, no coverage means depletion
            if self.use_gc and self._gc and c in self._gc:
                v = gc_correct(v, self._gc[c])
            v[~keep] = np.nan
            if self.segment_cn:                  # smooth in log space; CN is multiplicative
                lg = np.log2(np.clip(v, 1e-3, None))
                lg[~np.isfinite(v)] = np.nan
                v = np.exp2(segment(lg, beta=self.seg_beta))
                v[~keep] = np.nan
            per_chrom[c] = v
            pooled.append(v[np.isfinite(v)])
        if not pooled:
            return None
        allv = np.concatenate(pooled)
        if len(allv) < 1000 or np.isfinite(allv).mean() < self.min_finite_frac:
            return None
        hi = np.percentile(allv, self.winsor_hi) if self.winsor_hi < 100 else np.inf
        lo = np.percentile(allv, self.winsor_lo) if self.winsor_lo > 0 else -np.inf
        med = float(np.median(allv))
        if med <= 0:
            return None
        for c, v in per_chrom.items():
            st, en = self._grid[c]
            vv = np.clip(v, lo, hi) / med
            ok = np.isfinite(vv)
            segs.extend(zip([c] * int(ok.sum()), st[ok].tolist(), en[ok].tolist(), vv[ok].tolist()))
        return CNTrack(segs, ploidy=2.0, source=f"ChipInputInferred:{srx}") if segs else None

    def track(self, key):
        if key not in self._cache:
            cov = self.raw_bins(key)
            self._cache[key] = self.build_from_bins(cov, self.srx_for.get(key, "?")) if cov else None
        return self._cache[key]


# ----------------------------------------------------------------------------- admitted lines (FINDINGS §24)
HG38_AUTOSOMES = {
    "chr1": 248956422, "chr2": 242193529, "chr3": 198295559, "chr4": 190214555, "chr5": 181538259,
    "chr6": 170805979, "chr7": 159345973, "chr8": 145138636, "chr9": 138394717, "chr10": 133797422,
    "chr11": 135086622, "chr12": 133275309, "chr13": 114364328, "chr14": 107043718, "chr15": 101991189,
    "chr16": 90338345, "chr17": 83257441, "chr18": 80373285, "chr19": 58617616, "chr20": 64444167,
    "chr21": 46709983, "chr22": 50818468}


class BinnedInputCN(CNProvider):
    """Input-inferred CN from PRE-FETCHED 50 kb bins (phase1/scripts/20_cn_input_calibration.py fetch, one
    `<SRX>.npz` per input), so scoring needs no network. Same inference as the calibration (covered mean,
    zero = deletion, blacklist, no GC, no segmentation), then the global compression slope fitted against
    DepMap WGS: ratio -> ratio ** slope (log2 truth ~ 0.81 x log2 inferred, FINDINGS §24).

    `input_for` maps a scoring key -> the admitted input SRX (the line's most coherent input).

    Also reads measured CN from subsampled WGS (phase1/scripts/27_wgs_cn.py: reads per bin divided by a normal
    reference) with slope=1 and zero_is_deletion=False: there NaN marks a bin the reference cannot map, which must
    stay missing, while a true 0 (reference reads, none in the line) is still a deletion. segment=True smooths the
    50 kb ratios in log space first: on 9 calibration lines it raised agreement with DepMap WGS at SE loci from
    r 0.74 to 0.83 (best ChIP input 0.69) and amplicon recovery from 72% to 83% (input 5%) (FINDINGS §57).
    """

    def __init__(self, bins_dir, input_for, blacklist=None, slope=0.81, bin_size=50_000, zero_is_deletion=True,
                 segment=False):
        self.bins_dir = bins_dir
        self.input_for = dict(input_for)
        self.slope = slope
        grid = {c: (np.arange(int(L // bin_size), dtype=np.int64) * bin_size,
                    np.arange(1, int(L // bin_size) + 1, dtype=np.int64) * bin_size)
                for c, L in HG38_AUTOSOMES.items()}
        self._inf = ChipInputInferredCN({}, blacklist=blacklist or {}, bin_size=bin_size,
                                        zero_is_deletion=zero_is_deletion, segment_cn=segment)
        self._inf._grid, self._inf._bmask, self._inf._gc = grid, blacklist_mask(grid, blacklist or {}), None
        self._cache = {}

    def track(self, key):
        if key in self._cache:
            return self._cache[key]
        import os
        srx = self.input_for.get(key)
        fn = os.path.join(self.bins_dir, f"{srx}.npz") if srx else None
        tr = None
        if fn and os.path.exists(fn):
            z = np.load(fn)
            raw = self._inf.build_from_bins({c: z[c] for c in z.files if c.startswith("chr")}, srx)
            if raw is not None:
                segs = [(c, int(s), int(e), float(r) ** self.slope)
                        for c, (st, en, ra) in raw._chrom.items() for s, e, r in zip(st, en, ra)]
                tr = CNTrack(segs, ploidy=2.0, source=f"InputInferred:{srx}^{self.slope}")
        self._cache[key] = tr
        return tr
