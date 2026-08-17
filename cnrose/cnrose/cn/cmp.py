"""Cell Model Passports gene-level CN backend (DESIGN.md §4, tier 1).

Reads CMP `cnv_abs_copy_number_picnic_20191101.csv` (gene x model, **ABSOLUTE integer** copy number —
2.0 = diploid-neutral, unlike DepMap's *relative* ratio) and places each gene's CN at its locus
(Ensembl GTF) to form a positional CNTrack.

Promoted above Progenetix on 2026-08-15, overriding the reach-based order in DESIGN.md §4. Progenetix
turned out to carry only categorical CN (four EFO states, no continuous value), and mixing a coarse
categorical correction for some lines with DepMap's continuous correction for others would make
correction strength vary by CN provider — which correlates with how well-studied a line is, hence with
lineage. In a statistic whose whole purpose is detecting lineage specificity, that is the artifact the
correction exists to remove. CMP is continuous, so it keeps correction strength comparable.

Absolute -> relative: CNTrack wants a linear ratio vs neutral (1.0 = neutral). Default `neutral="median"`
divides by the line's own median gene CN, matching what DepMapGeneCN does with `recenter=True`, so the
two backends produce comparably-scaled tracks. `neutral="ploidy"` uses CMP's `ploidy_snp6` instead
(PICNIC is the SNP6 pipeline, so that is the matching ploidy estimate) — truer per-line biology, but not
scale-matched to the DepMap arm; prefer it only if every line in a scoring run comes from CMP.

Join: CMP `model_list` carries `RRID` (a CVCL accession) and `BROAD_ID` (an ACH-* DepMap ModelID)
directly, so no Cellosaurus DR hop is needed — `track()` accepts a CVCL, an ACH-* ModelID or a raw
SIDM and resolves internally.
"""
from __future__ import annotations

import csv
import gzip
import os

import numpy as np

from .base import CNTrack, CNProvider

# The PICNIC matrix carries three header rows before the data:
#   0: model_id,,SIDM00003,SIDM00023,...
#   1: model_name,,M14,TE-12,...
#   2: gene_id,symbol,,,...
# then gene rows: SIDG00001,A1BG,3.0,3.0,...
_N_HEADER_ROWS = 3
_COL_GENE_ID = 0
_COL_SYMBOL = 1


def load_model_index(model_list_csv):
    """{alias: SIDM} from CMP `model_list`, covering CVCL (RRID), ACH-* (BROAD_ID) and SIDM itself."""
    index = {}
    ploidy = {}
    with open(model_list_csv, newline="") as fh:
        for row in csv.DictReader(fh):
            sidm = (row.get("model_id") or "").strip()
            if not sidm:
                continue
            index[sidm] = sidm
            for col in ("RRID", "BROAD_ID"):
                v = (row.get(col) or "").strip()
                if v and v not in ("NA", "nan"):
                    index[v] = sidm
            p = (row.get("ploidy_snp6") or "").strip()
            if p and p not in ("NA", "nan"):
                try:
                    ploidy[sidm] = float(p)
                except ValueError:
                    pass
    return index, ploidy


class CellModelPassportsCN(CNProvider):
    """CNProvider backed by CMP PICNIC absolute CN + gene coordinates."""

    def __init__(self, cn_csv, model_list_csv, gene_coords, neutral="median"):
        if neutral not in ("median", "ploidy"):
            raise ValueError(f"neutral must be 'median' or 'ploidy', got {neutral!r}")
        self.cn_csv = cn_csv
        self.coords = gene_coords
        self.neutral = neutral
        self.alias, self.ploidy = load_model_index(model_list_csv)

        with open(cn_csv, newline="") as fh:
            r = csv.reader(fh)
            header = next(r)
        self.sidm_col = {s: i for i, s in enumerate(header) if s.startswith("SIDM")}
        self._cache = {}

    def resolve(self, key):
        """CVCL / ACH-* / SIDM -> a SIDM that actually has a CN column, or None."""
        sidm = self.alias.get(str(key).strip())
        return sidm if sidm in self.sidm_col else None

    def preload(self, keys):
        """One file pass for many lines — the matrix is gene x model, so a pass costs the same for 1 or N."""
        want = {}
        for k in keys:
            s = self.resolve(k)
            if s and s not in self._cache:
                want[s] = self.sidm_col[s]
        if not want:
            return
        cols = list(want.values())
        vals = {s: [] for s in want}
        genes = []
        with open(self.cn_csv, newline="") as fh:
            r = csv.reader(fh)
            for _ in range(_N_HEADER_ROWS):
                next(r)
            for row in r:
                sym = row[_COL_SYMBOL]
                if sym not in self.coords:
                    continue
                genes.append(sym)
                for s, ci in want.items():
                    vals[s].append(row[ci] if ci < len(row) else "")
        for s in want:
            self._cache[s] = self._build(genes, vals[s], s)

    def track(self, key):
        sidm = self.resolve(key)
        if sidm is None:
            return None
        if sidm not in self._cache:
            self.preload([sidm])
        return self._cache.get(sidm)

    def _build(self, genes, values, sidm):
        segs, ratios = [], []
        for sym, v in zip(genes, values):
            if v in ("", "NA", "nan"):
                continue
            try:
                cn = float(v)
            except ValueError:
                continue
            c, s, e = self.coords[sym]
            segs.append((c, s, e, cn))
            ratios.append(cn)
        if not segs:
            return None

        if self.neutral == "ploidy":
            denom = self.ploidy.get(sidm)
            if not denom or denom <= 0:                 # no ploidy call -> fall back rather than skip
                denom = float(np.median(ratios))
        else:
            denom = float(np.median(ratios))
        if not denom or denom <= 0:
            return None

        segs = [(c, s, e, cn / denom) for c, s, e, cn in segs]
        return CNTrack(segs, ploidy=self.ploidy.get(sidm, 2.0), source="CMPPicnicAbsolute")


class CellModelPassportsWesCN(CNProvider):
    """CMP WES pureCN gene-level CN (2025-02-07 release) — the preferred CMP backend.

    Chosen over the 2019 PICNIC matrix on 2026-08-15 for one reason above all: **PICNIC is integer-
    quantised**. Sampled over 231,876 values it emits only whole numbers 0-14, so after recentering on a
    line's median CN (typically 3) the finest representable step is ~0.33 — a 1.1x or 1.2x gain simply
    cannot be expressed. DepMap's relative CN is continuous, so PICNIC lines would receive a coarser
    correction than DepMap lines, and since CN provider correlates with how well-studied a line is, and
    hence with lineage, that is the same confound that got Progenetix demoted, just milder. This file's
    `gatk_mean_log2_copy_ratio` is continuous (1,012 distinct values per 60k rows sampled).

    Secondary wins: coordinates travel with the CN (no GTF join, so no genome-build mismatch between the
    CN source and the gene model), 1,253 models vs 986, and allele-specific `minor_copy_number`/`loh`.

    Caveats handled here:
      * `gatk_mean_log2_copy_ratio` reaches -29 in zero-coverage regions; left raw, one such gene would
        dominate a length-weighted mean. Clipped to +/-`log2_clip` (default 4 => ratio 0.0625..16).
      * ~4.8% of (model, gene) pairs appear more than once; they are averaged in log2 space rather than
        emitted as duplicate overlapping intervals, which would double-count in `region_cn`.
      * `total_copy_number` is missing on ~1.8% of rows; the log2 ratio never is, so it is the value used.

    The file is ~980 MB gzipped and long-format, so a scan is expensive. `preload()` does ONE pass for all
    requested lines and writes a per-model cache; later runs read the cache.
    """

    def __init__(self, cn_csv_gz, model_list_csv, cache_dir=None, log2_clip=4.0, recenter=True):
        self.cn_csv_gz = cn_csv_gz
        self.cache_dir = cache_dir
        self.log2_clip = float(log2_clip)
        self.recenter = recenter
        self.alias, self.ploidy = load_model_index(model_list_csv)
        self._cache = {}
        if cache_dir:
            os.makedirs(cache_dir, exist_ok=True)

    def resolve(self, key):
        return self.alias.get(str(key).strip())

    def _cache_file(self, sidm):
        return os.path.join(self.cache_dir, f"{sidm}.tsv.gz") if self.cache_dir else None

    def preload(self, keys):
        want = {s for s in (self.resolve(k) for k in keys) if s and s not in self._cache}
        if self.cache_dir:                                   # satisfy what we can from disk first
            for s in sorted(want):
                fn = self._cache_file(s)
                if os.path.exists(fn):
                    self._cache[s] = self._from_rows(_read_cache(fn), s)
            want = {s for s in want if s not in self._cache}
        if not want:
            return

        acc = {s: {} for s in want}                          # sidm -> (chrom,start,end,symbol) -> [log2...]
        with gzip.open(self.cn_csv_gz, "rt", newline="") as fh:
            for row in csv.DictReader(fh):
                s = row.get("model_id")
                if s not in acc:
                    continue
                lr = row.get("gatk_mean_log2_copy_ratio")
                chrom, a, b = row.get("chr_name"), row.get("chr_start"), row.get("chr_end")
                if not lr or not chrom or not a or not b:
                    continue
                try:
                    key = (chrom if chrom.startswith("chr") else "chr" + chrom, int(a), int(b))
                    acc[s].setdefault(key, []).append(float(lr))
                except ValueError:
                    continue

        for s in want:
            rows = [(c, a, b, float(np.mean(v))) for (c, a, b), v in acc[s].items()]   # dedup by averaging
            if self.cache_dir and rows:
                _write_cache(self._cache_file(s), rows)
            self._cache[s] = self._from_rows(rows, s)

    def track(self, key):
        sidm = self.resolve(key)
        if sidm is None:
            return None
        if sidm not in self._cache:
            self.preload([sidm])
        return self._cache.get(sidm)

    def _from_rows(self, rows, sidm):
        if not rows:
            return None
        clip = self.log2_clip
        segs = [(c, a, b, float(2.0 ** min(max(lr, -clip), clip))) for c, a, b, lr in rows]
        if self.recenter:
            med = float(np.median([r for *_, r in segs]))
            if med > 0:
                segs = [(c, a, b, r / med) for c, a, b, r in segs]
        return CNTrack(segs, ploidy=self.ploidy.get(sidm, 2.0), source="CMPWesPureCN")


def _write_cache(path, rows):
    with gzip.open(path, "wt", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerows(rows)


def _read_cache(path):
    with gzip.open(path, "rt", newline="") as fh:
        return [(c, int(a), int(b), float(lr)) for c, a, b, lr in csv.reader(fh, delimiter="\t")]
