"""DepMap gene-level CN backend (DESIGN.md §4, tier 1).

Reads DepMap `OmicsCNGeneWGS.csv` (gene x model, LINEAR relative CN, 1.0 = neutral — verified: OVCAR3
median 0.985, MCF7 1.064) and places each gene's CN at its locus (Ensembl GTF) to form a positional
CNTrack. Join key = DepMap ModelID (col 4; the default profile flagged IsDefaultEntryForModel == "Yes").

This is the first concrete backend; the segment-level `OmicsCNSegmentsWGS.csv` (denser, no intergenic
gaps) is the eventual preferred DepMap source and drops in behind the same CNProvider interface.
"""
from __future__ import annotations

import csv
import gzip
import os

import numpy as np

from .base import CNTrack, CNProvider

CANON = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}


def _norm_chrom(c):
    if c in ("MT", "M"):
        return "chrM"
    return c if c.startswith("chr") else "chr" + c


def load_gene_coords(gtf_path, cache_path=None):
    """{gene_symbol: (chrom, start0, end)} from an Ensembl GTF (feature==gene), canonical chroms, chr-prefixed.

    Coordinates are BED-style (0-based start). Caches to a small TSV for fast reuse.
    """
    if cache_path and os.path.exists(cache_path):
        coords = {}
        with open(cache_path) as fh:
            for line in fh:
                s, c, a, b = line.rstrip("\n").split("\t")
                coords[s] = (c, int(a), int(b))
        return coords

    op = gzip.open if gtf_path.endswith(".gz") else open
    coords = {}
    with op(gtf_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t")
            if len(f) < 9 or f[2] != "gene":
                continue
            chrom = _norm_chrom(f[0])
            if chrom not in CANON:
                continue
            attr = f[8]
            k = attr.find('gene_name "')
            if k < 0:
                continue
            name = attr[k + 11: attr.find('"', k + 11)]
            start, end = int(f[3]) - 1, int(f[4])
            if name in coords:                 # widen to span duplicate/patch entries on the same chrom
                c0, s0, e0 = coords[name]
                if c0 == chrom:
                    coords[name] = (chrom, min(s0, start), max(e0, end))
            else:
                coords[name] = (chrom, start, end)

    if cache_path:
        with open(cache_path, "w") as fh:
            for s, (c, a, b) in coords.items():
                fh.write(f"{s}\t{c}\t{a}\t{b}\n")
    return coords


class DepMapGeneCN(CNProvider):
    """CNProvider backed by DepMap gene-level relative CN + gene coordinates."""

    def __init__(self, cn_csv, gene_coords, recenter=True):
        self.cn_csv = cn_csv
        self.coords = gene_coords
        self.recenter = recenter
        with open(cn_csv) as fh:
            header = next(csv.reader(fh))
        self.ci_model = header.index("ModelID")
        self.ci_def = header.index("IsDefaultEntryForModel")
        self.gene_col = {}
        for i, h in enumerate(header):
            if i <= self.ci_def:               # skip the leading metadata columns
                continue
            self.gene_col[h.split(" (")[0]] = i
        self.usable = [g for g in self.gene_col if g in self.coords]
        self._cache = {}

    def preload(self, model_ids):
        """One file pass to parse the default-profile rows for a set of ModelIDs (efficient for many lines)."""
        need = {m for m in model_ids if m and m not in self._cache}
        if not need:
            return
        rows = {}
        with open(self.cn_csv) as fh:
            r = csv.reader(fh)
            next(r)
            for row in r:
                if row[self.ci_model] in need and row[self.ci_def] == "Yes":
                    rows[row[self.ci_model]] = row
                    if len(rows) == len(need):
                        break
        for m in need:
            self._cache[m] = self._build(rows.get(m))

    def track(self, model_id):
        if model_id not in self._cache:
            self.preload([model_id])
        return self._cache.get(model_id)

    def _build(self, row):
        if row is None:
            return None
        segs, ratios = [], []
        for g in self.usable:
            v = row[self.gene_col[g]]
            if v in ("", "NA"):
                continue
            r = float(v)
            c, s, e = self.coords[g]
            segs.append((c, s, e, r))
            ratios.append(r)
        if not segs:
            return None
        if self.recenter:
            med = float(np.median(ratios))
            if med > 0:
                segs = [(c, s, e, r / med) for c, s, e, r in segs]
        return CNTrack(segs, ploidy=2.0, source="DepMapGeneWGS")
