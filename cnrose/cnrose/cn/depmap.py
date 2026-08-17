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


class DepMapMcWesCN(DepMapGeneCN):
    """DepMap WES gene-level CN, keyed by ModelCondition (`OmicsCNGeneMC_WES.csv`, 26Q1+).

    DepMap DISCONTINUED the merged `OmicsCNGene.csv` after 24Q4 (GOTCHA 90); from 25Q2 copy number is
    split by modality. For 26Q1 the WES arm is this file, and unlike the WGS one it carries no `ModelID`
    column — rows are `ModelConditionID` + `IsDefaultEntryForMC`. Values are on the same scale as
    `OmicsCNGeneWGS.csv` (LINEAR relative CN, 1.0 = neutral; verified medians 0.999 / 1.045 / 1.020), so
    the parent's `_build` applies unchanged.

    Two wrinkles this class absorbs:
      * `IsDefaultEntryForMC` is not Yes/No — non-default rows carry `No_CDS-<id>` naming the CDS profile
        that won, so test `== "Yes"` rather than `!= "No"`.
      * 88 models have more than one default ModelCondition. Ties are broken by sorted ModelConditionID
        so a run is reproducible, rather than by file order.
    """

    def __init__(self, cn_csv, model_condition_csv, gene_coords, recenter=True):
        self.cn_csv = cn_csv
        self.coords = gene_coords
        self.recenter = recenter
        with open(cn_csv) as fh:
            header = next(csv.reader(fh))
        self.ci_mc = header.index("ModelConditionID")
        self.ci_def = header.index("IsDefaultEntryForMC")
        self.gene_col = {}
        for i, h in enumerate(header):
            if i <= self.ci_def:
                continue
            self.gene_col[h.split(" (")[0]] = i
        self.usable = [g for g in self.gene_col if g in self.coords]

        mc_to_model = {}
        with open(model_condition_csv, newline="") as fh:
            for row in csv.DictReader(fh):
                mc, mid = row.get("ModelConditionID"), row.get("ModelID")
                if mc and mid:
                    mc_to_model[mc] = mid
        self.mc_to_model = mc_to_model
        self._cache = {}

    def _mcs_for(self, model_ids):
        """ModelID -> the one ModelConditionID to read (default entry, ties broken by sorted MC id)."""
        want = set(model_ids)
        found = {}
        with open(self.cn_csv) as fh:
            r = csv.reader(fh)
            next(r)
            for row in r:
                mc = row[self.ci_mc]
                mid = self.mc_to_model.get(mc)
                if mid in want and row[self.ci_def] == "Yes":
                    found.setdefault(mid, []).append(mc)
        return {mid: sorted(mcs)[0] for mid, mcs in found.items()}

    def preload(self, model_ids):
        need = {m for m in model_ids if m and m not in self._cache}
        if not need:
            return
        chosen = self._mcs_for(need)
        target = {mc: mid for mid, mc in chosen.items()}
        rows = {}
        with open(self.cn_csv) as fh:
            r = csv.reader(fh)
            next(r)
            for row in r:
                mid = target.get(row[self.ci_mc])
                if mid:
                    rows[mid] = row
                    if len(rows) == len(target):
                        break
        for m in need:
            track = self._build(rows.get(m))
            if track is not None:
                track.source = "DepMapGeneWES"
            self._cache[m] = track
