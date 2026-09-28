"""CCLE 2019 SNP6 copy number as a CNTrack backend (FINDINGS 2026-09-28 §21).

Source: cBioPortal `ccle_broad_2019` — `data_cna_hg19.seg` (Affymetrix SNP6, 1,041 lines, continuous log2
segment means on hg19) and `data_clinical_sample.txt` (SAMPLE_ID -> DEPMAPID). Useful for lines whose
array CN predates DepMap's move to WGS/WES and therefore has no 26Q1 CN.

Each segment is lifted to hg38 by its two endpoints (pyliftover, UCSC hg19ToHg38 chain); segments whose ends
land on different chromosomes, the minus strand, or out of order are dropped and counted. Long segments
are split into <= `piece` bp so `CNTrack.region_cn`, which bounds its scan to segments starting within
3 Mb, still finds them. Values are linear ratios 2**seg.mean, re-centred to a length-weighted median of
1.0 per line (as DepMap's `recenter`). Keyed by DepMap ModelID.
"""
from __future__ import annotations

import numpy as np

from .base import CNTrack, CNProvider


class Ccle2019SnpCN(CNProvider):
    def __init__(self, seg_path, clinical_path, chain=None, recenter=True, piece=1_000_000):
        import pandas as pd
        cs = pd.read_csv(clinical_path, sep="\t", comment="#", usecols=["SAMPLE_ID", "DEPMAPID"])
        self.sample_for = {d: s for s, d in zip(cs.SAMPLE_ID, cs.DEPMAPID) if isinstance(d, str)}
        seg = pd.read_csv(seg_path, sep="\t", dtype={"chrom": str})
        self.seg = {k: g for k, g in seg.groupby("ID")}
        self.recenter = recenter
        self.piece = piece
        self._lo = None
        self.chain = chain          # hg19ToHg38.over.chain.gz path; pyliftover's own download can fail
        self._cache = {}
        self.dropped = {}

    def _lift(self):
        if self._lo is None:
            from pyliftover import LiftOver
            self._lo = LiftOver(self.chain) if self.chain else LiftOver("hg19", "hg38")
        return self._lo

    def track(self, model_id):
        if model_id in self._cache:
            return self._cache[model_id]
        sid = self.sample_for.get(model_id)
        g = self.seg.get(sid) if sid else None
        if g is None:
            self._cache[model_id] = None
            return None
        lo = self._lift()
        segs, drop = [], 0
        for c, s, e, v in zip(g["chrom"], g["loc.start"], g["loc.end"], g["seg.mean"]):
            ch = "chr" + str(c).replace("chr", "")
            if ch in ("chr23",):
                ch = "chrX"
            a = lo.convert_coordinate(ch, int(s) - 1)
            b = lo.convert_coordinate(ch, int(e) - 1)
            if not a or not b or a[0][0] != b[0][0] or a[0][2] != "+" or b[0][2] != "+" or b[0][1] <= a[0][1]:
                drop += 1
                continue
            c2, s2, e2 = a[0][0], int(a[0][1]), int(b[0][1]) + 1
            r = float(2.0 ** v) if np.isfinite(v) else np.nan
            if not np.isfinite(r):
                continue
            for p in range(s2, e2, self.piece):
                segs.append((c2, p, min(p + self.piece, e2), r))
        self.dropped[model_id] = drop
        if not segs:
            self._cache[model_id] = None
            return None
        if self.recenter:
            L = np.array([e - s for _, s, e, _ in segs], dtype=float)
            R = np.array([r for *_, r in segs])
            o = np.argsort(R)
            med = R[o][np.searchsorted(np.cumsum(L[o]), L.sum() / 2)]
            if med > 0:
                segs = [(c, s, e, r / med) for c, s, e, r in segs]
        self._cache[model_id] = CNTrack(segs, ploidy=2.0, source=f"CCLE2019_SNP6:{sid}")
        return self._cache[model_id]
