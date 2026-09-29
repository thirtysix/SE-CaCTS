"""CN segments from a plain table: key, chrom, start, end, ratio (linear, 1.0 = neutral).

For sources converted once, offline, so scoring needs neither the converter nor its dependencies, e.g. CCLE
2019 SNP6 lifted to hg38 (phase1/scripts/21_ccle_snp6_validate.py --export-keys; cnrose.cn.ccle).
"""
from __future__ import annotations

from .base import CNTrack, CNProvider


class SegmentFileCN(CNProvider):
    def __init__(self, path, source="segment-file"):
        import pandas as pd
        d = pd.read_csv(path, sep="\t", dtype={"key": str, "chrom": str})
        self._tracks = {k: CNTrack(list(zip(g.chrom, g.start, g.end, g.ratio)), ploidy=2.0, source=f"{source}:{k}")
                        for k, g in d.groupby("key")}

    def track(self, key):
        return self._tracks.get(key)
