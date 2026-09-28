#!/usr/bin/env python3
"""Remove the library-layout effect from the SE x experiment matrix (sensitivity arm, FINDINGS §16).

The per-SE effect is the one layout_batch.py estimates from WITHIN-line contrasts only (mean log1p signal
of a line's paired-end experiments minus its single-end ones, averaged over the 125 lines that have both),
so it carries no lineage biology. It is subtracted from every paired-end column in log1p space:

    S'[se, PE] = expm1(max(log1p(S[se, PE]) - effect[se], 0))

Single-end columns are unchanged. Streams the matrix row by row, so memory stays small.

  python layout_correct.py --signal atlas.s3.se_signal.tsv --effects layout_sensitive_ses.tsv.gz \
      --layout srx_layout.tsv --out atlas.s3.layoutcorr.se_signal.tsv
"""
from __future__ import annotations

import argparse
import gzip

import numpy as np


def opener(p, mode="rt"):
    return gzip.open(p, mode) if p.endswith(".gz") else open(p, mode)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--signal", required=True)
    ap.add_argument("--effects", required=True, help="layout_sensitive_ses.tsv.gz (column pe_minus_se)")
    ap.add_argument("--layout", required=True, help="srx_layout.tsv")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    eff = {}
    with opener(a.effects) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        j = hdr.index("pe_minus_se")
        for l in fh:
            f = l.rstrip("\n").split("\t")
            eff[f[0]] = float(f[j])
    paired = set()
    with open(a.layout) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        j = hdr.index("layout")
        for l in fh:
            f = l.rstrip("\n").split("\t")
            if f[j] == "PAIRED":
                paired.add(f[0])

    n = changed = 0
    with opener(a.signal) as fin, opener(a.out, "wt") as fout:
        head = fin.readline()
        fout.write(head)
        cols = head.rstrip("\n").split("\t")[1:]
        pe = np.array([c in paired for c in cols])
        for l in fin:
            f = l.rstrip("\n").split("\t")
            e = eff.get(f[0], 0.0)
            if e:
                v = np.asarray(f[1:], dtype=np.float64)
                v[pe] = np.expm1(np.maximum(np.log1p(v[pe]) - e, 0.0))
                f[1:] = [f"{x:.6g}" for x in v]
                changed += 1
            fout.write("\t".join(f) + "\n")
            n += 1
    print(f"[layout_correct] {n} SEs, {changed} corrected, {pe.sum()} paired-end columns of {len(cols)}")


if __name__ == "__main__":
    main()
