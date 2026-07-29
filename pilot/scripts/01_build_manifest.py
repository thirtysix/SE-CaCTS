#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Step B pilot — 01: build the H3K27ac experiment manifest from the LOCAL ChIP-Atlas AllCell bed.

No download. The file His.ALL.50.H3K27ac.AllCell.bed.gz annotates every peak's name field with
    ID=SRX...;Name=...;Title=GSM...;Cell group=...;cell line=...;cell type=...;chip antibody=...
(URL-encoded). We scan the file once and emit one row per unique SRX experiment.

Runs LOCALLY (needs the 3.1 GB bed). Output data/manifest.tsv is tiny and can be staged to HPC.
"""
import os
import argparse, gzip, re, sys
from urllib.parse import unquote

_ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
sys.path.insert(0, _ROOT)
from secacts_env import DATAROOT                      # noqa: E402  (local paths live in .env)

DEFAULT_BED = os.path.join(DATAROOT, "chip-atlas/00.data/Histone/His.ALL.50.H3K27ac.AllCell.bed.gz")

RE_SRX   = re.compile(r'ID=(SRX\d+)')
RE_GSM   = re.compile(r'Title=(GSM\d+)')
RE_NAME  = re.compile(r'Name=([^;]*)')
RE_GROUP = re.compile(r'Cell%20group=([^;<]*)')
RE_LINE  = re.compile(r'cell%20line=([^;<]*)')
RE_TYPE  = re.compile(r'cell%20type=([^;<]*)')
RE_AB    = re.compile(r'(?:chip%20antibody|chip%20epitope)=([^;<]*)')

COLS = ["srx", "gsm", "name", "cell_group", "cell_line", "cell_type", "antibody"]


def dec(x):
    return unquote(x).strip() if x else ""


def grab(rx, s):
    m = rx.search(s)
    return m.group(1) if m else ""


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--bed", default=DEFAULT_BED, help="ChIP-Atlas H3K27ac AllCell bed.gz")
    ap.add_argument("--out", default="../data/manifest.tsv")
    ap.add_argument("--limit", type=int, default=0, help="stop after N lines (0=all; for quick tests)")
    a = ap.parse_args()

    seen = {}
    n = 0
    with gzip.open(a.bed, "rt") as f:
        for line in f:
            if line.startswith("track"):
                continue
            n += 1
            if a.limit and n > a.limit:
                break
            parts = line.split("\t", 4)
            if len(parts) < 4:
                continue
            col4 = parts[3]
            m = RE_SRX.search(col4)
            if not m:
                continue
            srx = m.group(1)
            if srx in seen:
                continue
            seen[srx] = dict(
                srx=srx,
                gsm=grab(RE_GSM, col4),
                name=dec(grab(RE_NAME, col4)),
                cell_group=dec(grab(RE_GROUP, col4)),
                cell_line=dec(grab(RE_LINE, col4)),
                cell_type=dec(grab(RE_TYPE, col4)),
                antibody=dec(grab(RE_AB, col4)),
            )
            if n % 5_000_000 == 0:
                print(f"  ...{n:,} lines scanned, {len(seen):,} unique SRX", file=sys.stderr)

    import os
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    rows = sorted(seen.values(), key=lambda r: (r["cell_group"], r["cell_line"] or r["cell_type"], r["srx"]))
    with open(a.out, "w") as o:
        o.write("\t".join(COLS) + "\n")
        for r in rows:
            o.write("\t".join(r[c] for c in COLS) + "\n")
    print(f"[01] {len(rows):,} unique H3K27ac experiments  ->  {a.out}", file=sys.stderr)
    # quick cell-group census to stderr
    from collections import Counter
    c = Counter(r["cell_group"] for r in rows)
    print("[01] top cell groups: " + ", ".join(f"{k}={v}" for k, v in c.most_common(12)), file=sys.stderr)


if __name__ == "__main__":
    main()
