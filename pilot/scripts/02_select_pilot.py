#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Step B pilot — 02: choose the ~12 pilot H3K27ac experiments from data/manifest.tsv.

Design goal: a balanced lineage x cell-line grid so that LINEAGE (biology) and STUDY (batch) are NOT
collinear, plus a cross-study replicate probe (one cell line represented by >=2 experiments) as the
cleanest direct read on batch. See ../PILOT.md.

Runs LOCALLY, no network by default. `--resolve-study` (opt-in) best-effort resolves each chosen GSM's
GEO series (GSE) so the batch axis is labelled before download; it is not required.
"""
import argparse, re, sys
from collections import defaultdict

# Curated cancer cell lines by lineage (aliases as they appear across ChIP-Atlas free text).
# Ovary/HGSOC is pinned in (ties to the parent OVCAR3 manuscript). Lines chosen for abundant H3K27ac.
TARGETS = {
    "Ovary_HGSOC":      ["OVCAR-3", "NIHOVCAR3", "OVCAR3", "Kuramochi", "OVSAHO", "CAOV3", "COV362", "JHOS-2"],
    "Breast":           ["MCF-7", "MCF7", "T-47D", "T47D", "MDA-MB-231", "ZR-75-1", "HCC1954", "SK-BR-3"],
    "Colorectal":       ["HCT-116", "HCT116", "HT-29", "HT29", "Caco-2", "LoVo", "SW480", "DLD-1", "RKO"],
    "Blood_leukemia":   ["K-562", "K562", "HL-60", "HL60", "MV4-11", "Jurkat", "THP-1", "Kasumi-1", "NB4"],
    "Lung":             ["A-549", "A549", "NCI-H358", "NCI-H1299", "NCI-H23", "H358", "H1299"],
    "Prostate":         ["LNCaP", "VCaP", "PC-3", "PC3", "22Rv1", "DU145"],
    "Liver":            ["HepG2", "Hep-G2", "Huh-7", "Huh7", "SNU-475"],
}


def norm(s):
    return re.sub(r"[^a-z0-9]", "", (s or "").lower())


# canonical alias -> (lineage, canonical_name)
ALIAS = {}
for lin, names in TARGETS.items():
    for nm in names:
        ALIAS[norm(nm)] = (lin, nm)


def match_row(r):
    """Return (lineage, canonical_line) if a manifest row matches a curated target, else None.
    Longest-alias-first so 'NCI-H358' wins over a bare 'H358' substring collision."""
    hay = " ".join(norm(r[k]) for k in ("cell_line", "cell_type", "name"))
    for alias in sorted(ALIAS, key=len, reverse=True):
        if alias and alias in hay:
            return ALIAS[alias]
    return None


def read_manifest(path):
    rows = []
    with open(path) as f:
        hdr = f.readline().rstrip("\n").split("\t")
        for line in f:
            rows.append(dict(zip(hdr, line.rstrip("\n").split("\t"))))
    return rows


def resolve_gse(gsm, timeout=15):
    """Best-effort GSM -> GSE via GEO acc.cgi (text view). Returns '' on any failure/offline."""
    import urllib.request
    url = (f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gsm}"
           "&targ=self&form=text&view=quick")
    try:
        with urllib.request.urlopen(url, timeout=timeout) as resp:
            txt = resp.read().decode("utf-8", "replace")
        m = re.search(r"Sample_series_id\s*=\s*(GSE\d+)", txt)
        return m.group(1) if m else ""
    except Exception:
        return ""


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--manifest", default="../data/manifest.tsv")
    ap.add_argument("--out", default="../data/selection.tsv")
    ap.add_argument("--n-lineages", type=int, default=4)
    ap.add_argument("--per-lineage", type=int, default=3, help="distinct cell lines per lineage")
    ap.add_argument("--pin", default="Ovary_HGSOC", help="lineage always included if available")
    ap.add_argument("--batch-probe", default="", help="canonical line to over-sample (>=2 SRX) as batch probe; "
                                                      "'' = auto-pick the line with the most experiments")
    ap.add_argument("--resolve-study", action="store_true", help="opt-in: GSM->GSE via GEO (network)")
    a = ap.parse_args()

    rows = read_manifest(a.manifest)

    # bucket matching experiments: lineage -> canonical_line -> [rows]
    buckets = defaultdict(lambda: defaultdict(list))
    for r in rows:
        hit = match_row(r)
        if hit:
            lin, canon = hit
            buckets[lin][canon].append(r)

    if not buckets:
        sys.exit("[02] no curated cancer lines matched the manifest — check ../data/manifest.tsv")

    # rank lineages by how many distinct curated lines are actually present
    lin_order = sorted(buckets, key=lambda L: (-len(buckets[L]), L))
    chosen_lins = []
    if a.pin in buckets:
        chosen_lins.append(a.pin)
    for L in lin_order:
        if len(chosen_lins) >= a.n_lineages:
            break
        if L not in chosen_lins:
            chosen_lins.append(L)

    # auto batch-probe: the single (lineage,line) with the most experiments across all buckets
    probe_line = a.batch_probe
    if not probe_line:
        best = max(((len(v), L, c) for L in buckets for c, v in buckets[L].items()),
                   default=(0, "", ""))
        probe_line = best[2]

    selected = []  # dicts with role
    seen_srx = set()

    def take(r, lineage, canon, role):
        if r["srx"] in seen_srx:
            return
        seen_srx.add(r["srx"])
        selected.append(dict(srx=r["srx"], lineage=lineage, cell_line=canon, gsm=r["gsm"],
                             antibody=r["antibody"], name=r["name"], role=role))

    for L in chosen_lins:
        # distinct lines in this lineage, most-populated first
        lines_here = sorted(buckets[L], key=lambda c: (-len(buckets[L][c]), c))
        for canon in lines_here[:a.per_lineage]:
            take(buckets[L][canon][0], L, canon, "lineage")  # one representative experiment per line

    # cross-study batch probe: add up to 2 extra experiments of probe_line, chosen from DIFFERENT studies
    # than each other and than the representative already taken. Without a GSE lookup we proxy "study" by
    # (antibody vendor, GSM submission block) — GSMs within one GEO series are ~contiguous.
    def study_key(r):
        ab = norm(r.get("antibody", ""))
        digits = re.sub(r"\D", "", r.get("gsm", "") or "")
        blk = int(digits) // 500 if digits else 0
        return (ab, blk)

    def n_probe():
        return sum(1 for s in selected if s["role"] == "batch-probe")

    if probe_line:
        for L in buckets:
            if probe_line in buckets[L]:
                allexp = buckets[L][probe_line]
                used_keys = {study_key(r) for r in allexp if r["srx"] in seen_srx}
                pool = [r for r in allexp if r["srx"] not in seen_srx]
                for r in pool:                      # prefer different-study replicates
                    if n_probe() >= 2:
                        break
                    if study_key(r) not in used_keys:
                        used_keys.add(study_key(r))
                        take(r, L, probe_line, "batch-probe")
                for r in pool:                      # fallback: fill remaining if <2 distinct studies found
                    if n_probe() >= 2:
                        break
                    take(r, L, probe_line, "batch-probe")
                break

    if a.resolve_study:
        print("[02] resolving GSE for chosen experiments (network)...", file=sys.stderr)
        cache = {}
        for s in selected:
            gsm = s["gsm"]
            s["gse"] = cache.get(gsm) or resolve_gse(gsm)
            cache[gsm] = s["gse"]

    # write
    import os
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    cols = ["srx", "lineage", "cell_line", "role", "gsm"] + (["gse"] if a.resolve_study else []) + \
           ["antibody", "name"]
    with open(a.out, "w") as o:
        o.write("\t".join(cols) + "\n")
        for s in selected:
            o.write("\t".join(str(s.get(c, "")) for c in cols) + "\n")

    # summary to stderr
    print(f"[02] selected {len(selected)} experiments across {len(chosen_lins)} lineages "
          f"-> {a.out}", file=sys.stderr)
    by_lin = defaultdict(list)
    for s in selected:
        by_lin[s["lineage"]].append(s)
    for L in chosen_lins:
        items = ", ".join(f"{s['cell_line']}({s['srx']}{'*' if s['role']=='batch-probe' else ''})"
                          for s in by_lin[L])
        print(f"      {L:16s} {items}", file=sys.stderr)
    print(f"      batch-probe line = {probe_line}  (* = extra cross-study replicate)", file=sys.stderr)


if __name__ == "__main__":
    main()
