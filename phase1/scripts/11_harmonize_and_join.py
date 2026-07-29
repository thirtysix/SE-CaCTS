#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Phase 1 — 11: build the H3K27ac x DepMap manifest via Cellosaurus.

Pipeline:
  ChIP-Atlas experimentList.tab  --filter hg38 + H3K27ac-->  H3K27ac experiments (+ per-experiment QC)
  Cellosaurus flat file          --ID + synonyms----------->  free-text cell name -> CVCL accession
  DepMap Model.csv               --RRID(=CVCL)------------->  ModelID + Oncotree lineage/subtype
  DepMap CN / dependency         --ModelID---------------->  has_cn / has_dep flags

Output: data/phase1_manifest.tsv (one row per H3K27ac experiment, with the DepMap join + QC) and a coverage
report to stderr. No heavy compute — pure metadata harmonization. See ../../DATA_SOURCES.md.
"""
import os
import argparse, csv, re, sys
from collections import defaultdict

_ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
sys.path.insert(0, _ROOT)
from secacts_env import DATAROOT                      # noqa: E402  (local paths live in .env)
D = DATAROOT
# ChIP-Atlas free-text "not available" sentinels (and other non-specific fillers). These must NEVER be
# used as a name-resolution key: `cell_desc` is literally "NA" for thousands of under-annotated
# experiments, and Cellosaurus happens to list "NA" as a synonym of a Chinese-hamster line (CVCL_E4I6),
# so an unguarded cascade sinks all of them into that one bogus "line". Guard by sentinel + min length.
NAME_SENTINELS = {"na", "nd", "none", "null", "unknown", "notavailable", "notapplicable",
                  "other", "mixed", "control", "ctrl", "input", "wt"}


def norm(s):
    return re.sub(r"[^a-z0-9]", "", (s or "").lower())


def usable_key(cand):
    """A normalized name is safe to resolve on only if it is specific enough — not a not-available
    sentinel (the "NA"/cell_desc sink) and not a 1-char token. Legit 2-char lines (SR, KG) are kept."""
    k = norm(cand)
    return k if (len(k) >= 2 and k not in NAME_SENTINELS) else ""


# ---------- Cellosaurus: normalized name -> CVCL ----------
def load_cellosaurus(path):
    name2cvcl = {}
    ac = None
    names = []
    n_entries = 0
    with open(path, encoding="utf-8", errors="replace") as f:
        for line in f:
            tag = line[:2]
            if tag == "ID" or tag == "SY":
                val = line[5:].rstrip("\n")
                if tag == "SY":
                    names += [s.strip() for s in val.split(";") if s.strip()]
                else:
                    names += [val.strip()]
            elif tag == "AC":
                ac = line[5:].strip()
            elif line.startswith("//"):
                if ac:
                    n_entries += 1
                    for nm in names:
                        k = norm(nm)
                        # first writer wins, but a primary AC line's own ID should not be overwritten by a
                        # later entry's synonym: keep first occurrence (Cellosaurus lists primary entries once)
                        if k and k not in name2cvcl:
                            name2cvcl[k] = ac
                ac = None
                names = []
    print(f"[11] Cellosaurus: {n_entries:,} entries -> {len(name2cvcl):,} normalized name keys", file=sys.stderr)
    return name2cvcl


# ---------- DepMap ----------
def load_depmap():
    cvcl2model = {}          # CVCL -> (ModelID, lineage, subtype, stripped, rrid)
    stripped2model = {}      # normalized stripped/full name -> same tuple (fallback path)
    with open(f"{D}/DepMap/Model.csv") as f:
        for row in csv.DictReader(f):
            rrid = row.get("RRID", "")
            info = (row["ModelID"], row.get("OncotreeLineage", ""),
                    row.get("OncotreeSubtype", ""), row.get("StrippedCellLineName", ""), rrid)
            if rrid.startswith("CVCL"):
                cvcl2model.setdefault(rrid, info)
            for nm in (row.get("StrippedCellLineName", ""), row.get("CellLineName", "")):
                k = usable_key(nm)
                if k:
                    stripped2model.setdefault(k, info)
    cn = set()
    with open(f"{D}/DepMap/2026q1/OmicsCNGeneWGS.csv") as f:
        for row in csv.DictReader(f):
            if row.get("ModelID"):
                cn.add(row["ModelID"])
    dep = set()
    with open(f"{D}/DepMap/CRISPRGeneEffect.csv") as f:
        first = f.readline().split(",")[0]  # header
        for line in f:
            mid = line.split(",", 1)[0].strip().strip('"')
            if mid.startswith("ACH"):
                dep.add(mid)
    print(f"[11] DepMap: {len(cvcl2model):,} models w/ CVCL RRID; {len(cn):,} w/ CN; {len(dep):,} w/ dependency",
          file=sys.stderr)
    return cvcl2model, stripped2model, cn, dep


# ---------- ChIP-Atlas experimentList.tab ----------
def iter_h3k27ac(path):
    """Yield hg38 H3K27ac rows. experimentList.tab columns (tab-sep):
       0 SRX  1 genome  2 antigen_class  3 antigen  4 cell_class  5 cell  6 cell_desc
       7 processing_logs (comma: reads,%mapped,%dup,npeaks)  8 title  9+ attributes."""
    with open(path, encoding="utf-8", errors="replace") as f:
        for line in f:
            c = line.rstrip("\n").split("\t")
            if len(c) < 6:
                continue
            if c[1] != "hg38" or c[3] != "H3K27ac":
                continue
            qc = c[7].split(",") if len(c) > 7 else []
            def qn(i):
                try:
                    return qc[i]
                except IndexError:
                    return ""
            yield dict(srx=c[0], cell_class=c[4], cell=c[5], cell_desc=c[6] if len(c) > 6 else "",
                       reads=qn(0), pct_mapped=qn(1), pct_dup=qn(2), n_peaks=qn(3))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--experiment-list", default=f"{D}/chip-atlas/experimentList.tab")
    ap.add_argument("--cellosaurus", default=f"{D}/cellosaurus/cellosaurus.txt")
    ap.add_argument("--out", default="../data/phase1_manifest.tsv")
    ap.add_argument("--min-reads", type=int, default=5_000_000, help="QC: min mapped-ish reads")
    ap.add_argument("--min-peaks", type=int, default=1_000, help="QC: min called peaks")
    a = ap.parse_args()

    name2cvcl = load_cellosaurus(a.cellosaurus)
    cvcl2model, stripped2model, cn, dep = load_depmap()

    import os
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    cols = ["srx", "cell_class", "cell", "cvcl", "model_id", "lineage", "subtype",
            "has_cn", "has_dep", "qc_pass", "reads", "pct_mapped", "pct_dup", "n_peaks"]

    tot = via_cvcl = matched = with_cn = qc_pass = qc_pass_cn = 0
    lines = set(); cn_lines = set()
    lin_lines = defaultdict(set)
    out = open(a.out, "w")
    out.write("\t".join(cols) + "\n")

    for r in iter_h3k27ac(a.experiment_list):
        tot += 1
        # name -> CVCL  (try cell, then cell+class, then description tokens; skip sentinel/short keys)
        cvcl = ""
        for cand in (r["cell"], f"{r['cell']} {r['cell_class']}", r["cell_desc"]):
            k = usable_key(cand)
            if k and k in name2cvcl:
                cvcl = name2cvcl[k]; break
        model = None
        if cvcl:
            via_cvcl += 1
            model = cvcl2model.get(cvcl)
        if model is None:                                   # fallback: direct name match to DepMap
            k = usable_key(r["cell"])
            model = stripped2model.get(k) if k else None
        # DepMap RRID is the AUTHORITATIVE identity for any line DepMap contains: prefer it over the
        # Cellosaurus name-match, which can land on a same-named non-human homonym (HAP1, KG-1, ...).
        if model and model[4].startswith("CVCL"):
            cvcl = model[4]
        # QC
        def toi(x):
            try:
                return int(float(x))
            except (ValueError, TypeError):
                return 0
        qcok = toi(r["reads"]) >= a.min_reads and toi(r["n_peaks"]) >= a.min_peaks
        if qcok:
            qc_pass += 1
        mid = lineage = subtype = ""
        has_cn = has_dep = False
        if model:
            matched += 1
            mid, lineage, subtype, _, _rrid = model
            has_cn = mid in cn; has_dep = mid in dep
            lines.add(mid); lin_lines[lineage].add(mid)
            if has_cn:
                with_cn += 1; cn_lines.add(mid)
                if qcok:
                    qc_pass_cn += 1
        out.write("\t".join(str(x) for x in [
            r["srx"], r["cell_class"], r["cell"], cvcl, mid, lineage, subtype,
            int(has_cn), int(has_dep), int(qcok), r["reads"], r["pct_mapped"], r["pct_dup"], r["n_peaks"]]) + "\n")
    out.close()

    print(f"\n[11] PHASE-1 PICTURE (Cellosaurus-bridged)  -> {a.out}", file=sys.stderr)
    print(f"  H3K27ac hg38 experiments:              {tot:,}", file=sys.stderr)
    print(f"  mapped to a CVCL (Cellosaurus):        {via_cvcl:,}", file=sys.stderr)
    print(f"  joined to a DepMap model:              {matched:,} experiments on {len(lines):,} lines", file=sys.stderr)
    print(f"  ...on lines WITH copy number:          {with_cn:,} experiments on {len(cn_lines):,} lines", file=sys.stderr)
    print(f"  QC-pass (>= {a.min_reads:,} reads & >= {a.min_peaks:,} peaks): {qc_pass:,} experiments "
          f"({qc_pass_cn:,} also DepMap+CN)", file=sys.stderr)
    print("  top lineages (DepMap-joined lines / experiments):", file=sys.stderr)
    lin_exp = defaultdict(int)
    # recount experiments per lineage from the file would need a second pass; use line counts as proxy label
    for k in sorted(lin_lines, key=lambda x: -len(lin_lines[x]))[:12]:
        print(f"     {k:30s} {len(lin_lines[k]):3d} lines", file=sys.stderr)


if __name__ == "__main__":
    main()
