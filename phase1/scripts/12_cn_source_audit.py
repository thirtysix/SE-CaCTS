#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Phase 1 — 12: per-line CN-source audit (makes the DESIGN.md §4 tier split reproducible).

Scripts the one-off 2026-07-18 audit that produced the 254/119/25 net-new tier split. For every
QC-pass cell line in scope it emits a `cn_source` / tier by routing Cellosaurus `DR` cross-refs and
auditing matched ChIP-Atlas Input-control experiments — all from LOCAL files, no network.

NOTE: this now runs on the identity-CORRECTED manifest (script 11's sentinel-sink + RRID-canonicalization
fixes), so the split is **260 / 113 / 25**, not the original 254 / 119 / 25. The method reproduces §4; the
small shift is expected — several homonym/subline CVCLs were corrected to their canonical human accessions,
whose Cellosaurus `DR` blocks carry more measured-CN cross-refs (so lines moved C→B). See CN_COVERAGE.md.

Tiers (partition the full QC-pass line set exactly):
  A depmap_wgs      already in the atlas pull  (manifest has_cn==1: ModelID in OmicsCNGeneWGS.csv)
  B measured_other  net-new, CN obtainable from another DB via a Cellosaurus DR cross-ref:
                      DepMap (-> merged OmicsCNGene WES/SNP) | cancercelllines (Progenetix)
                      | Cell_Model_Passport (CMP) | Cosmic-CLP (COSMIC Cell Line Project)
  C input_inferred  net-new, no measured CN, but the line HAS a matched hg38 Input-control experiment
                      (=> CN inferrable from the input bigWig: HMCan/Control-FREEC/CNVkit)
  D none            net-new orphan: no measured CN, no matched input -> agnostic catalog only

Inputs (all local): the phase-1 manifest (line set + Tier-A membership), Cellosaurus (DR blocks +
the same name->CVCL map script 11 uses), ChIP-Atlas experimentList.tab (Input-control rows).
Output: data/cn_source.tsv (one row per CVCL line) + a reconciliation report vs DESIGN.md §4 to stderr.

Note on GDSC: the GDSC id == COSMIC ID_SAMPLE == Cosmic-CLP id (DESIGN §4), so GDSC adds ~no lines
beyond Cosmic-CLP; it is reported for reach but NOT part of the measured predicate (matches §4's
"DepMap-merged / Progenetix / CMP / COSMIC" tier-1 set).
"""
import argparse, csv, os, re, sys
from collections import defaultdict

_ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
sys.path.insert(0, _ROOT)
from secacts_env import DATAROOT                      # noqa: E402  (local paths live in .env)
D = DATAROOT
# measured-CN resources -> the Cellosaurus DR token that signals the line is present there
MEASURED_DR = {
    "depmap":     "DepMap",              # -> merged OmicsCNGene.csv (WES/SNP), recovers non-WGS lines
    "progenetix": "cancercelllines",     # Progenetix registry (cancercelllines.org), CVCL-native CNV
    "cmp":        "Cell_Model_Passport",  # Cell Model Passports, absolute integer CN
    "cosmicclp":  "Cosmic-CLP",          # COSMIC Cell Line Project (NOT the 18k `Cosmic;` mutation rows)
}
GDSC_DR = "GDSC"  # reported for reach only; == cosmicclp id, not in the measured predicate


def norm(s):
    return re.sub(r"[^a-z0-9]", "", (s or "").lower())


def load_cellosaurus(path):
    """One pass -> (name2cvcl: normalized name -> CVCL, cvcl2dr: CVCL -> set(DR resource tokens)).

    name2cvcl replicates script 11 exactly (ID/SY, first-writer-wins) so input-experiment cell names
    resolve identically to the manifest's H3K27ac cells."""
    name2cvcl, cvcl2dr = {}, defaultdict(set)
    ac, names, dr = None, [], set()
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
            elif tag == "DR":
                dr.add(line[5:].split(";", 1)[0].strip())
            elif line.startswith("//"):
                if ac:
                    n_entries += 1
                    for nm in names:
                        k = norm(nm)
                        if k and k not in name2cvcl:
                            name2cvcl[k] = ac
                    if dr:
                        cvcl2dr[ac] |= dr
                ac, names, dr = None, [], set()
    print(f"[12] Cellosaurus: {n_entries:,} entries; {len(name2cvcl):,} name keys; "
          f"{len(cvcl2dr):,} entries with DR cross-refs", file=sys.stderr)
    return name2cvcl, cvcl2dr


def resolve_cvcl(name2cvcl, cell, cell_class="", cell_desc=""):
    """Same name->CVCL cascade the manifest uses (cell, cell+class, cell_desc)."""
    for cand in (cell, f"{cell} {cell_class}", cell_desc):
        c = norm(cand)
        if c in name2cvcl:
            return name2cvcl[c]
    return ""


def input_line_sets(explist, name2cvcl):
    """hg38 Input-control experiments -> (set of CVCLs, set of normalized raw cell strings) that have a
    matched input. Cols per script 11: 1 genome, 2 antigen_class, 4 cell_class, 5 cell, 6 cell_desc."""
    cvcls, cells = set(), set()
    n = 0
    with open(explist, encoding="utf-8", errors="replace") as f:
        for line in f:
            c = line.rstrip("\n").split("\t")
            if len(c) < 6 or c[1] != "hg38" or c[2] != "Input control":
                continue
            n += 1
            cell = c[5]
            cells.add(norm(cell))
            v = resolve_cvcl(name2cvcl, cell, c[4], c[6] if len(c) > 6 else "")
            if v:
                cvcls.add(v)
    print(f"[12] Input audit: {n:,} hg38 Input-control experiments -> {len(cvcls):,} CVCL lines "
          f"(+{len(cells):,} raw cell strings) have a matched input", file=sys.stderr)
    return cvcls, cells


def load_manifest(path):
    """QC-pass lines from the manifest. Returns per-CVCL aggregate:
       cvcl -> dict(cell, model_id, lineage, subtype, n_exp, cells:set(norm), tierA:bool)."""
    lines = {}
    with open(path) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["qc_pass"] != "1" or not r["cvcl"]:
                continue
            cv = r["cvcl"]
            d = lines.setdefault(cv, dict(cell=r["cell"], model_id=r["model_id"], lineage=r["lineage"],
                                          subtype=r["subtype"], n_exp=0, cells=set(), tierA=False))
            d["n_exp"] += 1
            d["cells"].add(norm(r["cell"]))
            # backfill representative labels from any non-empty row
            for k in ("cell", "model_id", "lineage", "subtype"):
                if not d[k] and r[k]:
                    d[k] = r[k]
            if r["has_cn"] == "1":
                d["tierA"] = True
    return lines


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--manifest", default="../data/phase1_manifest.tsv")
    ap.add_argument("--experiment-list", default=f"{D}/chip-atlas/experimentList.tab")
    ap.add_argument("--cellosaurus", default=f"{D}/cellosaurus/cellosaurus.txt")
    ap.add_argument("--out", default="../data/cn_source.tsv")
    a = ap.parse_args()

    lines = load_manifest(a.manifest)
    name2cvcl, cvcl2dr = load_cellosaurus(a.cellosaurus)
    in_cvcls, in_cells = input_line_sets(a.experiment_list, name2cvcl)

    res_cols = list(MEASURED_DR)  # depmap, progenetix, cmp, cosmicclp
    cols = ["cvcl", "cell", "model_id", "lineage", "subtype", "n_h3k27ac_exp",
            "tier", "cn_source"] + [f"dr_{k}" for k in res_cols] + ["dr_gdsc", "has_matched_input"]

    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    out = open(a.out, "w")
    out.write("\t".join(cols) + "\n")

    # counters
    tier_ct = defaultdict(int)                 # over all 721
    nn_reach = defaultdict(int)                 # per-resource reach among net-new (398)
    nn_measured = nn_input = nn_orphan = nn = 0
    for cv in sorted(lines):
        d = lines[cv]
        dr = cvcl2dr.get(cv, set())
        flags = {k: (MEASURED_DR[k] in dr) for k in res_cols}
        gdsc = GDSC_DR in dr
        measured = any(flags.values())
        matched_input = (cv in in_cvcls) or bool(d["cells"] & in_cells)

        if d["tierA"]:
            tier, src = "A", "depmap_wgs"
        else:
            nn += 1
            for k in res_cols:
                if flags[k]:
                    nn_reach[k] += 1
            if gdsc:
                nn_reach["gdsc"] += 1
            if measured:
                tier, src = "B", "measured_other"; nn_measured += 1
            elif matched_input:
                tier, src = "C", "input_inferred"; nn_input += 1
            else:
                tier, src = "D", "none"; nn_orphan += 1
        tier_ct[tier] += 1

        out.write("\t".join(str(x) for x in [
            cv, d["cell"], d["model_id"], d["lineage"], d["subtype"], d["n_exp"],
            tier, src, *[int(flags[k]) for k in res_cols], int(gdsc), int(matched_input)]) + "\n")
    out.close()

    tot = len(lines)
    A, B, C, Dn = tier_ct["A"], tier_ct["B"], tier_ct["C"], tier_ct["D"]
    def pct(x, d): return f"{100*x/d:.1f}%" if d else "n/a"
    print(f"\n[12] CN-SOURCE AUDIT  -> {a.out}", file=sys.stderr)
    print(f"  QC-pass lines (CVCL):                 {tot:,}", file=sys.stderr)
    print(f"  A depmap_wgs     (in the atlas pull): {A:>4}  {pct(A,tot)}", file=sys.stderr)
    print(f"  B measured_other (DR lookup)        : {B:>4}  {pct(B,tot)}", file=sys.stderr)
    print(f"  C input_inferred (matched input)    : {C:>4}  {pct(C,tot)}", file=sys.stderr)
    print(f"  D none           (orphan)           : {Dn:>4}  {pct(Dn,tot)}", file=sys.stderr)
    print(f"  CN available (A+B+C):                 {A+B+C:>4}  {pct(A+B+C,tot)}", file=sys.stderr)
    print(f"\n  Net-new lines (B+C+D):                {nn:,}", file=sys.stderr)
    print(f"    measured {nn_measured} ({pct(nn_measured,nn)}) / input {nn_input} ({pct(nn_input,nn)}) "
          f"/ orphan {nn_orphan} ({pct(nn_orphan,nn)})", file=sys.stderr)
    print("    per-resource reach among net-new (overlapping):", file=sys.stderr)
    for k in res_cols + ["gdsc"]:
        print(f"      {k:12s} {nn_reach[k]:>4}", file=sys.stderr)
    print("\n  DESIGN.md §4 (2026-07-18, PRE-fix identities) said: net-new 398 = 254 / 119 / 25;",
          file=sys.stderr)
    print("    reach Progenetix 245, CMP 153, DepMap 140, COSMIC 95, GDSC 98. The shift above is EXPECTED:",
          file=sys.stderr)
    print("    the identity fix (script 11) corrected homonym/subline CVCLs, whose true human DR blocks add",
          file=sys.stderr)
    print("    measured-CN cross-refs (C->B). This corrected split supersedes §4.", file=sys.stderr)


if __name__ == "__main__":
    main()
