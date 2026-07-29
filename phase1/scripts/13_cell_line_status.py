#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Phase 1 — 13: human-readable per-line status (identity + lineage + capability).

Joins the CN-source audit (`cn_source.tsv`, script 12) with the recovered lineage
(`lineage_resolved.tsv`, script 14) into one browsable table — for every cell line: its canonical
identity and confidence, whether we hold its data, where its copy-number and its lineage come from,
and whether it is included in the cancer analysis (and if not, why). Narrative + summary:
`../CELL_LINE_IDENTITY.md`.

Reproduce (after 12 and 14):  cd phase1/scripts && python3 13_cell_line_status.py
"""
import csv, sys
from collections import Counter

CN = "../data/cn_source.tsv"
LIN = "../data/lineage_resolved.tsv"
OUT = "../data/cell_line_status.tsv"

CN_STATUS = {"A": "have (DepMap WGS)", "B": "lookup (other DB)",
             "C": "infer (from input)", "D": "none"}
RES_LABEL = {"dr_depmap": "DepMap-merged", "dr_progenetix": "Progenetix",
             "dr_cmp": "CMP", "dr_cosmicclp": "COSMIC"}


def main():
    cn = {r["cvcl"]: r for r in csv.DictReader(open(CN), delimiter="\t")}
    lin = {r["cvcl"]: r for r in csv.DictReader(open(LIN), delimiter="\t")}

    cols = ["cvcl", "cell_line", "category", "lineage", "subtype", "lineage_source",
            "n_experiments", "identity", "data_status", "cn_status", "cn_from",
            "tier", "include", "exclude_reason"]
    with open(OUT, "w") as out:
        out.write("\t".join(cols) + "\n")
        for cvcl in sorted(cn, key=lambda c: (lin[c]["include"] == "0", lin[c]["lineage"], cn[c]["cell"].lower())):
            r, L = cn[cvcl], lin[cvcl]
            t = r["tier"]
            identity = "confirmed (Cellosaurus+DepMap)" if r["model_id"] else "Cellosaurus-only"
            data_status = "in atlas" if t == "A" else "needs pull"
            if t == "A":
                cn_from = "DepMap WGS"
            elif t == "B":
                cn_from = ",".join(l for k, l in RES_LABEL.items() if r[k] == "1") or "?"
            elif t == "C":
                cn_from = "matched input"
            else:
                cn_from = "-"
            out.write("\t".join([
                cvcl, r["cell"], L["category"], L["lineage"], L["subtype"], L["lineage_source"],
                r["n_h3k27ac_exp"], identity, data_status, CN_STATUS[t], cn_from,
                t, L["include"], L["exclude_reason"]]) + "\n")

    # console summary (the "simple version")
    inc = [c for c in cn if lin[c]["include"] == "1"]
    exc = [c for c in cn if lin[c]["include"] == "0"]
    src = Counter(lin[c]["lineage_source"] for c in inc)
    tier = Counter(cn[c]["tier"] for c in inc)
    excr = Counter(lin[c]["exclude_reason"] for c in exc)
    n = len(cn)
    print(f"[13] wrote {OUT}: {n} lines", file=sys.stderr)
    print(f"  INCLUDED (human cancer, has lineage): {len(inc)}   EXCLUDED: {len(exc)} {dict(excr)}", file=sys.stderr)
    print(f"  lineage source: DepMap {src['depmap']} + recovered {src['cellosaurus_ncit']+src['cellosaurus_term']}"
          f" (crosswalk {src['cellosaurus_ncit']} / curated {src['cellosaurus_term']})", file=sys.stderr)
    print(f"  included by CN bucket: A {tier['A']} / B {tier['B']} / C {tier['C']} / D {tier['D']}", file=sys.stderr)


if __name__ == "__main__":
    main()
