#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Phase 1 — 14: recover Oncotree lineage from Cellosaurus for the non-DepMap lines.

The phase-1 join only read DepMap's OncotreeLineage, so the 259 lines without a DepMap record came out
lineage-less — even though Cellosaurus records their disease. This recovers it:

  1. Parse Cellosaurus per line: CA (category), OX (species), DI (NCIt disease code + term).
  2. Build an NCIt -> Oncotree crosswalk EMPIRICALLY from the ~1,860 DepMap models that carry both an
     OncotreeCode (Model.csv) and an NCIt disease (Cellosaurus) — no external Oncotree file needed.
     Majority vote per NCIt code; written to data/ncit_oncotree_crosswalk.tsv for audit.
  3. Assign each line a lineage by precedence:
       depmap            -> the DepMap OncotreeLineage/Subtype (authoritative; the 462 confirmed lines)
       cellosaurus_ncit  -> crosswalk hit on the line's NCIt disease code
       cellosaurus_term  -> curated term/keyword fallback for NCIt codes DepMap never used
       (none)            -> unresolved
  4. Gate inclusion (per user decision: cancer lines only) on THREE tests:
       - CA == "Cancer cell line"          (else exclude: non_cancer)
       - human (OX = Homo sapiens/9606)     (else exclude: non_human_homonym — misresolved name collision)
       - a lineage was resolved             (else exclude: no_lineage)
     DepMap-confirmed lines are authoritative and kept regardless of Cellosaurus CA/OX (their CVCL may be a
     homonym, but the RRID-matched DepMap lineage is correct); their species mismatches are FLAGGED, not
     dropped, in cvcl_species_flag for later CVCL correction.

Output: data/lineage_resolved.tsv (one row per line) + data/ncit_oncotree_crosswalk.tsv + a report.
Reproduce: cd phase1/scripts && python3 14_cellosaurus_lineage.py
"""
import argparse, csv, os, sys
from collections import Counter, defaultdict

_ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
sys.path.insert(0, _ROOT)
from secacts_env import DATAROOT                      # noqa: E402  (local paths live in .env)
D = DATAROOT
# Curated fallback for NCIt disease terms DepMap never annotated (so absent from the empirical crosswalk).
# Lineage names match DepMap OncotreeLineage vocabulary. Keep explicit + auditable; unmapped lines are
# reported, never silently guessed.
TERM_LINEAGE = {
    "nasopharyngeal carcinoma": "Head and Neck",
    "nasal type extranodal nk/t-cell lymphoma": "Lymphoid",
    "adult t-cell leukemia/lymphoma": "Lymphoid",
    "aids-related non-hodgkin lymphoma": "Lymphoid",
    "myxoid liposarcoma": "Soft Tissue",
    "childhood desmoplastic small round cell tumor": "Soft Tissue",
    "amelanotic cutaneous melanoma": "Skin",
    "clivus chordoma": "Bone",
    "renal cell carcinoma associated with xp11.2 translocations/tfe3 gene fusions": "Kidney",
}
# broad keyword backstop (applied only if the exact term isn't curated and the crosswalk missed)
KEYWORD_LINEAGE = [
    ("nasopharyng", "Head and Neck"), ("chordoma", "Bone"), ("osteosarcoma", "Bone"),
    ("ewing", "Bone"), ("liposarcoma", "Soft Tissue"), ("rhabdomyosarcoma", "Soft Tissue"),
    ("leiomyosarcoma", "Soft Tissue"), ("lymphoma", "Lymphoid"), ("nk/t-cell", "Lymphoid"),
    ("lymphoblastic", "Lymphoid"), ("hodgkin", "Lymphoid"), ("melanoma", "Skin"),
    ("renal cell carcinoma", "Kidney"),
]


def is_human(ox):
    return bool(ox) and ("Homo sapiens" in ox or "TaxID=9606" in ox)


def parse_di(v):
    p = [x.strip() for x in v.split(";")]
    return (p[1], p[2]) if len(p) >= 3 and p[0] == "NCIt" else None


def load_cellosaurus(path, keep):
    """CVCL -> dict(ca, ox, ncit=[(code,term)]) for the CVCLs in `keep`."""
    info = {}
    cur = ca = ox = None; ncit = []
    with open(path, encoding="utf-8", errors="replace") as f:
        for line in f:
            t = line[:2]
            if t == "AC": cur = line[5:].strip()
            elif t == "CA": ca = line[5:].strip()
            elif t == "OX" and ox is None: ox = line[5:].strip()
            elif t == "DI":
                d = parse_di(line[5:].strip())
                if d: ncit.append(d)
            elif line.startswith("//"):
                if cur in keep:
                    info[cur] = dict(ca=ca, ox=ox, ncit=ncit)
                cur = ca = ox = None; ncit = []
    return info


def load_models(path):
    """RRID(CVCL) -> dict(lineage, primary, subtype, ocode) for every DepMap model."""
    m = {}
    for row in csv.DictReader(open(path)):
        rrid = row.get("RRID", "")
        if rrid.startswith("CVCL"):
            m[rrid] = dict(lineage=row.get("OncotreeLineage", ""), primary=row.get("OncotreePrimaryDisease", ""),
                           subtype=row.get("OncotreeSubtype", ""), ocode=row.get("OncotreeCode", ""))
    return m


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cn-source", default="../data/cn_source.tsv")
    ap.add_argument("--model", default=f"{D}/DepMap/Model.csv")
    ap.add_argument("--cellosaurus", default=f"{D}/cellosaurus/cellosaurus.txt")
    ap.add_argument("--out", default="../data/lineage_resolved.tsv")
    ap.add_argument("--xwalk-out", default="../data/ncit_oncotree_crosswalk.tsv")
    a = ap.parse_args()

    lines = list(csv.DictReader(open(a.cn_source), delimiter="\t"))
    models = load_models(a.model)
    keep = {r["cvcl"] for r in lines} | set(models)          # targets + crosswalk training set
    cello = load_cellosaurus(a.cellosaurus, keep)

    # ---- build the empirical NCIt -> Oncotree crosswalk from co-annotated DepMap models ----
    # Vote LINEAGE first (robust), then take the modal (primary,subtype,ocode) WITHIN the winning lineage.
    # Reporting lineage_purity keeps the recovered lineage honest; subtype is best-effort (modal).
    xw_lin = defaultdict(Counter)     # code -> Counter(lineage)
    xw_tup = defaultdict(Counter)     # code -> Counter((lineage,primary,subtype,ocode))
    for cvcl, m in models.items():
        ci = cello.get(cvcl)
        if not m["lineage"] or not ci or not ci["ncit"] or not is_human(ci["ox"]):
            continue
        tup = (m["lineage"], m["primary"], m["subtype"], m["ocode"])
        for code, _t in ci["ncit"]:
            xw_lin[code][m["lineage"]] += 1
            xw_tup[code][tup] += 1
    xwalk_best = {}
    with open(a.xwalk_out, "w") as xf:
        xf.write("ncit_code\tn_support\tlineage_purity\tlineage\tprimary_disease\tsubtype\toncotree_code\n")
        for code in sorted(xw_lin):
            lin, linn = xw_lin[code].most_common(1)[0]
            tot = sum(xw_lin[code].values())
            tup = max((t for t in xw_tup[code] if t[0] == lin), key=lambda t: xw_tup[code][t])
            xwalk_best[code] = tup
            xf.write(f"{code}\t{tot}\t{linn/tot:.2f}\t" + "\t".join(tup) + "\n")

    def resolve_cello(cvcl):
        """(lineage, primary, subtype, source) from Cellosaurus for a non-DepMap line, or None."""
        ci = cello.get(cvcl)
        if not ci:
            return None
        for code, term in ci["ncit"]:
            if code in xwalk_best:
                lin, pri, sub, _oc = xwalk_best[code]
                return (lin, pri, sub, "cellosaurus_ncit")
        for _code, term in ci["ncit"]:            # curated/keyword fallback (lineage only)
            t = term.lower()
            if t in TERM_LINEAGE:
                return (TERM_LINEAGE[t], "", "", "cellosaurus_term")
            for kw, lin in KEYWORD_LINEAGE:
                if kw in t:
                    return (lin, "", "", "cellosaurus_term")
        return None

    cols = ["cvcl", "cell_line", "category", "species_human", "is_cancer",
            "ncit_code", "cello_disease", "lineage", "primary_disease", "subtype",
            "lineage_source", "include", "exclude_reason", "cvcl_species_flag"]
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    out = open(a.out, "w"); out.write("\t".join(cols) + "\n")

    ct = Counter(); src_ct = Counter(); excl_ct = Counter(); rec_lin = Counter()
    unmapped = []
    for r in lines:
        cvcl = r["cvcl"]; ci = cello.get(cvcl, {})
        ca = ci.get("ca", ""); ox = ci.get("ox", "")
        human = is_human(ox); is_cancer = (ca == "Cancer cell line")
        ncit0 = ci.get("ncit") or []
        ncit_code = ncit0[0][0] if ncit0 else ""
        cello_dis = ncit0[0][1] if ncit0 else ""
        depmap = bool(r["model_id"])
        species_flag = ""

        if depmap:
            # authoritative: use manifest DepMap lineage/subtype; keep regardless of homonym CA/OX
            lineage = r["lineage"]; subtype = r["subtype"]
            primary = models.get(cvcl, {}).get("primary", "")
            source = "depmap"
            dep_noncancer = (not lineage) or lineage.strip().lower() == "normal" \
                or primary.strip().lower() == "non-cancerous"
            if ci and not human:
                species_flag = ox.split("!")[-1].strip() if "!" in ox else ox  # e.g. "Mus musculus (Mouse)"
            if not dep_noncancer:
                include, reason = 1, ""
            else:
                include, reason = 0, "depmap_non_cancer"
        else:
            res = resolve_cello(cvcl)
            lineage = primary = subtype = ""; source = ""
            if res:
                lineage, primary, subtype, source = res
            # gates (species first: a non-human record means the name resolved to a homonym, so its
            # CA/DI are the wrong organism's and can't be trusted)
            if not human:
                include, reason = 0, "non_human_homonym"
            elif not is_cancer:
                include, reason = 0, "non_cancer"
            elif not lineage:
                include, reason = 0, "no_lineage"
                if ncit_code:
                    unmapped.append((cvcl, r.get("cell", ""), cello_dis))
            else:
                include, reason = 1, ""

        if include and not depmap:
            rec_lin[lineage] += 1
        ct[("include" if include else "exclude")] += 1
        if include: src_ct[source] += 1
        if not include and reason: excl_ct[reason] += 1

        out.write("\t".join(str(x) for x in [
            cvcl, r.get("cell", ""), ca, int(human), int(is_cancer), ncit_code, cello_dis,
            lineage, primary, subtype, source, include, reason, species_flag]) + "\n")
    out.close()

    n = len(lines)
    print(f"[14] crosswalk: {len(xw_lin)} NCIt codes from co-annotated DepMap models -> {a.xwalk_out}", file=sys.stderr)
    print(f"[14] lineage_resolved.tsv: {n} lines", file=sys.stderr)
    print(f"  INCLUDED {ct['include']}  /  EXCLUDED {ct['exclude']}", file=sys.stderr)
    print("  lineage source (included):", dict(src_ct), file=sys.stderr)
    print("  exclude reasons:", dict(excl_ct), file=sys.stderr)
    depmap_lin = src_ct["depmap"]
    recov = src_ct["cellosaurus_ncit"] + src_ct["cellosaurus_term"]
    print(f"  => lines with a lineage: {depmap_lin+recov}  (DepMap {depmap_lin} + RECOVERED {recov})", file=sys.stderr)
    print("  recovered lineage distribution:", file=sys.stderr)
    for k, v in rec_lin.most_common():
        print(f"     {v:3d}  {k}", file=sys.stderr)
    if unmapped:
        print(f"  !! {len(unmapped)} cancer lines had an NCIt disease but no lineage mapping (add to curation):",
              file=sys.stderr)
        for cv, nm, dis in unmapped:
            print(f"       {cv} {nm}  {dis}", file=sys.stderr)


if __name__ == "__main__":
    main()
