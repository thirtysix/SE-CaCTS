#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Phase 1 — 17: H3K27ac in GEO for copy-number-measured lines that have none in ChIP-Atlas (FINDINGS §19).

1. Targets: cancer lines with a copy-number record (DepMap WGS, DepMap MC_WES, or CMP WES pureCN) and no
   ChIP-Atlas H3K27ac (not among the 631 resolved lines, the v2 pull set or the v3 candidates). A record is
   not a CN track; the track check happens when a line is actually added (FINDINGS 2026-08-18 §4).
2. GEO: every human GSM whose indexed text mentions H3K27ac, from ONE E-utilities query; summaries (title,
   short summary, SRX, release date) are fetched in batches and cached.
3. Match titles + summaries to target-line names (DepMap, CMP, Cellosaurus ID + synonyms). A name shared
   by two lines, or by a target and a line we already have, is dropped as ambiguous.
4. Verify each hit on its full GEO record: a cell-line / source field must name the same line, and the
   antibody or title must say H3K27ac. Names of <= 4 characters count ONLY from a cell-line field (the H9
   trap: FINDINGS §17).
5. Classify each verified sample's SRX against the ChIP-Atlas experiment list: hg38 H3K27ac under another
   cell label, hg38 under another antigen, filed only under a non-human genome (spike-in), or absent.

Outputs (phase2/analysis/out/): geo_h3k27ac_scan.samples.tsv, geo_h3k27ac_scan.lines.tsv,
geo_h3k27ac_scan.all_gsm_status.tsv (every human H3K27ac GSM vs ChIP-Atlas, for the misfiling rate).

  python3 phase1/scripts/17_geo_h3k27ac_scan.py --experiment-list <current experimentList.tab>
"""
import argparse
import gzip
import json
import os
import re
import sys
import time
import urllib.parse
import urllib.request

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, ROOT)
from secacts_env import DATAROOT, cache_path                     # noqa: E402

P1 = os.path.join(ROOT, "phase1", "data")
P2 = os.path.join(ROOT, "phase2", "data")
OUT = os.path.join(ROOT, "phase2", "analysis", "out")
DEPMAP = os.path.join(DATAROOT, "DepMap", "2026q1")
CMP = os.path.join(DATAROOT, "CellModelPassports")
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
QUERY = ('(H3K27ac OR H3K27Ac OR "H3K27 acetylation" OR "H3K27-ac" OR "H3 K27ac" OR '
         '"acetyl-histone H3 (Lys27)") AND "Homo sapiens"[Organism] AND gsm[Entry Type]')
K27 = re.compile(r"h3\s*[-_ ]?\s*k\s*27\s*[-_ ]?\s*ac|acetyl[- ]?histone h3[^;]{0,12}(lys|k)\s*27", re.I)
CELL_KEYS = re.compile(r"^(cell[ _]?line( name)?|cell|cell type|cellline|cell_line_name|source|cells)$", re.I)
strip = lambda s: re.sub(r"[^A-Z0-9]", "", str(s).upper())


def get(url, tries=4):
    for i in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=120) as r:
                return r.read().decode("utf-8", "replace")
        except Exception as e:                                        # noqa: BLE001
            if i == tries - 1:
                raise
            time.sleep(2 * (i + 1))
    return ""


# ----------------------------------------------------------------------------- targets
def cn_lines():
    allm = pd.read_csv(os.path.join(DEPMAP, "Model.csv"),
                       usecols=["ModelID", "RRID", "CellLineName", "OncotreeLineage", "OncotreePrimaryDisease"])
    depmap_rrid = set(allm.RRID.dropna())
    m = allm[allm.RRID.notna() & (allm.OncotreePrimaryDisease != "Non-Cancerous")]
    wgs = set(pd.read_csv(os.path.join(DEPMAP, "OmicsCNGeneWGS.csv"), usecols=["ModelID"]).ModelID)
    mc = pd.read_csv(os.path.join(DEPMAP, "ModelCondition.csv"), usecols=["ModelConditionID", "ModelID"])
    mcw_mc = set(pd.read_csv(os.path.join(DEPMAP, "OmicsCNGeneMC_WES.csv"),
                             usecols=["ModelConditionID"]).ModelConditionID)
    mcw = set(mc[mc.ModelConditionID.isin(mcw_mc)].ModelID)
    cache = cache_path("cmp_wes_models.txt")
    if not os.path.exists(cache):
        ids = set()
        for ch in pd.read_csv(os.path.join(CMP, "WES_pureCN_CNV_genes_latest.csv.gz"), usecols=["model_id"],
                              chunksize=5_000_000):
            ids |= set(ch.model_id.unique())
        open(cache, "w").write("\n".join(sorted(ids)))
    cmp_ids = set(open(cache).read().split())
    ml = pd.read_csv(os.path.join(CMP, "model_list_20240110.csv"),
                     usecols=["model_id", "model_name", "RRID", "tissue", "cancer_type"])
    ml = ml[ml.model_id.isin(cmp_ids) & ml.RRID.notna()]
    ml = ml[~ml.cancer_type.fillna("").str.contains("Non-Cancerous", case=False)]
    dep = m[m.ModelID.isin(wgs | mcw)].assign(source=lambda d: d.ModelID.map(
        lambda x: "depmap_wgs" if x in wgs else "depmap_mc_wes"))
    rows = [(r.RRID, r.CellLineName, r.OncotreeLineage, r.OncotreePrimaryDisease, r.source, r.ModelID)
            for r in dep.itertuples()]
    seen = set(dep.RRID)
    cancer_rrid = set(m.RRID)
    for r in ml.itertuples():
        if r.RRID in seen:
            continue
        if r.RRID in depmap_rrid and r.RRID not in cancer_rrid:
            continue                                                   # DepMap calls it non-cancerous
        rows.append((r.RRID, r.model_name, r.tissue, r.cancer_type, "cmp_wes", ""))
        seen.add(r.RRID)
    return pd.DataFrame(rows, columns=["cvcl", "cell_line", "lineage", "disease", "cn_record", "model_id"]) \
        .drop_duplicates("cvcl")


def human_lexicon():
    """stripped name -> set of CVCLs, over EVERY human Cellosaurus entry (ID + synonyms)."""
    lex, cur, ac, human = {}, [], None, False
    with open(os.path.join(DATAROOT, "cellosaurus", "cellosaurus.txt"), encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("ID   "):
                cur, human = [line[5:].strip()], False
            elif line.startswith("AC   "):
                ac = line[5:].strip()
            elif line.startswith("SY   "):
                cur += [s.strip() for s in line[5:].split(";") if s.strip()]
            elif line.startswith("OX   ") and "NCBI_TaxID=9606" in line:
                human = True
            elif line.startswith("//"):
                if human and ac:
                    for n in cur:
                        k = strip(n)
                        if len(k) >= 2 and not k.isdigit():
                            lex.setdefault(k, set()).add(ac)
                cur, ac, human = [], None, False
    return lex


def variant_map(lines, lex):
    """stripped name -> cvcl for target lines, keeping only names that no other human line shares."""
    for c, n in zip(lines.cvcl, lines.cell_line):
        lex.setdefault(strip(n), set()).add(c)
    tset = set(lines.cvcl)
    return {k: next(iter(cs)) for k, cs in lex.items() if len(cs) == 1 and next(iter(cs)) in tset}


# ----------------------------------------------------------------------------- GEO
def geo_summaries():
    fn = cache_path("geo_h3k27ac_gsm.jsonl")
    if os.path.exists(fn):
        return [json.loads(l) for l in open(fn)]
    q = urllib.parse.urlencode({"db": "gds", "term": QUERY, "usehistory": "y", "retmax": 0, "retmode": "json"})
    es = json.loads(get(f"{EUTILS}/esearch.fcgi?{q}"))["esearchresult"]
    n, web, key = int(es["count"]), es["webenv"], es["querykey"]
    print(f"[17] GEO: {n:,} human H3K27ac GSMs", file=sys.stderr)
    recs = []
    for start in range(0, n, 400):
        u = f"{EUTILS}/esummary.fcgi?db=gds&query_key={key}&WebEnv={web}&retstart={start}&retmax=400&retmode=json"
        res = json.loads(get(u))["result"]
        for uid in res.get("uids", []):
            r = res[uid]
            srx = [x["targetobject"] for x in r.get("extrelations", []) if x.get("targetobject", "").startswith(("SRX", "ERX", "DRX"))]
            recs.append({"gsm": r["accession"], "title": r.get("title", ""), "summary": r.get("summary", ""),
                         "taxon": r.get("taxon", ""), "pdat": r.get("pdat", ""), "gse": r.get("gse", ""),
                         "srx": ";".join(srx)})
        time.sleep(0.4)
        if start % 4000 == 0:
            print(f"[17]   {start + 400:,}/{n:,}", file=sys.stderr)
    with open(fn, "w") as fh:
        for r in recs:
            fh.write(json.dumps(r) + "\n")
    return recs


def ngrams(text, k=4):
    t = [x for x in re.split(r"[\s_\-,;:()\[\]/|.+=]+", str(text)) if x]
    for i in range(len(t)):
        for j in range(i + 1, min(i + k, len(t)) + 1):
            yield strip("".join(t[i:j]))


LEX = {}


def match_text(text, vmap, k=5):
    """Longest match against the whole human lexicon, left to right: a name inside a longer name
    (SKN in SK-N-SH, HEC1 in HEC-1-B) never counts."""
    t = [x for x in re.split(r"[\s_\-,;:()\[\]/|.+=]+", str(text)) if x]
    hits, i = set(), 0
    while i < len(t):
        best = None
        for j in range(min(i + k, len(t)), i, -1):
            g = strip("".join(t[i:j]))
            if g in LEX or g in vmap:
                best = (j, g)
                break
        if best:
            j, g = best
            if g in vmap:
                hits.add(vmap[g])
            i = j
        else:
            i += 1
    return hits


def full_record(gsm):
    fn = cache_path(os.path.join("geo_gsm", f"{gsm}.txt"))
    os.makedirs(os.path.dirname(fn), exist_ok=True)
    if not os.path.exists(fn):
        txt = get(f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gsm}&targ=self&form=text&view=brief")
        open(fn, "w").write(txt)
        time.sleep(0.35)
    txt = open(fn).read()
    f = lambda k: re.findall(r"^!Sample_" + k + r" = (.*)$", txt, re.M)
    ch = f("characteristics_ch1")
    return {"title": " ".join(f("title")), "source": " ".join(f("source_name_ch1")), "chars": ch,
            "organism": " ".join(f("organism_ch1")),
            "srx": ";".join(sorted(set(re.findall(r"([SED]RX\d+)", " ".join(f("relation"))))))}


def verify(rec, cvcl, vmap):
    cell_vals = [rec["source"]] + [c.split(":", 1)[1] for c in rec["chars"]
                                   if ":" in c and CELL_KEYS.match(c.split(":", 1)[0].strip())]
    names = {k for k, c in vmap.items() if c == cvcl}
    norm = lambda v: strip(re.sub(r"\b(cells?|cell line|parental|wild[- ]?type|wt|human|line)\b", " ", v, flags=re.I))
    exact = any(norm(v) in names for v in cell_vals)
    longn = {k for k in names if len(k) > 4}
    in_cell_field = exact or any(cvcl in match_text(v, {k: vmap[k] for k in longn}) for v in cell_vals)
    k27 = bool(K27.search(rec["title"] + " " + " ".join(rec["chars"])))
    return in_cell_field, k27, "Homo sapiens" in rec["organism"]


# ----------------------------------------------------------------------------- ChIP-Atlas
def chip_atlas(path):
    cols = ["srx", "genome", "antigen_class", "antigen", "cell_class", "cell"]
    e = pd.read_csv(path, sep="\t", header=None, usecols=range(6), names=cols, dtype=str,
                    on_bad_lines="skip", quoting=3)
    out = {}
    for s, g, ag, c in zip(e.srx, e.genome, e.antigen, e.cell):
        out.setdefault(s, []).append((g, ag, c))
    return out


def status(srxs, ca):
    rows = [x for s in srxs.split(";") if s for x in ca.get(s, [])]
    if not rows:
        return "absent"
    human = [r for r in rows if r[0] in ("hg38", "hg19")]
    if not human:
        return "non-human genome only"
    if any(r[1] == "H3K27ac" for r in human):
        return "hg38 H3K27ac"
    return f"hg38 as {human[0][1]}"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--experiment-list", required=True, help="current ChIP-Atlas experimentList.tab")
    a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)

    have = set(pd.read_csv(os.path.join(P1, "cell_line_cn_status_2026-08-15.tsv"), sep="\t").cvcl) \
        | set(pd.read_csv(os.path.join(P2, "pull_set.v2.tsv"), sep="\t").cvcl.dropna()) \
        | set(pd.read_csv(os.path.join(P1, "expansion_v3_candidates.tsv"), sep="\t").cvcl.dropna())
    ps = pd.read_csv(os.path.join(P2, "pull_set.v2.tsv"), sep="\t")
    have_models = set(ps.key.dropna()) | set(ps.model_id.dropna())
    cn = cn_lines()
    tgt = cn[~cn.cvcl.isin(have) & ~cn.model_id.isin(have_models)].reset_index(drop=True)
    print(f"[17] CN-measured cancer lines {len(cn):,}; with ChIP-Atlas H3K27ac or candidates "
          f"{int(cn.cvcl.isin(have).sum()):,}; targets {len(tgt):,}", file=sys.stderr)
    global LEX
    LEX = human_lexicon()
    vmap = variant_map(tgt, LEX)
    short = {s for s in vmap if len(s) <= 4}
    print(f"[17] {len(vmap):,} unambiguous name variants ({len(short):,} short, cell-field only)", file=sys.stderr)

    recs = geo_summaries()
    ca = chip_atlas(a.experiment_list)
    allst = pd.DataFrame([{"gsm": r["gsm"], "pdat": r["pdat"], "srx": r["srx"],
                           "status": status(r["srx"], ca) if r["srx"] else "no SRA link"} for r in recs])
    allst.to_csv(os.path.join(OUT, "geo_h3k27ac_scan.all_gsm_status.tsv"), sep="\t", index=False)
    print("[17] all human H3K27ac GSMs vs ChIP-Atlas:", allst.status.value_counts().to_dict(), file=sys.stderr)

    long_map = {s: c for s, c in vmap.items() if s not in short}
    cand = []
    for r in recs:
        hits = match_text(r["title"] + " ; " + r["summary"], long_map) | \
            match_text(r["title"] + " ; " + r["summary"], {s: vmap[s] for s in short})
        for c in hits:
            cand.append((r, c))
    print(f"[17] first-pass matches: {len(cand):,} (GSM, line) pairs on {len({c for _, c in cand}):,} lines; "
          f"verifying on full records", file=sys.stderr)

    out = []
    for i, (r, c) in enumerate(cand):
        rec = full_record(r["gsm"])
        cell_ok, k27, hs = verify(rec, c, vmap)
        out.append({"gsm": r["gsm"], "gse": r["gse"], "pdat": r["pdat"], "cvcl": c, "title": rec["title"][:120],
                    "srx": rec["srx"] or r["srx"], "cell_field_match": cell_ok, "h3k27ac": k27, "human": hs,
                    "chip_atlas": status(rec["srx"] or r["srx"], ca)})
        if i % 200 == 0:
            print(f"[17]   verified {i:,}/{len(cand):,}", file=sys.stderr)
    s = pd.DataFrame(out).merge(tgt, on="cvcl")
    s["verified"] = s.cell_field_match & s.h3k27ac & s.human
    s.to_csv(os.path.join(OUT, "geo_h3k27ac_scan.samples.tsv"), sep="\t", index=False)
    v = s[s.verified].drop_duplicates(["gsm", "cvcl"])
    L = v.groupby(["cvcl", "cell_line", "lineage", "disease", "cn_record"]).agg(
        n_gsm=("gsm", "nunique"), n_gse=("gse", lambda x: len({g.split(";")[0] for g in x})),
        absent=("chip_atlas", lambda x: int((x == "absent").sum())),
        misfiled=("chip_atlas", lambda x: int((x == "non-human genome only").sum())),
        ca_other=("chip_atlas", lambda x: int(x.str.startswith("hg38").sum())),
        first=("pdat", "min"), last=("pdat", "max")).reset_index().sort_values("n_gsm", ascending=False)
    L.to_csv(os.path.join(OUT, "geo_h3k27ac_scan.lines.tsv"), sep="\t", index=False)
    print(f"[17] verified: {len(v):,} samples on {v.cvcl.nunique():,} target lines; lines with >= 2 GEO series "
          f"{int((L.n_gse >= 2).sum())}; ChIP-Atlas status {v.chip_atlas.value_counts().to_dict()}", file=sys.stderr)


if __name__ == "__main__":
    main()
