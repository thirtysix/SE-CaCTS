#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Phase 1 — 18: the H3K27ac experiment registry, one row per experiment x cell line (FINDINGS §19-20).

Every H3K27ac experiment this project has considered, from every route, with where it stands:

  resolved       ChIP-Atlas hg38 H3K27ac resolved to a line in phase 1 (phase1_manifest.tsv)
  rescue         ChIP-Atlas entries rescued from mislabelled antigen or GEO-title cell names (v3 candidates)
  identity-gap   lines ChIP-Atlas resolves whose copy number DepMap files under a parent or subclone
  geo-new        GEO samples on CN-measured lines with no ChIP-Atlas H3K27ac (17_geo_h3k27ac_scan.py)
  geo-have       GEO samples on lines we already have (17_geo_h3k27ac_scan.py --scope have)

Columns: identity (srx, gsm, gse, release date, library strategy), ChIP-Atlas status in the current list
(genome, antigen and cell label as filed), line (cvcl, name, lineage, route, scoring key, CN note),
progress (in the v2 pull set, in the scored v2 atlas, manifest it was pulled under, pull status), QC and
signal-check results, and a decision with its reason. No network access: it reads the cached/derived
tables only, so it can be rebuilt at any time; GEO is re-queried only by 17_geo_h3k27ac_scan.py --refresh.

  python3 phase1/scripts/18_h3k27ac_registry.py [--pulled phase2/data/pulled_srx.tsv]
Output: phase1/data/h3k27ac_registry.tsv (+ a summary on stderr)
"""
import argparse
import glob
import gzip
import os
import re
import sys

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, ROOT)
from secacts_env import cache_path                               # noqa: E402

P1 = os.path.join(ROOT, "phase1", "data")
P2 = os.path.join(ROOT, "phase2", "data")
OUT = os.path.join(ROOT, "phase2", "analysis", "out")


def rd(path, **kw):
    return pd.read_csv(path, sep="\t", **kw) if os.path.exists(path) else pd.DataFrame()


def explode_srx(df, col="srx"):
    df = df.copy()
    df[col] = df[col].fillna("").astype(str).str.split(";")
    df = df.explode(col)
    return df


def geo_strategy(gsm):
    fn = cache_path(os.path.join("geo_gsm", f"{gsm}.txt"))
    if not os.path.exists(fn):
        return ""
    m = re.findall(r"^!Sample_library_strategy = (.*)$", open(fn).read(), re.M)
    return m[0] if m else ""


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--experiment-list", default=cache_path("experimentList.2026-09.tab"))
    ap.add_argument("--pulled", default=os.path.join(P2, "pulled_srx.tsv"),
                    help="srx<TAB>status (done|failed) from Roihu out/ (see nextsession for the command)")
    ap.add_argument("--out", default=os.path.join(P1, "h3k27ac_registry.tsv"))
    a = ap.parse_args()

    # ---- line identity tables
    st = rd(os.path.join(P1, "cell_line_cn_status_2026-08-15.tsv"))
    names = dict(zip(st.cvcl, st.cell_line))
    cnstat = dict(zip(st.cvcl, st.cn_status))
    lineage = dict(zip(st.cvcl, st.lineage))

    rows = []
    # resolved (ChIP-Atlas, phase 1)
    man = rd(os.path.join(P1, "phase1_manifest.tsv"))
    for r in man.itertuples():
        rows.append(dict(srx=r.srx, cvcl=r.cvcl, route="resolved", qc_pass_meta=r.qc_pass, n_peaks_meta=r.n_peaks))
        names.setdefault(r.cvcl, r.cell)
    # rescue (v3 candidates)
    cand = rd(os.path.join(P1, "expansion_v3_candidates.tsv"))
    for r in cand.itertuples():
        rows.append(dict(srx=r.srx, cvcl=r.cvcl, route=f"rescue-{r.route}", qc_pass_meta=int(bool(r.qc_pass)),
                         n_peaks_meta=r.n_peaks))
    # identity gap
    idg = rd(os.path.join(OUT, "identity_gap_lines.tsv"))
    for r in idg.itertuples():
        for s in str(r.srx).split(";"):
            if s:
                rows.append(dict(srx=s, cvcl=r.cvcl, route="identity-gap", cn_note=f"CN from {r.relation} "
                                 f"{r.relative} ({r.model_id}); {r.handling}"))
                names.setdefault(r.cvcl, r.cell_line); lineage.setdefault(r.cvcl, r.lineage)
    # GEO scans
    for sfx, route in (("", "geo-new"), (".have", "geo-have")):
        g = rd(os.path.join(OUT, f"geo_h3k27ac_scan{sfx}.samples.tsv"))
        if g.empty:
            continue
        g = g[g.verified]
        for r in explode_srx(g).itertuples():
            rows.append(dict(srx=r.srx or f"(no SRA) {r.gsm}", gsm=r.gsm, gse=str(r.gse).split(";")[0],
                             geo_release=r.pdat, cvcl=r.cvcl, route=route))
            names.setdefault(r.cvcl, r.cell_line); lineage.setdefault(r.cvcl, r.lineage)
    R = pd.DataFrame(rows)
    for c in ("gsm", "gse", "geo_release", "cn_note", "qc_pass_meta", "n_peaks_meta"):
        if c not in R:
            R[c] = pd.NA

    # one row per (srx, cvcl): keep every route that found it
    agg = {c: "first" for c in R.columns if c not in ("srx", "cvcl", "route")}
    agg["route"] = lambda x: ",".join(sorted(set(x)))
    R = R.groupby(["srx", "cvcl"], as_index=False).agg(agg)
    R["cell_line"] = R.cvcl.map(names)
    R["lineage"] = R.cvcl.map(lineage)

    # ---- GEO identity for ChIP-Atlas rows (from the cached GEO summaries)
    gj = cache_path("geo_h3k27ac_gsm.jsonl")
    if os.path.exists(gj):
        geo = pd.read_json(gj, lines=True)
        geo = explode_srx(geo)[["srx", "gsm", "gse", "pdat"]].drop_duplicates("srx").set_index("srx")
        miss = R.gsm.isna() if "gsm" in R else pd.Series(True, index=R.index)
        R.loc[miss, "gsm"] = R.loc[miss, "srx"].map(geo.gsm)
        R.loc[miss, "gse"] = R.loc[miss, "srx"].map(geo.gse).astype(str).str.split(";").str[0]
        R.loc[miss, "geo_release"] = R.loc[miss, "srx"].map(geo.pdat)
    R["library_strategy"] = R.gsm.fillna("").map(lambda g: geo_strategy(g) if g else "")

    # ---- ChIP-Atlas as filed now
    cols = ["srx", "genome", "antigen_class", "antigen", "cell_class", "cell"]
    e = pd.read_csv(a.experiment_list, sep="\t", header=None, usecols=range(6), names=cols, dtype=str,
                    quoting=3, on_bad_lines="skip")
    hg = e[e.genome.isin(["hg38"])].drop_duplicates("srx").set_index("srx")
    anyg = e.groupby("srx").genome.agg(lambda x: ",".join(sorted(set(x))))
    R["chip_atlas_genomes"] = R.srx.map(anyg).fillna("")
    R["chip_atlas_antigen"] = R.srx.map(hg.antigen)
    R["chip_atlas_cell"] = R.srx.map(hg.cell)
    R["chip_atlas_status"] = [
        "absent" if not g else ("hg38" if "hg38" in g else "non-human genome only")
        for g in R.chip_atlas_genomes]

    # ---- progress
    ps = rd(os.path.join(P2, "pull_set.v2.tsv"))
    R["in_v2_pull"] = R.srx.isin(set(ps.srx))
    sig = os.path.join(ROOT, "phase2", "results_v2", "atlas.s3.se_signal.tsv.gz")
    atlas = set(gzip.open(sig, "rt").readline().rstrip("\n").split("\t")[1:]) if os.path.exists(sig) else set()
    R["in_v2_atlas"] = R.srx.isin(atlas)
    R["key"] = R.srx.map(ps.drop_duplicates("srx").set_index("srx").key)
    man_of = {}
    for f in sorted(glob.glob(os.path.join(P2, "pull_srx.v3*.txt"))):
        for s in open(f).read().split():
            man_of.setdefault(s, os.path.basename(f).replace("pull_srx.", "").replace(".txt", ""))
    R["v3_manifest"] = R.srx.map(man_of).fillna("")
    pulled = rd(a.pulled, header=None, names=["srx", "status"])
    R["pull_status"] = R.srx.map(pulled.drop_duplicates("srx").set_index("srx").status if len(pulled) else {}) \
        .fillna(R.in_v2_pull.map({True: "done (v2)", False: ""}))

    # ---- QC and signal check (every *candidate_check*.tsv in phase2/analysis/out)
    chk = pd.concat([rd(f) for f in glob.glob(os.path.join(OUT, "*candidate_check*.tsv"))], ignore_index=True)
    if len(chk):
        chk = chk.drop_duplicates("srx", keep="last").set_index("srx")
        for c in ("n_peaks", "dyn_range", "r_global", "own_line_rank", "best_line", "status", "flags"):
            if c in chk:
                R[f"check_{c}"] = R.srx.map(chk[c])

    # ---- decision (first pass; Jev labelling and review can overturn it)
    def decide(r):
        if r.in_v2_atlas:
            return "in atlas v2", ""
        if r.in_v2_pull:
            return "excluded v2", "QC gate (< 2,000 peaks) or failed pull"
        if isinstance(r.get("check_flags"), str) and r.check_flags:
            return "exclude", r.check_flags
        if r.get("check_status") == "pass":
            return "v3 include", "signal check passed; pending Jev labels"
        if r.chip_atlas_status == "absent":
            if r.library_strategy and r.library_strategy != "ChIP-Seq":
                return "exclude", f"{r.library_strategy} (CUT&Tag/CUT&RUN, different assay)"
            if str(r.srx).startswith("(no SRA)"):
                return "exclude", "no raw data in SRA"
            return "needs own processing", "not in ChIP-Atlas (Stage B)"
        if r.chip_atlas_status == "non-human genome only":
            return "needs own processing", "ChIP-Atlas filed it under a spike-in genome only"
        if r.v3_manifest and not r.pull_status:
            return "pulling", r.v3_manifest
        if r.pull_status == "done":
            return "pulled, check pending", r.v3_manifest
        if r.qc_pass_meta == 0:
            return "excluded", "ChIP-Atlas QC (reads/peaks) failed"
        if "resolved" in r.route and r.cvcl not in cnstat:
            return "out of scope", "line is not one of the 631 human cancer lines (non-cancer or excluded identity)"
        if cnstat.get(r.cvcl, "").startswith(("lookup", "infer", "none")) and "identity-gap" not in r.route:
            return "not pulled", f"line has no copy-number track ({cnstat[r.cvcl]})"
        return "not yet pulled", "in ChIP-Atlas hg38"
    R[["decision", "reason"]] = R.apply(lambda r: pd.Series(decide(r)), axis=1)

    order = ["srx", "gsm", "gse", "geo_release", "library_strategy", "cvcl", "cell_line", "lineage", "route",
             "key", "cn_note", "chip_atlas_status", "chip_atlas_genomes", "chip_atlas_antigen", "chip_atlas_cell",
             "qc_pass_meta", "n_peaks_meta", "in_v2_pull", "in_v2_atlas", "v3_manifest", "pull_status"] + \
            [c for c in R.columns if c.startswith("check_")] + ["decision", "reason"]
    R = R[[c for c in order if c in R.columns]].sort_values(["decision", "lineage", "cell_line", "srx"])
    R.to_csv(a.out, sep="\t", index=False)
    print(f"[18] registry: {len(R):,} experiment x line rows, {R.srx.nunique():,} experiments, "
          f"{R.cvcl.nunique():,} lines -> {os.path.relpath(a.out, ROOT)}", file=sys.stderr)
    print(R.groupby("decision").agg(rows=("srx", "size"), lines=("cvcl", "nunique")).to_string(), file=sys.stderr)


if __name__ == "__main__":
    main()
