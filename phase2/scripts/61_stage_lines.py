#!/usr/bin/env python3
"""Phase 2 — 61: per-line data for the dashboard's Genomic View (igv.js) and the SE Atlas line level.

For each line (scoring key = DepMap ModelID, or CVCL for lines DepMap does not hold):
  docs/data/lines/<key>.json           summary: labels, experiments (SRX, study, layout, ChIP-Atlas bigWig URL), and four
                                       comparisons, each with its comparison set, whether it could be tested (and why
                                       not), the number of specific SEs, and the top TOP of them (compact rows: `cols`)
  docs/data/lines/<key>.called.bed.gz  union SE loci called in any of the line's experiments
                                       (score = share of its experiments that called it, x1000)
  docs/data/lines/<key>.vs{all,sub,dis,lin}.bed.gz   specific SEs (FDR <= 0.1, up to TRACK_MAX by rank) vs all lines,
                                       vs the lines of the same subtype, primary disease, lineage
  docs/data/lines/<key>.cn.bedgraph.gz log2 CN ratio as scored (the line's CN source), runs of equal value merged
  docs/data/lines/index.json           the lines that have a page

Comparisons come from the `l*` arms of phase2/roihu/score_arm.slurm (consensus of studies, q25; the null shuffles
labels among all lines, or among the lines of the same subtype / disease / lineage; FINDINGS §37). A line is not
tested with one study, with no label at the level, or with fewer than 4 lines in its group.

Coverage is NOT copied: the page streams ChIP-Atlas's own per-experiment bigWigs (CORS + range requests,
verified 2026-09-29), raw RPM, unnormalised.

  python3 phase2/scripts/61_stage_lines.py                 # the prototype set
  python3 phase2/scripts/61_stage_lines.py --lines all     # every scored line
"""
import argparse
import gzip
import json
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "cnrose"))
from secacts_env import DATAROOT, cache_path                     # noqa: E402

P1, P2 = os.path.join(ROOT, "phase1", "data"), os.path.join(ROOT, "phase2", "data")
RES = os.environ.get("SECACTS_RES", os.path.join(ROOT, "phase2", "results_v2"))
SC = os.environ.get("SECACTS_SC", os.path.join(ROOT, "phase2", "scores_v2"))
DOCS = os.environ.get("SECACTS_DOCS", os.path.join(ROOT, "docs"))   # stage a copy, e.g. for screenshots
OUT = os.path.join(DOCS, "data", "lines")
PROTO = ["NIHOVCAR3", "MCF7", "K562", "A549", "KELLY", "JURKAT", "SKNBE2"]
BW = "https://chip-atlas.dbcls.jp/data/hg38/eachData/bw/{}.bw"
TOP, TRACK_MAX, STRATA_MIN = 200, 2000, 4
SC_LINES = os.environ.get("SECACTS_SC_LINES", os.path.join(SC, "out_lines"))
# comparison -> (arm prefix, track tag, level label, stratum column in the labels)
CMP = {"all": ("atlas.s3.lines.all", "vsall", "all lines", None),
       "lineage": ("atlas.s3.lines.lin", "vslin", "lineage", "lineage"),
       "disease": ("atlas.s3.lines.dis", "vsdis", "primary disease", "disease"),
       "subtype": ("atlas.s3.lines.sub", "vssub", "subtype", "subtype")}
FULL = {"vsall": "224,130,20", "vssub": "27,120,55", "vsdis": "192,57,43", "vslin": "123,50,148"}    # SE called in this line
PALE = {"vsall": "246,214,170", "vssub": "178,221,190", "vsdis": "240,190,184", "vslin": "215,190,228"}  # specific signal, not an SE here
# one table per line, rows = the "vs all" basis (every relatives call is also a vs-all call, checked 2026-09-30),
# with one FDR column per comparison: a number when called (<= 0.10), null when tested but not called
COLS = ["rank", "se", "chrom", "start", "end", "gene", "dist_kb", "jsd", "fdr_all", "fdr_lineage", "fdr_disease",
        "fdr_subtype", "cn", "called", "signal_rank", "flag", "pass"]
LETTER = {"all": "A", "lineage": "L", "disease": "D", "subtype": "S"}     # which comparisons call the SE, e.g. "ALD"


def clean(o):
    """NaN (a missing label) is not JSON: write null."""
    if isinstance(o, float) and o != o:
        return None
    if isinstance(o, dict):
        return {k: clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [clean(v) for v in o]
    return o


def flag_of(chrom, cn):
    """Known artifact classes, removed in v3.1 (ROADMAP): chrY presence depends on the line's sex, and dividing by
    a copy number near 0 inflates noise in deep deletions."""
    return "chrY" if chrom == "chrY" else ("CN<0.3" if cn is not None and cn < 0.3 else "")


def protein_coding_coords(gtf, cache):
    """{symbol: (chrom, start, end)} for protein-coding genes only: the nearest-gene label of an SE should
    name a gene (MYCN), not the lncRNA or pseudogene that happens to sit closer."""
    if os.path.exists(cache):
        return {g: (c, int(a), int(b)) for g, c, a, b in (l.rstrip("\n").split("\t") for l in open(cache))}
    out = {}
    with gzip.open(gtf, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t")
            if f[2] != "gene" or 'gene_biotype "protein_coding"' not in f[8]:
                continue
            m = [x for x in f[8].split(";") if x.strip().startswith("gene_name")]
            if not m:
                continue
            g = m[0].split('"')[1]
            c = f[0] if f[0].startswith("chr") else "chr" + f[0]
            out.setdefault(g, (c, int(f[3]) - 1, int(f[4])))
    with open(cache, "w") as fh:
        for g, (c, a, b) in out.items():
            fh.write(f"{g}\t{c}\t{a}\t{b}\n")
    return out


def nearest_genes(loci, coords):
    by = {}
    for g, (c, s, e) in coords.items():
        by.setdefault(c, []).append((s, e, g))
    idx = {c: (np.array([x[0] for x in v]), np.array([x[1] for x in v]), [x[2] for x in v]) for c, v in by.items()}
    out = []
    for c, s, e in loci:
        if c not in idx:
            out.append(("", -1)); continue
        gs, ge, names = idx[c]
        d = np.maximum(0, np.maximum(gs - e, s - ge))
        i = int(np.argmin(d))
        out.append((names[i], int(d[i] // 1000)))
    return out


def overlap_share(query, loci, called, share):
    """For each locus: the largest share of the line's experiments calling ANY union SE that overlaps it. The
    union catalog nests small loci inside large ones (37% of loci overlap another: the 25% reciprocal-overlap
    merge keeps a 2 kb SE apart from the 90 kb SE around it), so membership alone under-counts calls."""
    by = {}
    for se in called:
        c = loci.loc[se]
        by.setdefault(c.chrom, []).append((int(c.start), int(c.end), float(share[se])))
    idx = {c: tuple(np.array(x) for x in zip(*sorted(v))) for c, v in by.items()}
    out = {}
    for se, c in query.iterrows():
        if c.chrom not in idx:
            out[se] = 0.0; continue
        st, en, sh = idx[c.chrom]
        hit = (st < c.end) & (en > c.start)
        out[se] = float(sh[hit].max()) if hit.any() else 0.0
    return out


def cn_segments(track, tol=0.15, max_gap=1_000_000):
    """Gene-level CN -> segments: merge neighbouring genes whose log2 ratio stays within `tol` of the running
    length-weighted mean (DepMap's per-gene values jitter, so exact-equality merging keeps every gene)."""
    rows = []
    for c, (st, en, ra) in track._chrom.items():
        o = np.argsort(st, kind="stable")
        cur = None                                   # [chrom, start, end, sum(v*len), sum(len)]
        for s, e, r in zip(st[o], en[o], ra[o]):
            v, L = float(np.log2(max(r, 1e-3))), max(int(e) - int(s), 1)
            if cur and abs(v - cur[3] / cur[4]) < tol and s <= cur[2] + max_gap:
                cur[2] = max(cur[2], int(e)); cur[3] += v * L; cur[4] += L; continue
            if cur:
                end = int(s) if int(s) - cur[2] <= max_gap else cur[2]      # fill the gap to the next gene
                rows.append((cur[0], cur[1], max(min(end, int(s)), cur[1] + 1), round(cur[3] / cur[4], 3)))
            cur = [c, int(s), int(e), v * L, L]
        if cur:
            rows.append((cur[0], cur[1], cur[2], round(cur[3] / cur[4], 3)))
    return [r for r in rows if r[2] > r[1]]


_STAGE = None


def _stage_worker(key):
    return _STAGE(key)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--lines", default=",".join(PROTO), help="comma-separated line names (stripped) or 'all'")
    a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)

    ps = pd.read_csv(os.environ.get("SECACTS_PS", os.path.join(P2, "pull_set.v2.tsv")), sep="\t")
    cols = gzip.open(os.path.join(RES, "atlas.s3.se_signal.tsv.gz"), "rt").readline().rstrip("\n").split("\t")[1:]
    ps = ps[ps.srx.isin(cols)]
    model = pd.read_csv(os.path.join(DATAROOT, "DepMap", "2026q1", "Model.csv"),
                        usecols=["ModelID", "CellLineName", "StrippedCellLineName", "OncotreeLineage",
                                 "OncotreePrimaryDisease", "OncotreeSubtype"]).set_index("ModelID")
    lin = pd.read_csv(os.path.join(P1, "lineage_resolved.tsv"), sep="\t").drop_duplicates("cvcl").set_index("cvcl")
    man = json.load(open(os.path.join(DOCS, "data", "manifest.json")))
    line_groups = man["levels"]["line"]["groups"]
    study = pd.read_csv(os.path.join(P2, "srx_study.tsv"), sep="\t").set_index("srx").iloc[:, 0]
    layout = pd.read_csv(os.path.join(P2, "srx_layout.tsv"), sep="\t").set_index("srx").layout

    def labels(key):
        if key in model.index:
            m = model.loc[key]
            return m.StrippedCellLineName, m.CellLineName, m.OncotreeLineage, m.OncotreePrimaryDisease, m.OncotreeSubtype
        cv = ps[ps.key == key].cvcl.iloc[0]
        r = lin.loc[cv]
        name = str(r.cell_line)
        return "".join(ch for ch in name.upper() if ch.isalnum()), name, r.lineage, r.primary_disease, r.subtype
    keys = sorted(ps.key.unique())
    group_of = {k: labels(k)[0] for k in keys}
    want = keys if a.lines == "all" else [k for k in keys if group_of[k] in set(a.lines.split(","))]
    print(f"[61] {len(want)} line(s)", file=sys.stderr)

    cat = pd.read_csv(os.path.join(RES, "atlas.s3.union_catalog.bed.gz"), sep="\t", header=None,
                      usecols=[0, 1, 2, 3], names=["chrom", "start", "end", "se"]).set_index("se")
    pres = pd.read_csv(os.path.join(RES, "atlas.s3.se_presence.tsv.gz"), sep="\t", index_col=0,
                       usecols=lambda c: c == "se_id" or c in set(ps[ps.key.isin(want)].srx))
    want_srx = set(ps[ps.key.isin(want)].srx)
    sig = pd.read_csv(os.path.join(RES, "atlas.s3.se_signal.tsv.gz"), sep="\t", index_col=0,
                      usecols=lambda c: c == "se_id" or c in want_srx)
    from cnrose.cn.depmap import load_gene_coords, DepMapGeneCN, DepMapMcWesCN
    from cnrose.cn.cmp import CellModelPassportsWesCN
    coords = load_gene_coords(None, cache_path=cache_path("gene_coords.GRCh38.106.tsv"))
    pc = protein_coding_coords(os.path.join(DATAROOT, "0.human_genome", "Homo_sapiens.GRCh38.106.chr.gtf.gz"),
                               cache_path("gene_coords.GRCh38.106.protein_coding.tsv"))
    genes = dict(zip(cat.index, nearest_genes(cat[["chrom", "start", "end"]].itertuples(index=False), pc)))
    spec, top_any, strata, members = {}, {}, {}, {}
    for c, (pref, tag, _, _) in CMP.items():
        d = pd.read_csv(os.path.join(SC_LINES, f"{pref}.line.specific.tsv.gz"), sep="\t")
        spec[c] = {g: x.sort_values("rank") for g, x in d[d.fdr <= 0.1].groupby("group")}
        t = pd.read_csv(os.path.join(SC_LINES, f"{pref}.line.top_specific.tsv"), sep="\t")
        top_any[c] = {g: x.sort_values("rank") for g, x in t.groupby("group")}       # rankings, for untested lines
        sp = os.path.join(SC_LINES, f"{pref}.line.strata.tsv")
        if os.path.exists(sp):
            st = pd.read_csv(sp, sep="\t", keep_default_na=False)
            strata[c] = st.set_index("group")
            members[c] = st[st.stratum != ""].groupby("stratum").group.apply(frozenset).to_dict()
    H = pd.read_csv(os.path.join(SC_LINES, "atlas.s3.lines.all.hierarchy_summary.tsv"), sep="\t")
    n_units = H[H.level == "line"].set_index("group").n_lines.to_dict()                    # studies per line
    n_panel = int((H.level == "line").sum())

    src = ps.drop_duplicates("key").set_index("key")
    D = os.path.join(DATAROOT, "DepMap", "2026q1")
    wgs = DepMapGeneCN(os.path.join(D, "OmicsCNGeneWGS.csv"), coords)
    wgs.preload([k for k in want if src.cn_provider.get(k) == "depmap_wgs"])
    mcw = cmp = None
    if any(src.cn_provider.get(k) == "depmap_mc_wes" for k in want):
        mcw = DepMapMcWesCN(os.path.join(D, "OmicsCNGeneMC_WES.csv"), os.path.join(D, "ModelCondition.csv"), coords)
        mcw.preload([k for k in want if src.cn_provider.get(k) == "depmap_mc_wes"])
    if any(src.cn_provider.get(k) == "cmp_wes" for k in want):
        C = os.path.join(DATAROOT, "CellModelPassports")
        cmp = CellModelPassportsWesCN(os.path.join(C, "WES_pureCN_CNV_genes_latest.csv.gz"),
                                      os.path.join(C, "model_list_20240110.csv"), cache_dir=cache_path("cmp_wes"))

    def comparison(c, grp):
        """The comparison set of `grp` at level `c`, and whether it could be tested (and if not, why)."""
        nu = int(n_units.get(grp, 0))
        if c == "all":
            out = {"stratum": None, "n_lines": n_panel}
        else:
            st = strata[c].loc[grp] if grp in strata[c].index else None
            label = "" if st is None else str(st.stratum)
            out = {"stratum": label or None, "n_lines": 0 if st is None or not label else int(st.lines_in_stratum)}
        reason = None
        if nu < 2:
            reason = "one study only: calls need two independent studies"
        elif c != "all" and not out["stratum"]:
            reason = f"no {CMP[c][2]} label for this line"
        elif c != "all" and out["n_lines"] < STRATA_MIN:
            reason = f"only {out['n_lines']} line{'s' if out['n_lines'] != 1 else ''} in its {CMP[c][2]} (at least {STRATA_MIN} needed)"
        out.update(level=CMP[c][2], testable=reason is None, reason=reason, same_as=None)
        return out

    def stage_one(key):
        grp, name, lineage, disease, subtype = labels(key)
        exps = ps[ps.key == key].srx.tolist()
        called = pres[[s for s in exps if s in pres.columns]]
        share_all = called.mean(axis=1)
        share = share_all[share_all > 0]
        n_exp = called.shape[1]
        mean_sig = sig[[c for c in exps if c in sig.columns]].mean(axis=1)
        srank = mean_sig.loc[share.index].rank(ascending=False, method="first").astype(int)   # 1 = strongest called SE
        fn = lambda ext: os.path.join(OUT, f"{key}.{ext}")                                   # noqa: E731
        gz = lambda ext: gzip.open(fn(ext + ".gz"), "wt")                                    # noqa: E731
        with gz("called.bed") as fh:
            for se, v in share.items():
                c = cat.loc[se]
                # no spaces in the name: igv.js splits BED on any whitespace, which would shift the colour column
                fh.write(f"{c.chrom}\t{c.start}\t{c.end}\t{genes[se][0] or se}[{srank[se]}]\t{int(round(300 + 700 * v))}\n")
        spec_se = set()
        for c in CMP:
            if grp in spec[c]:
                spec_se |= set(spec[c][grp].se.head(TRACK_MAX))
            if grp in top_any[c]:
                spec_se |= set(top_any[c][grp].se.head(TOP))
        ov = overlap_share(cat.loc[sorted(spec_se)], cat, list(share.index), share) if spec_se else {}
        comps = {}
        infos = {c: comparison(c, grp) for c in CMP}
        called_in = {c: (spec[c][grp].set_index("se").fdr if infos[c]["testable"] and grp in spec[c] else pd.Series(dtype=float))
                     for c in CMP}
        pass_of = lambda se: "".join(LETTER[c] for c in CMP if se in called_in[c].index)    # noqa: E731
        for c, (_, tag, _, _) in CMP.items():                                           # tracks: every call, to TRACK_MAX
            d = spec[c].get(grp) if infos[c]["testable"] else None
            with gz(f"{tag}.bed") as fh:
                if d is not None:
                    for r in d.head(TRACK_MAX).itertuples():
                        cc = cat.loc[r.se]; g, _ = genes[r.se]
                        rgb = FULL[tag] if ov.get(r.se, 0.0) > 0 else PALE[tag]
                        fh.write(f"{cc.chrom}\t{cc.start}\t{cc.end}\t{g}[{int(r.rank)}]\t{int(1000 * (1 - min(r.fdr, 1)))}"
                                 f"\t.\t{cc.start}\t{cc.end}\t{rgb}\n")
            infos[c].update(n=0 if d is None else int(len(d)),
                            n_called_here=0 if d is None else int(sum(ov.get(x, 0) > 0 for x in d.se)))
            comps[c] = infos[c]
        # the table: the vs-all calls by rank (or, untested, the top of the vs-all ranking)
        base = spec["all"].get(grp) if infos["all"]["testable"] else None
        base = base if base is not None else top_any["all"].get(grp)
        rows = []
        if base is not None:
            for r in base.head(TOP).itertuples():
                cc = cat.loc[r.se]; g, dist = genes[r.se]
                cn = None if pd.isna(r.cn_mean) else round(float(r.cn_mean), 2)
                fdrs = [(round(float(called_in[c][r.se]), 4) if r.se in called_in[c].index else None) for c in CMP]
                rows.append([int(r.rank), r.se, cc.chrom, int(cc.start), int(cc.end), g, dist, round(float(r.jsd), 4),
                             *fdrs, cn, int(round(ov.get(r.se, 0.0) * n_exp)),
                             int(srank[r.se]) if r.se in srank.index else None, flag_of(cc.chrom, cn), pass_of(r.se)])
        # the same comparison set at two levels (a subtype that IS its disease): say so rather than repeat silently
        for lo, hi in (("subtype", "disease"), ("disease", "lineage")):
            a_, b_ = comps[lo], comps[hi]
            if a_["stratum"] and b_["stratum"] and members[lo].get(a_["stratum"]) == members[hi].get(b_["stratum"]):
                a_["same_as"] = hi
        prov = src.cn_provider.get(key, "depmap_wgs")
        tr = (mcw.track(key) if prov == "depmap_mc_wes" else
              cmp.track(src.cvcl.get(key)) if prov == "cmp_wes" else wgs.track(key))
        cn_max = 3
        if tr is not None:
            segs = cn_segments(tr)
            cn_max = int(min(8, max(3, np.ceil(max(v for *_, v in segs)))))
            with gz("cn.bedgraph") as fh:
                fh.write("track type=bedGraph\n")
                for c, st_, en, v in segs:
                    fh.write(f"{c}\t{st_}\t{en}\t{v}\n")
        for old in ("vsrel.bed.gz",):                                                          # pre-v3 file name
            if os.path.exists(fn(old)):
                os.remove(fn(old))
        summary = {
            "key": key, "group": grp, "name": name, "lineage": lineage, "disease": disease, "subtype": subtype,
            "cn_source": prov, "has_cn": tr is not None, "cn_max": cn_max,
            "n_called": int(len(share)), "n_experiments": int(n_exp), "n_studies": int(n_units.get(grp, 0)),
            "experiments": [{"srx": s_, "study": study.get(s_, ""), "layout": layout.get(s_, ""), "bw": BW.format(s_)}
                            for s_ in exps],
            "cols": COLS, "rows": rows, "comparisons": comps,
        }
        json.dump(clean(summary), open(fn("json"), "w"), separators=(",", ":"), allow_nan=False)
        return key, {"group": grp, "name": name, "lineage": lineage, "search": line_groups.get(grp, {}).get("search", name),
                     "n": {c: comps[c]["n"] for c in CMP}}

    index = {}
    workers = int(os.environ.get("STAGE_WORKERS", 8))
    global _STAGE
    _STAGE = stage_one                                   # forked workers inherit it; a closure cannot be pickled
    if workers > 1 and len(want) > 1:
        import multiprocessing as mp
        with mp.get_context("fork").Pool(workers) as pool:
            done = pool.map(_stage_worker, want, chunksize=4)
    else:
        done = [stage_one(k) for k in want]
    for k, v in done:
        index[k] = v
    json.dump(clean(index), open(os.path.join(OUT, "index.json"), "w"), separators=(",", ":"), sort_keys=True, allow_nan=False)
    print(f"[61] staged {len(done)} line(s) -> {OUT}", file=sys.stderr)

if __name__ == "__main__":
    main()
