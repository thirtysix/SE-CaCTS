#!/usr/bin/env python3
"""Phase 2 — 61: per-line data for the dashboard's Line view (igv.js), FINDINGS 2026-09-29 §27.

For each line (scoring key = DepMap ModelID, or CVCL for lines DepMap does not hold):
  docs/data/lines/<key>.json          summary: labels, experiments (SRX, study, layout, ChIP-Atlas bigWig URL),
                                      the comparison set for "vs relatives", and the top specific SEs in both
                                      senses, each with locus and nearest gene
  docs/data/lines/<key>.called.bed.gz union SE loci called in any of the line's experiments
                                      (score = share of its experiments that called it, x1000)
  docs/data/lines/<key>.vsall.bed.gz  specific vs the whole panel (consensus, studies as members; FDR <= 0.1)
  docs/data/lines/<key>.vsrel.bed.gz  specific vs relatives (within disease or lineage; FDR <= 0.1)
  docs/data/lines/<key>.cn.bedgraph.gz log2 CN ratio as scored (the line's CN source), runs of equal value merged
  docs/data/lines/index.json          the lines that have a page

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
RES, SC = os.path.join(ROOT, "phase2", "results_v2"), os.path.join(ROOT, "phase2", "scores_v2")
OUT = os.path.join(ROOT, "docs", "data", "lines")
PROTO = ["NIHOVCAR3", "MCF7", "K562", "A549", "KELLY", "JURKAT", "SKNBE2"]
BW = "https://chip-atlas.dbcls.jp/data/hg38/eachData/bw/{}.bw"
TOP = 50
FULL = {"vsall": "224,130,20", "vsrel": "192,57,43"}      # SE called in this line
PALE = {"vsall": "246,214,170", "vsrel": "240,190,184"}  # specific signal, but no experiment of the line calls an SE


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


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--lines", default=",".join(PROTO), help="comma-separated line names (stripped) or 'all'")
    a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)

    ps = pd.read_csv(os.path.join(P2, "pull_set.v2.tsv"), sep="\t")
    cols = gzip.open(os.path.join(RES, "atlas.s3.se_signal.tsv.gz"), "rt").readline().rstrip("\n").split("\t")[1:]
    ps = ps[ps.srx.isin(cols)]
    model = pd.read_csv(os.path.join(DATAROOT, "DepMap", "2026q1", "Model.csv"),
                        usecols=["ModelID", "CellLineName", "StrippedCellLineName", "OncotreeLineage",
                                 "OncotreePrimaryDisease", "OncotreeSubtype"]).set_index("ModelID")
    lin = pd.read_csv(os.path.join(P1, "lineage_resolved.tsv"), sep="\t").drop_duplicates("cvcl").set_index("cvcl")
    man = json.load(open(os.path.join(ROOT, "docs", "data", "manifest.json")))
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
    spec = {}
    for tag, f in (("vsall", "atlas.s3.conss.line.specific.tsv.gz"), ("vsrel", "atlas.s3.hybl.line.specific.tsv.gz")):
        d = pd.read_csv(os.path.join(SC, f), sep="\t")
        spec[tag] = {g: x.sort_values("rank") for g, x in d[d.fdr <= 0.1].groupby("group")}
    strata = pd.read_csv(os.path.join(SC, "atlas.s3.hybl.line.strata.tsv"), sep="\t").set_index("group")

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

    index = json.load(open(os.path.join(OUT, "index.json"))) if os.path.exists(os.path.join(OUT, "index.json")) else {}
    for key in want:
        grp, name, lineage, disease, subtype = labels(key)
        exps = ps[ps.key == key].srx.tolist()
        called = pres[[s for s in exps if s in pres.columns]]
        share_all = called.mean(axis=1)
        share = share_all[share_all > 0]
        n_exp = called.shape[1]
        mean_sig = sig[[c for c in exps if c in sig.columns]].mean(axis=1)
        srank = mean_sig.loc[share.index].rank(ascending=False, method="first").astype(int)   # 1 = strongest called SE
        fn = lambda ext: os.path.join(OUT, f"{key}.{ext}")
        gz = lambda ext: gzip.open(fn(ext + ".gz"), "wt")
        with gz("called.bed") as fh:
            for se, v in share.items():
                c = cat.loc[se]
                # no spaces in the name: igv.js splits BED on any whitespace, which would shift the colour column
                fh.write(f"{c.chrom}\t{c.start}\t{c.end}\t{genes[se][0] or se}[{srank[se]}]\t{int(round(300 + 700 * v))}\n")   # floor: a rarely called SE stays visible
        tops = {}
        spec_se = set().union(*[set(spec[t][grp].se) for t in ("vsall", "vsrel") if grp in spec[t]])
        ov = overlap_share(cat.loc[sorted(spec_se)], cat, list(share.index), share) if spec_se else {}
        for tag in ("vsall", "vsrel"):
            d = spec[tag].get(grp)
            rows = []
            with gz(f"{tag}.bed") as fh:
                if d is not None:
                    for r in d.itertuples():
                        c = cat.loc[r.se]
                        g, dist = genes[r.se]
                        sh = ov.get(r.se, 0.0)             # overlap-based: called here if any called SE covers it
                        rgb = FULL[tag] if sh > 0 else PALE[tag]       # pale = specific signal, but not an SE here
                        fh.write(f"{c.chrom}\t{c.start}\t{c.end}\t{g}[{int(r.rank)}]\t{int(1000 * (1 - min(r.fdr, 1)))}"
                                 f"\t.\t{c.start}\t{c.end}\t{rgb}\n")
                        if len(rows) < TOP:
                            rows.append({"rank": int(r.rank), "se": r.se, "chrom": c.chrom, "start": int(c.start),
                                         "end": int(c.end), "gene": g, "dist_kb": dist, "jsd": round(float(r.jsd), 4),
                                         "fdr": round(float(r.fdr), 4), "cn": round(float(r.cn_mean), 2),
                                         "called": int(round(sh * n_exp)), "signal_rank": int(srank[r.se]) if r.se in srank.index else None,
                                         "member": r.se in share.index})
            n_called = 0 if d is None else int(sum(ov.get(x, 0) > 0 for x in d.se))
            tops[tag] = {"n": 0 if d is None else len(d), "n_called_here": n_called, "top": rows}
        prov = src.cn_provider.get(key, "depmap_wgs")
        tr = (mcw.track(key) if prov == "depmap_mc_wes" else
              cmp.track(src.cvcl.get(key)) if prov == "cmp_wes" else wgs.track(key))
        cn_max = 3
        if tr is not None:
            segs = cn_segments(tr)
            cn_max = int(min(8, max(3, np.ceil(max(v for *_, v in segs)))))
            with gz("cn.bedgraph") as fh:
                fh.write("track type=bedGraph\n")
                for c, s, e, v in segs:
                    fh.write(f"{c}\t{s}\t{e}\t{v}\n")
        st = strata.loc[grp] if grp in strata.index else None
        summary = {
            "key": key, "group": grp, "name": name, "lineage": lineage, "disease": disease, "subtype": subtype,
            "cn_source": prov, "has_cn": tr is not None, "cn_max": cn_max,
            "comparison": None if st is None else {"stratum": st.stratum, "n_lines": int(st.lines_in_stratum)},
            "n_called": int(len(share)), "n_experiments": int(n_exp),
            "experiments": [{"srx": s, "study": study.get(s, ""), "layout": layout.get(s, ""), "bw": BW.format(s)}
                            for s in exps],
            "vsall": tops["vsall"], "vsrel": tops["vsrel"],
        }
        json.dump(summary, open(fn("json"), "w"), indent=0)
        index[key] = {"group": grp, "name": name, "lineage": lineage,
                      "search": line_groups.get(grp, {}).get("search", name)}
        print(f"[61] {name}: {len(exps)} experiments, {len(share)} called SEs, vs all {tops['vsall']['n']}, "
              f"vs relatives {tops['vsrel']['n']}, CN {prov}{'' if tr is not None else ' (missing)'}", file=sys.stderr)
    json.dump(index, open(os.path.join(OUT, "index.json"), "w"), indent=0, sort_keys=True)


if __name__ == "__main__":
    main()
