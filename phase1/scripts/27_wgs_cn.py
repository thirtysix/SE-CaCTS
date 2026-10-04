#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Phase 1 — 27: measured copy number from subsampled open WGS for lines scored on input-inferred CN (v3.2;
phase1/data/cn_gap_sources_2026-10.md, "cheapest route" 3).

One run per line, the first SPOTS spots only (the full 30-190 GB runs are never fetched), one read per spot,
Bowtie2 to GRCh38, reads per 50 kb bin (phase2/roihu/wgs_cn.slurm + phase2/wgs_bin.py). Read counts carry the
mappability of each bin, which ChIP-input covered-mean bins do not, so every line is divided by a NORMAL reference
run through the same pipeline (1000 Genomes high-coverage, 3 samples): bins the reference cannot map (< 30% of its
median) become NaN and stay missing (BinnedInputCN zero_is_deletion=False); a true 0 is still a deletion.
`noise` = median |first difference| of log2 ratio between neighbouring 50 kb bins: the truth-free quality measure
a target is gated on, against its range on the calibration lines.
Calibration lines (role calib) have DepMap WGS, the same data's CCLE reads, so `evaluate` measures the estimator
itself: agreement at the atlas SE loci and on 5 Mb bins, against the input-inferred track where one exists.

  python3 phase1/scripts/27_wgs_cn.py manifest --out phase2/data/wgs_cn_manifest.tsv
  # Roihu: sbatch --array=0-<n-1> --export=ALL,PROJ=<project>,MANIFEST=<tsv> wgs_cn.slurm; fetch wgs_cn/bins/
  python3 phase1/scripts/27_wgs_cn.py ratio --manifest phase2/data/wgs_cn_manifest.tsv --bins <bins> --out <ratio dir>
  # sorted-alignment runs (FINDINGS §59): Roihu csra_cov.slurm, then
  python3 phase1/scripts/27_wgs_cn.py csra --dumps <dir of .ref.tsv.gz> --out <ratio dir>
  python3 phase1/scripts/27_wgs_cn.py evaluate --manifest phase2/data/wgs_cn_manifest.tsv --ratio <ratio dir>
"""
import argparse
import io
import os
import sys
import urllib.parse
import urllib.request

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "cnrose"))
from secacts_env import DATAROOT, cache_path                              # noqa: E402

# (cell as in the pull set, run, mode) from cn_gap_sources_2026-10.tsv, open short-read WGS (low-pass for GSC11/23),
# 10x linked reads (read 2), ONT (cut to 100 bp). RERF-LC-OK is out of scope (identity drop, v3.1).
TARGETS = [("MIC", "SRR11235318", None), ("POE", "SRR11235333", None), ("GSC11", "SRR4009298", None),
           ("GSC23", "SRR4009300", None), ("SNU-484", "DRR878403", None), ("HepaRG", "ERR2004597", None),
           ("II-18", "DRR016885", None), ("BCBL1", "SRR18286498", None),
           ("TL-Om1", "DRR248583", None), ("LA-N-5", "ERR10075187", None),
           ("SMS-KCNR", "ERR10075184", None),
           ("C4-2B", "SRR15368035", "10x"), ("LAPC-4", "SRR15368039", "10x"), ("C666-1", "SRR37324005", "ont")]
# Runs submitted as COORDINATE-SORTED alignments: their first spots all come from the start of chr1 (FINDINGS §57),
# so "the first N spots" is not a sample. The two low-pass ones are fetched whole (SPOTS_ALL); the others need a
# different sampler and stay on input inference for now: OCI-LY1 SRR1236466, MUTZ-3 SRR30002267, NBL-S SRR34067880,
# ASPS1 SRR33168805. CCLE 2019 WGS (PRJNA523380) is sorted too, so it cannot calibrate this estimator.
SPOTS_ALL = {"GSC11": 40_000_000, "GSC23": 40_000_000}
# Calibration: WGS in sequencing order of lines with DepMap WGS, from other labs (so DepMap is independent truth):
# Princess Maxima neuroblastoma panel (PRJEB54725; MYCN-amplified IMR-32, CHP-134, NGP, NB-1643; non-amplified
# SK-N-AS, SH-SY5Y) and the Institut Curie Ewing panel (PRJNA610192; A673, SK-N-MC, TC-71, RD-ES)
CALIB = [("IMR32", "ACH-000310", "ERR10075183"), ("CHP134", "ACH-001338", "ERR10075181"),
         ("NGP", "ACH-001366", "ERR10075191"), ("NB1643", "ACH-001303", "ERR10075189"),
         ("SKNAS", "ACH-000260", "ERR10075198"), ("SHSY5Y", "ACH-001188", "ERR10075193"),
         ("A673", "ACH-000052", "SRR11235335"), ("SKNMC", "ACH-000039", "SRR11235329"),
         ("TC71", "ACH-000424", "SRR11235327"), ("RDES", "ACH-000041", "SRR25068618")]
REF = [("ref_HG03814", "ERR3534515"), ("ref_HG00097", "ERR3535646"), ("ref_HG00106", "ERR3535780")]
AUTO = [f"chr{i}" for i in range(1, 23)]


def ena_layout(runs):
    q = " OR ".join(f'run_accession="{r}"' for r in runs)
    u = ("https://www.ebi.ac.uk/ena/portal/api/search?result=read_run&fields=run_accession,library_layout,"
         "instrument_model,read_count,base_count&format=tsv&limit=0&query=" + urllib.parse.quote(q))
    with urllib.request.urlopen(u, timeout=120) as r:
        return pd.read_csv(io.StringIO(r.read().decode()), sep="\t").set_index("run_accession")


def cmd_manifest(a):
    ps = pd.read_csv(os.path.join(ROOT, "phase2/data/pull_set.v31.tsv"), sep="\t")
    key_of = ps.drop_duplicates("cell").set_index("cell")["key"].to_dict()
    lay = ena_layout([r for _, r, _ in TARGETS] + [r for *_, r in CALIB] + [r for _, r in REF])
    rows = []
    for cell, run, mode in TARGETS:
        # not in ENA (DRR878403): fastq-dump reads NCBI SRA; the gap table lists it as a paired WGS library
        m = mode or ("pe" if lay["library_layout"].get(run, "PAIRED") == "PAIRED" else "se")
        rows.append((cell, key_of[cell], "target", run, m,
                     SPOTS_ALL.get(cell, a.ont_spots if m == "ont" else a.spots)))
    for cell, key, run in CALIB:
        rows.append((cell, key, "calib", run, "pe", a.spots))
    for name, run in REF:
        rows.append((name, name, "ref", run, "pe", a.ref_spots))
    out = pd.DataFrame(rows, columns=["name", "key", "role", "run", "mode", "spots"])
    out.to_csv(a.out, sep="\t", index=False)
    print(out.merge(lay, left_on="run", right_index=True, how="left").to_string(index=False))


def load_bins(fn):
    z = np.load(fn)
    return {c: z[c].astype(np.float64) for c in z.files if c.startswith("chr")}, int(z["_reads"][0])


def cmd_ratio(a):
    man = pd.read_csv(a.manifest, sep="\t")
    os.makedirs(a.out, exist_ok=True)
    for f in os.listdir(a.out):           # a line skipped now must not leave an older ratio file behind
        if f.endswith(".npz"):
            os.remove(os.path.join(a.out, f))
    refs = []
    for n in man.loc[man.role == "ref", "name"]:
        b, nr = load_bins(os.path.join(a.bins, f"{n}.npz"))
        tot = sum(b[c].sum() for c in AUTO)
        refs.append({c: v / tot for c, v in b.items()})
        print(f"[ratio] reference {n}: {nr:,} reads", file=sys.stderr)
    ref = {c: np.mean([r[c] for r in refs], axis=0) for c in refs[0]}
    med = np.median(np.concatenate([ref[c][ref[c] > 0] for c in AUTO]))
    rows = []
    for r in man[man.role != "ref"].itertuples():
        fn = os.path.join(a.bins, f"{r.name}.npz")
        if not os.path.exists(fn):
            print(f"[ratio] missing {r.name}", file=sys.stderr)
            continue
        b, nr = load_bins(fn)
        tot = sum(b[c].sum() for c in AUTO)
        empty = np.mean(np.concatenate([b[c][ref[c] >= a.min_ref * med] == 0 for c in AUTO]))
        if empty > 0.05:                  # a coordinate-sorted run: the reads cover only the start of the genome
            print(f"[ratio] {r.name}: {empty:.0%} of mappable bins empty, not a genome-wide sample; skipped",
                  file=sys.stderr)
            continue
        out = {}
        for c in ref:
            ok = ref[c] >= a.min_ref * med
            v = np.full(len(ref[c]), np.nan)
            v[ok] = (b[c][ok] / tot) / ref[c][ok]
            out[c] = v.astype(np.float32)
        np.savez_compressed(os.path.join(a.out, f"{r.name}.npz"), **out, _reads=np.array([nr]))
        rows.append((r.name, r.key, r.role, nr, "wgs_fastq"))
        print(f"[ratio] {r.name} ({r.role}): {nr:,} reads", file=sys.stderr)
    pd.DataFrame(rows, columns=["name", "key", "role", "reads", "source"]).to_csv(os.path.join(a.out, "index.tsv"),
                                                                                   sep="\t", index=False)


ASSEMBLY_CHR1 = {247_249_719: "hg18", 249_250_621: "hg19", 248_956_422: "hg38"}
CHAINS = {"hg18": os.path.join(DATAROOT, "0.human_genome/hg18ToHg38.over.chain"),
          "hg19": os.path.expanduser("~/.pyliftover/hg19ToHg38.over.chain.gz")}
HG38_LEN = {"chr1": 248956422, "chr2": 242193529, "chr3": 198295559, "chr4": 190214555, "chr5": 181538259,
            "chr6": 170805979, "chr7": 159345973, "chr8": 145138636, "chr9": 138394717, "chr10": 133797422,
            "chr11": 135086622, "chr12": 133275309, "chr13": 114364328, "chr14": 107043718, "chr15": 101991189,
            "chr16": 90338345, "chr17": 83257441, "chr18": 80373285, "chr19": 58617616, "chr20": 64444167,
            "chr21": 46709983, "chr22": 50818468, "chrX": 156040895}


def csra_bins(fn, bin_size=50_000):
    """5 kb chunk first-alignment ids of a sorted cSRA run (csra_cov.slurm) -> reads per hg38 50 kb bin."""
    d = pd.read_csv(fn, sep="\t", header=None, names=["name", "start", "len", "first"], dtype={"name": str})
    d["first"] = pd.to_numeric(d["first"], errors="coerce")
    f = d["first"].to_numpy()
    nonempty = np.where(np.isfinite(f))[0]
    cnt = np.zeros(len(d))
    cnt[nonempty[:-1]] = np.diff(f[nonempty])                  # ids are contiguous in sorted order
    cnt[nonempty[-1]] = np.nan                                  # the last chunk's count is unknown
    d["count"] = cnt
    d["chrom"] = "chr" + d["name"].str.replace("chr", "", regex=False)
    c1 = d[d.chrom == "chr1"]
    L1 = int((c1.start + c1.len - 1).max())
    near = min(ASSEMBLY_CHR1, key=lambda x: abs(x - L1))          # a dump can end a few chunks short (OCI-LY1 hg18)
    asm = ASSEMBLY_CHR1[near] if abs(near - L1) < 200_000 else None
    if asm is None:
        raise SystemExit(f"[csra] {fn}: chr1 length {L1} matches no known assembly")
    d = d[d.chrom.isin(HG38_LEN) & d["count"].notna()]
    out = {c: np.zeros(L // bin_size) for c, L in HG38_LEN.items()}
    lo = None
    if asm != "hg38":
        from pyliftover import LiftOver
        lo = LiftOver(CHAINS[asm])
    lost = 0
    for ch, s, ln, n in zip(d.chrom, d.start, d["len"], d["count"]):
        mid = int(s) - 1 + int(ln) // 2
        if lo is not None:
            r = lo.convert_coordinate(ch, mid)
            if not r or r[0][0] != ch:
                lost += n
                continue
            mid = int(r[0][1])
        i = mid // bin_size
        if i < len(out[ch]):
            out[ch][i] += n
    return out, asm, float(lost / max(d["count"].sum(), 1))


def cmd_csra(a):
    """Sorted-alignment runs: bins from the cSRA REFERENCE table, divided by a pooled median of the CCLE runs
    (leave-one-out for the calibration lines); written next to the WGS ratios and appended to index.tsv."""
    man = pd.read_csv(a.manifest, sep="\t")
    bins, info = {}, []
    for r in man.itertuples():
        fn = os.path.join(a.dumps, f"{r.name}.ref.tsv.gz")
        if not os.path.exists(fn):
            print(f"[csra] missing {r.name}", file=sys.stderr)
            continue
        b, asm, lost = csra_bins(fn)
        tot = sum(b[c].sum() for c in AUTO)
        bins[r.name] = {c: v / tot for c, v in b.items()}
        info.append((r.name, r.key, r.role, int(tot), asm, round(lost, 4)))
        print(f"[csra] {r.name}: {int(tot):,} reads, {asm}, {lost:.2%} lost in liftover", file=sys.stderr)
    pool = [n for n, _, role, *_ in info if role in ("calib", "pool")]
    idx_f = os.path.join(a.out, "index.tsv")
    idx = pd.read_csv(idx_f, sep="\t") if os.path.exists(idx_f) else pd.DataFrame(columns=["name", "key", "role", "reads"])
    idx = idx[~idx.name.isin([i[0] for i in info])]
    rows = []
    for name, key, role, tot, asm, lost in info:
        use = [p for p in pool if p != name]
        ref = {c: np.median([bins[p][c] for p in use], axis=0) for c in HG38_LEN}
        med = np.median(np.concatenate([ref[c][ref[c] > 0] for c in AUTO]))
        out = {}
        for c in HG38_LEN:
            ok = ref[c] >= a.min_ref * med
            v = np.full(len(ref[c]), np.nan)
            v[ok] = bins[name][c][ok] / ref[c][ok]
            out[c] = v.astype(np.float32)
        np.savez_compressed(os.path.join(a.out, f"{name}.npz"), **out, _reads=np.array([tot]))
        rows.append({"name": name, "key": key, "role": "calib" if role == "pool" else role, "reads": tot,
                     "source": "wgs_csra"})
    pd.concat([idx, pd.DataFrame(rows)], ignore_index=True).to_csv(idx_f, sep="\t", index=False)
    pd.DataFrame(info, columns=["name", "key", "role", "reads", "assembly", "lost_in_liftover"]).to_csv(
        os.path.join(a.out, "csra_info.tsv"), sep="\t", index=False)


def cmd_evaluate(a):
    from cnrose.cn.inferred import BinnedInputCN, load_blacklist
    from cnrose.cn.depmap import load_gene_coords, DepMapGeneCN
    man = pd.read_csv(os.path.join(a.ratio, "index.tsv"), sep="\t")   # what `ratio` and `csra` wrote
    man = man[man.role != "ref"]
    bl = load_blacklist(a.blacklist)
    wgs = BinnedInputCN(a.ratio, dict(zip(man.key, man.name)), blacklist=bl, slope=1.0, zero_is_deletion=False,
                        segment=True)
    im = pd.read_csv(os.path.join(ROOT, "phase2/analysis/out/cn_admissions.tsv"), sep="\t").dropna(subset=["cn_input"])
    inf = BinnedInputCN(a.input_bins, dict(zip(im["key"], im["cn_input"])), blacklist=bl)
    cat = pd.read_csv(a.catalog, sep="\t", header=None, usecols=[0, 1, 2], names=["c", "s", "e"])
    cat = cat[cat.c.isin(AUTO)]
    tiles = [(c, s, s + 5_000_000) for c in AUTO for s in range(0, 250_000_000, 5_000_000)]
    truth = None
    if (man.role == "calib").any():
        truth = DepMapGeneCN(a.depmap_wgs, load_gene_coords(None, cache_path=a.gene_cache))
        truth.preload(list(man.loc[man.role == "calib", "key"]))

    def at(tr, regs):
        return np.array([tr.region_cn(c, s, e) if tr is not None else np.nan for c, s, e in regs])

    def cmp(x, y):
        ok = np.isfinite(x) & np.isfinite(y)
        lx, ly = np.log2(np.clip(x[ok], 0.05, None)), np.log2(np.clip(y[ok], 0.05, None))
        return np.corrcoef(lx, ly)[0, 1], np.median(np.abs(lx - ly))

    rows = []
    for r in man.itertuples():
        w = wgs.track(r.key)
        if w is None:
            print(f"[eval] {r.name}: no WGS track", file=sys.stderr)
            continue
        i = inf.track(r.key)
        regs = list(zip(cat.c, cat.s, cat.e))
        xw = at(w, regs)
        z = np.load(os.path.join(a.ratio, f"{r.name}.npz"))
        lg = [np.log2(np.clip(z[c][np.isfinite(z[c]) & (z[c] > 0)], 1e-3, None)) for c in AUTO]
        noise = float(np.median(np.concatenate([np.abs(np.diff(v)) for v in lg if len(v) > 1])))
        row = {"name": r.name, "key": r.key, "role": r.role, "noise": round(noise, 3), "se_amp_wgs": int(np.sum(xw >= 2)),
               "se_altered_wgs": round(float(np.mean(np.abs(np.log2(np.clip(xw, 0.05, None))) > 0.3)), 3)}
        if i is not None:
            row["r_se_wgs_vs_inferred"], _ = cmp(xw, at(i, regs))
        if r.role == "calib":
            t = truth.track(r.key)
            xt = at(t, regs)
            row["r_se_wgs_vs_truth"], row["mad_se_wgs_vs_truth"] = cmp(xw, xt)
            row["r_5mb_wgs_vs_truth"], _ = cmp(at(w, tiles), at(t, tiles))
            amp = xt >= 2
            row["amp_sens"] = round(float(np.mean(xw[amp] >= 2)), 3) if amp.any() else np.nan
            row["false_amp"] = round(float(np.mean(xw[np.abs(np.log2(np.clip(xt, 0.05, None))) < 0.3] >= 2)), 4)
            if i is not None:
                row["r_se_inferred_vs_truth"], _ = cmp(at(i, regs), xt)
        rows.append(row)
        print(f"[eval] {row}", file=sys.stderr)
    out = pd.DataFrame(rows)
    out.to_csv(a.out, sep="\t", index=False, float_format="%.3f")
    # gate, per source (FASTQ WGS, sorted-run WGS, each array platform): the source passes when its calibration lines
    # (DepMap WGS truth) reach median r >= MIN_R at SE loci (input inference: 0.69) and median false amplification
    # <= MAX_FA; a target passes when its source does and its bin noise is within that source's calibration (+10%)
    src = man.set_index("name")["source"] if "source" in man else pd.Series(dtype=str)
    out["source"] = out.name.map(src).fillna("wgs_fastq")
    cal = out[out.role == "calib"].groupby("source").agg(r=("r_se_wgs_vs_truth", "median"),
                                                         fa=("false_amp", "median"), noise=("noise", "max"))
    cal["ok"] = (cal.r >= a.min_r) & (cal.fa <= a.max_fa)
    print(cal.round(3).to_string(), file=sys.stderr)
    g = out[out.role == "target"][["key", "name", "noise", "source"]].copy()
    g["limit"] = g.source.map((cal.noise * 1.1).round(3))
    g["source_ok"] = g.source.map(cal.ok).fillna(False)
    g["accept"] = np.where(g.source_ok & (g.noise <= g.limit), "yes", "no")
    g["provider"] = np.where(g.source.str.startswith("array"), "array_cgh", "wgs_reads")
    g.to_csv(os.path.join(ROOT, "phase2/data/wgs_cn_gate.tsv"), sep="\t", index=False)
    # the track each accepted line is scored on, one per key (WGS before arrays; then the first written): the readers
    # use this, not index.tsv, where a line with two tracks (OCI-LY1 WGS + array) would resolve to the LAST row
    use = g[g.accept == "yes"].assign(o=lambda d: (d.provider != "wgs_reads").astype(int)).sort_values("o", kind="stable")
    use.drop_duplicates("key")[["key", "name", "provider"]].to_csv(os.path.join(a.ratio, "use.tsv"), sep="\t", index=False)
    print(f"[gate] accepted {', '.join(g.loc[g.accept == 'yes', 'name'])}; "
          f"rejected {', '.join(g.loc[g.accept == 'no', 'name'])}", file=sys.stderr)
    with pd.option_context("display.width", 250):
        print(out.round(3).to_string(index=False))


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    m = sub.add_parser("manifest")
    m.add_argument("--out", required=True)
    m.add_argument("--spots", type=int, default=5_000_000)
    m.add_argument("--ref-spots", type=int, default=10_000_000)
    m.add_argument("--ont-spots", type=int, default=150_000)
    r = sub.add_parser("ratio")
    r.add_argument("--manifest", required=True)
    r.add_argument("--bins", required=True)
    r.add_argument("--out", required=True)
    r.add_argument("--min-ref", type=float, default=0.3, help="reference bins below this share of its median -> NaN")
    c = sub.add_parser("csra")
    c.add_argument("--manifest", default=os.path.join(ROOT, "phase2/data/csra_cn_manifest.tsv"))
    c.add_argument("--dumps", required=True, help="dir of <name>.ref.tsv.gz from csra_cov.slurm")
    c.add_argument("--out", required=True, help="the ratio dir of `ratio` (run that first: it clears the dir)")
    c.add_argument("--min-ref", type=float, default=0.3)
    e = sub.add_parser("evaluate")
    e.add_argument("--manifest", required=True)
    e.add_argument("--ratio", required=True)
    e.add_argument("--input-bins", default=os.path.join(ROOT, ".cache/cncal_bins"))
    e.add_argument("--blacklist", default=os.path.join(ROOT, ".cache/hg38-blacklist.v2.bed.gz"))
    e.add_argument("--catalog", default=os.path.join(ROOT, "phase2/results_v31f/atlas.s3.union_catalog.bed.gz"))
    e.add_argument("--depmap-wgs", default=os.path.join(DATAROOT, "DepMap/2026q1/OmicsCNGeneWGS.csv"))
    e.add_argument("--gene-cache", default=cache_path("gene_coords.GRCh38.106.tsv"))
    e.add_argument("--out", default=os.path.join(ROOT, "phase2/analysis/out/wgs_cn_eval.tsv"))
    e.add_argument("--min-r", type=float, default=0.69, help="calibration median r a source needs (input: 0.69)")
    e.add_argument("--max-fa", type=float, default=0.01, help="calibration median false-amplification rate")
    a = ap.parse_args()
    {"manifest": cmd_manifest, "ratio": cmd_ratio, "csra": cmd_csra, "evaluate": cmd_evaluate}[a.cmd](a)


if __name__ == "__main__":
    main()
