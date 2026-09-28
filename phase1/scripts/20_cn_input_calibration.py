#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Phase 1 — 20: copy number from ChIP input, and a truth-free test of WHEN it can be trusted (FINDINGS §21).

The August pilot (11 lines) found per-line agreement with DepMap WGS anywhere from r = 0.15 to 0.96, and
nothing it tried narrowed that. This scales it to every line with both measured CN and a ChIP-Atlas input
('calib'), and computes the same inferred tracks for the lines held back without CN ('target'), so a gate
learnt on calib can be applied to target directly.

  manifest   lines x input SRX (<= 3 deepest per line), with ChIP-Atlas input QC     [laptop, no network]
  fetch      50 kb bins per input from bigWig zoom levels over HTTP (no download); one .npz per SRX,
             kept so every later method variant reruns without re-querying            [Roihu array]
  evaluate   per input: inferred vs DepMap WGS at 5 Mb (calib only), and truth-free features for all:
             lag-1 autocorrelation and first-difference noise of the 50 kb track, zero-coverage fraction,
             top-1% coverage share (enrichment leak), depth/duplication, agreement between a line's inputs
Inference = cnrose ChipInputInferredCN.build_from_bins with the August best settings (covered mean,
zero = deletion, ENCODE blacklist, no GC, no segmentation).

  python3 phase1/scripts/20_cn_input_calibration.py manifest --out <tsv>
  python3 phase1/scripts/20_cn_input_calibration.py fetch --manifest <tsv> --bins-dir <dir> --task I --ntasks N
  python3 phase1/scripts/20_cn_input_calibration.py evaluate --manifest <tsv> --bins-dir <dir> \
        --depmap-wgs OmicsCNGeneWGS.csv --gene-coords gene_coords.GRCh38.106.tsv --blacklist <bed.gz> --out <prefix>
"""
import argparse
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
CHIP_ATLAS_BW = "https://chip-atlas.dbcls.jp/data/hg38/eachData/bw/{srx}.bw"
HG38 = {"chr1": 248956422, "chr2": 242193529, "chr3": 198295559, "chr4": 190214555, "chr5": 181538259,
        "chr6": 170805979, "chr7": 159345973, "chr8": 145138636, "chr9": 138394717, "chr10": 133797422,
        "chr11": 135086622, "chr12": 133275309, "chr13": 114364328, "chr14": 107043718, "chr15": 101991189,
        "chr16": 90338345, "chr17": 83257441, "chr18": 80373285, "chr19": 58617616, "chr20": 64444167,
        "chr21": 46709983, "chr22": 50818468}
BIN, COARSE = 50_000, 5_000_000


# ----------------------------------------------------------------------------- manifest (laptop)
def cmd_manifest(a):
    import pandas as pd
    P1, P2, OUT = (os.path.join(ROOT, *p) for p in (("phase1", "data"), ("phase2", "data"), ("phase2", "analysis", "out")))
    cols = ["srx", "genome", "antigen_class", "antigen", "cell_class", "cell", "desc", "qc"]
    e = pd.read_csv(a.experiment_list, sep="\t", header=None, usecols=range(8), names=cols, dtype=str,
                    quoting=3, on_bad_lines="skip")
    inp = e[(e.genome == "hg38") & (e.antigen_class == "Input control")].copy()
    q = inp.qc.fillna("").str.split(",", expand=True)
    inp["reads"] = pd.to_numeric(q[0], errors="coerce")
    inp["pct_mapped"] = pd.to_numeric(q[1], errors="coerce")
    inp["pct_dup"] = pd.to_numeric(q[2], errors="coerce")
    by_label = {k: g.sort_values("reads", ascending=False) for k, g in inp.groupby("cell")}

    rows = []
    ps = pd.read_csv(os.path.join(P2, "pull_set.v2.tsv"), sep="\t")
    calib = ps[ps.cn_provider == "depmap_wgs"]
    for key, g in calib.groupby("key"):
        rows.append(("calib", key, g.cvcl.iloc[0], g.cell.iloc[0], set(g.cell)))
    reg = pd.read_csv(os.path.join(P1, "h3k27ac_registry.tsv"), sep="\t", low_memory=False)
    held = reg[reg.decision == "not pulled"]
    for cv, g in held.groupby("cvcl"):
        rows.append(("target", cv, cv, g.cell_line.iloc[0], set(g.chip_atlas_cell.dropna())))
    out = []
    for role, key, cv, name, labels in rows:
        cand = pd.concat([by_label[l] for l in labels if l in by_label]) if any(l in by_label for l in labels) \
            else pd.DataFrame()
        if cand.empty:
            continue
        cand = cand.drop_duplicates("srx").sort_values("reads", ascending=False).head(a.max_inputs)
        for r in cand.itertuples():
            out.append(dict(role=role, key=key, cvcl=cv, cell_line=name, input_srx=r.srx, input_label=r.cell,
                            reads=r.reads, pct_mapped=r.pct_mapped, pct_dup=r.pct_dup))
    m = pd.DataFrame(out)
    m.to_csv(a.out, sep="\t", index=False)
    print(f"[20] manifest: {m.groupby('role').key.nunique().to_dict()} lines, {m.input_srx.nunique()} inputs "
          f"-> {a.out}", file=sys.stderr)


# ----------------------------------------------------------------------------- fetch (Roihu array)
def grid():
    return {c: np.arange(int(L // BIN), dtype=np.int64) * BIN for c, L in HG38.items()}


def cmd_fetch(a):
    import pyBigWig
    srxs = sorted({l.split("\t")[4] for l in open(a.manifest).read().splitlines()[1:]})
    mine = srxs[a.task::a.ntasks]
    os.makedirs(a.bins_dir, exist_ok=True)
    g = grid()
    ok = skip = fail = 0
    for srx in mine:
        fn = os.path.join(a.bins_dir, f"{srx}.npz")
        if os.path.exists(fn):
            skip += 1
            continue
        bw, local = None, None
        url = CHIP_ATLAS_BW.format(srx=srx)
        if not pyBigWig.remote:                   # wheel built without libcurl (Roihu venv): download, bin, delete
            import subprocess
            local = os.path.join(a.bins_dir, ".tmp", f"{srx}.bw")
            os.makedirs(os.path.dirname(local), exist_ok=True)
            if subprocess.call(["curl", "-sfL", "--retry", "3", "-o", local, url]) != 0:
                local = None
            url = local
        for _ in range(3):
            try:
                bw = pyBigWig.open(url) if url else None
                break
            except Exception:                                      # noqa: BLE001
                bw = None
        if bw is None:
            if local and os.path.exists(local):
                os.remove(local)
            fail += 1
            print(f"[20] {srx} open failed", file=sys.stderr)
            continue
        try:
            have = bw.chroms()
            hdr = bw.header()
            arrs = {}
            for c, st in g.items():
                if c in have:
                    raw = bw.stats(c, 0, int(st[-1] + BIN), nBins=len(st), type="mean")
                    arrs[c] = np.array([x if x is not None else np.nan for x in raw], dtype=np.float32)
            tmp = fn + ".part.npz"
            np.savez_compressed(tmp, **arrs, _covered=np.array([hdr.get("nBasesCovered", 0)]),
                                _sum=np.array([hdr.get("sumData", 0)]))
            os.replace(tmp, fn)                                        # never leave a truncated file
            ok += 1
        except Exception as e:                                         # noqa: BLE001
            fail += 1
            print(f"[20] {srx} failed: {e}", file=sys.stderr)
        finally:
            bw.close()
            if local and os.path.exists(local):
                os.remove(local)
    print(f"[20] task {a.task}: ok={ok} skipped={skip} failed={fail} of {len(mine)}", file=sys.stderr)


# ----------------------------------------------------------------------------- evaluate
def coarse_log2(track, chroms=HG38):
    out = []
    for c, L in chroms.items():
        for s in range(0, L - COARSE, COARSE):
            v = track.region_cn(c, s, s + COARSE)
            out.append(np.log2(v) if v and v > 0 else np.nan)
    return np.array(out, dtype=float)


def features(cov, bmask):
    zero, lag1, dnoise, logs = [], [], [], []
    for c, v in cov.items():
        keep = ~bmask[c]
        x = np.nan_to_num(v[keep].astype(float), nan=0.0)
        zero.append(x == 0)
        lg = np.log2(x[x > 0])
        if lg.size > 10:
            logs.append(lg)
            d = lg - np.median(lg)
            lag1.append(np.corrcoef(d[:-1], d[1:])[0, 1])
            dnoise.append(np.median(np.abs(np.diff(lg))))
    allx = np.concatenate([np.nan_to_num(v[~bmask[c]].astype(float)) for c, v in cov.items()])
    srt = np.sort(allx)[::-1]
    top = srt[: max(1, len(srt) // 100)].sum() / max(srt.sum(), 1e-9)
    lg = np.concatenate(logs) if logs else np.array([np.nan])
    return {"zero_frac": float(np.concatenate(zero).mean()), "lag1_autocorr": float(np.nanmedian(lag1)),
            "diff_noise": float(np.nanmedian(dnoise)), "log2_iqr": float(np.subtract(*np.percentile(lg, [75, 25]))),
            "top1pct_share": float(top)}


def cmd_evaluate(a):
    import pandas as pd
    from scipy.stats import spearmanr
    sys.path.insert(0, a.cnrose)
    from cnrose.cn.inferred import ChipInputInferredCN, load_blacklist, blacklist_mask
    from cnrose.cn.depmap import DepMapGeneCN, load_gene_coords

    m = pd.read_csv(a.manifest, sep="\t")
    g = {c: (st, st + BIN) for c, st in grid().items()}
    bl = load_blacklist(a.blacklist)
    prov = ChipInputInferredCN({}, blacklist=bl)
    prov._grid, prov._bmask, prov._gc = g, blacklist_mask(g, bl), None           # no network for the grid
    coords = load_gene_coords(None, cache_path=a.gene_coords)
    truth_p = DepMapGeneCN(a.depmap_wgs, coords)
    calib_keys = sorted(m[m.role == "calib"].key.unique())
    truth_p.preload(calib_keys)
    truth = {k: coarse_log2(truth_p.track(k)) for k in calib_keys if truth_p.track(k) is not None}
    print(f"[20] truth tracks: {len(truth)} of {len(calib_keys)} calib lines", file=sys.stderr)

    rows, prof = [], {}
    for r in m.itertuples():
        fn = os.path.join(a.bins_dir, f"{r.input_srx}.npz")
        row = r._asdict()
        row.pop("Index", None)
        if not os.path.exists(fn):
            row["status"] = "no bins"; rows.append(row); continue
        z = np.load(fn)
        cov = {c: z[c] for c in z.files if c.startswith("chr")}
        tr = prov.build_from_bins(cov, r.input_srx)
        if tr is None:
            row["status"] = "track failed"; rows.append(row); continue
        row.update(features(cov, prov._bmask))
        x = coarse_log2(tr)
        prof[(r.key, r.input_srx)] = x
        if r.key in truth:
            y = truth[r.key]
            ok = np.isfinite(x) & np.isfinite(y)
            if ok.sum() > 50:
                row["r_pearson"] = float(np.corrcoef(x[ok], y[ok])[0, 1])
                row["r_spearman"] = float(spearmanr(x[ok], y[ok]).correlation)
                row["n_bins"] = int(ok.sum())
        row["status"] = "ok"
        rows.append(row)
    d = pd.DataFrame(rows)
    # truth-free consistency: agreement between a line's own inputs
    agree = {}
    for key, gg in d[d.status == "ok"].groupby("key"):
        ss = list(gg.input_srx)
        if len(ss) >= 2:
            rs = []
            for i in range(len(ss)):
                for j in range(i + 1, len(ss)):
                    x, y = prof[(key, ss[i])], prof[(key, ss[j])]
                    ok = np.isfinite(x) & np.isfinite(y)
                    if ok.sum() > 50:
                        rs.append(np.corrcoef(x[ok], y[ok])[0, 1])
            agree[key] = float(np.median(rs)) if rs else np.nan
    d["inter_input_r"] = d.key.map(agree)
    d.to_csv(a.out + ".inputs.tsv", sep="\t", index=False)
    c = d[(d.role == "calib") & d.r_pearson.notna()]
    print(f"[20] calib inputs with truth: {len(c)} on {c.key.nunique()} lines; median r {c.r_pearson.median():.3f}", file=sys.stderr)
    feats = ["lag1_autocorr", "diff_noise", "zero_frac", "top1pct_share", "log2_iqr", "reads", "pct_dup", "inter_input_r"]
    for f in feats:
        if f in c and c[f].notna().sum() > 20:
            print(f"   Spearman(r, {f:14s}) = {spearmanr(c[f], c.r_pearson, nan_policy='omit').correlation:+.2f}", file=sys.stderr)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("manifest")
    s.add_argument("--experiment-list", default=os.path.join(ROOT, ".cache", "experimentList.2026-09.tab"))
    s.add_argument("--max-inputs", type=int, default=3)
    s.add_argument("--out", required=True)
    s = sub.add_parser("fetch")
    s.add_argument("--manifest", required=True); s.add_argument("--bins-dir", required=True)
    s.add_argument("--task", type=int, default=0); s.add_argument("--ntasks", type=int, default=1)
    s = sub.add_parser("evaluate")
    s.add_argument("--manifest", required=True); s.add_argument("--bins-dir", required=True)
    s.add_argument("--depmap-wgs", required=True); s.add_argument("--gene-coords", required=True)
    s.add_argument("--blacklist", required=True)
    s.add_argument("--cnrose", default=os.path.join(ROOT, "cnrose"))
    s.add_argument("--out", required=True)
    a = ap.parse_args()
    {"manifest": cmd_manifest, "fetch": cmd_fetch, "evaluate": cmd_evaluate}[a.cmd](a)


if __name__ == "__main__":
    main()
