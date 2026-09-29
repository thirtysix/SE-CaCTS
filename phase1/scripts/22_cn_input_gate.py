#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Phase 1 — 22: does input-inferred CN help, per input, and can that be predicted without truth? (FINDINGS §24)

Correlation (script 20) is the wrong yardstick for a correction: a near-diploid genome has little variance, so
r is low even when the inferred track is flat and harmless. What matters is whether APPLYING the inferred CN
beats applying nothing. For every calib input, on 5 Mb log2 profiles (both median-centred):

  gain       1 - MSE(b * inferred, truth) / MSE(0, truth); > 0 means the correction helps. b is ONE global
             compression slope fitted on all calib bins (inferred amplitudes are compressed, FINDINGS
             2026-08-18 §6), i.e. the correction we would actually apply
  amp_sens   truth-amplified bins (log2 >= 1) called amplified (b * inferred >= 0.58)
  false_amp  truth-neutral bins (|log2| < 0.3) called amplified: the harmful error for SE calling

Then: which truth-free features predict gain > 0, with grouped cross-validation so no line's own inputs train
its prediction; a gate at the probability giving >= 90% precision; applied to the target lines.

  python3 phase1/scripts/22_cn_input_gate.py --bins-dir .cache/cncal_bins \
      --depmap-wgs $DATAROOT/DepMap/2026q1/OmicsCNGeneWGS.csv
"""
import argparse
import importlib.util
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
OUT = os.path.join(ROOT, "phase2", "analysis", "out")
FEATS = ["lag1_autocorr", "diff_noise", "zero_frac", "top1pct_share", "log2_iqr", "log_reads", "pct_dup",
         "inter_input_r"]


def load20():
    spec = importlib.util.spec_from_file_location("cal20", os.path.join(HERE, "20_cn_input_calibration.py"))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bins-dir", required=True)
    ap.add_argument("--depmap-wgs", required=True)
    ap.add_argument("--inputs", default=os.path.join(OUT, "cn_input_calibration.inputs.tsv"))
    ap.add_argument("--gene-coords", default=os.path.join(ROOT, ".cache", "gene_coords.GRCh38.106.tsv"))
    ap.add_argument("--blacklist", default=os.path.join(ROOT, ".cache", "hg38-blacklist.v2.bed.gz"))
    ap.add_argument("--precision", type=float, default=0.90)
    a = ap.parse_args()
    sys.path.insert(0, os.path.join(ROOT, "cnrose"))
    from cnrose.cn.inferred import ChipInputInferredCN, load_blacklist, blacklist_mask
    from cnrose.cn.depmap import DepMapGeneCN, load_gene_coords
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import roc_auc_score
    from sklearn.model_selection import GroupKFold
    c20 = load20()

    d = pd.read_csv(a.inputs, sep="\t")
    d = d[d.status == "ok"].reset_index(drop=True)
    g = {c: (st, st + c20.BIN) for c, st in c20.grid().items()}
    bl = load_blacklist(a.blacklist)
    prov = ChipInputInferredCN({}, blacklist=bl)
    prov._grid, prov._bmask, prov._gc = g, blacklist_mask(g, bl), None
    prof = {}
    for s in d.input_srx.unique():
        z = np.load(os.path.join(a.bins_dir, f"{s}.npz"))
        tr = prov.build_from_bins({c: z[c] for c in z.files if c.startswith("chr")}, s)
        if tr is not None:
            x = c20.coarse_log2(tr)
            prof[s] = x - np.nanmedian(x)
    keys = sorted(d[d.role == "calib"].key.unique())
    wgs = DepMapGeneCN(a.depmap_wgs, load_gene_coords(None, cache_path=a.gene_coords))
    wgs.preload(keys)
    truth = {}
    for k in keys:
        t = wgs.track(k)
        if t is not None:
            y = c20.coarse_log2(t)
            truth[k] = y - np.nanmedian(y)
    print(f"[22] profiles: {len(prof)} inputs; truth: {len(truth)} calib lines", file=sys.stderr)

    cal = d[(d.role == "calib") & d.key.isin(truth) & d.input_srx.isin(prof)].copy()
    X = np.concatenate([prof[s] for s in cal.input_srx])
    Y = np.concatenate([truth[k] for k in cal.key])
    ok = np.isfinite(X) & np.isfinite(Y)
    b = float((X[ok] * Y[ok]).sum() / (X[ok] ** 2).sum())
    print(f"[22] global compression slope b = {b:.2f} (truth ~ b x inferred)", file=sys.stderr)

    def metrics(s, k):
        x, y = b * prof[s], truth[k]
        ok = np.isfinite(x) & np.isfinite(y)
        x, y = x[ok], y[ok]
        m0 = np.mean(y ** 2)
        amp, neu = y >= 1, np.abs(y) < 0.3
        return {"gain": 1 - np.mean((x - y) ** 2) / m0 if m0 > 0 else np.nan, "truth_sd": float(np.std(y)),
                "amp_bins": int(amp.sum()), "amp_sens": float((x[amp] >= 0.58).mean()) if amp.any() else np.nan,
                "false_amp": float((x[neu] >= 0.58).mean()) if neu.any() else np.nan}
    cal = pd.concat([cal.reset_index(drop=True),
                     pd.DataFrame([metrics(s, k) for s, k in zip(cal.input_srx, cal.key)])], axis=1)
    d["log_reads"] = np.log10(d.reads.clip(lower=1e5))
    cal["log_reads"] = np.log10(cal.reads.clip(lower=1e5))
    print(f"[22] calib inputs {len(cal)} on {cal.key.nunique()} lines: gain > 0 for {(cal.gain > 0).mean():.0%} "
          f"(median gain {cal.gain.median():+.2f}); median false_amp {cal.false_amp.median():.3f}; "
          f"amp_sens median {cal.amp_sens.median():.2f} over {int(cal.amp_bins.gt(0).sum())} inputs with amplified bins",
          file=sys.stderr)
    from scipy.stats import spearmanr
    for f in FEATS:
        print(f"   Spearman(gain, {f:14s}) {spearmanr(cal[f], cal.gain, nan_policy='omit').correlation:+.2f}   "
              f"(r vs truth_sd {spearmanr(cal.r_pearson, cal.truth_sd).correlation:+.2f})" if f == FEATS[0] else
              f"   Spearman(gain, {f:14s}) {spearmanr(cal[f], cal.gain, nan_policy='omit').correlation:+.2f}", file=sys.stderr)

    # grouped CV logistic model: P(gain > 0) from truth-free features
    feats = [f for f in FEATS if f != "inter_input_r"]                 # most lines have one input: keep it out
    Xf = cal[feats].fillna(cal[feats].median()).values
    yb = (cal.gain > 0).values.astype(int)
    pred = np.zeros(len(cal))
    for tr, te in GroupKFold(n_splits=10).split(Xf, yb, groups=cal.key):
        mu, sd = Xf[tr].mean(0), Xf[tr].std(0) + 1e-9
        m = LogisticRegression(max_iter=1000).fit((Xf[tr] - mu) / sd, yb[tr])
        pred[te] = m.predict_proba((Xf[te] - mu) / sd)[:, 1]
    cal["p_help"] = pred
    auc = roc_auc_score(yb, pred)
    ths = np.linspace(0.3, 0.95, 66)
    prec = [(yb[pred >= t].mean() if (pred >= t).any() else np.nan) for t in ths]
    ok_t = [t for t, p in zip(ths, prec) if p == p and p >= a.precision]
    t_star = min(ok_t) if ok_t else None
    print(f"[22] CV AUROC {auc:.2f}; gate p >= {t_star} for precision >= {a.precision:.0%}", file=sys.stderr)
    # final model on all calib, applied to everything
    mu, sd = Xf.mean(0), Xf.std(0) + 1e-9
    m = LogisticRegression(max_iter=1000).fit((Xf - mu) / sd, yb)
    Xa = d[feats].fillna(cal[feats].median()).values
    d["p_help"] = m.predict_proba((Xa - mu) / sd)[:, 1]
    d = d.merge(cal[["input_srx", "gain", "truth_sd", "amp_sens", "false_amp"]], on="input_srx", how="left")
    d.loc[d.role == "calib", "p_help"] = d.loc[d.role == "calib", "input_srx"].map(cal.set_index("input_srx").p_help)
    best = d.sort_values("p_help", ascending=False).drop_duplicates("key")
    if t_star is not None:
        best["admit"] = best.p_help >= t_star
        cb = best[best.role == "calib"]
        tb = best[best.role == "target"]
        print(f"[22] per line (best input): calib admitted {int(cb.admit.sum())}/{len(cb)}, of which gain > 0 "
              f"{(cb[cb.admit].gain > 0).mean():.0%} (median gain {cb[cb.admit].gain.median():+.2f}, false_amp "
              f"{cb[cb.admit].false_amp.median():.3f}); rejected lines gain > 0 {(cb[~cb.admit].gain > 0).mean():.0%}",
              file=sys.stderr)
        reg = pd.read_csv(os.path.join(ROOT, "phase1", "data", "h3k27ac_registry.tsv"), sep="\t", low_memory=False)
        held = reg[reg.decision == "not pulled"].groupby("cvcl").srx.nunique()
        tb = tb.assign(n_exp=tb.key.map(held).fillna(0).astype(int))
        print(f"[22] targets admitted: {int(tb.admit.sum())}/{len(tb)} lines, "
              f"{int(tb[tb.admit].n_exp.sum())}/{int(tb.n_exp.sum())} held-back experiments", file=sys.stderr)
    d.to_csv(os.path.join(OUT, "cn_input_gate.inputs.tsv"), sep="\t", index=False)
    best.to_csv(os.path.join(OUT, "cn_input_gate.lines.tsv"), sep="\t", index=False)


if __name__ == "__main__":
    main()
