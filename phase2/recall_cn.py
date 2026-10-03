#!/usr/bin/env python3
"""Calling-time CN correction from the retained per-sample outputs, without a re-pull (GOTCHAS 59, v3.1 item 13).

The pull ran cnrose with no `--cn`, so only the agnostic `<SRX>.se.bed` exists. The calling signal is kept in
`<SRX>.enhancers.tsv` (stitched regions + SIGNAL), and calling-time correction acts on exactly that vector
(cnrose.pipeline.call_sample), so re-calling from it reproduces `cnrose call --cn ... --correct-at calling`:
region CN (length-weighted mean), log2-offset, beta 1, floor 1.0 (amplify-only: a deleted enhancer must never
become an SE; cnrose DESIGN §6.3), then the ROSE tangent cutoff. Writes `<SRX>.cn.enhancers.tsv` (with
REGION_CN) and `<SRX>.cn.se.bed` next to the originals and one summary row per sample.

`--mode fused` (user proposal 2026-10-03, FINDINGS §54) writes `<SRX>.fu.se.bed` instead: the cutoff is computed on the
CN-corrected signal (amplicons cannot inflate it) and applied to the UNCORRECTED signal, so amplified SEs stay and
SEs the inflated cutoff hid come back. Each call is labelled in column 7: `core` (an SE corrected and uncorrected),
`amplified` (passes only because its region is amplified: corrected signal below the cutoff), `unmasked` (an SE only
once the cutoff is not inflated: not an agnostic SE). Corrected SEs are a subset of the fused set (amplify-only
correction never raises a signal); an agnostic SE can fall out only where the corrected cutoff is higher.

Each line reads the CN source it is scored on (manifest `cn_provider`, as score_pilot.py): DepMap WGS, CMP WES
(by `cn_cvcl` or `cvcl`), DepMap MC_WES, CCLE SNP6, or input-inferred. A line with no track is skipped.

    python3 phase2/recall_cn.py --manifest recall_manifest.tsv --out-dir $W/out --task 0 --ntasks 16 \\
        --inferred-bins $W/../cncal/bins --blacklist $W/../cncal/hg38-blacklist.v2.bed.gz --summary sum.0.tsv
"""
import argparse, os, sys, time
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
SECACTS = os.path.dirname(HERE)
sys.path.insert(0, SECACTS)
from secacts_env import DATAROOT, cache_path                      # noqa: E402
sys.path.insert(0, os.path.join(SECACTS, "cnrose"))
from cnrose.callsuper import call_super                            # noqa: E402
from cnrose.cn.base import correct                                 # noqa: E402
from cnrose.cn.depmap import load_gene_coords, DepMapGeneCN, DepMapMcWesCN   # noqa: E402
from cnrose.cn.cmp import CellModelPassportsWesCN                  # noqa: E402
from cnrose.cn.inferred import BinnedInputCN, load_blacklist       # noqa: E402
from cnrose.cn.segfile import SegmentFileCN                        # noqa: E402
from cnrose.pipeline import _write_catalog                         # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--manifest", required=True, help="TSV: srx, key, cn_provider, cvcl, cn_cvcl")
    ap.add_argument("--out-dir", required=True, help="the pull's out/ (sharded out/<last 2 chars>/)")
    ap.add_argument("--task", type=int, default=0)
    ap.add_argument("--ntasks", type=int, default=1)
    ap.add_argument("--summary", required=True)
    ap.add_argument("--cn-floor", type=float, default=1.0, help="1.0 = amplify-only (calling time)")
    ap.add_argument("--cn-gene-csv", default=os.path.join(DATAROOT, "DepMap/2026q1/OmicsCNGeneWGS.csv"))
    ap.add_argument("--cmp-wes", default=os.path.join(DATAROOT, "CellModelPassports/WES_pureCN_CNV_genes_latest.csv.gz"))
    ap.add_argument("--cmp-model-list", default=os.path.join(DATAROOT, "CellModelPassports/model_list_20240110.csv"))
    ap.add_argument("--mc-wes", default=os.path.join(DATAROOT, "DepMap/2026q1/OmicsCNGeneMC_WES.csv"))
    ap.add_argument("--model-condition", default=os.path.join(DATAROOT, "DepMap/2026q1/ModelCondition.csv"))
    ap.add_argument("--ccle-segments", default=os.path.join(SECACTS, "phase2/data/cn_ccle_snp6.hg38.tsv.gz"))
    ap.add_argument("--inferred-bins")
    ap.add_argument("--inferred-map", default=os.path.join(SECACTS, "phase2/analysis/out/cn_admissions.tsv"))
    ap.add_argument("--inferred-slope", type=float, default=0.81)
    ap.add_argument("--blacklist")
    ap.add_argument("--gtf", default=os.path.join(DATAROOT, "0.human_genome/Homo_sapiens.GRCh38.106.chr.gtf.gz"))
    ap.add_argument("--gene-cache", default=cache_path("gene_coords.GRCh38.106.tsv"))
    ap.add_argument("--force", action="store_true", help="rewrite existing outputs")
    ap.add_argument("--mode", choices=["cn", "fused"], default="cn",
                    help="cn: <SRX>.cn.{enhancers.tsv,se.bed}; fused: <SRX>.fu.se.bed (labelled)")
    a = ap.parse_args()

    man = pd.read_csv(a.manifest, sep="\t").drop_duplicates("srx")
    keys = sorted(man["key"].dropna().unique())
    mine = set(keys[a.task::a.ntasks])                              # split by LINE so each task builds a track once
    man = man[man["key"].isin(mine)].sort_values(["key", "srx"])
    lines = man.drop_duplicates("key").set_index("key")
    by_src = {p: [k for k in lines.index if lines.at[k, "cn_provider"] == p]
              for p in ("depmap_wgs", "cmp_wes", "depmap_mc_wes", "ccle_snp6", "input_inferred")}
    print(f"[recall] task {a.task}/{a.ntasks}: {len(man)} samples on {len(lines)} lines; "
          + ", ".join(f"{p}={len(v)}" for p, v in by_src.items()), file=sys.stderr, flush=True)

    gc = load_gene_coords(a.gtf, cache_path=a.gene_cache)
    prov = DepMapGeneCN(a.cn_gene_csv, gc)
    if by_src["depmap_wgs"]:
        prov.preload(by_src["depmap_wgs"])
    cvcl_of = {k: (lines.at[k, "cn_cvcl"] if isinstance(lines.at[k, "cn_cvcl"], str) else lines.at[k, "cvcl"])
               for k in lines.index}
    cmp_prov = mcw_prov = ccle_prov = inf_prov = None
    if by_src["cmp_wes"]:
        cmp_prov = CellModelPassportsWesCN(a.cmp_wes, a.cmp_model_list, cache_dir=cache_path("cmp_wes"))
        cmp_prov.preload([cvcl_of[k] for k in by_src["cmp_wes"]])
    if by_src["depmap_mc_wes"]:
        mcw_prov = DepMapMcWesCN(a.mc_wes, a.model_condition, gc)
        mcw_prov.preload(by_src["depmap_mc_wes"])
    if by_src["ccle_snp6"]:
        ccle_prov = SegmentFileCN(a.ccle_segments, "ccle_snp6")
    if by_src["input_inferred"]:
        if not (a.inferred_bins and a.blacklist):
            sys.exit("[recall] input_inferred lines need --inferred-bins and --blacklist")
        im = pd.read_csv(a.inferred_map, sep="\t").dropna(subset=["cn_input"])
        inf_prov = BinnedInputCN(a.inferred_bins, dict(zip(im["key"], im["cn_input"])),
                                 blacklist=load_blacklist(a.blacklist), slope=a.inferred_slope)

    def track_for(key):
        src = lines.at[key, "cn_provider"]
        if src == "ccle_snp6":
            return ccle_prov.track(key)
        if src == "input_inferred":
            return inf_prov.track(key)
        if src == "cmp_wes":
            return cmp_prov.track(cvcl_of[key])
        if src == "depmap_mc_wes":
            return mcw_prov.track(key)
        return prov.track(key)

    rows, t0, track, cur = [], time.time(), None, None
    for r in man.itertuples(index=False):
        if r.key != cur:
            cur, track = r.key, track_for(r.key)
        base = os.path.join(a.out_dir, r.srx[-2:], r.srx)
        row = dict(srx=r.srx, key=r.key, cn_provider=lines.at[r.key, "cn_provider"])
        if not os.path.exists(base + ".enhancers.tsv"):
            rows.append({**row, "status": "no_enhancers"}); continue
        if track is None:
            rows.append({**row, "status": "no_cn_track"}); continue
        done_file = base + (".fu.se.bed" if a.mode == "fused" else ".cn.se.bed")
        if os.path.exists(done_file) and not a.force:
            rows.append({**row, "status": "exists"}); continue
        E = pd.read_csv(base + ".enhancers.tsv", sep="\t")
        sig = E["SIGNAL"].to_numpy(float)
        cn = np.array([track.region_cn(c, s, e, agg="wlen") for c, s, e in
                       zip(E["CHROM"], E["START"], E["STOP"])])
        csig = correct(sig, cn, model="log2offset", beta=1.0, floor=a.cn_floor)
        cut, sup = call_super(csig)
        ag = E["isSuper"].to_numpy(int).astype(bool)
        if a.mode == "fused":
            sup = np.asarray(sup, bool)
            fu = sig > cut                                          # corrected cutoff, uncorrected signal
            lab = np.where(~sup, "amplified", np.where(ag, "core", "unmasked"))
            order = np.argsort(-sig, kind="stable")
            rank = np.empty(len(sig), int); rank[order] = np.arange(1, len(sig) + 1)
            with open(base + ".fu.se.bed", "w") as fh:
                for i in np.flatnonzero(fu):
                    fh.write(f"{E.CHROM.iat[i]}\t{E.START.iat[i]}\t{E.STOP.iat[i]}\tSE_{i}\t{sig[i]:.6g}\t{rank[i]}\t{lab[i]}\n")
            rows.append({**row, "status": "ok", "n_regions": len(E), "n_super": int(ag.sum()), "n_super_cn": int(sup.sum()),
                         "n_fused": int(fu.sum()), "n_core": int((fu & sup & ag).sum()),
                         "n_amplified": int((fu & ~sup).sum()), "n_unmasked": int((sup & ~ag).sum()),
                         "n_agnostic_lost": int((ag & ~fu).sum()), "cn_cutoff": cut,
                         "max_region_cn": float(cn.max()) if len(cn) else np.nan})
            continue
        regions = [dict(chrom=c, start=s, end=e, num_loci=n, constituent_size=z) for c, s, e, n, z in
                   zip(E["CHROM"], E["START"], E["STOP"], E["NUM_LOCI"], E["CONSTITUENT_SIZE"])]
        n_cn = _write_catalog(base, ".cn", regions, csig, sup, extra_cols=[("REGION_CN", cn)])
        _, sup0 = call_super(sig)                                   # the agnostic call must reproduce exactly
        rows.append({**row, "status": "ok", "agnostic_reproduced": bool((np.asarray(sup0, bool) == ag).all()), "n_regions": len(E), "n_super": int(ag.sum()), "n_super_cn": n_cn,
                     "cn_cutoff": cut, "only_agnostic": int((ag & ~sup).sum()), "only_corrected": int((~ag & sup).sum()),
                     "max_region_cn": float(cn.max()) if len(cn) else np.nan})
        if len(rows) % 100 == 0:
            print(f"[recall]   {len(rows)}/{len(man)} ({time.time() - t0:.0f}s)", file=sys.stderr, flush=True)
    S = pd.DataFrame(rows)
    S.to_csv(a.summary, sep="\t", index=False)
    print(f"[recall] task {a.task} done in {time.time() - t0:.0f}s: {S['status'].value_counts().to_dict()}",
          file=sys.stderr)


if __name__ == "__main__":
    main()
