"""cnrose CLI.

    cnrose call --bw SRX.bw --peaks SRX.20.bed --out SRX [--grid grid.20.bed] [--window 12500]

Run with the atac_hdac env (pyBigWig + numpy + scipy):
    ~/miniconda3/envs/atac_hdac/bin/python -m cnrose.cli call ...
"""
from __future__ import annotations

import argparse
import sys

from .pipeline import call_sample


def main(argv=None):
    p = argparse.ArgumentParser(prog="cnrose", description="bigWig-native ROSE-style super-enhancer caller")
    sub = p.add_subparsers(dest="cmd", required=True)

    c = sub.add_parser("call", help="call super-enhancers for one sample")
    c.add_argument("--bw", required=True, help="per-sample H3K27ac bigWig")
    c.add_argument("--peaks", required=True, help="constituent peak BED (e.g. ChIP-Atlas bed20)")
    c.add_argument("--out", required=True, help="output prefix")
    c.add_argument("--grid", action="append", default=None,
                   help="fixed quantification grid BED (writes {out}.signal.tsv). REPEATABLE: every grid is "
                        "quantified in the same pass while the bigWig is open — free against a "
                        "download-bound pull, and it keeps a looser-threshold future exact "
                        "(PULL_DESIGN.md §8.5.2). Later grids write {out}.<gridlabel>.signal.tsv/.f32.")
    c.add_argument("--signal-format", choices=["tsv", "f32", "both"], default="tsv",
                   help="grid column format. 'f32' = raw little-endian float32 in grid row order: exact and "
                        "7.5x smaller than TSV (15 MB -> 2.0 MB/sample), which is what the pull should use.")
    c.add_argument("--bins", type=int, default=None, metavar="BP",
                   help="also write a genome-wide binned archive {out}.bin<BP>.f16.gz (+ .json layout). "
                        "100 recommended. This is INSURANCE that retires the irreversible grid-threshold "
                        "choice — any future region set is recoverable to ~1%% median without re-downloading "
                        "— but NOT a replacement for exact quantification (see "
                        "phase2/analysis/binned_storage_eval.py).")
    c.add_argument("--chunk-mb", type=float, default=5.0,
                   help="bigWig interval-read window in Mb (default 5). Peak RSS grows with per-window "
                        "interval density, so on very dense/large bigWigs (~1 GB+) a smaller window keeps the "
                        "1-CPU/2-GiB task under its memory limit — see PULL_DESIGN.md §4.1. No BU cost.")
    c.add_argument("--qc", action="store_true",
                   help="write {out}.qc.json: signal percentiles, dynamic range (p99/median — the statistic "
                        "S3norm's exponent tracks), FRiP proxy, tangent cutoff, peak/SE counts. Only the "
                        "bigWig can provide these, so emit them before it is deleted.")
    c.add_argument("--window", type=int, default=12500, help="stitch distance in bp (default 12500)")
    c.add_argument("--agg", choices=["sum", "mean"], default="sum", help="region signal aggregation")
    c.add_argument("--all-chroms", action="store_true", help="keep non-canonical chroms (default: chr1..22,X,Y)")
    # CN correction (DESIGN.md §4): optional; produces a second, CN-corrected catalog ({out}.cn.*)
    c.add_argument("--cn", choices=["depmap"], default=None, help="CN backend (calling-time correction)")
    c.add_argument("--cn-key", default=None, help="cell-line key for the CN backend (DepMap ModelID)")
    c.add_argument("--cn-gene-csv", default=None, help="DepMap OmicsCNGeneWGS.csv")
    c.add_argument("--cn-gtf", default=None, help="Ensembl GRCh38 GTF for gene coordinates")
    c.add_argument("--cn-gene-cache", default=None, help="gene-coord TSV cache (written/reused)")
    c.add_argument("--cn-model", choices=["log2offset", "divide"], default="log2offset")
    c.add_argument("--cn-agg", choices=["wlen", "max", "median"], default="wlen")
    c.add_argument("--cn-beta", type=float, default=1.0)
    c.add_argument("--cn-floor", type=float, default=1.0,
                   help="clamp CN before correcting; 1.0=amplify-only (calling default), 0.1=symmetric")

    a = p.parse_args(argv)
    if a.cmd == "call":
        cn_provider = _build_cn_provider(a) if a.cn else None
        correct_at = "calling" if cn_provider else "none"
        s = call_sample(a.bw, a.peaks, a.out, grid_bed=a.grid, window=a.window,
                        agg=a.agg, canonical_only=not a.all_chroms,
                        cn_provider=cn_provider, cn_key=a.cn_key, correct_at=correct_at,
                        cn_model=a.cn_model, cn_agg=a.cn_agg, cn_beta=a.cn_beta, cn_floor=a.cn_floor,
                        signal_format=a.signal_format, bins=a.bins, qc=a.qc, chunk=int(a.chunk_mb*1_000_000))
        msg = (f"[cnrose] {a.out}: {s['n_super']} SEs / {s['n_regions']} stitched regions "
               f"(cutoff={s['cutoff']:.6g})")
        if s.get("grids"):
            msg += " | grids: " + ", ".join(f"{k}={v}" for k, v in s["grids"].items())
        if "archive_bins" in s:
            msg += f" | archive {s['archive_bins']} bins"
        if s.get("qc"):
            dr = s["qc"].get("dynamic_range_p99_over_median")
            fr = s["qc"].get("frip_proxy")
            msg += f" | QC dyn.range={dr:.1f}" if dr else ""
            msg += f" FRiP~{fr:.3f}" if fr else ""
        if "n_super_cn" in s:
            msg += (f" | CN-corrected: {s['n_super_cn']} SEs (cutoff={s['cn_cutoff']:.6g}); "
                    f"agnostic-only={s['cn_only_agnostic']}, corrected-only={s['cn_only_corrected']}")
        elif "cn_status" in s:
            msg += f" | {s['cn_status']}"
        print(msg, file=sys.stderr)


def _build_cn_provider(a):
    if a.cn == "depmap":
        if not (a.cn_gene_csv and a.cn_gtf):
            sys.exit("[cnrose] --cn depmap needs --cn-gene-csv and --cn-gtf")
        from .cn.depmap import load_gene_coords, DepMapGeneCN
        coords = load_gene_coords(a.cn_gtf, cache_path=a.cn_gene_cache)
        return DepMapGeneCN(a.cn_gene_csv, coords)
    return None


if __name__ == "__main__":
    main()
