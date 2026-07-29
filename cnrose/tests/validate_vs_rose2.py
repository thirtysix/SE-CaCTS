"""Validation contract (DESIGN.md §7): cnrose --cn none must reproduce ROSE2 on the pilot bigWigs.

Two independent checks per SRX, using the real ROSE2 as oracle (no BAM needed):
  1. STITCH: cnrose.stitch(constituents)  vs  rose2.utils stitchCollection (se_rose env).  -> exact match?
  2. CUTOFF: feed a BYTE-IDENTICAL signal table (built once over the ROSE-stitched regions) to both
     cnrose.call_super and the real ROSE2_callSuper.R.  -> SE-set Jaccard, isSuper agreement.

Run under atac_hdac (pyBigWig+numpy+scipy+Rscript):
  ~/miniconda3/envs/atac_hdac/bin/python tests/validate_vs_rose2.py \
      --bw-dir <pilot bigwigs> --bed20-dir <rehearse/bed20> --srx-list <pilot_srx.txt> --workdir <scratch>
"""
from __future__ import annotations

import argparse
import os
import subprocess
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from cnrose.io import quantify, read_bed3, CANONICAL
from cnrose.stitch import stitch
from cnrose.callsuper import call_super

HERE = os.path.dirname(os.path.abspath(__file__))
SE_ROSE_PY = os.path.expanduser("~/miniconda3/envs/se_rose/bin/python")
ROSE_CALLSUPER = os.path.expanduser("~/miniconda3/envs/se_rose/bin/ROSE2_callSuper.R")


def jaccard(a, b):
    a, b = set(a), set(b)
    return 1.0 if not a and not b else len(a & b) / len(a | b)


def rose_stitch(constituents_bed, window):
    out = subprocess.run([SE_ROSE_PY, os.path.join(HERE, "_rose_stitch_helper.py"),
                          constituents_bed, str(window)], capture_output=True, text=True)
    if out.returncode != 0:
        raise RuntimeError(f"rose stitch helper failed:\n{out.stderr}")
    regs = []
    for line in out.stdout.splitlines():
        c, s, e = line.split("\t")
        regs.append((c, int(s), int(e)))
    return regs


def run_rose_callsuper(regions, signal, workdir, name):
    """Write the ROSE ENHANCER_REGION_MAP and run the real ROSE2_callSuper.R; return set of super REGION_IDs."""
    mapfile = os.path.join(workdir, f"{name}_REGION_MAP.txt")
    with open(mapfile, "w") as fh:
        fh.write("REGION_ID\tCHROM\tSTART\tSTOP\tNUM_LOCI\tCONSTITUENT_SIZE\tSIGNAL\n")
        for i, ((c, s, e), sig) in enumerate(zip(regions, signal)):
            fh.write(f"reg_{i}\t{c}\t{s}\t{e}\t1\t{e - s}\t{sig:.6f}\n")
    env = dict(os.environ, R_DEFAULT_DEVICE="png")
    out = subprocess.run([ROSE_CALLSUPER, workdir + os.sep, mapfile, name, "NONE"],
                         capture_output=True, text=True, env=env)
    super_tbl = os.path.join(workdir, f"{name}_SuperEnhancers.table.txt")
    if not os.path.exists(super_tbl):
        raise RuntimeError(f"ROSE2_callSuper.R produced no SuperEnhancers table (rc={out.returncode}).\n"
                           f"STDOUT tail:\n{out.stdout[-800:]}\nSTDERR tail:\n{out.stderr[-800:]}")
    supers = set()
    with open(super_tbl) as fh:
        for line in fh:
            if line.startswith("#") or line.startswith("REGION_ID"):
                continue
            f = line.split("\t")
            if f and f[0].startswith("reg_"):
                supers.add(f[0])
    return supers


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--bw-dir", required=True)
    ap.add_argument("--bed20-dir", required=True)
    ap.add_argument("--srx-list", required=True)
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--window", type=int, default=12500)
    a = ap.parse_args()
    os.makedirs(a.workdir, exist_ok=True)

    with open(a.srx_list) as fh:
        srxs = [x.strip() for x in fh if x.strip()]

    print(f"{'SRX':<14}{'peaks':>8}{'stitch_cn':>10}{'stitch_rose':>12}{'stitch_ok':>10}"
          f"{'SE_cn':>7}{'SE_rose':>8}{'SE_jacc':>9}")
    rows = []
    for srx in srxs:
        bw = os.path.join(a.bw_dir, f"{srx}.bw")
        bed20 = os.path.join(a.bed20_dir, f"{srx}.20.bed")
        if not (os.path.exists(bw) and os.path.exists(bed20)):
            print(f"{srx:<14}  MISSING bw or bed20 -> skip"); continue

        # constituents: cols 1-3, canonical chroms
        peaks = [(c, s, e) for c, s, e in read_bed3(bed20) if c in CANONICAL]
        constituents = os.path.join(a.workdir, f"{srx}.constituents.bed")
        with open(constituents, "w") as fh:
            for c, s, e in sorted(peaks):
                fh.write(f"{c}\t{s}\t{e}\n")

        # 1) stitch
        cn_regs = [(r["chrom"], r["start"], r["end"]) for r in stitch(peaks, window=a.window)]
        rose_regs = rose_stitch(constituents, a.window)
        stitch_ok = set(cn_regs) == set(rose_regs)

        # 2) cutoff — identical regions + signal to both callers (use the ROSE regions)
        signal = quantify(bw, rose_regs, agg="sum")
        _, is_super_cn = call_super(signal)
        cn_super = {f"reg_{i}" for i, v in enumerate(is_super_cn) if v}
        rose_super = run_rose_callsuper(rose_regs, signal, a.workdir, srx)
        se_jacc = jaccard(cn_super, rose_super)

        print(f"{srx:<14}{len(peaks):>8}{len(cn_regs):>10}{len(rose_regs):>12}"
              f"{('YES' if stitch_ok else 'NO'):>10}{len(cn_super):>7}{len(rose_super):>8}{se_jacc:>9.4f}")
        rows.append((srx, stitch_ok, len(cn_super), len(rose_super), se_jacc))

    print("-" * 78)
    n = len(rows)
    stitch_pass = sum(1 for r in rows if r[1])
    mean_jacc = np.mean([r[4] for r in rows]) if rows else 0.0
    min_jacc = min((r[4] for r in rows), default=0.0)
    print(f"stitch exact-match: {stitch_pass}/{n}   SE-set Jaccard mean={mean_jacc:.4f} min={min_jacc:.4f}")
    ok = stitch_pass == n and min_jacc >= 0.98
    print(f"VALIDATION: {'PASS' if ok else 'REVIEW'} (contract: stitch all-exact AND min Jaccard >= 0.98)")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
