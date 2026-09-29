#!/usr/bin/env python3
"""Cut the DepMap CN tables down to the models a pull set scores, for staging on the HPC (score_arm.slurm).

The full OmicsCNGeneMC_WES.csv is 1.3 GB; the scorer only reads rows of the models it scores. Keeps rows of
every ModelID in the pull set's `model_id` and `key` columns (WGS by ModelID; MC_WES through ModelCondition).
Usage:  71_filter_depmap_cn.py <pull_set.tsv> <out_dir>
"""
import csv
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..", "..")))
from secacts_env import DATAROOT                                  # noqa: E402

DEPMAP = os.path.join(DATAROOT, "DepMap", "2026q1")
csv.field_size_limit(10 ** 9)


def main(pull_set, out):
    os.makedirs(out, exist_ok=True)
    with open(pull_set) as fh:
        models = {v for r in csv.DictReader(fh, delimiter="\t") for v in (r["model_id"], r["key"])
                  if v and v.startswith("ACH-")}
    with open(os.path.join(DEPMAP, "ModelCondition.csv"), newline="") as fh:
        mc_ok = {r["ModelConditionID"] for r in csv.DictReader(fh) if r.get("ModelID") in models}
    for name, col, ok in (("OmicsCNGeneWGS.csv", "ModelID", models),
                          ("OmicsCNGeneMC_WES.csv", "ModelConditionID", mc_ok)):
        n = kept = 0
        with open(os.path.join(DEPMAP, name), newline="") as src, open(os.path.join(out, name), "w", newline="") as dst:
            header = src.readline(); dst.write(header)
            ci = next(csv.reader([header])).index(col)
            for line in src:
                n += 1
                if next(csv.reader([line]))[ci] in ok:
                    dst.write(line); kept += 1
        print(f"[71] {name}: kept {kept}/{n} rows")
    print(f"[71] models in pull set: {len(models)}")


if __name__ == "__main__":
    main(*sys.argv[1:3])
