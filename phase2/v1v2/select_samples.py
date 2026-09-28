#!/usr/bin/env python3
"""Pick the samples for the ChIP-Atlas v1 vs v2 processing comparison (ROADMAP Phase 3c, FINDINGS §13).

From the v2 atlas (S3norm, QC-gated): ChIP-Atlas v1 QC of 10-40M reads (bounded compute) and >= 80% mapped,
library layout and FASTQ URLs from ENA; half single-end, half paired-end; one per cell line; spread over
lineages, with lymphoid lines over-represented because HLA / alt-haplotype loci are the suspected v2 bias.
Deterministic (fixed seed). Output: phase2/v1v2/samples.tsv.

  ~/miniconda3/envs/atac_hdac/bin/python phase2/v1v2/select_samples.py [--n 20]
"""
import argparse, io, os, subprocess, sys, time
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
PRIORITY = {"Lymphoid": 3, "Myeloid": 2, "Breast": 2, "Lung": 2, "Bowel": 2, "Peripheral Nervous System": 2,
            "Ovary/Fallopian Tube": 1, "Skin": 1, "Bone": 1, "CNS/Brain": 1, "Esophagus/Stomach": 1, "Prostate": 1,
            "Kidney": 1, "Pancreas": 1, "Liver": 1}

def ena(srx):
    out = []
    for i in range(0, len(srx), 150):
        q = " OR ".join(f'experiment_accession="{s}"' for s in srx[i:i + 150])
        r = subprocess.run(["curl", "-s", "-X", "POST", "https://www.ebi.ac.uk/ena/portal/api/search",
                            "--data-urlencode", "result=read_run", "--data-urlencode",
                            "fields=experiment_accession,run_accession,library_layout,read_count,fastq_ftp,fastq_bytes",
                            "--data-urlencode", f"query={q}", "--data-urlencode", "format=tsv"],
                           capture_output=True, text=True, timeout=180)
        if r.stdout.strip():
            out.append(pd.read_csv(io.StringIO(r.stdout), sep="\t"))
        time.sleep(0.3)
    return pd.concat(out)

def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--n", type=int, default=20); a = ap.parse_args()
    man = pd.read_csv(f"{ROOT}/phase1/data/phase1_manifest.tsv", sep="\t").set_index("srx")
    ps = pd.read_csv(f"{ROOT}/phase2/data/pull_set.v2.tsv", sep="\t")
    kept = set(pd.read_csv(f"{ROOT}/phase2/results_v2/atlas.s3.s3norm_params.tsv.gz", sep="\t")["sample"])
    d = ps[ps.srx.isin(kept)].copy()
    d["reads"] = d.srx.map(man.reads); d["pct_mapped"] = d.srx.map(man.pct_mapped)
    d = d[(d.reads.between(10e6, 40e6)) & (d.pct_mapped >= 80)]
    runs = ena(list(d.srx))
    g = runs.groupby("experiment_accession").agg(layout=("library_layout", "first"), n_runs=("run_accession", "size"),
                                                  ftp=("fastq_ftp", lambda x: ";".join(map(str, x))))
    d = d.join(g, on="srx").dropna(subset=["layout"])
    d = d[d.ftp.str.contains("fastq", na=False)]
    rng = np.random.default_rng(20260928)
    d = d.sample(frac=1, random_state=20260928).drop_duplicates("key")         # one experiment per line
    per_layout = a.n // 2
    pick = []
    for layout in ("SINGLE", "PAIRED"):
        pool = d[d.layout == layout]
        quota = {k: max(1, round(v * per_layout / sum(PRIORITY.values()))) for k, v in PRIORITY.items()}
        chosen = []
        for lin, q in sorted(quota.items(), key=lambda kv: -kv[1]):
            chosen += list(pool[pool.lineage == lin].head(q).srx)
        rest = pool[~pool.srx.isin(chosen)]
        chosen = (chosen + list(rest.srx))[:per_layout]
        pick += chosen
    out = d[d.srx.isin(pick)][["srx", "cell", "key", "lineage", "subtype", "layout", "n_runs", "reads", "pct_mapped", "ftp"]]
    out.sort_values(["layout", "lineage"]).to_csv(f"{HERE}/samples.tsv", sep="\t", index=False)
    print(out.drop(columns="ftp").sort_values(["layout", "lineage"]).to_string(index=False))
    print(f"\n{len(out)} samples: {out.layout.value_counts().to_dict()}", file=sys.stderr)

if __name__ == "__main__":
    main()
