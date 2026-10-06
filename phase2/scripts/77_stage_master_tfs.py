#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Known lineage master regulators against the atlas (the poster's Fig. 5, redrawn per release), for the dashboard
and the README.

The 22 transcription factors chosen from the literature before scoring (poster_figures.IDENTITY). For each, the
super-enhancer within 100 kb that is best in the gene's own lineage (lowest FDR, ties by JSD; poster_figures.
identity_table says why not the lowest JSD), with its JSD rank, permutation FDR and CN status in every lineage.

  <docs>/data/master_tfs.json        {release, n_ses, lineages, own: [called, of], other: [called, of],
                                      rows: [{gene, lineage, se, locus, cells: [[rank, fdr, cn], ...]}]}
                                      cn = "r" CN-robust / "u" CN-unmasked, on called cells (FDR <= 0.10) only
  <fig-out>/master-regulators.{png,svg,pdf}  the static figure; without --fig-out the JSON only

  ~/miniconda3/envs/atac_hdac/bin/python phase2/scripts/77_stage_master_tfs.py --docs docs \
      --scores phase2/scores_v31f --results phase2/results_v31f --fig-out assets/figures
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import os

SECACTS = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
FDR = 0.10


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--docs", required=True)
    ap.add_argument("--scores", required=True)
    ap.add_argument("--results", required=True)
    ap.add_argument("--labels", help="fused_labels.groups table (default <scores>/fused_labels.groups.tsv.gz if present)")
    ap.add_argument("--fig-out")
    a = ap.parse_args()
    labels = a.labels or os.path.join(a.scores, "fused_labels.groups.tsv.gz")
    labels = labels if os.path.exists(labels) else None

    spec = importlib.util.spec_from_file_location("poster_figures", os.path.join(SECACTS, "phase2", "figures", "poster_figures.py"))
    pf = importlib.util.module_from_spec(spec); spec.loader.exec_module(pf)
    cols, recs = pf.identity_table(a.scores, a.results, select="fdr", labels=labels)

    import pandas as pd
    n_ses = len(pd.read_csv(os.path.join(a.scores, "atlas.s3.perm.OncotreeLineage.fdr.tsv.gz"), sep="\t", usecols=[0]))
    code = {"robust": "r", "unmasked": "u"}
    own, other, rows = [0, 0], [0, 0], []
    for r in recs:
        cells = []
        for j, c in enumerate(cols):
            f = float(r["fdr"][j])
            cells.append([r["rank"][j], round(f, 4), code.get(r["cn"][j]) if f <= FDR else None])
            k = own if c == r["lineage"] else other
            k[0] += f <= FDR; k[1] += 1
        rows.append(dict(gene=r["gene"], lineage=r["lineage"], se=r["se"],
                         locus=f"{r['chrom']}:{r['start'] + 1}-{r['end']}", cells=cells))
    meta = json.load(open(os.path.join(a.docs, "data", "meta.json")))
    out = dict(release=meta.get("release", {}).get("version"), n_ses=n_ses, lineages=cols, own=own, other=other, rows=rows)
    with open(os.path.join(a.docs, "data", "master_tfs.json"), "w") as fh:
        json.dump(out, fh, separators=(",", ":"))
    print(f"[77] {own[0]}/{own[1]} called in their own lineage; elsewhere {other[0]}/{other[1]} "
          f"({100 * other[0] / max(other[1], 1):.1f}%); missed {[r['gene'] for r in rows if r['cells'][cols.index(r['lineage'])][1] > FDR]}")

    if a.fig_out:
        pf.OUT = a.fig_out
        pf.fig_identity(table=(cols, recs), name="master-regulators", cn_mark=labels is not None, contrast_text=True)


if __name__ == "__main__":
    main()
