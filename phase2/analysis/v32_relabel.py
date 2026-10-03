#!/usr/bin/env python3
"""Baseline labels for the experiments the 2026-09 ChIP-Atlas metadata adds (v3.2, FINDINGS §57).

`phase1/data/v32_new_h3k27ac.tsv` lists 121 QC-pass in-scope hg38 H3K27ac experiments. The atlas is ChIP-seq only
(CUT&Tag / CUT&RUN are a different assay, not planned: ROADMAP, FINDINGS §13), and the library strategy closes
ChIP-Atlas's title field ("...; Homo sapiens; <strategy>"): `OTHER` and every CUT&RUN / CUT&Tag strategy or
attribute is `other_assay`. That leaves two ChIP-seq experiments, labelled here by the v3.1 rules (§47 prompt)
instead of a model run:
  SRX6712584  MCF-7 E2                                      perturbed (hormone)
  SRX6712585  MCF-7 vehicle                                 unclear: title says H3K27ac, the GEO record names both
              an H3K4me3 and an H3K27ac antibody, and ChIP-Atlas listed the run as H3K4me3 on hg38 until 2026-09
              (still on hg19)

    python3 phase2/analysis/v32_relabel.py
"""
import os
import re

import pandas as pd

SECACTS = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
MANUAL = {"SRX6712584": ("perturbed", "-", "E2 treatment arm"),
          "SRX6712585": ("unclear", "-", "title-antibody conflict: GEO antibody fields H3K4me3 and H3K27ac; "
                                         "ChIP-Atlas hg19 entry H3K4me3")}


def main():
    new = pd.read_csv(os.path.join(SECACTS, "phase1/data/v32_new_h3k27ac.tsv"), sep="\t")
    want = set(new.srx)
    meta = {}
    with open(os.path.join(SECACTS, ".cache/experimentList.2026-09.tab"), encoding="utf-8", errors="replace") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[0] in want and f[1] == "hg38":
                meta[f[0]] = (f[8], " ".join(f[9:]))
    rows = []
    for s in sorted(want):
        title, attrs = meta[s]
        strategy = title.split("; ")[-1]
        cutx = re.search(r"CUT&(amp;)?(RUN|Tag)", f"{title} {attrs}", re.I)
        if s in MANUAL:
            c, t, note = MANUAL[s]
        elif strategy == "OTHER" or cutx:
            c, t, note = "other_assay", "-", ("CUT&RUN/CUT&Tag" if cutx else "library strategy OTHER")
        else:
            c, t, note = "unclear", "-", "unlabelled ChIP-seq"
        rows.append({"srx": s, "final_class": c, "final_tier": t, "final_note": note, "strategy": strategy[:60]})
    out = pd.DataFrame(rows)
    d = os.path.join(SECACTS, "phase2/analysis/out/v32_relabel")
    os.makedirs(d, exist_ok=True)
    out.to_csv(os.path.join(d, "final_labels.tsv"), sep="\t", index=False)
    print(out.final_class.value_counts().to_string())


if __name__ == "__main__":
    main()
