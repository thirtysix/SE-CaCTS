#!/usr/bin/env python3
"""v3.1 item 2: relabel v3 "control" experiments that sit in studies with treated arms (FINDINGS §36).

Writes chunked inputs for Sonnet labellers (`out/v31_relabel/chunk_NN.txt`, one TARGET block per experiment with
the rest of its study as context, like the §36 adjudication rows) and the shared prompt (`PROMPT.md`), which
carries the tier policy the user adopted on 2026-09-30 and the rule fix from §36. Labellers write
`labels_NN.tsv`; `v31_relabel_merge.py` compares them with the v3 labels and lists disagreements for Opus.

    python3 phase2/analysis/v31_relabel_build.py --chunk 95
"""
import argparse, glob, os
import pandas as pd

SECACTS = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
O = os.path.join(SECACTS, "phase2/analysis/out")
MAX_CONTEXT = 14

PROMPT = """# v3.1 treatment relabel: instructions

You label H3K27ac ChIP-seq experiments on cancer cell lines for a BASELINE-ONLY atlas: we want each line's
unperturbed enhancer landscape. Each block in your chunk file is one TARGET experiment (its ChIP-Atlas / GEO
metadata: title | attribute string) followed by the other experiments of the same study, for context. Every
TARGET was previously labelled "control" by a model that made one systematic error (below). Label the TARGET only.

Read the WHOLE attribute string: treatment fields often sit at the end of long ENCODE-style strings
(e.g. "... treatment=5-Ph-IAA"), and the study's other arms usually reveal what the design is.

## Classes
- `untreated`  no treatment, standard culture, no genetic modification that is the tested variable.
- `control`    the control arm of a design, and the cells are otherwise baseline: vehicle (DMSO, EtOH, PBS),
               non-targeting / scrambled si/sh/sg RNA, empty vector, mock, 0 h, uninduced.
- `perturbed`  any drug, inhibitor, degrader, hormone, ligand, cytokine stimulation, irradiation, hypoxia, or a
               genetic knockdown / knockout / overexpression / induced degron that is the tested variable.
               **A control arm counts as control only if the cells are otherwise baseline. A GENETIC control arm
               that also receives a drug, hormone or stimulation is `perturbed`** (siCTRL + DHT, siNT + E2,
               sgCtrl + PMA/ionomycin, shScr + doxorubicin): the treatment decides. This was the previous
               labeller's main error, so check every TARGET for it.
               Always perturbed too: synchronised populations; GSI "mock washout" / reference arms that stay on
               drug; FACS-sorted subpopulations; resistant, in-vivo-selected or otherwise derived sublines of a
               parent line; cell-fusion hybrids; engineered insertions that are the tested variable.
- `not_h3k27ac` an input, IgG, or another mark/protein, whatever the title says. If the title says H3K27ac but
               the antibody field names another protein and the study does not resolve it, use `unclear` with
               reason "title-antibody conflict".
- `out_of_scope` no cancer cell line (hESC/iPSC and their derivatives, primary tissue, normal cells).
- `unclear`    the metadata cannot decide; say why.

## Tier (only for `untreated` and `control`; else `-`)
- `1` keep: untreated; vehicle / non-targeting / empty-vector controls; the routine maintenance cytokine of a
  factor-dependent line (F-36P GM-CSF, UT-7 EPO, TF-1 GM-CSF); the native virus of a virus-positive line (BCBL1
  KSHV, uninduced RTA); fixation or method variants (FA vs DSG crosslinking, crosslink titration).
- `2` keep only when the line has no tier-1 experiment: dox-induced reporter / GFP / BFP control arms;
  empty-vector single-cell clones; defined-medium control arms; dCas9-only, dCas9-KRAB/VP64 with a non-targeting
  guide, or cutting-control guides at expressed genes; inactive-compound controls (inactive enantiomer, inactive
  PROTAC).

## Also flag in `note`
- `derivative:<name>` when the TARGET is a catalogued derivative with its own identity (e.g. NCI/ADR-RES, not
  OVCAR-8): it is its own line, never pooled with the parent.

## Output
Write a TSV with header `srx	class	tier	confidence	reason	note` (tab-separated, one row per TARGET, every
TARGET exactly once; confidence high/medium/low; reason under 15 words; note may be `-`). Do not label the
context experiments.
"""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--chunk", type=int, default=95)
    ap.add_argument("--out", default=os.path.join(O, "v31_relabel"))
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    inp = pd.concat([pd.read_csv(p, sep="\t") for p in sorted(glob.glob(os.path.join(O, "atlas_labels/input_*.tsv")))])
    inp = inp.drop_duplicates("srx")
    lab = pd.read_csv(os.path.join(O, "atlas_treatment_labels.tsv"), sep="\t", usecols=["srx", "class"])
    extra = [pd.read_csv(p, sep="\t", usecols=["srx", "class"])
             for p in sorted(glob.glob(os.path.join(O, "atlas_labels/labels_1[0-9]*.tsv")))]
    L = pd.concat([lab, *extra]).drop_duplicates("srx", keep="last").set_index("srx")["class"]
    inp["class"] = inp.srx.map(L)
    ps = pd.read_csv(os.path.join(SECACTS, "phase2/data/pull_set.v3.tsv"), sep="\t")
    treated = set(inp.loc[inp["class"] == "perturbed", "group"])
    T = inp[(inp["class"] == "control") & inp.group.isin(treated) & inp.srx.isin(ps.srx)]
    T = T.sort_values(["group", "srx"])
    T[["srx", "cell_line", "group", "class"]].to_csv(os.path.join(a.out, "targets.tsv"), sep="\t", index=False)

    by_group = {g: d for g, d in inp.groupby("group")}
    blocks = []
    for r in T.itertuples():
        ctx = by_group[r.group]
        ctx = ctx[ctx.srx != r.srx]
        lines = [f"=== {r.srx}", f"TARGET: {r.cell_line} | antibody class: {r.atlas_antigen} | {r.title}",
                 f"STUDY {r.group}: {len(ctx)} other experiments"
                 + (f" (first {MAX_CONTEXT} shown)" if len(ctx) > MAX_CONTEXT else "")]
        lines += [f"  - {c.cell_line} | {c.title}" for c in ctx.head(MAX_CONTEXT).itertuples()]
        blocks.append("\n".join(lines))
    n = 0
    for i in range(0, len(blocks), a.chunk):
        with open(os.path.join(a.out, f"chunk_{n:02d}.txt"), "w") as fh:
            fh.write("\n\n".join(blocks[i:i + a.chunk]) + "\n")
        n += 1
    with open(os.path.join(a.out, "PROMPT.md"), "w") as fh:
        fh.write(PROMPT)
    print(f"{len(T)} targets in {T.group.nunique()} studies -> {n} chunks of <= {a.chunk} in {a.out}")


if __name__ == "__main__":
    main()
