#!/usr/bin/env python3
"""Merge the v3.1 Sonnet relabels and pick the rows Opus adjudicates (v3.1 item 2).

Every target was labelled "control" before. Opus reads:
  - every row the relabel moved out of keep-tier-1 (class != untreated/control, or tier 2, or low confidence);
  - rows still kept whose TARGET text names a treatment (a keyword net: some Sonnet labellers skimmed study
    context, so a missed "genetic control + drug" arm can only be caught from the target's own text).
Writes `merged.tsv` and `opus_rows_{a,b}.txt` (blocks copied from the chunks, split in two).

    python3 phase2/analysis/v31_relabel_merge.py
"""
import glob, os, re
import pandas as pd

SECACTS = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
D = os.path.join(SECACTS, "phase2/analysis/out/v31_relabel")
TREAT = re.compile(
    r"\b(DHT|R1881|E2|estradiol|oestradiol|dex|dexamethasone|PMA|TPA|ionomycin|LPS|TNF|IFN|IL-?\d+|TGF|"
    r"dox|doxycycline|tamoxifen|4-?OHT|fulvestrant|enzalutamide|JQ1|THZ1|SAHA|vorinostat|decitabine|azacytidine|"
    r"inhibitor|degrader|PROTAC|dTAG|IAA|auxin|drug|treated|treatment=(?!none|untreated|no|control|vehicle|dmso)|"
    r"stimulat\w*|induced|irradiat\w*|hypoxi\w*|nM|uM|µM|mg/ml|ug/ml|hrs?\b|\d+\s?h\b|\d+\s?hr|\d+\s?min)",
    re.I)


def blocks():
    out = {}
    for p in sorted(glob.glob(os.path.join(D, "chunk_*.txt"))):
        for b in open(p).read().split("\n\n"):
            m = re.match(r"=== (\S+)", b.strip())
            if m:
                out[m.group(1)] = b.strip()
    return out


def main():
    L = pd.concat([pd.read_csv(p, sep="\t", dtype=str) for p in sorted(glob.glob(os.path.join(D, "labels_*.tsv")))])
    T = pd.read_csv(os.path.join(D, "targets.tsv"), sep="\t", dtype=str).rename(columns={"class": "prev_class"})
    assert set(L.srx) == set(T.srx) and not L.srx.duplicated().any(), "labels do not cover the targets 1:1"
    M = T.merge(L, on="srx")
    B = blocks()
    target_line = {s: b.split("\n")[1] for s, b in B.items()}
    keep = M["class"].isin(["untreated", "control"])
    moved = ~keep | (M["tier"].astype(str) == "2") | (M["confidence"] == "low")
    net = keep & ~moved & M.srx.map(lambda s: bool(TREAT.search(target_line[s].split("|", 2)[-1])))
    M["opus"] = (moved | net).map({True: "yes", False: "no"})
    M["why_opus"] = ""
    M.loc[moved, "why_opus"] = "moved"
    M.loc[net, "why_opus"] = "keyword"
    M.to_csv(os.path.join(D, "merged.tsv"), sep="\t", index=False)
    print("Sonnet classes:", M["class"].value_counts().to_dict())
    print("tier among kept:", M.loc[keep, "tier"].value_counts().to_dict())
    print(f"to Opus: {int(moved.sum())} moved + {int(net.sum())} keyword-flagged = {int((moved | net).sum())}")
    sel = M[M.opus == "yes"].sort_values(["group", "srx"])
    half = (len(sel) + 1) // 2
    for tag, part in (("a", sel.iloc[:half]), ("b", sel.iloc[half:])):
        with open(os.path.join(D, f"opus_rows_{tag}.txt"), "w") as fh:
            fh.write("\n\n".join(B[s] for s in part.srx) + "\n")
        print(f"opus_rows_{tag}.txt: {len(part)} rows, {sum(len(B[s]) for s in part.srx) // 4:,} tokens approx")


if __name__ == "__main__":
    main()
