#!/usr/bin/env python3
"""Add the experiments that phase2/data/srx_layout.tsv lacks, from ENA run metadata (one row per experiment).

The table was built for v2; v3/v3.1 experiments missing from it read an empty layout on the dashboard and are left
out of the layout arm (FINDINGS §16). Same columns and conventions as the existing rows: runs summed per
experiment, read_len = bases / reads (halved for PAIRED), instrument of the experiment's first run.

    python3 phase2/analysis/srx_layout_fill.py --signal phase2/results_v31f/atlas.s3.se_signal.tsv.gz
"""
import argparse, gzip, io, os, sys, time, urllib.parse, urllib.request
import pandas as pd

SECACTS = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
ENA = "https://www.ebi.ac.uk/ena/portal/api/search"
FIELDS = "run_accession,experiment_accession,library_layout,read_count,base_count,instrument_model,study_accession"


def ena(accs):
    q = " OR ".join(f'experiment_accession="{a}"' for a in accs)
    body = urllib.parse.urlencode({"result": "read_run", "query": q, "fields": FIELDS, "format": "tsv",
                                   "limit": 0}).encode()
    for k in range(4):
        try:
            with urllib.request.urlopen(urllib.request.Request(ENA, data=body), timeout=120) as r:
                txt = r.read().decode()
            return pd.read_csv(io.StringIO(txt), sep="\t") if txt.strip() else pd.DataFrame(columns=FIELDS.split(","))
        except Exception as e:                                    # ENA drops connections now and then
            print(f"[layout] retry {k + 1}: {e}", file=sys.stderr)
            time.sleep(5 * (k + 1))
    sys.exit("[layout] ENA unreachable")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--signal", required=True, help="SE x experiment matrix whose columns must all be covered")
    ap.add_argument("--table", default=os.path.join(SECACTS, "phase2/data/srx_layout.tsv"))
    ap.add_argument("--batch", type=int, default=100)
    a = ap.parse_args()
    with gzip.open(a.signal, "rt") if a.signal.endswith(".gz") else open(a.signal) as fh:
        cols = fh.readline().rstrip("\n").split("\t")[1:]
    lay = pd.read_csv(a.table, sep="\t")
    miss = [c for c in cols if c not in set(lay["srx"])]
    print(f"[layout] {len(miss)} of {len(cols)} experiments missing from {os.path.basename(a.table)}")
    if not miss:
        return
    runs = pd.concat([ena(miss[i:i + a.batch]) for i in range(0, len(miss), a.batch)], ignore_index=True)
    rows = []
    for srx, d in runs.groupby("experiment_accession", sort=False):
        reads, bases = int(d["read_count"].fillna(0).sum()), int(d["base_count"].fillna(0).sum())
        layout = d["library_layout"].iloc[0]
        rl = round(bases / reads / (2 if layout == "PAIRED" else 1), 1) if reads else None
        rows.append({"srx": srx, "layout": layout, "n_runs": len(d), "reads": reads, "bases": bases,
                     "instrument": d["instrument_model"].iloc[0], "read_len": rl, "study": d["study_accession"].iloc[0]})
    new = pd.DataFrame(rows)
    left = sorted(set(miss) - set(new["srx"]))
    print(f"[layout] ENA returned {len(new)} experiments ({new['layout'].value_counts().to_dict()}); "
          f"not found: {len(left)} {left[:10]}")
    lay = pd.concat([lay, new], ignore_index=True).sort_values("srx")
    lay.to_csv(a.table, sep="\t", index=False)
    print(f"[layout] wrote {len(lay)} rows")


if __name__ == "__main__":
    main()
