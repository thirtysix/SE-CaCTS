#!/usr/bin/env bash
# Phase 2 — 20: build the fixed quantification GRID (see phase2/PULL_DESIGN.md §3).
#
# Grid = merged union of all ChIP-Atlas H3K27ac peaks belonging to our analysis-set SRX
# (qc_pass & has_cn = 2,917 experiments), on canonical chromosomes. This is the immutable,
# read-only reference every per-sample pull task quantifies against; union SEs are later
# recovered as sums of grid rows, so each bigWig is fetched only once (single-pass design).
#
# Streams the .gz directly (pigz) — no 19 GB uncompressed copy on disk. ~2-3 min, low RAM.
#
# Outputs (phase2/data/): grid.bed, pull_set.tsv, pull_srx.txt, grid_seen_srx.txt
set -euo pipefail

REPO="$(cd "$(dirname "$0")/../.." && pwd)"
HERE="$(cd "$(dirname "$0")" && pwd)"
SECACTS="$(cd "$HERE/../.." && pwd)"
. "$SECACTS/secacts_env.sh"
need_dataroot
GZ="${GZ:-$SECACTS_DATAROOT/chip-atlas/00.data/Histone/His.ALL.50.H3K27ac.AllCell.bed.gz}"
MAN="$REPO/phase1/data/phase1_manifest.tsv"
OUT="$REPO/phase2/data"
TMP="${TMP:-${TMPDIR:-/tmp}/secacts_grid.$$}"
mkdir -p "$OUT" "$TMP"
trap 'rm -rf "$TMP"' EXIT

[[ -f "$GZ" ]]  || { echo "missing AllCell bed: $GZ"; exit 1; }
[[ -f "$MAN" ]] || { echo "missing manifest: $MAN"; exit 1; }

# 1) analysis set = qc_pass(10) & has_cn(8); emit the SRX key list + the richer pull_set.tsv
awk -F'\t' 'NR>1 && $8==1 && $10==1 {print $1}' "$MAN" | sort -u > "$OUT/pull_srx.txt"
awk -F'\t' 'BEGIN{OFS="\t"} NR==1{print "srx","cell","cvcl","model_id","lineage","subtype"; next}
            $8==1 && $10==1 {print $1,$3,$4,$5,$6,$7}' "$MAN" > "$OUT/pull_set.tsv"
NSRX=$(wc -l < "$OUT/pull_srx.txt")
echo "[20] analysis-set SRX: $NSRX"

# 2) stream the .gz: keep canonical-chrom peaks whose col-4 ID= accession is in our set;
#    emit 3-col intervals (-> sort -> bedtools merge) and record which SRX were actually seen.
echo "[20] streaming peaks -> grid (compressed: ~1 min; raw 100GB+: several min) ..."
# accept either a gzipped (.gz, streamed via pigz) or a plain uncompressed .bed (ChIP-Atlas serves raw)
{ case "$GZ" in *.gz) pigz -dc "$GZ";; *) cat "$GZ";; esac; } \
 | awk -F'\t' -v seen="$OUT/grid_seen_srx.txt" '
     NR==FNR { keep[$1]=1; next }
     {
       if ($1 !~ /^chr([1-9][0-9]?|X|Y)$/) next          # drop chrM / alt / random contigs
       s=$4; sub(/^.*ID=/,"",s); sub(/[^A-Za-z0-9].*$/,"",s)
       if (s in keep) { seenarr[s]=1; print $1"\t"$2"\t"$3 }
     }
     END { for (k in seenarr) print k > seen }
   ' "$OUT/pull_srx.txt" - \
 | sort -k1,1 -k2,2n --parallel=8 -S 4G -T "$TMP" \
 | bedtools merge -i - \
 > "$OUT/grid.bed"

# 3) report
REGIONS=$(wc -l < "$OUT/grid.bed")
SEEN=$(sort -u "$OUT/grid_seen_srx.txt" | wc -l)
BP=$(awk '{s+=$3-$2} END{print s}' "$OUT/grid.bed")
echo "-------------------------------------------------------------"
echo "[20] grid.bed regions : $REGIONS"
printf  "[20] grid bp covered  : %d  (%.2f%% of ~3.1 Gb genome)\n" "$BP" "$(awk -v b=$BP 'BEGIN{print 100*b/3.1e9}')"
echo "[20] SRX with >=1 peak : $SEEN / $NSRX  (missing $((NSRX-SEEN)))"
echo "[20] chrom distribution:"; cut -f1 "$OUT/grid.bed" | sort | uniq -c | sort -k2,2V | awk '{printf "        %-6s %s\n",$2,$1}'
echo "[20] outputs -> $OUT/{grid.bed,pull_set.tsv,pull_srx.txt,grid_seen_srx.txt}"