#!/bin/bash
#SBATCH --partition=small
#SBATCH --job-name=wgs-cn
#SBATCH --time=03:00:00
#SBATCH --mem=12G
#SBATCH --cpus-per-task=4
#SBATCH --output=logs/%x-%A-%a.out
#SBATCH --error=logs/%x-%A-%a.err
# Measured copy number from subsampled open WGS (v3.2; phase1/data/cn_gap_sources_2026-10.md "cheapest route" 3):
# the first SPOTS spots of one run, one read per spot, Bowtie2 to GRCh38 (the Stage B binary and index), reads per
# 50 kb bin (phase2/wgs_bin.py). Ratios to a normal reference and the calibration against DepMap WGS are local
# (phase1/scripts/27_wgs_cn.py). Manifest columns: name key role run mode spots (27_wgs_cn.py manifest).
#   mode pe: read 1 · se: the read · 10x: read 2 (read 1 starts with the 16 bp barcode + 7 bp) ·
#   ont: long reads cut into 100 bp pieces, aligned --local
# Submit from $SB (score_v3/), one task per manifest row:
#   sbatch -A <project> --array=0-<n-1> --export=ALL,PROJ=<project>,MANIFEST=<tsv> wgs_cn.slurm
source /appl/profile/zz-csc-env.sh 2>/dev/null || true
set -euo pipefail
: "${PROJ:?}" "${MANIFEST:?}"
R=/scratch/$PROJ/se-cacts; S=$R/stageB; SB=${SB:-$R/score_v3}
BT2="$S/bin/bowtie2-2.2.2/bowtie2"; IDX="$S/ref/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index"
SRA="apptainer exec -B /scratch $S/sif/sra-tools.sif"
PY=/projappl/$PROJ/secacts_venv/bin/python; REPO=$SB/dataroot/002.AI_projects/SE-CaCTS
OUT=$SB/wgs_cn/bins; mkdir -p "$OUT"
IFS=$'\t' read -r name key role run mode spots < <(awk -F'\t' -v i=$((SLURM_ARRAY_TASK_ID + 2)) 'NR==i' "$MANIFEST")
[ -s "$OUT/$name.npz" ] && { echo "[wgs] $name done, skipped"; exit 0; }
T=$(mktemp -d -p "$R/tmp" "wgs.$name.XXXX" 2>/dev/null || mktemp -d -p /scratch/$PROJ "wgs.$name.XXXX")
trap 'rm -rf "$T"' EXIT
echo "[wgs] $name ($key, $role) $run mode=$mode spots=$spots on $(hostname) $(date -Is)"
for k in 1 2 3; do
  $SRA fastq-dump -X "$spots" --split-files --skip-technical -O "$T" "$run" && break
  echo "[wgs] fastq-dump attempt $k failed"; rm -f "$T"/*.fastq; sleep $((60 * k))
done
ls -la "$T"
case "$mode" in
  pe)  fq="$T/${run}_1.fastq" ;;
  10x) fq="$T/${run}_2.fastq" ;;
  se)  fq=$(ls "$T"/${run}*.fastq | head -1) ;;
  ont) src=$(ls "$T"/${run}*.fastq | head -1); fq="$T/chunks.fastq"
       awk 'NR%4==1{n=substr($1,2)} NR%4==2{s=$0} NR%4==0{for(i=1;i+99<=length(s);i+=100)
            printf "@%s_%d\n%s\n+\n%s\n", n, i, substr(s,i,100), substr($0,i,100)}' "$src" > "$fq"; rm -f "$src" ;;
  *) echo "[wgs] unknown mode $mode"; exit 1 ;;
esac
[ -s "$fq" ] || { echo "[wgs] no reads in $fq"; exit 1; }
echo "[wgs] reads: $(( $(wc -l < "$fq") / 4 ))"
loc=""; [ "$mode" = ont ] && loc="--local"
"$BT2" -p "${SLURM_CPUS_PER_TASK:-4}" $loc --no-unal -x "$IDX" -U "$fq" 2> "$OUT/$name.bowtie2.log" \
  | "$PY" "$REPO/phase2/wgs_bin.py" --out "$T/$name.npz"
mv "$T/$name.npz" "$OUT/$name.npz"
tail -n 6 "$OUT/$name.bowtie2.log"
echo "[wgs] $name done $(date -Is)"
