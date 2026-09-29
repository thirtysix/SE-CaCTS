#!/bin/bash
#SBATCH --partition=small
#SBATCH --job-name=sB-post
#SBATCH --time=12:00:00
#SBATCH --mem=6G
#SBATCH --cpus-per-task=2
#SBATCH --output=logs/%x-%A_%a.out
# The single-threaded tail of sraTailor.sh: rmdup (-s for single-end), RPM bigWig (genomecov -scale 1e6/N
# after dedup) and MACS2 2.1.1 at q=1e-20 run side by side on 2 cores, bedClip -> bed20; then the atlas's own
# cnrose pass (grids + archive + QC), exactly as 40_pull_one.sh does for ChIP-Atlas files.
# mode=validate keeps everything under $S/validate/<SRX>; mode=run publishes to the pull tree ($W/out, bed20).
set -euo pipefail
source ./common.sh
# Roihu caps a user at 200 queued / 100 running jobs on `small`, so each array task takes SAMPLES_PER_TASK
# consecutive manifest rows, one after another; a failed row never stops the rest of the task.
process() {
  r=$(row "$1"); srx=$(cut -f1 <<<"$r"); layout=$(cut -f2 <<<"$r"); mode=$(cut -f4 <<<"$r")
  d="$S/work/$srx"
  if [ "$mode" = validate ]; then [ -s "$S/validate/$srx/$srx.se.bed" ] && { echo "[post] $srx done"; exit 0; }
  else [ -f "$W/out/${srx: -2}/$srx.done" ] && { echo "[post] $srx done"; exit 0; }; fi
  [ -f "$d/.aligned" ] || { echo "[post] $srx not aligned"; exit 1; }
  cd "$d"; export TMPDIR="$d"
  t0=$(date +%s)
  opt=""; [ "$layout" = PAIRED ] || opt="-s"
  $SAMT rmdup $opt dup.bam final.bam
  nb=$($SAMT view -c dup.bam); na=$($SAMT view -c final.bam); rm -f dup.bam
  scale=$(awk -v n="$na" 'BEGIN{printf "%.10f", 1000000/n}')
  ( $MACS2 callpeak -t final.bam -f BAM -g hs -n "$srx.20" -q 1e-20 --outdir . > macs.log 2>&1 ) & mp=$!
  $BEDT genomecov -scale "$scale" -ibam final.bam -bg -g "$SIZES" | LC_ALL=C sort -S 1G -k1,1 -k2,2n > "$srx.bg"
  $BG2BW "$srx.bg" "$SIZES" "$srx.bw"; rm -f "$srx.bg"
  wait $mp
  LC_ALL=C sort -k1,1 -k2,2n "${srx}.20_peaks.narrowPeak" > "$srx.20.unclip"
  $BEDCLIP "$srx.20.unclip" "$SIZES" "$srx.20.bed"
  rm -f "$srx.20.unclip" "${srx}.20_peaks.narrowPeak" "${srx}.20_peaks.xls" "${srx}.20_summits.bed" "${srx}.20_model.r" final.bam
  ns=$(cat nspots); np=$(wc -l < "$srx.20.bed")
  printf '{"srx":"%s","layout":"%s","reads":%s,"n_before_rmdup":%s,"n_after_rmdup":%s,"pct_mapped":%.2f,"pct_dup":%.2f,"n_peaks":%s,"pipeline":"stageB v1-emulation (bowtie2 2.2.2, samtools 0.1.19, bedtools 2.17.0, macs2 2.1.1, NCBI GRCh38 no_alt)"}\n' \
    "$srx" "$layout" "$ns" "$nb" "$na" "$(awk -v a="$nb" -v b="$ns" -v m="$([ "$layout" = PAIRED ] && echo 2 || echo 1)" 'BEGIN{print 100*a/(b*m)}')" "$(awk -v a="$na" -v b="$nb" 'BEGIN{print 100*(b-a)/b}')" "$np" > "$srx.stageB_qc.json"
  t1=$(date +%s)
  if [ "$mode" = validate ]; then out="$S/validate/$srx"; else out="$W/out/${srx: -2}/.staging.$srx.$SLURM_JOB_ID"; fi
  mkdir -p "$out"
  PYTHONPATH="$CNROSE" "$PY" -m cnrose.cli call --bw "$srx.bw" --peaks "$srx.20.bed" \
    --grid "$W/data/grid.20.bed" --grid "$W/data/grid.10.bed" --signal-format f32 --bins 100 --qc --out "$out/$srx"
  cp "$srx.stageB_qc.json" "$out/"
  echo "post_s=$((t1-t0)) cnrose_s=$(( $(date +%s) - t1 ))" | tee timing.post
  if [ "$mode" = validate ]; then
    cp "$srx.bw" "$srx.20.bed" timing.align timing.post bowtieReport.txt "$out/"
  else
    cp "$srx.20.bed" "$W/data/bed20/$srx.20.bed"
    dest="$W/out/${srx: -2}"
    for f in "$out/$srx".*; do mv -f "$f" "$dest/$(basename "$f")"; done
    rmdir "$out"; touch "$dest/$srx.done"
  fi
  cd "$S/work"; rm -rf "$d"
}
K="${SAMPLES_PER_TASK:-1}"; N=$(( $(wc -l < "$MANIFEST") - 1 ))
for (( i = SLURM_ARRAY_TASK_ID * K; i < (SLURM_ARRAY_TASK_ID + 1) * K && i < N; i++ )); do
  ( process "$i" ) || echo "[task $SLURM_ARRAY_TASK_ID] row $i failed"
done
