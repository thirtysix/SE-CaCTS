# Stage B: ChIP-Atlas v1, emulated (sraTailor.sh, FINDINGS 2026-09-28 §13, §25). Sourced by every job.
# Needs PROJ in the environment (sbatch --export=ALL,PROJ=...). No project id is hard-coded.
: "${PROJ:?set PROJ}"
S="${S:-/scratch/$PROJ/se-cacts/stageB}"            # this stage's tree
W="${W:-/scratch/$PROJ/se-cacts/phase2}"             # the atlas pull (grids, out/, bed20)
BT2="$S/bin/bowtie2-2.2.2/bowtie2"
IDX="$S/ref/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index"
SIZES="$S/ref/hg38.chrom.sizes"
AX="apptainer exec -B /scratch"
SAMT="$AX $S/sif/samtools_0.1.19.sif samtools"
BEDT="$AX $S/sif/bedtools_2.17.0.sif bedtools"
MACS2="$AX $S/sif/macs2_2.1.1.sif macs2"
BG2BW="$AX $S/sif/ucsc_bg2bw.sif bedGraphToBigWig"
BEDCLIP="$AX $S/sif/ucsc_bedclip.sif bedClip"
PY="/projappl/$PROJ/secacts_venv/bin/python"
CNROSE="$W/cnrose"
MANIFEST="${MANIFEST:-$S/manifest.tsv}"             # srx, layout, fastq (;-separated ENA ftp paths), mode
mkdir -p "$S"/{logs,work,validate}
row() { awk -F'\t' -v i=$(( $1 + 2 )) 'NR==i' "$MANIFEST"; }

# The three steps, one manifest row each ($1 = row index), always run through `step`: they `exit` on
# completion or failure. The one-step arrays (fetch/align/post.slurm) and the worker pool (worker.slurm)
# share them, so production runs exactly the code that was validated.
#
# `step` runs one in a fresh bash, so `set -e` holds inside it. Called from `if` or the left of `||`, a
# function or ( subshell ) runs with errexit silently off and a failed command is followed by the next one:
# that is how a fetch killed mid-download was still aligned as single-end on 2026-09-29 (GOTCHAS 112, 115).
step() { bash -c 'set -euo pipefail; source "$0"; "$@"' "$S/common.sh" "$@"; }

# Per experiment: runs concatenated per mate (sraTailor.sh concatenates runs per SRX). Download-bound.
fetch_one() {
  r=$(row "$1"); srx=$(cut -f1 <<<"$r"); layout=$(cut -f2 <<<"$r"); ftp=$(cut -f3 <<<"$r")
  runs=$(cut -f9 <<<"$r")
  d="$S/work/$srx"; mkdir -p "$d/runs"
  [ -f "$d/.fetched" ] && { echo "[fetch] $srx done"; exit 0; }
  t0=$(date +%s)
  # SRA first: fasterq-dump reproduces the read names ChIP-Atlas v1 aligned, which Bowtie2 uses to break ties
  # between equally good placements, so the output matches v1 exactly (validation 2026-09-29); ENA FASTQs carry
  # other names and differ slightly in repeats. --check-rs no: do not pre-download the references of
  # reference-compressed runs (84 chromosomes, one timed out); fasterq-dump resolves them as needed.
  # FETCH_THREADS > 1 is safe only if its output is byte-identical to --threads 1 (dumptest.slurm).
  fetch_sra() {
    local SRA="$AX $S/sif/sra-tools.sif" run
    rm -rf "$d/runs"; mkdir -p "$d/runs"
    for run in ${runs//;/ }; do
      # 3 tries: a transient "cannot resolve remote location" (2026-09-29) otherwise sends the row to inexact ENA
      for a in 1 2 3; do
        $SRA prefetch --max-size 500G --check-rs no -O "$d/runs" "$run" && break
        [ "$a" = 3 ] && return 1; sleep 60
      done
      $SRA fasterq-dump --split-files --threads "${FETCH_THREADS:-1}" -t "$d/runs" -O "$d/runs" "$d/runs/$run/$run.sra" || return 1
      rm -rf "$d/runs/$run"
    done
    # one run: move, not copy, so the FASTQ is never on disk twice (same bytes as the cat)
    mates() { local f=( $(ls "$@" | sort) ); if [ ${#f[@]} -eq 1 ]; then mv "${f[0]}" "$OUT"; else cat "${f[@]}" > "$OUT"; fi; }
    if [ "$layout" = PAIRED ]; then
      OUT="$d/r1.fq" mates "$d"/runs/*_1.fastq && OUT="$d/r2.fq" mates "$d"/runs/*_2.fastq || return 1
    else
      OUT="$d/r1.fq" mates "$d"/runs/*.fastq || return 1
    fi
    rm -rf "$d/runs"
    [ -s "$d/r1.fq" ] && { [ "$layout" != PAIRED ] || [ -s "$d/r2.fq" ]; } || return 1
    # gzip while the aligner is idle anyway: the reads then hold ~0.3x the disk for the whole alignment, and
    # Bowtie2 reads the same records from .gz (outputs byte-identical, valtest.slurm 2026-09-29)
    pigz -1 -p "${FETCH_THREADS:-1}" "$d"/r?.fq
  }
  if [ -n "$runs" ]; then
    if fetch_sra; then
      echo sra > "$d/route"; echo "fetch_s=$(( $(date +%s) - t0 ))" > "$d/timing.fetch"
      touch "$d/.fetched"; echo "[fetch] $srx via SRA"; ls -la "$d"; exit 0
    fi
    echo "[fetch] $srx: SRA route failed; falling back to ENA"; rm -f "$d"/r1.fq "$d"/r2.fq; rm -rf "$d/runs"; mkdir -p "$d/runs"
    [ -n "$ftp" ] || { echo "[fetch] $srx: no ENA FASTQ either"; exit 1; }
  fi
  for u in ${ftp//;/ }; do
    f=$(basename "$u")
    case "$layout:$f" in PAIRED:*_1.fastq.gz|PAIRED:*_2.fastq.gz|SINGLE:*) ;; *) continue ;; esac
    [ -s "$d/runs/$f" ] || { curl -sfL --retry 5 --retry-delay 15 -o "$d/runs/$f.part" "https://$u" && mv "$d/runs/$f.part" "$d/runs/$f"; }
  done
  if [ "$layout" = PAIRED ]; then
    cat $(ls "$d"/runs/*_1.fastq.gz | sort) > "$d/r1.fq.gz"; cat $(ls "$d"/runs/*_2.fastq.gz | sort) > "$d/r2.fq.gz"
  else
    cat $(ls "$d"/runs/*.fastq.gz | sort) > "$d/r1.fq.gz"
  fi
  rm -rf "$d/runs"
  # never mark a partial fetch complete: an empty mate file would be aligned as single-end (validation, 2026-09-29)
  [ -s "$d/r1.fq.gz" ] && { [ "$layout" != PAIRED ] || [ -s "$d/r2.fq.gz" ]; } || { echo "[fetch] $srx incomplete"; exit 1; }
  echo ena > "$d/route"; echo "fetch_s=$(( $(date +%s) - t0 ))" > "$d/timing.fetch"
  touch "$d/.fetched"; ls -la "$d"
}

# Bowtie2 2.2.2 defaults + --no-unal (pairs as pairs), samtools 0.1.19 sort: the multi-threaded part.
align_one() {
  r=$(row "$1"); srx=$(cut -f1 <<<"$r"); layout=$(cut -f2 <<<"$r")
  d="$S/work/$srx"; cd "$d"
  local p="${ALIGN_THREADS:-$SLURM_CPUS_PER_TASK}"
  [ -f .aligned ] && { echo "[align] $srx done"; exit 0; }
  [ -f .fetched ] || { echo "[align] $srx not fetched"; exit 1; }
  z=.gz; [ -s r1.fq.gz ] || z=""                        # both routes give .fq.gz; plain .fq from older fetches
  if [ "$layout" = PAIRED ] && [ -s "r2.fq$z" ]; then in=(-1 "r1.fq$z" -2 "r2.fq$z"); else in=(-U "r1.fq$z"); fi
  t0=$(date +%s)
  # -1: fast-compressed, not -u, so the unsorted BAM is ~1/4 the size; sort reads the same records either way
  "$BT2" -p "$p" -t --no-unal -x "$IDX" -q "${in[@]}" 2> bowtieReport.txt | $SAMT view -S -1 - > unsrt.bam
  t1=$(date +%s)
  rm -f .fetched r1.fq.gz r2.fq.gz r1.fq r2.fq          # free the reads before the sort; a rerun refetches
  $SAMT sort -@ "$p" -m 800M unsrt.bam dup
  rm -f unsrt.bam
  awk '/reads; of these/{print $1}' bowtieReport.txt > nspots
  echo "align_s=$((t1-t0)) sort_s=$(( $(date +%s) - t1 ))" > timing.align
  cat bowtieReport.txt timing.align; touch .aligned
}

# The single-threaded tail of sraTailor.sh: rmdup (-s for single-end), RPM bigWig (genomecov -scale 1e6/N
# after dedup) and MACS2 2.1.1 at q=1e-20 run side by side on 2 cores, bedClip -> bed20; then the atlas's own
# cnrose pass (grids + archive + QC), exactly as 40_pull_one.sh does for ChIP-Atlas files.
# mode=validate keeps everything under $VALIDATE_DIR/<SRX> ($S/validate); mode=run publishes to $W/out, bed20.
post_one() {
  r=$(row "$1"); srx=$(cut -f1 <<<"$r"); layout=$(cut -f2 <<<"$r"); mode=$(cut -f4 <<<"$r")
  d="$S/work/$srx"
  V="${VALIDATE_DIR:-$S/validate}"
  if [ "$mode" = validate ]; then [ -s "$V/$srx/$srx.se.bed" ] && { echo "[post] $srx done"; exit 0; }
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
  ns=$(cat nspots); np=$(wc -l < "$srx.20.bed"); route=$(cat route 2>/dev/null || echo unknown)
  printf '{"srx":"%s","layout":"%s","reads":%s,"n_before_rmdup":%s,"n_after_rmdup":%s,"pct_mapped":%.2f,"pct_dup":%.2f,"n_peaks":%s,"fetch_route":"%s","pipeline":"stageB v1-emulation (bowtie2 2.2.2, samtools 0.1.19, bedtools 2.17.0, macs2 2.1.1, NCBI GRCh38 no_alt)"}\n' \
    "$srx" "$layout" "$ns" "$nb" "$na" "$(awk -v a="$nb" -v b="$ns" -v m="$([ "$layout" = PAIRED ] && echo 2 || echo 1)" 'BEGIN{print 100*a/(b*m)}')" "$(awk -v a="$na" -v b="$nb" 'BEGIN{print 100*(b-a)/b}')" "$np" "$route" > "$srx.stageB_qc.json"
  t1=$(date +%s)
  if [ "$mode" = validate ]; then out="$V/$srx"; else out="$W/out/${srx: -2}/.staging.$srx.$SLURM_JOB_ID"; fi
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
