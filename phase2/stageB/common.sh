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
