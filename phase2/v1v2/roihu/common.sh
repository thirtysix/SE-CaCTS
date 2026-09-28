# Shared settings for the ChIP-Atlas v1 vs v2 comparison (FINDINGS 2026-09-28 §13). Sourced by every job.
# Needs PROJ in the environment (sbatch --export=ALL,PROJ=...). No project id is hard-coded.
: "${PROJ:?set PROJ}"
V2W="${V2W:-/scratch/$PROJ/se-cacts/v1v2}"            # this comparison's work tree
W="${W:-/scratch/$PROJ/se-cacts/phase2}"               # the atlas pull (v1 SE calls, grid, bed20)
SIF="$V2W/pipeline-v2.sif"
V2IMG="docker://ghcr.io/inutano/chip-atlas-pipeline-v2:v1.1.0"
REPO="$V2W/chip-atlas-pipeline-v2"
REF="$V2W/ref/hg38"                                    # hg38.fa (UCSC, as v2 uses), chrom.sizes, bwa-mem2 index
SAMPLES="$V2W/samples.tsv"
VENV="/projappl/$PROJ/secacts_venv"
mkdir -p "$V2W"/{logs,fastq,v1,v2,cmp,tmp} "$REF"
