#!/bin/bash
#SBATCH --partition=small
#SBATCH --job-name=sB-setup
#SBATCH --time=03:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/%x-%j.out
# Exact v1 tools: Bowtie2 2.2.2 (official binary), samtools 0.1.19, bedtools 2.17.0, MACS2 2.1.1 (the image
# sraTailor.sh names), UCSC bedGraphToBigWig/bedClip; NCBI GRCh38 no_alt analysis set with its prebuilt
# Bowtie2 index (no 16 GB index build).
set -euo pipefail
source ./common.sh
mkdir -p "$S"/{bin,sif,ref}
cd "$S/bin"
[ -x "$BT2" ] || { curl -sfL -o bt2.zip "https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.2/bowtie2-2.2.2-linux-x86_64.zip/download"; unzip -q bt2.zip; rm bt2.zip; }
"$BT2" --version | head -1
cd "$S/sif"
pull() { [ -s "$1" ] || apptainer -q pull "$1" "$2"; }
pull samtools_0.1.19.sif docker://quay.io/biocontainers/samtools:0.1.19--h94a8ba4_6
pull bedtools_2.17.0.sif docker://quay.io/biocontainers/bedtools:2.17.0--0
pull macs2_2.1.1.sif docker://quay.io/biocontainers/macs2:2.1.1--r3.2.2_0
pull ucsc_bg2bw.sif docker://quay.io/biocontainers/ucsc-bedgraphtobigwig:357--1
pull ucsc_bedclip.sif docker://quay.io/biocontainers/ucsc-bedclip:357--1
pull sra-tools.sif docker://quay.io/biocontainers/sra-tools:3.2.1--h4304569_1
cd "$S/ref"
N=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids
[ -s "$IDX.1.bt2" ] || { curl -sfL "$N/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index.tar.gz" | tar xz; }
[ -s "$SIZES" ] || { curl -sfL "$N/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai" | cut -f1,2 > "$SIZES"; }
ls -la "$S/ref" | head; wc -l "$SIZES"
($SAMT 2>&1 || true) | grep -m1 Version; $BEDT --version; $MACS2 --version 2>&1 | tail -1
