#!/bin/bash
# NOTE: --account is NOT hardcoded. secacts_env.sh exports SBATCH_ACCOUNT from
#       SECACTS_CSC_PROJECT (.env), which sbatch honours. Override per-run with
#       `sbatch --account=...`.
#SBATCH --partition=small            # verify with `sinfo -s` on Roihu; single-node CPU partition
#SBATCH --job-name=secacts-pilot
#SBATCH --time=01:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=8
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err
# SE-CaCTS Step-B pilot compute on Roihu: quantify (deepTools) + PCA (numpy).
# Bigwigs are expected already staged under $WORK/data/bigwigs (rsynced from the laptop).

source /appl/profile/zz-csc-env.sh    # before set -u (the CSC env script isn't -u clean)
set -uo pipefail
# Machine-local config: walk up to the repo root and source .env via secacts_env.sh.
# Depth-independent so this works from a staged copy on the cluster too; if the file is
# absent (partial stage), a real SECACTS_CSC_PROJECT in the environment still satisfies the
# `:?` checks below.
_d="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
while [ "$_d" != "/" ] && [ ! -f "$_d/secacts_env.sh" ]; do _d="$(dirname "$_d")"; done
[ -f "$_d/secacts_env.sh" ] && . "$_d/secacts_env.sh"

PROJ="${SECACTS_CSC_PROJECT:?set SECACTS_CSC_PROJECT in .env (sbatch --export carries it to the job)}"
WORK=/scratch/$PROJ/se-cacts/pilot
TYKKY=/projappl/$PROJ/secacts_tykky

# Toolchain: prefer a maintained module, else the tykky container built by 00_build_env.sh.
if module -t spider 2>/dev/null | grep -qi '^deeptools'; then
  module load deeptools
  module load python-data 2>/dev/null || true   # numpy/pandas/scipy/matplotlib for step 05
else
  export PATH=$TYKKY/bin:$PATH
fi

cd "$WORK/scripts"
echo "[pilot] python: $(command -v python); multiBigwigSummary: $(command -v multiBigwigSummary)"

python -u 04_quantify.py --mode bins --bin-size 10000 \
       --bw-dir "$WORK/data/bigwigs" --selection "$WORK/data/selection.tsv" \
       --out-prefix "$WORK/results/signal" --threads "$SLURM_CPUS_PER_TASK"

python -u 05_normalize_pca.py \
       --signal "$WORK/results/signal.tab" --selection "$WORK/data/selection.tsv" \
       --outdir "$WORK/results"

echo "[pilot] done -> $WORK/results/{pca_batch_vs_lineage.png,pca_summary.tsv}"
