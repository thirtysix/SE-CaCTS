#!/usr/bin/env bash
# Roihu — driver run FROM THE LAPTOP: stage the pilot to scratch and submit the SLURM job.
# Honors login-node discipline: bigWigs are downloaded LOCALLY (03_download_bigwigs.sh --go) and
# rsynced up, rather than pulled on the login node. Assumes an `ssh roihu` alias (see README.md).
set -euo pipefail
# Machine-local config: walk up to the repo root and source .env via secacts_env.sh.
# Depth-independent so this works from a staged copy on the cluster too; if the file is
# absent (partial stage), a real SECACTS_CSC_PROJECT in the environment still satisfies the
# `:?` checks below.
_d="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
while [ "$_d" != "/" ] && [ ! -f "$_d/secacts_env.sh" ]; do _d="$(dirname "$_d")"; done
[ -f "$_d/secacts_env.sh" ] && . "$_d/secacts_env.sh"
HOST="${HOST:-roihu}"
PROJ="${PROJ:-${SECACTS_CSC_PROJECT:?set SECACTS_CSC_PROJECT in .env, or pass PROJ=}}"
WORK="/scratch/${PROJ}/se-cacts/pilot"
HERE="$(cd "$(dirname "$0")/.." && pwd)"          # the pilot/ dir

echo "[stage] target: ${HOST}:${WORK}"
ssh "$HOST" "mkdir -p ${WORK}/{scripts,data,results}"

# 1) scripts + the small metadata (manifest/selection); NOT the big local bed or results
rsync -az "$HERE/scripts/" "${HOST}:${WORK}/scripts/"
rsync -az "$HERE/data/manifest.tsv" "$HERE/data/selection.tsv" "${HOST}:${WORK}/data/"

# 2) bigWigs — must exist locally first (bash ../scripts/03_download_bigwigs.sh --go)
if compgen -G "$HERE/data/bigwigs/SRX*.bw" > /dev/null; then
  echo "[stage] rsyncing $(ls "$HERE"/data/bigwigs/SRX*.bw | wc -l) bigWigs..."
  rsync -az --info=progress2 "$HERE/data/bigwigs/" "${HOST}:${WORK}/data/bigwigs/"
else
  echo "[stage] NO local bigWigs found. Run first:  (cd $HERE/scripts && bash 03_download_bigwigs.sh --go)"
  echo "        (or download on the Roihu login node if you prefer — compute nodes have no internet)"
  exit 1
fi

# 3) build the container once (idempotent-ish), then submit
echo "[stage] ensuring tykky env exists (login node build if missing)..."
ssh "$HOST" "source /appl/profile/zz-csc-env.sh; [ -x /projappl/${PROJ}/secacts_tykky/bin/multiBigwigSummary ] || bash ${WORK}/../pilot/roihu/00_build_env.sh || bash -s" < "$HERE/roihu/00_build_env.sh" || true

echo "[submit] sbatch pilot.slurm"
rsync -az "$HERE/roihu/pilot.slurm" "${HOST}:${WORK}/pilot.slurm"
JID=$(ssh "$HOST" "source /appl/profile/zz-csc-env.sh; cd ${WORK}; sbatch --parsable pilot.slurm")
echo "[submit] job $JID — monitor:  ssh $HOST 'squeue -j $JID -o \"%.15i %.8T %.10M %R\"'"
echo "[submit] pull results:  rsync -az ${HOST}:${WORK}/results/ $HERE/results/"
