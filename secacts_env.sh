# shellcheck shell=bash
# SE-CaCTS — machine-local configuration for the shell drivers.
#
# Source it, don't execute it. Callers already resolve the repo root from $0:
#
#   HERE="$(cd "$(dirname "$0")" && pwd)"
#   SECACTS="$(cd "$HERE/../.." && pwd)"
#   . "$SECACTS/secacts_env.sh"
#
# Exports SECACTS_DATAROOT / SECACTS_CSC_PROJECT / SECACTS_CACHE_DIR / SECACTS_PY from
# `.env` (gitignored), and SBATCH_ACCOUNT so SLURM scripts need no hardcoded
# `#SBATCH --account=`. Values already in the environment always win, so a one-off
# override works: SECACTS_DATAROOT=/data/mirror ./phase2/scripts/30_rehearse_local.sh
#
# Then call `need_dataroot` in any script that actually reads the shared datasets — it
# fails loudly with the fix rather than expanding to an empty string and building paths
# like "/DepMap/..." that fail confusingly ten lines later.

_secacts_env_root="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
SECACTS_ENV_FILE="${SECACTS_ENV_FILE:-$_secacts_env_root/.env}"

if [ -f "$SECACTS_ENV_FILE" ]; then
  # Read KEY=VALUE without executing the file, and let the real environment win.
  while IFS= read -r _line || [ -n "$_line" ]; do
    case "$_line" in ''|'#'*) continue ;; esac
    _line="${_line#export }"
    case "$_line" in *=*) ;; *) continue ;; esac
    _k="${_line%%=*}"; _v="${_line#*=}"
    _k="$(printf '%s' "$_k" | tr -d '[:space:]')"
    case "$_k" in ''|*[!A-Za-z0-9_]*) continue ;; esac       # skip anything not a plain var name
    _v="${_v%\"}"; _v="${_v#\"}"; _v="${_v%\'}"; _v="${_v#\'}"
    if [ -z "$(eval "printf '%s' \"\${$_k:-}\"")" ]; then
      export "$_k=$(eval "printf '%s' \"$_v\"" 2>/dev/null || printf '%s' "$_v")"
    fi
  done < "$SECACTS_ENV_FILE"
  unset _line _k _v
fi

: "${SECACTS_CACHE_DIR:=$_secacts_env_root/.cache}"
: "${SECACTS_PY:=$HOME/miniconda3/envs/atac_hdac/bin/python}"
export SECACTS_CACHE_DIR SECACTS_PY

# SLURM reads --account from this, so no .slurm file needs the allocation baked in.
[ -n "${SECACTS_CSC_PROJECT:-}" ] && export SBATCH_ACCOUNT="$SECACTS_CSC_PROJECT"

# Default cluster workdir, derived unless explicitly set.
[ -z "${SECACTS_SCRATCH:-}" ] && [ -n "${SECACTS_CSC_PROJECT:-}" ] \
  && export SECACTS_SCRATCH="/scratch/$SECACTS_CSC_PROJECT/se-cacts"

need_dataroot() {
  if [ -z "${SECACTS_DATAROOT:-}" ]; then
    echo "ERROR: SECACTS_DATAROOT is not set." >&2
    echo "  It is the parent directory holding DepMap/, chip-atlas/, cellosaurus/," >&2
    echo "  0.human_genome/ and the sibling pyCaCTS checkout." >&2
    echo "  Fix: cp '$_secacts_env_root/sample.env' '$_secacts_env_root/.env' and set it," >&2
    echo "       or export SECACTS_DATAROOT=... for this run." >&2
    exit 2
  fi
}

need_csc_project() {
  if [ -z "${SECACTS_CSC_PROJECT:-}" ]; then
    echo "ERROR: SECACTS_CSC_PROJECT is not set (needed only for the HPC pull)." >&2
    echo "  Fix: set it in '$_secacts_env_root/.env', or export SECACTS_CSC_PROJECT=..." >&2
    exit 2
  fi
}
