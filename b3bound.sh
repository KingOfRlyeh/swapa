#!/bin/bash
# ============================================================
#  b3bound.job -- Oscar (Brown CCV) batch job script for one b3bound run.
#  Submitted by b3submit.sh, which passes B3_* via --export.
#  Follows https://docs.ccv.brown.edu/oscar/submitting-jobs/batch
#  (resource flags: -n cores, -c cpus/task, -t HH:MM:SS, --mem per node).
# ============================================================
#SBATCH -J b3bound                 # job name
#SBATCH -n 1                       # one task
#SBATCH -c 1                       # one core per task (b3bound is single-threaded)
#SBATCH -t 04:00:00                # walltime HH:MM:SS
#SBATCH --mem=4G                   # memory per node
#SBATCH -o logs/b3-%j.out          # stdout  (%j = job id)
#SBATCH -e logs/b3-%j.err          # stderr
# #SBATCH -p batch                 # partition (uncomment/override to target one)
# #SBATCH --mail-type=END,FAIL
# #SBATCH --mail-user=you@brown.edu

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
mkdir -p logs results

# --- environment: load a Python with networkx/numpy/sympy ---------
# Oscar example (uncomment and adjust to an available module or your venv):
# module load python/3.11.0s
# source ~/envs/b3bound/bin/activate
# If b3submit.sh was given MODULES='...', run it here:
[[ -n "${MODULES:-}" ]] && eval "$MODULES"
PYTHON="${PYTHON:-python3}"

# --- inputs passed via --export from b3submit.sh ------------------
: "${B3_KIND:?B3_KIND not set}"          # FAM or EDGES
B3_EXTRA="${B3_EXTRA:-}"

if [[ "$B3_KIND" == "EDGES" ]]; then
    : "${B3_LABEL:?}"; : "${B3_EDGES:?}"
    OUT="results/${B3_LABEL}.txt"
    echo "host=$(hostname) job=${SLURM_JOB_ID:-NA} custom=${B3_LABEL} edges='${B3_EDGES}' $(date -Is)"
    "$PYTHON" b3bound.py --edges "$B3_EDGES" $B3_EXTRA | tee "$OUT"
else
    : "${B3_FAMILY:?}"; : "${B3_N:?}"
    OUT="results/${B3_FAMILY}_${B3_N}.txt"
    echo "host=$(hostname) job=${SLURM_JOB_ID:-NA} family=${B3_FAMILY} n=${B3_N} $(date -Is)"
    "$PYTHON" b3bound.py --family "$B3_FAMILY" --range "$B3_N" $B3_EXTRA | tee "$OUT"
fi
echo "done $(date -Is) -> $OUT"
