#!/bin/bash
# ============================================================
#  b3bound.sh -- SLURM array job for b3bound.py on Oscar (Brown CCV)
#  https://docs.ccv.brown.edu/oscar/submitting-jobs/array
#
#  The array index IS the parameter n, so submit a family over a range with:
#
#      sbatch --array=<range> b3bound.sh <family> [extra b3bound.py flags]
#
#  RANGE (standard SLURM --array syntax):
#      3-8        n = 3,4,5,6,7,8
#      3-9:2      n = 3,5,7,9        (step 2)
#      4          n = 4              (single)
#      3-6,9,12   comma list + ranges
#
#  EXAMPLES
#      sbatch --array=3-8   b3bound.sh cycle
#      sbatch --array=3-9:2 b3bound.sh sun
#      sbatch --array=4     b3bound.sh mobius
#      sbatch --array=3-8   b3bound.sh cycle --color-cap 5000000 --verbose
#
#  Resource / partition overrides go on the sbatch line, e.g.
#      sbatch --array=3-10 -t 12:00:00 --mem=8G -p batch b3bound.sh sun
# ============================================================

#SBATCH -J b3bound                     # job name
#SBATCH -c 1                           # 1 core (b3bound is single-threaded)
#SBATCH -t 04:00:00                    # walltime HH:MM:SS
#SBATCH --mem=4G                       # memory per node
#SBATCH -o logs/b3-%A_%a.out           # %A = array job id, %a = task id (= n)
#SBATCH -e logs/b3-%A_%a.err
# #SBATCH -p batch                     # partition (or pass -p on the sbatch line)
# #SBATCH --mail-type=END,FAIL
# #SBATCH --mail-user=you@brown.edu
#
# You may also bake the range in instead of passing --array on the command line:
# #SBATCH --array=3-8

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
mkdir -p logs results

# --- Python environment: needs networkx, numpy, sympy --------------
# Uncomment / edit for Oscar (check 'module avail python' for a real version):
# module load python/3.11.0s
# source ~/envs/b3bound/bin/activate
PYTHON="${PYTHON:-python3}"

# --- inputs --------------------------------------------------------
FAMILY="${1:?usage: sbatch --array=<range> b3bound.sh <family> [extra flags]}"
shift || true                          # remaining args pass through to b3bound.py
N="${SLURM_ARRAY_TASK_ID:?run me as an array job: sbatch --array=<range> ...}"

# Batch runs suppress pot output (keeps result files to the bound table).
# Drop --pots, --pot-cap <n>, and --render-pots [prefix] from forwarded flags.
PASS=()
while (( $# )); do
    case "$1" in
        --pots)         shift ;;
        --pot-cap)      shift 2 ;;
        --render-pots)  shift
                        # consume an optional prefix argument (not another flag)
                        if (( $# )) && [[ "$1" != -* ]]; then shift; fi ;;
        *)              PASS+=("$1"); shift ;;
    esac
done

OUT="results/${FAMILY}_${N}.txt"
echo "task n=$N  family=$FAMILY  host=$HOSTNAME  $(date -Is)"

"$PYTHON" b3bound.py --family "$FAMILY" --range "$N" "${PASS[@]}" | tee "$OUT"

echo "done $(date -Is) -> $OUT"