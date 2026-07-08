#!/bin/bash
# ============================================================
#  b3submit.sh -- submit b3bound jobs on Oscar by family + range.
#  Mirrors the loop-and-export pattern in the CCV batch-jobs docs:
#  https://docs.ccv.brown.edu/oscar/submitting-jobs/batch
#  One sbatch job is submitted per value of n.
#
#  USAGE
#    ./b3submit.sh <family> <start> [end] [step]
#    ./b3submit.sh <family> <start-end[:step]>
#    ./b3submit.sh --edges "0-1 1-2 2-0" [label]
#    Anything after '--' is forwarded to b3bound.py.
#
#  EXAMPLES
#    ./b3submit.sh cycle 3 8
#    ./b3submit.sh cycle 3-8
#    ./b3submit.sh sun 3 9 2
#    ./b3submit.sh mobius 4
#    ./b3submit.sh --edges "0-1 1-2 2-0 2-3" tri
#    ./b3submit.sh cycle 3 8 -- --color-cap 5000000 --verbose
#
#  Oscar resource overrides (env vars, per the CCV docs' sbatch flags):
#    TIME=12:00:00 MEM=8G CORES=1 PART=batch \
#      MODULES='module load python/3.11.0s' ./b3submit.sh sun 3 10
# ============================================================
set -euo pipefail

JOB="b3bound.sh"
[[ -f "$JOB" ]] || { echo "cannot find $JOB in $PWD" >&2; exit 1; }

# ---- resource overrides forwarded to sbatch (all optional) -------
SB=()
[[ -n "${TIME:-}"  ]] && SB+=(-t "$TIME")
[[ -n "${MEM:-}"   ]] && SB+=(--mem "$MEM")
[[ -n "${CORES:-}" ]] && SB+=(-c "$CORES")
[[ -n "${PART:-}"  ]] && SB+=(-p "$PART")
export PYTHON="${PYTHON:-python3}"
export MODULES="${MODULES:-}"          # e.g. 'module load python/3.11.0s'

usage() { sed -n '2,30p' "$0"; exit 1; }
[[ $# -lt 1 ]] && usage

# ---- split pass-through flags after '--' -------------------------
ARGS=(); EXTRA=(); PASS=false
for a in "$@"; do
    if $PASS; then EXTRA+=("$a")
    elif [[ "$a" == "--" ]]; then PASS=true
    else ARGS+=("$a"); fi
done
set -- "${ARGS[@]}"
EXTRA_STR="${EXTRA[*]:-}"

mkdir -p logs results

submit_one() { # args: NAME=VALUE pairs to export into the job environment
    # Export vars into the environment, then let the job inherit them
    # (--export=ALL). This avoids comma/space parsing issues in --export values,
    # per the loop-and-export pattern in the CCV docs.
    export PYTHON MODULES
    export B3_EXTRA="$EXTRA_STR"
    local kv
    for kv in "$@"; do export "${kv?}"; done
    sbatch "${SB[@]}" --export=ALL "$JOB"
}

if [[ "${1:-}" == "--edges" ]]; then
    EDGES="${2:?--edges needs an edge list}"; LABEL="${3:-custom}"
    echo "submitting custom graph '$LABEL'"
    submit_one "B3_KIND=EDGES" "B3_LABEL=${LABEL}" "B3_EDGES=${EDGES}"
    exit 0
fi

FAMILY="${1:?need a graph family}"
SPEC="${2:?need a range: start [end] [step]  or  start-end[:step]}"
STEP="${4:-1}"

if [[ "$SPEC" == *-* ]]; then
    BASE="${SPEC%%:*}"; [[ "$SPEC" == *:* ]] && STEP="${SPEC##*:}"
    START="${BASE%%-*}"; END="${BASE##*-}"
else
    START="$SPEC"; END="${3:-$START}"
fi

echo "submitting $FAMILY for n = $START..$END step $STEP"
for ((n=START; n<=END; n+=STEP)); do
    submit_one "B3_KIND=FAM" "B3_FAMILY=${FAMILY}" "B3_N=${n}"
done