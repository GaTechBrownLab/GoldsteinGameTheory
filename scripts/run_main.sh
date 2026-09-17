#!/usr/bin/env bash
# =============================================================================
# run_main.sh — main experiments: 2 fitness models x 4 conditions x 4 replicates
# -> results/{minimal,acute}/ (Fig 2, Fig 3, Fig 4A-B, S1, S3-S6, S10-S11)
#
# 100K recorded substitutions behind a 10K burn-in, one row every 10
# substitutions. Seeds are deterministic (seed_base + rep - 1).
#
#   scripts/run_main.sh                    # everything
#   REPS=8 scripts/run_main.sh             # more replicates
#   JOBS=6 scripts/run_main.sh             # fewer parallel jobs
#   MODELS=minimal scripts/run_main.sh     # one model only
#
# Progress:  ls results/{minimal,acute}/*/simulation.csv | wc -l
#            grep -l 'Done!' logs/main/*.log | wc -l
# ER/ER runs are the slowest (days at sigma = 0.01); ET/ET finishes in minutes.
# =============================================================================
set -uo pipefail
cd "$(dirname "$0")/.."

OUT=${OUT:-results}
REPS=${REPS:-4}
GENS=${GENS:-100000}
WRITE=${WRITE:-10}
JOBS=${JOBS:-12}
PY=${PYTHON:-python3}
MODELS=${MODELS:-"minimal acute"}
# Slowest first, so the long ER/ER runs start immediately
CONDS=${CONDS:-"ERhost_ERpath ERhost_ETpath EThost_ERpath EThost_ETpath"}

mkdir -p logs/main
CMDS=$(mktemp)
for c in $CONDS; do
  for m in $MODELS; do
    for r in $(seq 1 "$REPS"); do
      echo "$PY run_experiments.py -f $m -c $c --diploid --step-sizes 0.01 --rep $r --write-every $WRITE --gens $GENS -o $OUT > logs/main/${m}_${c}_rep${r}.log 2>&1" >> "$CMDS"
    done
  done
done

echo "$(date '+%F %T')  launching $(wc -l < "$CMDS" | tr -d ' ') runs, $JOBS at a time -> $OUT"
# NUL-separated whole lines: BSD xargs -I caps a substituted command at 255 bytes
tr '\n' '\0' < "$CMDS" | xargs -0 -n 1 -P "$JOBS" bash -c
echo "$(date '+%F %T')  all runs finished"
rm -f "$CMDS"
