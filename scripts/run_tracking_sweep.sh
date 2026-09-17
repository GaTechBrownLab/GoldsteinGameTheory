#!/usr/bin/env bash
# =============================================================================
# run_tracking_sweep.sh — tracking-strength sweep: k = 0-4 x 4 conditions x
# 8 replicates -> results/tracking/ (S8, S9)
#
# 10K recorded substitutions behind the FULL 10K burn-in. Checked against 100K
# runs: 10K windows reproduce SD(v) and corr(W_H,W_P) without bias, and with
# 8 reps the standard error is ~10x smaller than the k effect. Burn-in is
# passed explicitly because --gens alone would shrink it to gens/10.
#
#   scripts/run_tracking_sweep.sh                 # k = 0 .. 4, 4 conditions, 8 reps
#   KS="1.5 3" scripts/run_tracking_sweep.sh      # selected k values
#   REPS=4 JOBS=6 scripts/run_tracking_sweep.sh
#
# Resume-safe: runs whose log reports Done! (or that are still running) are skipped.
# =============================================================================
set -uo pipefail
cd "$(dirname "$0")/.."

OUT=${OUT:-results}
REPS=${REPS:-8}
GENS=${GENS:-10000}
BURN=${BURN:-10000}
WRITE=${WRITE:-10}
JOBS=${JOBS:-10}
PY=${PYTHON:-python3}
KS=${KS:-"0 0.5 1 1.5 2 3 4"}
# Slowest condition first so the long runs start immediately
CONDS=${CONDS:-"ERhost_ERpath ERhost_ETpath EThost_ERpath EThost_ETpath"}

mkdir -p logs/tracking
CMDS=$(mktemp)
for c in $CONDS; do
  for k in $KS; do
    for r in $(seq 1 "$REPS"); do
      log="logs/tracking/k${k}_${c}_rep${r}.log"
      # Resume-safe: skip runs that already finished or are still running
      grep -q 'Done!' "$log" 2>/dev/null && continue
      ps -eww -o command | grep -F -- "--k $k -c $c --diploid --step-sizes 0.01 --rep $r --write-every" \
        | grep -v -e grep -e 'bash -c' | grep -q . && continue
      echo "$PY run_experiments.py -f tracking --k $k -c $c --diploid --step-sizes 0.01 --rep $r --write-every $WRITE --gens $GENS --burn-in $BURN -o $OUT > $log 2>&1" >> "$CMDS"
    done
  done
done

echo "$(date '+%F %T')  tracking sweep: $(wc -l < "$CMDS" | tr -d ' ') runs, $JOBS at a time -> $OUT"
# NUL-separated whole lines: BSD xargs -I caps a substituted command at 255 bytes
tr '\n' '\0' < "$CMDS" | xargs -0 -n 1 -P "$JOBS" bash -c
# Wait for runs an earlier invocation left running (skipped above), so
# "finished" only prints once every tracking run is really done
while ps -eww -o command | grep -v grep | grep -q 'run_experiments.py -f tracking'; do sleep 60; done
echo "$(date '+%F %T')  tracking sweep finished"
rm -f "$CMDS"
