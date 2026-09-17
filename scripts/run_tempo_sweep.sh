#!/usr/bin/env bash
# =============================================================================
# run_tempo_sweep.sh — does mutational tempo control host boundary trapping?
# minimal model, 3 gamma x 3 conditions x 4 replicates -> results_gamma/
# (Fig 4D-E, S7)
#
# Equal population sizes (N_H = N_P = 1e4), so the tempo ratio
#     R = (1 - gamma) N_P / (gamma N_H)
# is set by gamma alone: gamma = 1e-4, 1e-2, 0.5  ->  R ~ 1e4, 1e2, 1.
# gamma > 0.5 is unnecessary with equal N: gamma = 0.99 mirrors gamma = 0.01
# with the mixed scenarios swapped. At R = 1 the two mixed scenarios must be
# mirror images with v and c exchanged (the minimal model is symmetric), which
# doubles as a check for asymmetries in the implementation.
#
# Written to results_gamma/ so these runs never pool with the main results
# (the loaders filter on gamma but not on population size).
# 10K recorded substitutions behind a full 10K burn-in, as for the tracking sweep.
#
#   scripts/run_tempo_sweep.sh
#   GAMMAS="0.5" REPS=8 scripts/run_tempo_sweep.sh
# =============================================================================
set -uo pipefail
cd "$(dirname "$0")/.."

OUT=${OUT:-results_gamma}
REPS=${REPS:-4}
GENS=${GENS:-10000}
BURN=${BURN:-10000}
WRITE=${WRITE:-10}
JOBS=${JOBS:-10}
N=${N:-10000}
PY=${PYTHON:-python3}
GAMMAS=${GAMMAS:-"0.0001 0.01 0.5"}
# Slowest first
CONDS=${CONDS:-"ERhost_ERpath ERhost_ETpath EThost_ERpath"}

mkdir -p logs/tempo
CMDS=$(mktemp)
for c in $CONDS; do
  for g in $GAMMAS; do
    for r in $(seq 1 "$REPS"); do
      echo "$PY run_experiments.py -f minimal -c $c --diploid --step-sizes 0.01 --gamma $g --pop-size $N --path-pop-size $N --rep $r --write-every $WRITE --gens $GENS --burn-in $BURN -o $OUT > logs/tempo/g${g}_${c}_rep${r}.log 2>&1" >> "$CMDS"
    done
  done
done

echo "$(date '+%F %T')  tempo sweep: $(wc -l < "$CMDS" | tr -d ' ') runs, $JOBS at a time -> $OUT"
# One whole line per argument, NUL-separated: BSD xargs -I caps the substituted
# command at 255 bytes, which these commands exceed (that silently ran 0 runs)
tr '\n' '\0' < "$CMDS" | xargs -0 -n 1 -P "$JOBS" bash -c
n_done=$(grep -l 'Done!' logs/tempo/g*_rep*.log 2>/dev/null | wc -l | tr -d ' ')
echo "$(date '+%F %T')  tempo sweep finished: $n_done of $(wc -l < "$CMDS" | tr -d ' ') runs completed"
rm -f "$CMDS"
