#!/usr/bin/env bash
# =============================================================================
# run_all.sh — the whole manuscript pipeline, in dependency order.
#
# Defaults: 3 replicates everywhere except the zoom runs, which stay at 8
# because each replicate is a lineage in the time-shift assay's allopatric
# control (8 reps = 64 pairings, 3 reps = 9).
#
# Every run uses the current model defaults, recorded in each config.json:
#   EQ_SELECTION = anchor     realised equilibrium follows the prior phenotype
#   STEP_GRID    = quantile   equal-probability Gaussian mutation steps
#
#   scripts/run_all.sh                 # everything
#   SKIP="main" scripts/run_all.sh     # skip a stage already done
#   REPS=4 scripts/run_all.sh
#
# Progress:  tail -f logs/run_all.log ; ls data/main/*/*/simulation.csv | wc -l
# =============================================================================
set -uo pipefail
cd "$(dirname "$0")/.."

REPS=${REPS:-3}
ZOOM_REPS=${ZOOM_REPS:-8}
SKIP=${SKIP:-""}
stage_wanted() { [[ " $SKIP " != *" $1 "* ]]; }
say() { echo "$(date '+%F %T')  $*"; }

if stage_wanted main; then
  say "stage 1/6: main runs (2 models x 4 scenarios x $REPS reps, 100K substitutions)"
  REPS=$REPS JOBS=8 scripts/run_main.sh
fi

if stage_wanted zoom; then
  say "stage 2/6: zoom runs (4 scenarios x $ZOOM_REPS reps, every substitution)"
  scripts/run_zoom.sh -r "$ZOOM_REPS" -j 8
fi

if stage_wanted tempo; then
  say "stage 3/6: tempo sweep (3 gammas x 3 scenarios x $REPS reps)"
  REPS=$REPS JOBS=${JOBS:-10} scripts/run_tempo_sweep.sh
fi

if stage_wanted tracking; then
  say "stage 4/6: tracking sweep (7 k x 4 scenarios x $REPS reps)"
  REPS=$REPS JOBS=${JOBS:-10} scripts/run_tracking_sweep.sh
fi

if stage_wanted oldrule; then
  say "stage 5/6: old-rule ER/ER arm for the SI comparison (minimal, $REPS reps)"
  CMDS=$(mktemp)
  for r in $(seq 1 "$REPS"); do
    echo "python3 src/run_experiments.py -f minimal -c ERhost_ERpath --diploid --step-sizes 0.01 \
--rep $r --write-every 10 --gens 100000 --burn-in 10000 --eq-selection mutant \
-o data/rule_mutant > logs/main/oldrule_rep${r}.log 2>&1" >> "$CMDS"
  done
  tr '\n' '\0' < "$CMDS" | xargs -0 -n 1 -P "$REPS" bash -c
  rm -f "$CMDS"
fi

if stage_wanted timeshift; then
  say "stage 6/6: time-shift summaries"
  scripts/run_timeshift.sh main
  scripts/run_timeshift.sh zoom
fi

say "pipeline finished — rebuild figures with: Rscript manuscript_figures.R"
