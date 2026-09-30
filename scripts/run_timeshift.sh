#!/usr/bin/env bash
# =============================================================================
# run_timeshift.sh — time-shift (cross-inoculation) summaries for Fig 5,
# written to data/timeshift/. Run after scripts/run_main.sh and
# scripts/run_zoom.sh have finished.
#
#   scripts/run_timeshift.sh main   # main runs, results/ (panels D-E use delta = 0)
#   scripts/run_timeshift.sh zoom   # write_every = 1 zoom runs, data/zoom/ (panels A-C)
#
# The main grid is in recorded steps, which are 10 substitutions apart.
# 'flow-contemporary' is the assay counterpart of the simulation's anchored
# rule, so the sympatric delta = 0 cell reproduces the recorded state; 'flow'
# starts each pairing from the uninduced levels, as a fresh infection would.
# =============================================================================
set -euo pipefail
cd "$(dirname "$0")/.."
OUT=data/timeshift
mkdir -p "$OUT"
PY=${PYTHON:-python3}

case "${1:-}" in
  main)
    "$PY" src/timeshift.py --roots data/main -f minimal --deltas="-20:20:1" \
        --axis sub --n-focal 120 --allopatric --settle flow flow-contemporary \
        --pairs "$OUT/timeshift_main_pairs.csv" \
        --diagnostics "$OUT/diagnostics_main.csv" \
        -o "$OUT/timeshift_main.csv" ;;
  zoom)
    "$PY" src/timeshift.py --roots data/zoom --deltas="-200:200:5" \
        --axis sub --n-focal 500 --allopatric --settle flow-contemporary --n-block 33 \
        --pairs "$OUT/timeshift_zoom_pairs.csv" \
        --diagnostics "$OUT/diagnostics_zoom.csv" \
        -o "$OUT/timeshift_zoom.csv" ;;
  *) echo "usage: $0 main|zoom" >&2; exit 2 ;;
esac
echo "$(date '+%F %T')  timeshift ${1} finished"
