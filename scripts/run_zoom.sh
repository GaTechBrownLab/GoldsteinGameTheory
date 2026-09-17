#!/usr/bin/env bash
# =============================================================================
# run_zoom.sh — "zoom" runs that record every substitution (write_every = 1):
# minimal model, 4 conditions x 8 replicates -> results_zoom/ (Fig 4C, S2, and
# the Fig 5 time-shift assay, which needs >1 lineage for allopatric controls)
#
# sigma 0.01, diploid, 10 000 substitutions, burn-in 1 000 (gens/10).
# Seeds are deterministic, seed_base + (rep - 1), so runs are reproducible.
#
#   scripts/run_zoom.sh                 # minimal, reps 1-8, 8 at a time
#   scripts/run_zoom.sh -m acute        # acute instead
#   scripts/run_zoom.sh -m "minimal acute" -r 8 -j 6
#   scripts/run_zoom.sh -n              # dry run: list what would be done
#
# Safe to re-run: any replicate whose simulation.csv is already complete is
# skipped, so an interrupted batch resumes where it stopped.
# =============================================================================
set -uo pipefail
cd "$(dirname "$0")/.."

MODELS="minimal"
MAX_REP=8
JOBS=8
DRY=0
GENS=10000
OUT="results_zoom"
PY="${PYTHON:-python3}"

while getopts "m:r:j:g:o:np:h" opt; do
  case $opt in
    m) MODELS="$OPTARG" ;;
    r) MAX_REP="$OPTARG" ;;
    j) JOBS="$OPTARG" ;;
    g) GENS="$OPTARG" ;;      # shorter runs, for smoke-testing the batch
    o) OUT="$OPTARG" ;;       # output root (default results_zoom)
    n) DRY=1 ;;
    p) PY="$OPTARG" ;;
    h) sed -n '2,17p' "$0"; exit 0 ;;
    *) exit 2 ;;
  esac
done

CONDITIONS="EThost_ETpath EThost_ERpath ERhost_ETpath ERhost_ERpath"
EXPECTED_ROWS=$((GENS * 2 + 1))   # pre+post per generation, plus header

LOGDIR="logs/zoom"
mkdir -p "$LOGDIR"

# Build the job list, skipping anything already complete.
JOBFILE="$(mktemp)"
trap 'rm -f "$JOBFILE"' EXIT
skipped=0
for model in $MODELS; do
  for cond in $CONDITIONS; do
    for rep in $(seq 1 "$MAX_REP"); do
      out="$OUT"
      dir="${out}/${model}/${cond}_sigma0.01_diploid_rep${rep}_zoom"
      csv="${dir}/simulation.csv"
      if [[ -f "$csv" ]] && [[ "$(wc -l < "$csv")" -eq "$EXPECTED_ROWS" ]]; then
        skipped=$((skipped + 1)); continue
      fi
      printf '%s\t%s\t%s\t%s\n' "$model" "$cond" "$rep" "$out" >> "$JOBFILE"
    done
  done
done

total=$(wc -l < "$JOBFILE" | tr -d ' ')
echo "models:      $MODELS"
echo "conditions:  $CONDITIONS"
echo "reps:        1..$MAX_REP   (already complete: $skipped)"
echo "to run:      $total        (parallel: $JOBS)"
echo "logs:        $LOGDIR/"
if [[ "$total" -eq 0 ]]; then echo "Nothing to do."; exit 0; fi
echo "est. wall time: ~$(( (total + JOBS - 1) / JOBS )) waves x ~50-70 min"
echo

if [[ "$DRY" -eq 1 ]]; then
  echo "--- dry run, nothing launched ---"
  awk -F'\t' '{printf "  %-9s %-14s rep%s -> %s\n", $1, $2, $3, $4}' "$JOBFILE"
  exit 0
fi

start=$(date +%s)

# Process pool. Not xargs -P: BSD xargs on macOS cannot assemble a -I{}
# command this long and aborts the whole batch. bash here is 3.2, which has no
# `wait -n`, so the pool is drained by polling the live PIDs.
run_one() {
  local model="$1" cond="$2" rep="$3" out="$4"
  local log="${LOGDIR}/${model}_${cond}_rep${rep}.log"
  if "$PY" run_experiments.py -f "$model" -c "$cond" \
        --step-sizes 0.01 --diploid --gens "$GENS" --write-every 1 \
        --rep "$rep" --tag zoom -o "$out" > "$log" 2>&1; then
    echo "  [$(date +%H:%M)] done   ${model}/${cond} rep${rep}"
  else
    echo "  [$(date +%H:%M)] FAILED ${model}/${cond} rep${rep}  (see $log)"
  fi
}

pids=""
launched=0
while IFS=$'\t' read -r model cond rep out; do
  [ -z "$model" ] && continue
  run_one "$model" "$cond" "$rep" "$out" &
  pids="$pids $!"
  launched=$((launched + 1))
  # Block until a slot frees up.
  while :; do
    alive=""
    for p in $pids; do
      if kill -0 "$p" 2>/dev/null; then alive="$alive $p"; fi
    done
    pids="$alive"
    n=$(echo $pids | wc -w | tr -d ' ')
    [ "$n" -lt "$JOBS" ] && break
    sleep 5
  done
done < "$JOBFILE"
wait

end=$(date +%s)
echo
echo "batch finished in $(( (end - start) / 60 )) min"
echo "Now regenerate the zoom time-shift summary: scripts/run_timeshift.sh zoom"
