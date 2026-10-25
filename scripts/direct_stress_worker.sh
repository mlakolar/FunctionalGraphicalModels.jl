#!/bin/bash
set -euo pipefail
out="$DFGM_OUTPUT_DIR/stress-${1:?configuration required}"
mkdir -p "$out"
printf 'job=%s step=%s\n' "$SLURM_JOB_ID" "$SLURM_STEP_ID" >"$out/step.txt"
exec /home1/mkolar/juliaPkg/julia-1.13.2/bin/julia --project="$DFGM_SOURCE_DIR" \
  "$DFGM_SOURCE_DIR/experiments/direct_stress.jl" --setting "$1" --replicates "${DFGM_STRESS_REPLICATES:-3}" \
  --out "$out" >"$out/worker.out" 2>&1
