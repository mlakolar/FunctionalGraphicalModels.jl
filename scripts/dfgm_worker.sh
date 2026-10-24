#!/bin/bash
set -euo pipefail
setting=${1:?}; replicate=${2:?}; mode=${3:?}; n=${4:?}; p=${5:?}
integration_points=${6:?}; nlambda=${7:?}; work_points=${8:?}
out="$DFGM_OUTPUT_DIR/setting-$setting-rep-$replicate"
mkdir -p "$out"
exec /home1/mkolar/juliaPkg/julia-1.13.2/bin/julia --project="$DFGM_SOURCE_DIR" \
    "$DFGM_SOURCE_DIR/experiments/dfgm2020.jl" --setting "$setting" --replicate "$replicate" \
    --mode "$mode" --n "$n" --p "$p" --integration-points "$integration_points" \
    --nlambda "$nlambda" --work-points "$work_points" --out "$out" >"$out/worker.out" 2>&1
