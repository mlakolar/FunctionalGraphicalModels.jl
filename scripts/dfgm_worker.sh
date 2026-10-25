#!/bin/bash
set -euo pipefail
setting=${1:?}; replicate=${2:?}; mode=${3:?}; n=${4:?}; p=${5:?}
integration_points=${6:?}; nlambda=${7:?}; work_points=${8:?}
out="$DFGM_OUTPUT_DIR/setting-$setting-rep-$replicate"
mkdir -p "$out"
printf 'job=%s step=%s\n' "$SLURM_JOB_ID" "$SLURM_STEP_ID" >"$out/step.txt"
extra=()
if [[ ${DFGM_PATH_REUSE:-0} == 1 ]]; then
    extra+=(--path-reuse 1)
fi
cache_folder="${DFGM_CACHE_ROOT:-}/setting-$setting-rep-$replicate"
if [[ -n ${DFGM_CACHE_ROOT:-} && -f "$cache_folder/covariances.bin" ]]; then
    extra+=(--covariance-cache "$cache_folder")
fi
if [[ -n ${DFGM_BIAS_CONSTANTS:-} ]]; then
    extra+=(--bias-constants "$DFGM_BIAS_CONSTANTS")
fi
exec /home1/mkolar/juliaPkg/julia-1.13.2/bin/julia --project="$DFGM_SOURCE_DIR" \
    "$DFGM_SOURCE_DIR/experiments/dfgm2020.jl" --setting "$setting" --replicate "$replicate" \
    --mode "$mode" --n "$n" --p "$p" --integration-points "$integration_points" \
    --nlambda "$nlambda" --work-points "$work_points" --bias-constant "${DFGM_BIAS_CONSTANT:-1}" \
    --lp-solver "${DFGM_LP_SOLVER:-highs}" --ntasks "${DFGM_THREADS_PER_WORKER:-1}" \
    "${extra[@]}" --out "$out" >"$out/worker.out" 2>&1
