#!/bin/bash
set -euo pipefail
source_dir=${1:?source directory required}
shift
export JULIA_DEPOT_PATH=/home1/mkolar/juliaPkg/depot:
export JULIA_NUM_THREADS=1 JULIA_NUM_PRECOMPILE_TASKS=1 JULIA_PKG_PRECOMPILE_AUTO=0
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
cd "$source_dir"
exec /home1/mkolar/juliaPkg/julia-1.13.2/bin/julia --project "$@"
