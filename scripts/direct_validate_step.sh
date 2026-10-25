#!/bin/bash
# One exclusive CPU and 4 GB inside an already authorized allocation.
set -euo pipefail
export JULIA_DEPOT_PATH=/home1/mkolar/juliaPkg/depot:
export JULIA_NUM_THREADS=1 JULIA_NUM_PRECOMPILE_TASKS=1 JULIA_PKG_PRECOMPILE_AUTO=0
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
cd "${1:?source required}"
JULIA=/home1/mkolar/juliaPkg/julia-1.13.2/bin/julia
echo "validation job=$SLURM_JOB_ID step=$SLURM_STEP_ID"
"$JULIA" --project -e 'using Pkg; Pkg.resolve(); Pkg.instantiate(); Pkg.precompile(); Pkg.test(;allow_reresolve=false)'
"$JULIA" --project experiments/benchmark_direct_lp.jl --out "${2:?output required}/backend"
"$JULIA" --project experiments/benchmark_dfgm_common.jl --out "$2/covariance"
