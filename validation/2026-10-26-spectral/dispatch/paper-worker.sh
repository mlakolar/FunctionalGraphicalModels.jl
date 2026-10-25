#!/bin/bash
set -euo pipefail
root=/project2/mkolar_1314/FunctionalGraphicalModels/20261025-direct
export DFGM_SPECTRAL_SOURCE=$root/spectral-20261026-v5-source
export JULIA_DEPOT_PATH=/home1/mkolar/juliaPkg/depot:
export JULIA_NUM_THREADS=4 JULIA_NUM_GC_THREADS=1 JULIA_PKG_PRECOMPILE_AUTO=0 JULIA_HEAP_SIZE_HINT=4000M
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
exec /home1/mkolar/juliaPkg/julia-1.13.2/bin/julia --project="$root/common-fast-v2-source" "$root/spectral-dispatch-v1/spectral_batch.jl" "$1" "$root/outputs/spectral-v5" 4
