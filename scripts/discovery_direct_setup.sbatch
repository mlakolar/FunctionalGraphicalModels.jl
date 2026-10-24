#!/bin/bash
#SBATCH --account=mkolar_1314
#SBATCH --partition=main
#SBATCH --constraint=epyc-7513
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=00:20:00
set -euo pipefail
export JULIA_DEPOT_PATH=/home1/mkolar/juliaPkg/depot:
export JULIA_NUM_THREADS=1 JULIA_NUM_PRECOMPILE_TASKS=4 JULIA_PKG_PRECOMPILE_AUTO=0
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
cd "${1:?source directory required}"
JULIA=/home1/mkolar/juliaPkg/julia-1.13.2/bin/julia
"$JULIA" scripts/instantiate.jl
"$JULIA" --project -e 'using Pkg; Pkg.precompile()'
"$JULIA" --project -e 'using Pkg; Pkg.test(;allow_reresolve=false)'
