#!/bin/bash
#SBATCH --account=mkolar_1314
#SBATCH --partition=main
#SBATCH --constraint=epyc-7513
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=01:15:00
# Pass source, tasks.tsv and output root. No implicit array or full campaign.
set -euo pipefail
module load parallel/20240522
export JULIA_DEPOT_PATH=/home1/mkolar/juliaPkg/depot:
export DFGM_THREADS_PER_WORKER=${DFGM_THREADS_PER_WORKER:-1}
[[ $DFGM_THREADS_PER_WORKER =~ ^[0-9]+$ ]] && (( DFGM_THREADS_PER_WORKER > 0 )) || exit 2
(( DFGM_THREADS_PER_WORKER <= SLURM_CPUS_PER_TASK )) || {
    echo 'Worker threads exceed CPUs per task' >&2; exit 2;
}
export JULIA_NUM_THREADS=$DFGM_THREADS_PER_WORKER JULIA_NUM_GC_THREADS=1 JULIA_NUM_PRECOMPILE_TASKS=1 JULIA_PKG_PRECOMPILE_AUTO=0
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
source_dir=${1:?source directory required}
task_file=${2:?task file required}
output_dir=${3:?output directory required}
mkdir -p "$output_dir"
cd "$source_dir"
export DFGM_SOURCE_DIR="$source_dir" DFGM_OUTPUT_DIR="$output_dir"
export DFGM_WORKER_MEM_MB=${DFGM_WORKER_MEM_MB:-1800}
[[ $DFGM_WORKER_MEM_MB =~ ^[0-9]+$ ]] && (( DFGM_WORKER_MEM_MB > 0 )) || exit 2
if [[ -n ${SLURM_MEM_PER_NODE:-} ]]; then
    (( SLURM_NTASKS * DFGM_WORKER_MEM_MB <= SLURM_MEM_PER_NODE )) || {
        echo 'Worker memory requests exceed allocation' >&2; exit 2;
    }
fi
cp "$task_file" "$output_dir/tasks.tsv"
if [[ ${DFGM_RUN_TESTS:-0} == 1 ]]; then
    /home1/mkolar/juliaPkg/julia-1.13.2/bin/julia --project \
      -e 'using FunctionalGraphicalModels, LinearAlgebra, Random, Test; import FunctionalDataAnalysis as FDA; include("test/direct.jl")'
fi
if [[ ${DFGM_CANARY:-0} == 1 ]]; then
    /home1/mkolar/juliaPkg/julia-1.13.2/bin/julia --project experiments/dfgm2020.jl \
      --setting 1 --replicate 1 --n 30 --p 10 --integration-points 31 \
      --nlambda 3 --work-points 15 --out "$output_dir/canary-refit"
fi
worker_status=0
parallel --colsep '\t' --jobs "$SLURM_NTASKS" --joblog "$output_dir/parallel.tsv" \
    --halt never \
    'srun --exclusive --exact -N1 -n1 -c"${DFGM_THREADS_PER_WORKER}" --mem="${DFGM_WORKER_MEM_MB}M" bash scripts/dfgm_worker.sh {1} {2} {3} {4} {5} {6} {7} {8}' \
    :::: "$task_file" || worker_status=$?
python3 experiments/summarize_dfgm.py "$output_dir"
exit "$worker_status"
