#!/bin/bash
#SBATCH --account=mkolar_1314
#SBATCH --partition=main
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:05:00
set -euo pipefail
root=/project2/mkolar_1314/FunctionalGraphicalModels/20261025-direct
/home1/mkolar/juliaPkg/runs/20261024-fgm-development/plot-env/bin/python "$root/final-paired-figure-v1/plot.py" "$root/outputs/paper-v4/reports/report-recovery-13447246-0" "$root/outputs/final-paired-figure-$SLURM_JOB_ID"
