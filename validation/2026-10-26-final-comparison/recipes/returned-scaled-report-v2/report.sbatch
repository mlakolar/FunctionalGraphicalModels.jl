#!/bin/bash
#SBATCH --account=mkolar_1314
#SBATCH --partition=main
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
set -euo pipefail
root=/project2/mkolar_1314/FunctionalGraphicalModels/20261025-direct
python3 "$root/returned-scaled-report-v2/report.py" "$root" "$root/outputs/returned-scaled-report-$SLURM_JOB_ID"
