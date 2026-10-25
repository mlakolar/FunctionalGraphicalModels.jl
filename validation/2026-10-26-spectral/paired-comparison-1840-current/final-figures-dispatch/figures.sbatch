#!/bin/bash
#SBATCH --account=mkolar_1314
#SBATCH --partition=main
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:10:00
set -euo pipefail
root=/project2/mkolar_1314/FunctionalGraphicalModels/20261025-direct
/home1/mkolar/juliaPkg/runs/20261024-fgm-development/plot-env/bin/python "$root/spectral-figures-final-v1/plot.py" "$root/outputs/spectral-full-final-v3-13448698" "$root/outputs/spectral-controlled-timing-13447505" "$root/outputs/spectral-final-figures-$SLURM_JOB_ID"
