#!/bin/bash
#SBATCH --account=mkolar_1314
#SBATCH --partition=main
#SBATCH --constraint=epyc-7513|epyc-7542
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=00:10:00
set -euo pipefail
python3 /project2/mkolar_1314/FunctionalGraphicalModels/20261025-direct/spectral-dispatch-v1/submit_paper100.py
