#!/usr/bin/env bash
#SBATCH --job-name=figure5d_actual
#SBATCH --partition=shared
#SBATCH --cpus-per-task=8
#SBATCH --mem=192G
#SBATCH --time=04:00:00
set -euo pipefail
export OPENBLAS_NUM_THREADS=8 OMP_NUM_THREADS=8 NUMBA_NUM_THREADS=8 MPLCONFIGDIR=/tmp/fig5d_actual PYTHONDONTWRITEBYTECODE=1
/dcs04/hicks/data/Jan/sim_project/sim_paper/env/albis-tutorial/bin/python -u /dcs04/hicks/data/Jan/sim_project/sim_paper/code/comparison_methods/figure_5D/actual_capture.py
