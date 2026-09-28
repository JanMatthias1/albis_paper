#!/usr/bin/env bash
#SBATCH --partition=shared
#SBATCH --constraint=sapphirerapids
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=1-00:00:00
set -euo pipefail
project=/dcs04/hicks/data/Jan/sim_project
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export PYTHONDONTWRITEBYTECODE=1
export PYTHONPATH="$project/albis"
export MPLCONFIGDIR="${OUTPUT}_cache/matplotlib" XDG_CACHE_HOME="${OUTPUT}_cache" NUMBA_CACHE_DIR="${OUTPUT}_cache/numba"
# The Python script creates OUTPUT exclusively; caches are lazy-created inside it.
exec "$project/comparison_methods/env/spider/bin/python" -u "$project/sim_paper/code/comparison_methods/compute/spider_native.py"
