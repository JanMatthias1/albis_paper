#!/usr/bin/env bash
#SBATCH --partition=shared
#SBATCH --constraint=sapphirerapids
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=02:00:00
set -euo pipefail
project=/dcs04/hicks/data/Jan/sim_project
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export PYTHONDONTWRITEBYTECODE=1 PYTHONHASHSEED=2025
export MPLCONFIGDIR="${OUTPUT}_cache/matplotlib" XDG_CACHE_HOME="${OUTPUT}_cache" NUMBA_CACHE_DIR="${OUTPUT}_cache/numba"
case "$METHOD" in
    albis) python_bin="$project/comparison_methods/env/our_method/bin/python" ;;
    spider) python_bin="$project/comparison_methods/env/spider/bin/python" ;;
    *) exit 2 ;;
esac
exec "$python_bin" -u "$RUN_ROOT/code/${METHOD}_native.py"
