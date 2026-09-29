#!/usr/bin/env bash
#SBATCH --partition=shared
#SBATCH --constraint=sapphirerapids
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
set -euo pipefail
project=/dcs04/hicks/data/Jan/sim_project
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export PYTHONDONTWRITEBYTECODE=1 PYTHONHASHSEED=2025
export MPLCONFIGDIR="$RUN_ROOT/cache/$MODE/matplotlib" XDG_CACHE_HOME="$RUN_ROOT/cache/$MODE"
export NUMBA_CACHE_DIR="$RUN_ROOT/cache/$MODE/numba"
exec "$project/comparison_methods/env/sccube/bin/python" -u "$RUN_ROOT/code/sccube_split.py" --root "$RUN_ROOT" --source "$SOURCE_ROOT" --mode "$MODE"
