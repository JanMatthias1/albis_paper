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
export MPLCONFIGDIR="$RUN_ROOT/cache/split/matplotlib" XDG_CACHE_HOME="$RUN_ROOT/cache/split"
export NUMBA_CACHE_DIR="$RUN_ROOT/cache/split/numba"
exec "$project/comparison_methods/env/sccube/bin/python" -u "$RUN_ROOT/code/sccube_split.py" --root "$RUN_ROOT" --source "$SOURCE_ROOT" --mode split --sizes 200000 400000 600000 --settings-dir "$RUN_ROOT/settings" --skip-output-hashes
