#!/usr/bin/env bash
#SBATCH --partition=shared
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:10:00
set -euo pipefail
project=/dcs04/hicks/data/Jan/sim_project
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
export MPLCONFIGDIR="$RUN_ROOT/cache/report/matplotlib" XDG_CACHE_HOME="$RUN_ROOT/cache/report"
exec "$project/comparison_methods/env/analysis/bin/python" "$RUN_ROOT/code/report_native_extension.py" --root "$RUN_ROOT"
