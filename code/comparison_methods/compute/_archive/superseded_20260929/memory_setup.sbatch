#!/usr/bin/env bash
#SBATCH --partition=shared
#SBATCH --constraint=sapphirerapids
#SBATCH --cpus-per-task=2
#SBATCH --mem=96G
#SBATCH --time=04:00:00
#SBATCH --array=0-23%4
set -euo pipefail
project=/dcs04/hicks/data/Jan/sim_project
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export PYTHONDONTWRITEBYTECODE=1 PYTHONHASHSEED=2025 PYTHONPATH="$project/albis"
"$project/comparison_methods/env/analysis/bin/python" "$RUN_ROOT/code/measure_memory_setup.py" --root "$RUN_ROOT" --task "$SLURM_ARRAY_TASK_ID"
"$project/comparison_methods/env/analysis/bin/python" "$RUN_ROOT/code/plot_compute_four_panel.py" --root "$RUN_ROOT" --task-only
