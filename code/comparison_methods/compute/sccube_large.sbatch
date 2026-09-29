#!/usr/bin/env bash
#SBATCH --partition=shared
#SBATCH --constraint=sapphirerapids
#SBATCH --cpus-per-task=2
#SBATCH --mem=96G
#SBATCH --time=04:00:00
set -euo pipefail
project=/dcs04/hicks/data/Jan/sim_project
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export PYTHONDONTWRITEBYTECODE=1 PYTHONHASHSEED=2025
export MPLCONFIGDIR="$RUN_ROOT/cache/split/matplotlib" XDG_CACHE_HOME="$RUN_ROOT/cache/split"
export NUMBA_CACHE_DIR="$RUN_ROOT/cache/split/numba"
phase_file="$RUN_ROOT/rss_phase.txt"
done_file="$RUN_ROOT/rss_done.txt"
"$project/comparison_methods/env/sccube/bin/python" -u "$RUN_ROOT/code/sccube_split.py" --root "$RUN_ROOT" --source "$SOURCE_ROOT" --mode split --sizes 1000000 2000000 5000000 --settings-dir "$RUN_ROOT/settings" --skip-output-hashes --phase-file "$phase_file" --memory-budget-gib 96 &
worker_pid=$!
"$project/comparison_methods/env/sccube/bin/python" "$RUN_ROOT/code/sample_rss.py" --pid "$worker_pid" --phase-file "$phase_file" --done-file "$done_file" --output "$RUN_ROOT/rss_samples.csv" &
sampler_pid=$!
worker_status=0
wait "$worker_pid" || worker_status=$?
touch "$done_file"
wait "$sampler_pid"
exit "$worker_status"
