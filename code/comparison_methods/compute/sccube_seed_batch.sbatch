#!/usr/bin/env bash
#SBATCH --partition=shared
#SBATCH --constraint=sapphirerapids
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=04:00:00
# scCube for one seed and one size (SIZES): train a VAE, then generate that size, with
# 100 ms RSS sampling.
# Env: RUN_ROOT, SEED, BATCH (folder name), SIZES (space-separated cell counts).
set -euo pipefail
project=/dcs04/hicks/data/Jan/sim_project
code="$project/sim_paper/code/comparison_methods/compute"
out="$RUN_ROOT/sccube/seed${SEED}_${BATCH}"
mkdir -p "$RUN_ROOT/sccube" "$RUN_ROOT/cache/sccube_seed${SEED}_${BATCH}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export PYTHONDONTWRITEBYTECODE=1 PYTHONHASHSEED=2025
export MPLCONFIGDIR="$RUN_ROOT/cache/sccube_seed${SEED}_${BATCH}/matplotlib" XDG_CACHE_HOME="$RUN_ROOT/cache/sccube_seed${SEED}_${BATCH}"
export NUMBA_CACHE_DIR="$RUN_ROOT/cache/sccube_seed${SEED}_${BATCH}/numba"
phase_file="$RUN_ROOT/sccube/seed${SEED}_${BATCH}_rss_phase.txt"
done_file="$RUN_ROOT/sccube/seed${SEED}_${BATCH}_rss_done.txt"
# shellcheck disable=SC2086
"$project/comparison_methods/env/sccube/bin/python" -u "$code/sccube_split.py" --root "$RUN_ROOT" --source "$RUN_ROOT" \
    --mode split --seed "$SEED" --sizes $SIZES --settings-dir "$RUN_ROOT/settings" --out "$out" \
    --skip-output-hashes --phase-file "$phase_file" --memory-budget-gib 96 &
worker_pid=$!
"$project/comparison_methods/env/sccube/bin/python" "$code/sample_rss.py" --pid "$worker_pid" --phase-file "$phase_file" \
    --done-file "$done_file" --output "$RUN_ROOT/sccube/seed${SEED}_${BATCH}_rss_samples.csv" &
sampler_pid=$!
worker_status=0
wait "$worker_pid" || worker_status=$?
touch "$done_file"
wait "$sampler_pid"
exit "$worker_status"
