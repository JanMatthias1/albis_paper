#!/usr/bin/env bash
#SBATCH --job-name=banksy_domains
#SBATCH --partition=shared
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --time=12:00:00
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/comparison_methods/clustering/domain/logs/pilot_%j.out
set -euo pipefail
project=/dcs04/hicks/data/Jan/sim_project
code="$project/sim_paper/code/comparison_methods/clustering/domain"
root="$project/sim_paper/data/comparison_methods/clustering/domain"
method="${1:?Supply albis or spider}"
case "$method" in albis|spider) ;; *) exit 2 ;; esac
run="$root/$method/pilot_${SLURM_JOB_ID:-local}"
mkdir -p "$run"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
export MPLCONFIGDIR="$run/cache/matplotlib" XDG_CACHE_HOME="$run/cache" NUMBA_CACHE_DIR="$run/cache/numba"
mkdir -p "$MPLCONFIGDIR" "$NUMBA_CACHE_DIR"
if [[ "$method" == spider ]]; then
    "$project/comparison_methods/env/spider/bin/python" -u "$code/generate_spider.py" \
        --out "$run/input" --pool "$project/sim_paper/data/comparison_methods/overview_native556_20260926/splatter" \
        > "$run/generation.log" 2>&1
    input="$run/input/cell.h5ad"
else
    input="$project/sim_paper/data/figure_4/cross_modality_alignment/strong_domain_mix_shift3x/data/cell/simulation_cell_z.h5ad"
fi
"$project/sim_paper/env/sim-app-banksy/bin/python" -u "$code/run_banksy.py" \
    --input "$input" --out "$run/banksy" > "$run/banksy.log" 2>&1
printf 'Finished: %s\n' "$run/banksy/metrics.json"
