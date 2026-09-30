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
root="$project/sim_paper/data/figure_5/clustering/domain"
method="${1:?Supply albis}"
case "$method" in albis) ;; *) exit 2 ;; esac
run="$root/$method/pilot_${SLURM_JOB_ID:-local}"
mkdir -p "$run"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
export MPLCONFIGDIR="$run/cache/matplotlib" XDG_CACHE_HOME="$run/cache" NUMBA_CACHE_DIR="$run/cache/numba"
mkdir -p "$MPLCONFIGDIR" "$NUMBA_CACHE_DIR"
input="$project/sim_paper/data/figure_4/cross_modality_alignment/strong_domain_mix_shift3x/data/cell/simulation_cell_z.h5ad"
"$project/sim_paper/env/sim-app-banksy/bin/python" -u "$code/run_banksy.py" \
    --input "$input" --out "$run/banksy" > "$run/banksy.log" 2>&1
printf 'Finished: %s\n' "$run/banksy/metrics.json"
