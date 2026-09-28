#!/bin/bash
#SBATCH --job-name=compute_stress
#SBATCH --array=0-71%3
#SBATCH --time=2-00:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --partition=shared
set -euo pipefail
project=/dcs04/hicks/data/Jan/sim_project
root="$1"
counts=(10000 20000 40000 50000 100000 200000 500000 1000000)
technologies=(cell bin spot)
methods=(albis spider sccube)
i="$SLURM_ARRAY_TASK_ID"
n="${counts[$((i / 9))]}"
tech="${technologies[$(((i / 3) % 3))]}"
method="${methods[$((i % 3))]}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
export MPLCONFIGDIR="$root/cache/${i}/matplotlib" NUMBA_CACHE_DIR="$root/cache/${i}/numba" XDG_CACHE_HOME="$root/cache/${i}"
mkdir -p "$MPLCONFIGDIR" "$NUMBA_CACHE_DIR"
exec "$project/comparison_methods/env/analysis/bin/python" "$project/sim_paper/code/comparison_methods/compute/run_compute.py" --n-cells "$n" --technology "$tech" --method "$method" --out-dir "$root/raw/n${n}_${tech}_${method}"
