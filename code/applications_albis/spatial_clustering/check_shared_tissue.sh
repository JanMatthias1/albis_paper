#!/usr/bin/env bash
#SBATCH --job-name=stagate_shared_check
#SBATCH --partition=gpu
#SBATCH --gres=gpu:tesh100:1
#SBATCH --mem=150G
#SBATCH --cpus-per-task=4
#SBATCH --time=00:45:00
set -euo pipefail
project=/dcs04/hicks/data/Jan/sim_project
out="${1:?new feasibility output directory}"
code="$project/sim_paper/code/applications_albis/spatial_clustering"
stagate_env="$project/sim_paper/env/stagate-pyg"
export LD_LIBRARY_PATH="$stagate_env/lib:${LD_LIBRARY_PATH:-}" R_HOME="$stagate_env/lib/R"
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
export PYTHONDONTWRITEBYTECODE=1 MPLCONFIGDIR="$out/cache/matplotlib" XDG_CACHE_HOME="$out/cache"
mkdir -p "$MPLCONFIGDIR"
status=0
for arm in 3d 2d; do
    "$stagate_env/bin/python" -u "$code/check_shared_tissue.py" --out "$out" --arm "$arm" --epochs 5 > "$out/${arm}.log" 2>&1 || status=1
done
exit "$status"
