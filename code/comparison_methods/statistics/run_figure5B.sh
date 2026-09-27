#!/usr/bin/env bash
#SBATCH --job-name=figure5B
#SBATCH --partition=shared
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/comparison_methods/statistics/logs/figure5B_%j.out
set -euo pipefail
project=/dcs04/hicks/data/Jan/sim_project
export MPLCONFIGDIR=/tmp/figure5B_mpl_${SLURM_JOB_ID:-local}
export OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 NUMBA_NUM_THREADS=2
for modality in cell bin spot; do
    "$project/comparison_methods/env/analysis/bin/python" -u "$project/sim_paper/code/comparison_methods/statistics/plot_distributions.py" --modality "$modality"
    "$project/comparison_methods/env/analysis/bin/python" -u "$project/sim_paper/code/comparison_methods/statistics/plot_distributions.py" --modality "$modality" --hvg-match
done
"$project/comparison_methods/env/analysis/bin/python" -u "$project/sim_paper/code/comparison_methods/statistics/plot_log_normalized.py"
