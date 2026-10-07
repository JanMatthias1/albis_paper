#!/usr/bin/env bash
#SBATCH --job-name=fig3_condensed
#SBATCH --partition=shared
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:15:00
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/albis_paper/data/figure_3/ari_recovery_summary/condensed_%j.log
# Run with bash for local plotting, or sbatch on the cluster. No clustering is rerun.
set -euo pipefail
project="${SIM_PROJECT_ROOT:-/dcs04/hicks/data/Jan/sim_project}"
export MPLCONFIGDIR="${TMPDIR:-/tmp}/fig3_condensed_mpl_${SLURM_JOB_ID:-$$}"
export XDG_CACHE_HOME="${TMPDIR:-/tmp}/fig3_condensed_cache_${SLURM_JOB_ID:-$$}"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
mkdir -p "$MPLCONFIGDIR" "$XDG_CACHE_HOME"
exec "$project/albis_paper/env/albis-tutorial/bin/python" -u \
    "$project/albis_paper/code/clustering/ari_recovery_summary/plot_banksy_vs_pca_condensed.py"
