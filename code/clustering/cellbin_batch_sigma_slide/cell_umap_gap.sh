#!/bin/bash
#SBATCH --job-name=cell_umap_gap
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/logs/cell_umap_gap_%A_%a.out
#SBATCH --array=0-3%4
#SBATCH --time=04:00:00
#SBATCH --mem=120G
#SBATCH --cpus-per-task=8
#SBATCH --partition=shared
#
# Before/after-Harmony UMAPs for the cell modality's remaining intermediate
# batch_sigma-slide points (cellbin_batch_sigma_slide/cell/bs{0.1,0.2,0.3,0.4}/)
# -- these currently only have PCA-only diagnostics. cell's bs0.5/bs1.0
# already got UMAPs from run_batch_sigma_slide_umap.sh, and bs0/bs1.5 have
# full diagnostics via banksy_batch_compare/cell/{prebatch,tuned}. Same call
# pattern/resources as run_batch_sigma_slide_umap.sh.
set -euo pipefail

mkdir -p /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/logs
source /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/_env.sh
cd /dcs04/hicks/data/Jan/sim_project

LABELS=(bs0.1 bs0.2 bs0.3 bs0.4)

i="${SLURM_ARRAY_TASK_ID:-0}"
LABEL="${LABELS[$i]}"

H5AD="sim_paper/data/figure_3/cellbin_batch_sigma_slide/cell/${LABEL}/banksy_matrix/simulation_cell_z_banksy_pca_harmony_qc.h5ad"

echo "[task ${i}] label=${LABEL}"
"${PYTHON_BIN}" sim_paper/code/clustering/plot_banksy_results.py \
    --input "${H5AD}"
echo "[done] task ${i}"
