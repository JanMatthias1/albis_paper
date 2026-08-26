#!/bin/bash
#SBATCH --job-name=replot_pca_harmony_qc
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/logs/replot_pca_harmony_qc_%A_%a.out
#SBATCH --array=0-2
#SBATCH --time=04:00:00
#SBATCH --mem=320G
#SBATCH --cpus-per-task=4
#SBATCH --partition=shared
#
# One-off: regenerate the pca_harmony.py plots (PCA + before/after UMAP) for
# the final Figure 3 cell-type-panel outputs after a plot-formatting change
# (bold title, "slice_id" -> "Slice Id"), without rerunning PCA/Harmony.
# Uses --plots-only against the already-computed figure_3/pca_harmony_single_cell
# h5ad for each modality.

set -euo pipefail

mkdir -p /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/logs
source /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/_env.sh
cd /dcs04/hicks/data/Jan/sim_project

MODALITIES=(spot bin cell)
MODALITY="${MODALITIES[${SLURM_ARRAY_TASK_ID:-0}]}"
H5AD="sim_paper/data/figure_3/pca_harmony_single_cell/${MODALITY}/pca_harmony_qc/simulation_${MODALITY}_z_pca_harmony_qc.h5ad"

echo "Modality: ${MODALITY}"
"${PYTHON_BIN}" sim_paper/code/clustering/pca_harmony.py \
    --modality "${MODALITY}" \
    --input "${H5AD}" \
    --output "${H5AD}" \
    --plots-only
