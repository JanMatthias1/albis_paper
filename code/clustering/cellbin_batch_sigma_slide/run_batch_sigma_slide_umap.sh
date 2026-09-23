#!/bin/bash
#SBATCH --job-name=batch_slide_umap
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/logs/batch_slide_umap_%A_%a.out
#SBATCH --array=0-5%6
#SBATCH --time=04:00:00
#SBATCH --mem=120G
#SBATCH --cpus-per-task=8
#SBATCH --partition=shared
#
# Before/after-Harmony UMAPs at the batch_sigma-slide cliff points. Cell's
# bs=0/1.5 already have UMAPs from job 35638000 (k_geom=200, unchanged) --
# only its 2 new intermediate points need building. Bin16um's k_geom changed
# (200 -> 100) so ALL 4 of its points need fresh UMAPs.
set -euo pipefail

mkdir -p /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/logs
source /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/_env.sh
cd /dcs04/hicks/data/Jan/sim_project

SLIDE_ROOT="sim_paper/data/figure_3/cellbin_batch_sigma_slide"
MODS=(cell cell bin bin bin bin)
LABELS=(bs0.5 bs1.0 bs0 bs0.25 bs0.45 bs0.7)

i="${SLURM_ARRAY_TASK_ID:-0}"
MOD="${MODS[$i]}"
LABEL="${LABELS[$i]}"

FOLDER="${MOD}"
if [[ "$MOD" == bin ]]; then FOLDER=bin16um; fi

H5AD="${SLIDE_ROOT}/${FOLDER}/${LABEL}/banksy_matrix/simulation_${MOD}_z_banksy_pca_harmony_qc.h5ad"

echo "[task ${i}] modality=${MOD} label=${LABEL}"
"${PYTHON_BIN}" sim_paper/code/clustering/plot_banksy_results.py \
    --input "${H5AD}"
echo "[done] task ${i}"
