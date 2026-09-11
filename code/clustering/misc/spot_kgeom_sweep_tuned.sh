#!/bin/bash
#SBATCH --job-name=spot_kgeom_sweep
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/logs/spot_kgeom_sweep_%A_%a.out
#SBATCH --array=0-11%12
#SBATCH --time=01:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --partition=shared
#
# Spot-specific lambda x k_geom sweep at the REAL tuned batch effect, closing
# the gap left by banksy_batch_compare.sh (which only tested spot at the
# cell-tuned k_geom=200 -- confirmed an oversmoothing artifact, see
# check_spot_kgeom15_ari.py). k_geom values here are scaled to spot's own
# density (~1000 obs/slice), NOT cell/bin's sweep range (15/60/200, tuned for
# ~40-60k obs/slice). A run counts as a genuine hit only if domain ARI > 0.10
# AND slice_leak_ARI < 0.30 (same bar as the original lambda_kgeom sweep).
set -euo pipefail

mkdir -p /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/logs
source /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/_env_banksy.sh
BANKSY_PYTHON="${PYTHON_BIN}"
source /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/_env.sh
TUTORIAL_PYTHON="${PYTHON_BIN}"
cd /dcs04/hicks/data/Jan/sim_project

SIM_QC="sim_paper/data/figure_4/spatial_clustering/sim_data/packing_pf0p04_log_mu_-2.5_bsigma03_strongmix/simulation_spot_z_qc.h5ad"
OUT_ROOT="sim_paper/data/figure_3/spot_kgeom_sweep_tuned"

LAMBDAS=(0.2 0.2 0.2 0.2 0.5 0.5 0.5 0.5 0.8 0.8 0.8 0.8)
KGEOMS=(6 15 30 60 6 15 30 60 6 15 30 60)

i="${SLURM_ARRAY_TASK_ID:-0}"
LAM="${LAMBDAS[$i]}"
KG="${KGEOMS[$i]}"
TAG="lambda${LAM}_kgeom${KG}"
RUN="${OUT_ROOT}/${TAG}"

echo "[task ${i}] lambda=${LAM} k_geom=${KG} (tuned/real batch)"
mkdir -p "${RUN}"

BANKSY_H5AD="${RUN}/banksy_matrix/simulation_spot_z_banksy_pca_harmony_qc.h5ad"
if [[ ! -f "${BANKSY_H5AD}" ]]; then
    "${BANKSY_PYTHON}" sim_paper/code/clustering/01_build_banksy_matrix.py \
        --modality spot --input "${SIM_QC}" \
        --lambda "${LAM}" --k-geom "${KG}" --max-m 1 \
        --stagger-scale 5 --skip-umap \
        --output-dir "${RUN}/banksy_matrix"
fi

"${TUTORIAL_PYTHON}" sim_paper/code/clustering/ari_vs_ground_truth.py \
    --modality spot --packing-tag "spot_kgeom_sweep_${TAG}" \
    --input "${BANKSY_H5AD}" --output-dir "${RUN}/ari"

"${TUTORIAL_PYTHON}" sim_paper/code/clustering/composition_recovery.py \
    --h5ad "${RUN}/ari/simulation_spot_z_ari_recovery.h5ad" --tag "spot_${TAG}" \
    --out-dir "${OUT_ROOT}/scores"

echo "[done] task ${i}  ->  ${RUN}"
