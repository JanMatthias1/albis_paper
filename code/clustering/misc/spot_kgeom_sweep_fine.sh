#!/bin/bash
#SBATCH --job-name=spot_kgeom_fine
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/logs/spot_kgeom_fine_%A_%a.out
#SBATCH --array=0-19%12
#SBATCH --time=01:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --partition=shared
#
# Fine sweep around spot's winning point from spot_kgeom_sweep_tuned.sh
# (lambda=0.2, k_geom=6 -> domain ARI 0.271, leak 0.236, real tuned batch).
# umap_pca_harmony_before_after_by_slice_id.png at that point shows Harmony
# merges most (~7/10) slices into one manifold but 1-2 slices (esp. slice 4)
# resist -- checking whether nearby lambda/k_geom values close that gap
# further, or whether 6 is already the local optimum.
set -euo pipefail

mkdir -p /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/logs
source /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/_env_banksy.sh
BANKSY_PYTHON="${PYTHON_BIN}"
source /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/_env.sh
TUTORIAL_PYTHON="${PYTHON_BIN}"
cd /dcs04/hicks/data/Jan/sim_project

SIM_QC="sim_paper/data/figure_4/spatial_clustering/sim_data/packing_pf0p04_log_mu_-2.5_bsigma03_strongmix/simulation_spot_z_qc.h5ad"
OUT_ROOT="sim_paper/data/figure_3/spot_kgeom_sweep_fine"

# 4 lambdas x 5 k_geoms = 20 combos, all TUNED (real) batch
LAMBDAS=(0.1 0.1 0.1 0.1 0.1 0.15 0.15 0.15 0.15 0.15 0.2 0.2 0.2 0.2 0.2 0.3 0.3 0.3 0.3 0.3)
KGEOMS=(4 5 6 8 10 4 5 6 8 10 4 5 6 8 10 4 5 6 8 10)

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
    --modality spot --packing-tag "spot_kgeom_fine_${TAG}" \
    --input "${BANKSY_H5AD}" --output-dir "${RUN}/ari"

"${TUTORIAL_PYTHON}" sim_paper/code/clustering/composition_recovery.py \
    --h5ad "${RUN}/ari/simulation_spot_z_ari_recovery.h5ad" --tag "spot_${TAG}" \
    --out-dir "${OUT_ROOT}/scores"

echo "[done] task ${i}  ->  ${RUN}"
