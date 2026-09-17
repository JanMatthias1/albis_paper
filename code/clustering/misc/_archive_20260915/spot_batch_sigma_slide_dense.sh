#!/bin/bash
#SBATCH --job-name=spot_batch_slide_dense
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/logs/spot_batch_slide_dense_%A_%a.out
#SBATCH --array=0-5%6
#SBATCH --time=02:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --partition=shared
#
# Brackets the two transition edges found in the coarser spot batch_sigma
# slide (batch_sigma_slide_fine.sh + spot_batch_sigma_slide.sh): ARI is
# clean ~0.42-0.53 through bs<=0.12, collapses to ~0 for bs in [0.15,0.25],
# then unexpectedly RECOVERS to ARI~0.36 at bs=0.3 (spot's own canonical
# value) -- a dip, not a cliff. This fills both edges of the dip:
#   down-edge: 0.13, 0.14   (between clean 0.12 and collapsed 0.15)
#   up-edge:   0.26, 0.27, 0.28, 0.29  (between collapsed 0.25 and recovered 0.30)
# Same pipeline/config as spot_batch_sigma_slide.sh (lambda=0.1/k_geom=8,
# log_mu=-2.5, theta=2.0, jitter=1.0, sphere_r_um=2050) so points slot
# straight into the existing cellbin_batch_sigma_slide/spot/ curve --
# plot_batch_sigma_slide.py auto-discovers every bs<value> dir, no edits needed.
set -euo pipefail

mkdir -p /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/logs
source /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/_env_banksy.sh
BANKSY_PYTHON="${PYTHON_BIN}"
source /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/_env.sh
TUTORIAL_PYTHON="${PYTHON_BIN}"
cd /dcs04/hicks/data/Jan/sim_project

OUT_ROOT="sim_paper/data/figure_3/cellbin_batch_sigma_slide"
LAM=0.1
KG=8

BATCH_SIGMAS=(0.13 0.14 0.26 0.27 0.28 0.29)
BS="${BATCH_SIGMAS[${SLURM_ARRAY_TASK_ID:-0}]}"
TAG="spot_batch_slide_bs${BS}"

SIM_RAW="sim_paper/data/noisy/${TAG}/simulation_spot_z.h5ad"
SIM_QC="sim_paper/data/noisy/${TAG}/simulation_spot_z_qc.h5ad"
RUN="${OUT_ROOT}/spot/bs${BS}"

echo "[task ${SLURM_ARRAY_TASK_ID:-0}] spot batch_sigma=${BS} lambda=${LAM} k_geom=${KG}"

if [[ ! -f "${SIM_RAW}" ]]; then
    "${BANKSY_PYTHON}" sim_paper/code/data/generate_simulation_noisy.py \
        --modality spot --sphere-r-um 2050 \
        --base-gene-lognormal -2.5 0.7 --theta 2.0 --theta-jitter 1.0 \
        --strong-domain-mix --batch-sigma "${BS}" \
        --out-tag "${TAG}"
fi
if [[ ! -f "${SIM_QC}" ]]; then
    "${TUTORIAL_PYTHON}" sim_paper/code/clustering/00_qc_filter.py \
        --modality spot --input "${SIM_RAW}" --output "${SIM_QC}"
fi

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
    --modality spot --packing-tag "cellbin_batch_slide_spot_bs${BS}" \
    --input "${BANKSY_H5AD}" --output-dir "${RUN}/ari"

"${TUTORIAL_PYTHON}" sim_paper/code/clustering/composition_recovery.py \
    --h5ad "${RUN}/ari/simulation_spot_z_ari_recovery.h5ad" --tag "spot_bs${BS}" \
    --out-dir "${OUT_ROOT}/scores"

echo "[done] bs${BS}  ->  ${RUN}"
