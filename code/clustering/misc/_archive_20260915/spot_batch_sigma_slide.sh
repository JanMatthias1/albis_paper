#!/bin/bash
#SBATCH --job-name=spot_batch_slide
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/logs/spot_batch_slide_%A_%a.out
#SBATCH --array=0-1%2
#SBATCH --time=02:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --partition=shared
#
# Extends the cell/bin16um batch_sigma slide to spot, at spot's own Phase-1
# winning config (lambda=0.1/k_geom=8, from spot_kgeom_sweep_fine.sh). Spot
# already has endpoints at bs=0 (prebatch, ARI~0.40-0.53) and bs=0.3
# (canonical, ARI~0.36-0.40) -- unlike cell/bin, spot does NOT fully collapse
# at its canonical value, so this checks where its cliff (if any) actually
# sits relative to 0.3. Params match spot's strongmix tag exactly
# (generate_strongmix.sh): log_mu=-2.5, theta=2.0, jitter=1.0, sphere_r_um=2050.
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

BATCH_SIGMAS=(0.1 0.2)
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
