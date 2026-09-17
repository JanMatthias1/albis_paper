#!/bin/bash
#SBATCH --job-name=strong_mix_spot
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/strong_mix/logs/strong_mix_spot_%A_%a.out
#SBATCH --array=0-15%8
#SBATCH --time=02:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --partition=shared
#
# Spot's strong-domain-mix batch_sigma slide (no batch effect -> canonical
# 0.3), consolidated 2026-09-17 from code/clustering/misc/
# spot_batch_sigma_slide_corrected.sh into data/strong_mix/ (was
# data/figure_3/cellbin_batch_sigma_slide/). Logic unchanged -- spot's config
# (BANKSY lambda=0.1/k_geom=8, dispersion log_mu=-2.0/theta=0.25/jitter=0.10)
# was already a single clean array script covering the full 16-point grid,
# just relocated.
#
# Full pipeline per point: generate (strong-domain-mix) -> QC ->
# BANKSY+Harmony -> 02_leiden_resolution_sweep.py -> composition_recovery.py
# (ARI + slice_id leakage).
set -euo pipefail

mkdir -p /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/strong_mix/logs
source /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/_env_banksy.sh
BANKSY_PYTHON="${PYTHON_BIN}"
source /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/_env.sh
TUTORIAL_PYTHON="${PYTHON_BIN}"
cd /dcs04/hicks/data/Jan/sim_project

OUT_ROOT="sim_paper/data/figure_3/cellbin_batch_sigma_slide"
LAM=0.1
KG=8

BATCH_SIGMAS=(0.0 0.05 0.1 0.12 0.13 0.14 0.15 0.18 0.2 0.22 0.25 0.26 0.27 0.28 0.29 0.3)
BS="${BATCH_SIGMAS[${SLURM_ARRAY_TASK_ID:-0}]}"
TAG="spot_strong_mix_bs${BS}"

SIM_RAW="sim_paper/data/noisy/${TAG}/simulation_spot_z.h5ad"
SIM_QC="sim_paper/data/noisy/${TAG}/simulation_spot_z_qc.h5ad"
RUN="${OUT_ROOT}/spot/bs${BS}"

echo "[task ${SLURM_ARRAY_TASK_ID:-0}] spot batch_sigma=${BS} lambda=${LAM} k_geom=${KG}"

if [[ ! -f "${SIM_RAW}" ]]; then
    "${BANKSY_PYTHON}" sim_paper/code/data/generate_simulation_noisy.py \
        --modality spot --sphere-r-um 2050 \
        --base-gene-lognormal -2.0 0.7 --theta 0.25 --theta-jitter 0.10 \
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

"${TUTORIAL_PYTHON}" sim_paper/code/clustering/02_leiden_resolution_sweep.py \
    --modality spot --packing-tag "strong_mix_spot_bs${BS}" \
    --input "${BANKSY_H5AD}" --output-dir "${RUN}/ari"

"${TUTORIAL_PYTHON}" sim_paper/code/clustering/composition_recovery.py \
    --h5ad "${RUN}/ari/simulation_spot_z_ari_recovery.h5ad" --tag "spot_bs${BS}" \
    --out-dir "${OUT_ROOT}/scores"

echo "[done] task ${SLURM_ARRAY_TASK_ID:-0} -> ${RUN}"
