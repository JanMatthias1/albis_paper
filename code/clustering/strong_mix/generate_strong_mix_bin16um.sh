#!/bin/bash
#SBATCH --job-name=strong_mix_bin16um
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/strong_mix/logs/strong_mix_bin16um_%A_%a.out
#SBATCH --array=0-6%4
#SBATCH --time=18:00:00
#SBATCH --mem=320G
#SBATCH --cpus-per-task=8
#SBATCH --partition=shared
#
# Bin16um's strong-domain-mix batch_sigma slide, consolidated 2026-09-17
# from the scattered points built across code/clustering/misc/
# bin16um_kgeom100_endpoints.sh (bs0/bs0.7) + cellbin_batch_sigma_slide.sh
# (bs0.25/0.45) + batch_sigma_slide_fine.sh (bs0.30/0.35/0.40) into one
# self-contained array script, output moved from
# data/figure_3/cellbin_batch_sigma_slide/bin/ to data/strong_mix/bin16um/.
# Config unchanged: BANKSY lambda=0.5/k_geom=100 (bin16um's own Phase-1
# winner, beats the older k_geom=200 -- 0.806 vs 0.784 domain ARI at
# bs=0), dispersion log_mu=-2.5 (defaults for theta/jitter), batch_sigma
# canonical=0.7, sphere_r_um=2050 (smaller_sphere disc), bin_size_um=16.
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
LAM=0.5
KG=100

BATCH_SIGMAS=(0 0.25 0.30 0.35 0.40 0.45 0.7)
BS="${BATCH_SIGMAS[${SLURM_ARRAY_TASK_ID:-0}]}"
TAG="bin16um_strong_mix_bs${BS}"

SIM_RAW="sim_paper/data/noisy/${TAG}/simulation_bin_z.h5ad"
SIM_QC="sim_paper/data/noisy/${TAG}/simulation_bin_z_qc.h5ad"
RUN="${OUT_ROOT}/bin16um/bs${BS}"

echo "[task ${SLURM_ARRAY_TASK_ID:-0}] bin16um batch_sigma=${BS} lambda=${LAM} k_geom=${KG}"

if [[ ! -f "${SIM_RAW}" ]]; then
    "${BANKSY_PYTHON}" sim_paper/code/data/generate_simulation_noisy.py \
        --modality bin --sphere-r-um 2050 --bin-size-um 16 \
        --base-gene-lognormal -2.5 0.7 \
        --strong-domain-mix --batch-sigma "${BS}" \
        --out-tag "${TAG}"
fi
if [[ ! -f "${SIM_QC}" ]]; then
    "${TUTORIAL_PYTHON}" sim_paper/code/clustering/00_qc_filter.py \
        --modality bin --input "${SIM_RAW}" --output "${SIM_QC}"
fi

mkdir -p "${RUN}"
BANKSY_H5AD="${RUN}/banksy_matrix/simulation_bin_z_banksy_pca_harmony_qc.h5ad"
if [[ ! -f "${BANKSY_H5AD}" ]]; then
    "${BANKSY_PYTHON}" sim_paper/code/clustering/01_build_banksy_matrix.py \
        --modality bin --input "${SIM_QC}" \
        --lambda "${LAM}" --k-geom "${KG}" --max-m 1 \
        --stagger-scale 5 --skip-umap \
        --output-dir "${RUN}/banksy_matrix"
fi

"${TUTORIAL_PYTHON}" sim_paper/code/clustering/02_leiden_resolution_sweep.py \
    --modality bin --packing-tag "strong_mix_bin16um_bs${BS}" \
    --input "${BANKSY_H5AD}" --output-dir "${RUN}/ari"

"${TUTORIAL_PYTHON}" sim_paper/code/clustering/composition_recovery.py \
    --h5ad "${RUN}/ari/simulation_bin_z_ari_recovery.h5ad" --tag "bin16um_bs${BS}" \
    --out-dir "${OUT_ROOT}/scores"

echo "[done] task ${SLURM_ARRAY_TASK_ID:-0} -> ${RUN}"
