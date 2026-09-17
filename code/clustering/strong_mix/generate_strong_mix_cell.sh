#!/bin/bash
#SBATCH --job-name=strong_mix_cell
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/strong_mix/logs/strong_mix_cell_%A_%a.out
#SBATCH --array=0-5%4
#SBATCH --time=18:00:00
#SBATCH --mem=320G
#SBATCH --cpus-per-task=8
#SBATCH --partition=shared
#
# Cell's strong-domain-mix batch_sigma slide, consolidated 2026-09-17 from
# the scattered points built across code/clustering/misc/
# cellbin_batch_sigma_slide.sh (bs0.5/1.0) + batch_sigma_slide_fine.sh
# (bs0.1/0.2/0.3/0.4) into one self-contained array script, output moved
# from data/figure_3/cellbin_batch_sigma_slide/cell/ to data/strong_mix/cell/.
# Config unchanged: BANKSY lambda=0.5/k_geom=200 (cell's own Phase-1 ceiling,
# near-optimal and unchanged from the domain-panel tuning), dispersion
# log_mu=-2.3/theta=0.40/jitter=0.15 (matching cell's Figure 2/3 shared tag),
# sphere_r_um=6000 (cell's native geometry, NOT the smaller_sphere disc).
#
# Grid: the 6 points already generated (0.1, 0.2, 0.3, 0.4, 0.5, 1.0) --
# bs=0 and the canonical bs=1.5 live separately under
# data/figure_3/banksy_batch_compare/cell/{prebatch,tuned}/ (built by a
# different script at the same lambda/k_geom); not reproduced here to avoid
# a duplicate, expensive 600k-cell regeneration -- add them as array indices
# 6/7 if this script needs to become the single source of truth for cell's
# whole curve.
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
KG=200

BATCH_SIGMAS=(0.1 0.2 0.3 0.4 0.5 1.0)
BS="${BATCH_SIGMAS[${SLURM_ARRAY_TASK_ID:-0}]}"
TAG="cell_strong_mix_bs${BS}"

SIM_RAW="sim_paper/data/noisy/${TAG}/simulation_cell_z.h5ad"
SIM_QC="sim_paper/data/noisy/${TAG}/simulation_cell_z_qc.h5ad"
RUN="${OUT_ROOT}/cell/bs${BS}"

echo "[task ${SLURM_ARRAY_TASK_ID:-0}] cell batch_sigma=${BS} lambda=${LAM} k_geom=${KG}"

if [[ ! -f "${SIM_RAW}" ]]; then
    "${BANKSY_PYTHON}" sim_paper/code/data/generate_simulation_noisy.py \
        --modality cell --sphere-r-um 6000 \
        --base-gene-lognormal -2.3 0.7 --theta 0.40 --theta-jitter 0.15 \
        --strong-domain-mix --batch-sigma "${BS}" \
        --out-tag "${TAG}"
fi
if [[ ! -f "${SIM_QC}" ]]; then
    "${TUTORIAL_PYTHON}" sim_paper/code/clustering/00_qc_filter.py \
        --modality cell --input "${SIM_RAW}" --output "${SIM_QC}"
fi

mkdir -p "${RUN}"
BANKSY_H5AD="${RUN}/banksy_matrix/simulation_cell_z_banksy_pca_harmony_qc.h5ad"
if [[ ! -f "${BANKSY_H5AD}" ]]; then
    "${BANKSY_PYTHON}" sim_paper/code/clustering/01_build_banksy_matrix.py \
        --modality cell --input "${SIM_QC}" \
        --lambda "${LAM}" --k-geom "${KG}" --max-m 1 \
        --stagger-scale 5 --skip-umap \
        --output-dir "${RUN}/banksy_matrix"
fi

"${TUTORIAL_PYTHON}" sim_paper/code/clustering/02_leiden_resolution_sweep.py \
    --modality cell --packing-tag "strong_mix_cell_bs${BS}" \
    --input "${BANKSY_H5AD}" --output-dir "${RUN}/ari"

"${TUTORIAL_PYTHON}" sim_paper/code/clustering/composition_recovery.py \
    --h5ad "${RUN}/ari/simulation_cell_z_ari_recovery.h5ad" --tag "cell_bs${BS}" \
    --out-dir "${OUT_ROOT}/scores"

echo "[done] task ${SLURM_ARRAY_TASK_ID:-0} -> ${RUN}"
