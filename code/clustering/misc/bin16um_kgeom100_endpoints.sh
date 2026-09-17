#!/bin/bash
#SBATCH --job-name=bin16um_kg100_endpoints
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/logs/bin16um_kg100_endpoints_%A_%a.out
#SBATCH --array=0-1%2
#SBATCH --time=18:00:00
#SBATCH --mem=320G
#SBATCH --cpus-per-task=8
#SBATCH --partition=shared
#
# Completes the bin16um batch_sigma curve at the Phase-1 winning k_geom=100
# (see cellbin_batch_sigma_slide.sh for the 2 intermediate points at
# bs=0.25/0.45). The 0 and 0.7 endpoints already exist at the OLD k_geom=200
# (data/figure_3/banksy_batch_compare/bin/{prebatch,tuned}/) -- no new
# simulation needed here, just a fresh BANKSY build at k_geom=100 on the
# SAME already-QC'd strongmix h5ad, so this is much cheaper than the
# intermediate points.
set -euo pipefail

mkdir -p /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/logs
source /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/_env_banksy.sh
BANKSY_PYTHON="${PYTHON_BIN}"
source /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/_env.sh
TUTORIAL_PYTHON="${PYTHON_BIN}"
cd /dcs04/hicks/data/Jan/sim_project

SIM_QC="sim_paper/data/figure_4/spatial_clustering/sim_data/packing_pf0p04_bin16um_log_mu_-2.5_bsigma07_strongmix/simulation_bin_z_qc.h5ad"
OUT_ROOT="sim_paper/data/figure_3/cellbin_batch_sigma_slide"
LAM=0.5
KG=100

BATCHES=(prebatch tuned)
LABELS=(bs0 bs0.7)
BATCH="${BATCHES[${SLURM_ARRAY_TASK_ID:-0}]}"
LABEL="${LABELS[${SLURM_ARRAY_TASK_ID:-0}]}"

RUN="${OUT_ROOT}/bin/${LABEL}"
mkdir -p "${RUN}"

PREBATCH_FLAG=()
[[ "${BATCH}" == "prebatch" ]] && PREBATCH_FLAG=(--use-pre-batch)

echo "[task ${SLURM_ARRAY_TASK_ID:-0}] bin16um ${BATCH} (${LABEL}) lambda=${LAM} k_geom=${KG}"

BANKSY_H5AD="${RUN}/banksy_matrix/simulation_bin_z_banksy_pca_harmony_qc.h5ad"
if [[ ! -f "${BANKSY_H5AD}" ]]; then
    "${BANKSY_PYTHON}" sim_paper/code/clustering/01_build_banksy_matrix.py \
        --modality bin --input "${SIM_QC}" "${PREBATCH_FLAG[@]}" \
        --lambda "${LAM}" --k-geom "${KG}" --max-m 1 \
        --stagger-scale 5 --skip-umap \
        --output-dir "${RUN}/banksy_matrix"
fi

"${TUTORIAL_PYTHON}" sim_paper/code/clustering/02_leiden_resolution_sweep.py \
    --modality bin --packing-tag "cellbin_batch_slide_bin_${LABEL}" \
    --input "${BANKSY_H5AD}" --output-dir "${RUN}/ari"

"${TUTORIAL_PYTHON}" sim_paper/code/clustering/misc/composition_recovery.py \
    --h5ad "${RUN}/ari/simulation_bin_z_ari_recovery.h5ad" --tag "bin_${LABEL}" \
    --out-dir "${OUT_ROOT}/scores"

echo "[done] ${LABEL}  ->  ${RUN}"
