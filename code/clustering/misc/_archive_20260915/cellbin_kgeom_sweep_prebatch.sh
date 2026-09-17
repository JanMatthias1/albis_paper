#!/bin/bash
#SBATCH --job-name=cellbin_kgeom_prebatch
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/logs/cellbin_kgeom_prebatch_%A_%a.out
#SBATCH --array=0-23%8
#SBATCH --time=18:00:00
#SBATCH --mem=320G
#SBATCH --cpus-per-task=8
#SBATCH --partition=shared
#
# Phase 1 (no Harmony, no batch): find cell/bin16um's own ceiling lambda x
# k_geom config on domain_true, decoupled from the batch/Harmony question.
# We already have ONE point here (lambda=0.5/k_geom=200 -> cell 0.777,
# bin16um 0.790, from job 35638436 reusing job 35608415's prebatch build) --
# this explores around and above it, since job 35638963's k_geom<=15 sweep
# (tuned batch, wrong regime for this question) already showed small k_geom
# is bad for these dense modalities. Feeds Phase 2: sliding batch_sigma at
# each modality's own best config here, with Harmony on, to find the
# recoverability threshold (see figure.md / project_figure3_banksy_domain_sweep
# memory, 2026-09-11).
set -euo pipefail

mkdir -p /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/logs
source /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/_env_banksy.sh
BANKSY_PYTHON="${PYTHON_BIN}"
source /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/_env.sh
TUTORIAL_PYTHON="${PYTHON_BIN}"
cd /dcs04/hicks/data/Jan/sim_project

SIM_DATA_DIR="sim_paper/data/figure_4/spatial_clustering/sim_data"
OUT_ROOT="sim_paper/data/figure_3/cellbin_kgeom_sweep_prebatch"

# 12 (lambda, k_geom) combos x 2 modalities = 24 tasks
LAMBDAS=(0.3 0.3 0.3 0.3 0.5 0.5 0.5 0.5 0.7 0.7 0.7 0.7)
KGEOMS=(100 200 300 400 100 200 300 400 100 200 300 400)
MODS=(cell bin)
TAGS=(
  "log_mu_-2.3_theta_0.40_jitter0.15_bsigma15_strongmix"
  "packing_pf0p04_bin16um_log_mu_-2.5_bsigma07_strongmix"
)

i="${SLURM_ARRAY_TASK_ID:-0}"
COMBO_IDX=$(( i % 12 ))
MOD_IDX=$(( i / 12 ))
LAM="${LAMBDAS[$COMBO_IDX]}"
KG="${KGEOMS[$COMBO_IDX]}"
MOD="${MODS[$MOD_IDX]}"
TAG_DATA="${TAGS[$MOD_IDX]}"

SIM_QC="${SIM_DATA_DIR}/${TAG_DATA}/simulation_${MOD}_z_qc.h5ad"
TAG="lambda${LAM}_kgeom${KG}"
RUN="${OUT_ROOT}/${MOD}/${TAG}"

echo "[task ${i}] modality=${MOD} lambda=${LAM} k_geom=${KG} (PREBATCH, no Harmony scoring)"
mkdir -p "${RUN}"

BANKSY_H5AD="${RUN}/banksy_matrix/simulation_${MOD}_z_banksy_pca_harmony_qc.h5ad"
if [[ ! -f "${BANKSY_H5AD}" ]]; then
    "${BANKSY_PYTHON}" sim_paper/code/clustering/01_build_banksy_matrix.py \
        --modality "${MOD}" --input "${SIM_QC}" --use-pre-batch \
        --lambda "${LAM}" --k-geom "${KG}" --max-m 1 \
        --stagger-scale 5 --skip-umap \
        --output-dir "${RUN}/banksy_matrix"
fi

"${TUTORIAL_PYTHON}" sim_paper/code/clustering/misc/check_banksy_no_harmony_ari.py \
    --modality "${MOD}" --input "${BANKSY_H5AD}" \
    > "${RUN}/no_harmony_ari.log" 2>&1
cat "${RUN}/no_harmony_ari.log"

echo "[done] task ${i}  ->  ${RUN}"
