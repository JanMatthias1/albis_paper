#!/bin/bash
#SBATCH --job-name=spot_batch_cmp_kg8
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/logs/spot_batch_cmp_kg8_%A_%a.out
#SBATCH --array=0-1%2
#SBATCH --time=01:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --partition=shared
#
# Rebuild of banksy_batch_compare/spot/{prebatch,tuned} at lambda=0.1/
# k_geom=8 -- the winner of spot_kgeom_sweep_fine.sh's finer sweep around
# the earlier lambda=0.2/k_geom=6 result (domain ARI 0.396, slice_leak_ARI
# 0.013 at real tuned batch, vs 0.271/0.236 for the coarser winner -- both
# better ARI and much cleaner leakage). Old lambda=0.2/k_geom=6 run archived
# at banksy_batch_compare/spot_pre_lambda0.1_fix_20260911/; the original
# borrowed lambda=0.5/k_geom=200 run before that is at
# banksy_batch_compare/spot_pre_kgeom_fix_20260911/.
set -euo pipefail

mkdir -p /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/logs
source /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/_env_banksy.sh
BANKSY_PYTHON="${PYTHON_BIN}"
source /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/_env.sh
TUTORIAL_PYTHON="${PYTHON_BIN}"
cd /dcs04/hicks/data/Jan/sim_project

SIM_QC="sim_paper/data/figure_4/spatial_clustering/sim_data/packing_pf0p04_log_mu_-2.5_bsigma03_strongmix/simulation_spot_z_qc.h5ad"
LAM=0.1
KG=8
BATCHES=(prebatch tuned)
BATCH="${BATCHES[${SLURM_ARRAY_TASK_ID:-0}]}"

SWEEP="sim_paper/data/figure_3/banksy_batch_compare"
RUN="${SWEEP}/spot/${BATCH}"
BANKSY_H5AD="${RUN}/banksy_matrix/simulation_spot_z_banksy_pca_harmony_qc.h5ad"
ARI_H5AD="${RUN}/ari/simulation_spot_z_ari_recovery.h5ad"

PREBATCH_FLAG=()
[[ "${BATCH}" == "prebatch" ]] && PREBATCH_FLAG=(--use-pre-batch)

echo "[task] MOD=spot batch=${BATCH} lambda=${LAM} k_geom=${KG}"
mkdir -p "${RUN}" "${SWEEP}/scores"

"${BANKSY_PYTHON}" sim_paper/code/clustering/01_build_banksy_matrix.py \
    --modality spot --input "${SIM_QC}" "${PREBATCH_FLAG[@]}" \
    --lambda "${LAM}" --k-geom "${KG}" --max-m 1 \
    --stagger-scale 5 --skip-umap \
    --output-dir "${RUN}/banksy_matrix"

"${TUTORIAL_PYTHON}" sim_paper/code/clustering/ari_vs_ground_truth.py \
    --modality spot --packing-tag "spot_kgeom6_${BATCH}" \
    --input "${BANKSY_H5AD}" --output-dir "${RUN}/ari"

"${TUTORIAL_PYTHON}" sim_paper/code/clustering/composition_recovery.py \
    --h5ad "${ARI_H5AD}" --tag "spot_${BATCH}" \
    --out-dir "${SWEEP}/scores"

echo "[done] ${BATCH}  ->  ${RUN}"
