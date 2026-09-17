#!/bin/bash
#SBATCH --job-name=spot_cmp_corrected
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/logs/spot_cmp_corrected_%A_%a.out
#SBATCH --array=0-1%2
#SBATCH --time=01:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --partition=shared
#
# Rebuild of the OFFICIAL banksy_batch_compare/spot/{prebatch,tuned} pair
# using spot's corrected dispersion (log_mu=-2.0, theta=0.25, jitter=0.10 --
# see spot_batch_sigma_slide_corrected.sh for why). Same lambda=0.1/k_geom=8
# as before. Reuses the bs=0.3 (canonical) simulation that
# spot_batch_sigma_slide_corrected.sh's array task 15 already builds --
# submit this with --dependency=afterok:<that job> so the QC'd file is
# guaranteed to exist; if it's somehow missing this will generate it fresh
# (same tag, so no duplicate work either way).
# Old stale-param official run archived at
# data/figure_3/_archive_20260914/banksy_batch_compare_spot_staleparams_20260914/
# (+ its score jsons in the _scores/ sibling dir) before this ran.
set -euo pipefail

mkdir -p /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/logs
source /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/_env_banksy.sh
BANKSY_PYTHON="${PYTHON_BIN}"
source /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/_env.sh
TUTORIAL_PYTHON="${PYTHON_BIN}"
cd /dcs04/hicks/data/Jan/sim_project

TAG="spot_batch_slide_corrected_bs0.3"
SIM_RAW="sim_paper/data/noisy/${TAG}/simulation_spot_z.h5ad"
SIM_QC="sim_paper/data/noisy/${TAG}/simulation_spot_z_qc.h5ad"

if [[ ! -f "${SIM_RAW}" ]]; then
    "${BANKSY_PYTHON}" sim_paper/code/data/generate_simulation_noisy.py \
        --modality spot --sphere-r-um 2050 \
        --base-gene-lognormal -2.0 0.7 --theta 0.25 --theta-jitter 0.10 \
        --strong-domain-mix --batch-sigma 0.3 \
        --out-tag "${TAG}"
fi
if [[ ! -f "${SIM_QC}" ]]; then
    "${TUTORIAL_PYTHON}" sim_paper/code/clustering/00_qc_filter.py \
        --modality spot --input "${SIM_RAW}" --output "${SIM_QC}"
fi

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

echo "[task] MOD=spot batch=${BATCH} lambda=${LAM} k_geom=${KG} (corrected dispersion)"
mkdir -p "${RUN}" "${SWEEP}/scores"

"${BANKSY_PYTHON}" sim_paper/code/clustering/01_build_banksy_matrix.py \
    --modality spot --input "${SIM_QC}" "${PREBATCH_FLAG[@]}" \
    --lambda "${LAM}" --k-geom "${KG}" --max-m 1 \
    --stagger-scale 5 --skip-umap \
    --output-dir "${RUN}/banksy_matrix"

"${TUTORIAL_PYTHON}" sim_paper/code/clustering/02_leiden_resolution_sweep.py \
    --modality spot --packing-tag "spot_corrected_${BATCH}" \
    --input "${BANKSY_H5AD}" --output-dir "${RUN}/ari"

"${TUTORIAL_PYTHON}" sim_paper/code/clustering/composition_recovery.py \
    --h5ad "${ARI_H5AD}" --tag "spot_${BATCH}" \
    --out-dir "${SWEEP}/scores"

echo "[done] ${BATCH}  ->  ${RUN}"
