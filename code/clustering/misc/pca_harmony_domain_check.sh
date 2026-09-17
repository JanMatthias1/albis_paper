#!/bin/bash
#SBATCH --job-name=pca_harmony_domain_check
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/logs/pca_harmony_domain_check_%A_%a.out
#SBATCH --array=0-5%6
#SBATCH --time=06:00:00
#SBATCH --mem=250G
#SBATCH --cpus-per-task=8
#SBATCH --partition=shared
#
# Domain-recovery sanity check requested 2026-09-14: run the SAME
# BANKSY(lambda,k_geom)->Harmony->Leiden domain-recovery pipeline used for
# the whole cellbin_batch_sigma_slide/banksy_batch_compare investigation, but
# directly on the EXACT QC'd h5ad files that back the real Figure 3
# pca_harmony_single_cell panel (cell_celltype_panel.sh /
# bin16um_celltype_panel.sh / spot_celltype_panel.sh) -- no new simulation.
# These are the WEAK (manuscript-baseline) domain_type_mix datasets, at each
# modality's real calibrated dispersion/mean and canonical batch_sigma, NOT
# the separately-generated --strong-domain-mix testbed used everywhere else
# in this investigation. --use-pre-batch reads layers['counts_pre_batch']
# from the same file, so prebatch+tuned both come from one h5ad, no fresh
# simulation needed either way.
#
# Motivation: spot's strong-domain-mix generation config (generate_strongmix.sh)
# uses log_mu=-2.5, theta=2.0, jitter=1.0 (generate_simulation_noisy.py
# defaults) -- this predates and was NEVER updated after Figure 2's spot
# config was retuned (2026-09-05/06) to log_mu=-2.0, theta=0.25, jitter=0.10.
# Cell and bin16um's strongmix configs still match their celltype-panel
# configs exactly; spot's does not. This check runs on spot's ACTUAL current
# calibrated data instead, to see whether the batch_sigma dip/cliff findings
# so far still hold.
set -euo pipefail

mkdir -p /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/logs
source /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/_env_banksy.sh
BANKSY_PYTHON="${PYTHON_BIN}"
source /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/_env.sh
TUTORIAL_PYTHON="${PYTHON_BIN}"
cd /dcs04/hicks/data/Jan/sim_project

OUT_ROOT="sim_paper/data/figure_3/pca_harmony_domain_figure2_data"

MODS=(cell cell bin16um bin16um spot spot)
BATCHES=(prebatch tuned prebatch tuned prebatch tuned)

i="${SLURM_ARRAY_TASK_ID:-0}"
MOD="${MODS[$i]}"
BATCH="${BATCHES[$i]}"

case "${MOD}" in
  cell)
    SIM_QC="sim_paper/data/figure_2/smaller_sphere/data/log_mu_-2.3_theta_0.40_jitter0.15_bsigma15/simulation_cell_z_qc.h5ad"
    BANKSY_MOD="cell"
    LAM=0.5; KG=200
    ;;
  bin16um)
    SIM_QC="sim_paper/data/figure_2/smaller_sphere/data/packing_pf0p04_bin16um_log_mu_-2.5_bsigma07/simulation_bin_z_qc.h5ad"
    BANKSY_MOD="bin"
    LAM=0.5; KG=100
    ;;
  spot)
    SIM_QC="sim_paper/data/figure_2/smaller_sphere/data/packing_pf0p04_log_mu_-2.0_theta_0.25_jitter0.10_bsigma03/simulation_spot_z_qc.h5ad"
    BANKSY_MOD="spot"
    LAM=0.1; KG=8
    ;;
esac

RUN="${OUT_ROOT}/${MOD}/${BATCH}"
BANKSY_H5AD="${RUN}/banksy_matrix/simulation_${BANKSY_MOD}_z_banksy_pca_harmony_qc.h5ad"

PREBATCH_FLAG=()
[[ "${BATCH}" == "prebatch" ]] && PREBATCH_FLAG=(--use-pre-batch)

echo "[task ${i}] modality=${MOD} batch=${BATCH} lambda=${LAM} k_geom=${KG} input=${SIM_QC}"
mkdir -p "${RUN}" "${OUT_ROOT}/scores"

"${BANKSY_PYTHON}" sim_paper/code/clustering/01_build_banksy_matrix.py \
    --modality "${BANKSY_MOD}" --input "${SIM_QC}" "${PREBATCH_FLAG[@]}" \
    --lambda "${LAM}" --k-geom "${KG}" --max-m 1 \
    --stagger-scale 5 --skip-umap \
    --output-dir "${RUN}/banksy_matrix"

"${TUTORIAL_PYTHON}" sim_paper/code/clustering/02_leiden_resolution_sweep.py \
    --modality "${BANKSY_MOD}" --packing-tag "pca_harmony_domain_check_${MOD}_${BATCH}" \
    --input "${BANKSY_H5AD}" --output-dir "${RUN}/ari"

"${TUTORIAL_PYTHON}" sim_paper/code/clustering/composition_recovery.py \
    --h5ad "${RUN}/ari/simulation_${BANKSY_MOD}_z_ari_recovery.h5ad" --tag "${MOD}_${BATCH}" \
    --out-dir "${OUT_ROOT}/scores"

echo "[done] ${MOD}/${BATCH}  ->  ${RUN}"
