#!/bin/bash
#SBATCH --job-name=cell_ari_lm23_probe
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/logs/cell_ari_lm23_probe_%A_%a.out
#SBATCH --array=0-2
#SBATCH --time=06:00:00
#SBATCH --mem=250G
#SBATCH --cpus-per-task=4
#SBATCH --partition=shared
#
# One-off probe (2026-08-27): measure cell_type_true ARI for the round-2
# batch_sigma_opt cell grid points that were only ever scored on
# count_distribution, never clustered. Question: does pre-compensating cell's
# dispersion (log_mu -2.3, theta 0.40) recover the cell_type_true ARI that
# batch_sigma=1.5 costs (0.47 at the adopted lm-2.5/th0.25 config vs 0.54
# pre-bump), while also giving the better count-distribution match those
# grid points already showed (theta_hat 0.082 vs 0.056, genes/cell 42 vs 30)?
#
# Reuses the QC'd h5ads already on disk under data/noisy/batch_opt2_cell_*/.
# pca_harmony.py -> ari_vs_ground_truth.py only (no leiden qualitative plots).
# Baselines already known: lm-2.5/th0.25/bs1.5 = 0.47, archived bs0.22 = 0.54.

set -euo pipefail
source /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/_env.sh
cd /dcs04/hicks/data/Jan/sim_project

TAGS=(batch_opt2_cell_lm-2.3_th0.40_bs1.0 batch_opt2_cell_lm-2.3_th0.40_bs1.5 batch_opt2_cell_lm-2.3_th0.25_bs1.5)
TAG="${TAGS[${SLURM_ARRAY_TASK_ID:-0}]}"
MOD=cell
SIM_QC="sim_paper/data/noisy/${TAG}/simulation_${MOD}_z_qc.h5ad"
OUT="sim_paper/data/figure_3/batch_sigma_sweep/${TAG}"
PCA_OUT="${OUT}/pca_harmony_qc/simulation_${MOD}_z_pca_harmony_qc.h5ad"

echo "[probe] tag=${TAG}  qc=${SIM_QC}"
[[ -f "${SIM_QC}" ]] || { echo "missing ${SIM_QC}" >&2; exit 1; }

if [[ ! -f "${PCA_OUT}" ]]; then
    echo "[pca_harmony]"
    "${PYTHON_BIN}" sim_paper/code/clustering/pca_harmony.py \
        --modality "${MOD}" --input "${SIM_QC}" --output "${PCA_OUT}"
fi

echo "[ari]"
"${PYTHON_BIN}" sim_paper/code/clustering/ari_vs_ground_truth.py \
    --modality "${MOD}" --packing-tag "${TAG}" \
    --input "${PCA_OUT}" --output-dir "${OUT}/ari_recovery_qc"

echo "[done] ${OUT}/ari_recovery_qc/ari_summary_${MOD}.json"
