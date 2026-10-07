#!/bin/bash
#SBATCH --job-name=fig2_replot_qc_modes
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/albis_paper/code/count_distribution/smaller_sphere/logs/fig2_replot_qc_modes_%j.out
#SBATCH --time=04:00:00
#SBATCH --mem=96G
#SBATCH --cpus-per-task=4
#SBATCH --partition=shared
#
# Replot only the manuscript modes (qc_filtered, qc_and_hvg_matched) of all 8 smaller_sphere
# pairs from existing QC'd h5ads -- no regeneration, no full_panel/hvg_matched (those need the
# pre-QC matrices and 250G). Same commands as the per-pair scripts; use after plot-only changes
# to count_distribution.py.
set -euo pipefail
source /dcs04/hicks/data/Jan/sim_project/albis_paper/code/count_distribution/_env.sh
cd /dcs04/hicks/data/Jan/sim_project

DATA="albis_paper/data/figure_2/smaller_sphere/data"
PLOTS="albis_paper/data/figure_2/smaller_sphere/plots"
REAL="albis_paper/data/real_data_qc"
BIN16="packing_pf0p04_bin16um_log_mu_-2.5_jitter0.6_bsigma07"
BIN8="packing_pf0p04_bin8um_from16umcfg_log_mu_-2.5_jitter0.6_bsigma07"
CELL="log_mu_-2.5_theta_0.40_jitter0.15_bsigma15"
SPOT="packing_pf0p04_log_mu_-2.25_sigma1.0_theta_0.25_jitter0.10_bsigma03"

# modality | sim tag | real label | real h5ad (under $REAL) | plot folder
PAIRS=(
    "bin|${BIN16}|breast_cancer_visium_hd_16um|breast_cancer_visium_hd_16um/breast_cancer_visium_hd_qc.h5ad|bin_vs_breast_cancer_visium_hd_16um"
    "bin|${BIN16}|human_pancreas_visium_hd_16um|human_pancreas_visium_hd_16um/human_pancreas_visium_hd_qc.h5ad|bin_vs_human_pancreas_visium_hd_16um"
    "bin|${BIN8}|breast_cancer_visium_hd|breast_cancer_visium_hd/breast_cancer_visium_hd_qc.h5ad|bin8um_same_tissue_as_16um_vs_breast_cancer_visium_hd"
    "bin|${BIN8}|human_pancreas_visium_hd|human_pancreas_visium_hd/human_pancreas_visium_hd_qc.h5ad|bin8um_same_tissue_as_16um_vs_human_pancreas_visium_hd"
    "cell|${CELL}|lung_cancer|lung_cancer/lung_cancer_qc.h5ad|cell_vs_lung_cancer"
    "cell|${CELL}|non_diseased_lung|non_diseased_lung/non_diseased_lung_qc.h5ad|cell_vs_non_diseased_lung"
    "spot|${SPOT}|lymph_node_visium|lymph_node_visium/lymph_node_visium_qc.h5ad|spot_vs_lymph_node_visium"
    "spot|${SPOT}|tonsil_visium|tonsil_visium/tonsil_visium_qc.h5ad|spot_vs_tonsil_visium"
)

for pair in "${PAIRS[@]}"; do
    IFS="|" read -r MODALITY TAG REAL_LABEL REAL_FILE OUT <<< "${pair}"
    SIM_QC="${DATA}/${TAG}/simulation_${MODALITY}_z_qc.h5ad"
    "${PYTHON_BIN}" albis_paper/code/count_distribution/count_distribution.py \
        --modality "${MODALITY}" --input "${SIM_QC}" \
        --compare-input "${REAL}/${REAL_FILE}" --compare-label "${REAL_LABEL}" \
        --output-dir "${PLOTS}/${OUT}/qc_filtered"
    "${PYTHON_BIN}" albis_paper/code/count_distribution/count_distribution.py \
        --modality "${MODALITY}" --input "${SIM_QC}" \
        --compare-input "${REAL}/${REAL_FILE}" --compare-label "${REAL_LABEL}" \
        --match-panel-size \
        --output-dir "${PLOTS}/${OUT}/qc_and_hvg_matched"
    echo "[done] ${OUT}"
done
