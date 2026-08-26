#!/bin/bash
#SBATCH --job-name=sweep_bin16um_logmu_hvg
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/data/logs/sweep_bin16um_logmu_hvg_%A_%a.out
#SBATCH --array=0-4
#SBATCH --time=02:00:00
#SBATCH --mem=250G
#SBATCH --cpus-per-task=4
#SBATCH --partition=shared
#
# Backfill: add the QC + HVG-matched (--match-panel-size) comparison mode on
# top of the qc_filtered-only sweep already run (sweep_bin16um_logmu.sh +
# _round2.sh, all 5 log_mu points landed). Real Visium HD is whole-transcriptome
# (18,085 genes) vs. sim's 556-gene panel -- total_counts/genes_per_cell need
# panel-matching to be a fair comparison (established convention, see
# bin_vs_breast_cancer_visium_hd.sh and figure.md: only total_counts_compare.png/
# genes_per_cell_compare.png should ever be cited from this mode -- mean_variance
# stays HVG-biased, qc_filtered remains the source for theta_hat/mean_variance).
# Reuses the already-generated simulation_bin_z_qc.h5ad for each log_mu -- no
# regeneration, just count_distribution.py --match-panel-size.

set -euo pipefail

mkdir -p /dcs04/hicks/data/Jan/sim_project/sim_paper/code/data/logs
source /dcs04/hicks/data/Jan/sim_project/sim_paper/code/count_distribution/_env.sh
cd /dcs04/hicks/data/Jan/sim_project

LOG_MUS=(-2.5 -1.5 -0.5 0.0 0.5)
LOG_MU="${LOG_MUS[${SLURM_ARRAY_TASK_ID:-0}]}"
TAG="packing_pf0p04_bin16um_log_mu_${LOG_MU}"
MODALITY="bin"

SIM_QC="sim_paper/data/noisy/${TAG}/simulation_${MODALITY}_z_qc.h5ad"
OUT_ROOT="sim_paper/data/count_distribution/sweeps/bin16um_logmu_sweep"

if [[ ! -f "${SIM_QC}" ]]; then
    echo "ERROR: ${SIM_QC} not found -- expected round 1/round 2 to have already generated this." >&2
    exit 1
fi

for real in breast_cancer_visium_hd human_pancreas_visium_hd; do
    echo "[compare] ${TAG} vs ${real}_16um (qc_and_hvg_matched)"
    "${PYTHON_BIN}" sim_paper/code/count_distribution/count_distribution.py \
        --modality "${MODALITY}" --input "${SIM_QC}" \
        --compare-input "sim_paper/data/real_data_qc/${real}_16um/${real}_qc.h5ad" \
        --compare-label "${real}_16um" \
        --match-panel-size \
        --output-dir "${OUT_ROOT}/${TAG}_qc_and_hvg_matched_vs_${real}_16um"
done

echo "[done] ${TAG}"
