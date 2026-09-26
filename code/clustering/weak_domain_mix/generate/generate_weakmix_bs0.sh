#!/bin/bash
#SBATCH --job-name=weakmix_bs0_gen
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/data/figure_3/weak_domain_mix/logs/weakmix_bs0_gen_%A_%a.out
#SBATCH --array=0-2
#SBATCH --time=06:00:00
#SBATCH --mem=96G
#SBATCH --cpus-per-task=4
#SBATCH --partition=shared
#
# Weak (Figure 2) domain-mix data with NO batch effect, for cell / bin16um / spot:
# the exact Figure 2 generator flags with --batch-sigma 0 (no --strong-domain-mix),
# so vs. Figure 2 only batch differs, and vs. the strong-mix data only the domain
# mix differs. Only batch_sigma 0 lives here; the batched weak-mix data is
# Figure 2's own. Generates + QC-filters only -- ../pca_harmony/bs0_celltype_panel.sh
# and ../banksy/ cluster it.
#
# History: written 2026-09-24 in data/banksy_cell/ to test whether BANKSY's
# slice collapse on Figure 2 data is entirely the batch effect (it is). Moved
# here 2026-09-25; its outputs were moved, not regenerated (jobs 35904914 +
# 35910892). It used to run PCA+Harmony too; that step now lives in
# ../pca_harmony/bs0_celltype_panel.sh.
#
# Outputs: data/figure_3/weak_domain_mix/data/<modality>/{<stem>.h5ad, <stem>_qc.h5ad, plots/}
set -euo pipefail

ROOT=/dcs04/hicks/data/Jan/sim_project/sim_paper/data/figure_3/weak_domain_mix
mkdir -p "${ROOT}/logs"
source /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/_env.sh
cd /dcs04/hicks/data/Jan/sim_project

case "${SLURM_ARRAY_TASK_ID:?run as an array task}" in
    0) NAME=cell; MOD=cell
       STEM="log_mu_-2.5_theta_0.40_jitter0.15_bsigma0"
       FLAGS=(--sphere-r-um 2050 --n-cells 24207 --base-gene-lognormal -2.5 0.7
              --theta 0.40 --theta-jitter 0.15) ;;
    1) NAME=bin16um; MOD=bin
       STEM="packing_pf0p04_bin16um_log_mu_-2.5_jitter0.6_bsigma0"
       FLAGS=(--sphere-r-um 2050 --bin-size-um 16 --base-gene-lognormal -2.5 0.7
              --theta 2.0 --theta-jitter 0.6) ;;
    2) NAME=spot; MOD=spot
       STEM="packing_pf0p04_log_mu_-2.25_sigma1.0_theta_0.25_jitter0.10_bsigma0"
       FLAGS=(--sphere-r-um 2050 --base-gene-lognormal -2.25 1.0
              --theta 0.25 --theta-jitter 0.10) ;;
esac
RUN="${ROOT}/data/${NAME}"
SIM_RAW="${RUN}/${STEM}.h5ad"
SIM_QC="${RUN}/${STEM}_qc.h5ad"
mkdir -p "${RUN}"
echo "[task ${SLURM_ARRAY_TASK_ID}] ${NAME} weak mix, batch_sigma=0"

if [[ ! -f "${SIM_RAW}" ]]; then
    "${PYTHON_BIN}" sim_paper/code/data/generate_simulation_noisy.py \
        --modality "${MOD}" "${FLAGS[@]}" --batch-sigma 0 --sync-unaligned-seed \
        --output-dir "${RUN}" --output-stem "${STEM}"
fi
if [[ ! -f "${SIM_QC}" ]]; then
    "${PYTHON_BIN}" sim_paper/code/clustering/00_qc_filter.py \
        --modality "${MOD}" --input "${SIM_RAW}" --output "${SIM_QC}"
fi
echo "[done] ${RUN}"
