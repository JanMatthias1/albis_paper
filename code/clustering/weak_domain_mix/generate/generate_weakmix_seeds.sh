#!/bin/bash
#SBATCH --job-name=weakmix_seed_gen
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/data/figure_3/weak_domain_mix/logs/weakmix_seed_gen_%A_%a.out
#SBATCH --array=0-11
#SBATCH --time=08:00:00
#SBATCH --mem=120G
#SBATCH --cpus-per-task=4
#SBATCH --partition=shared
#
# Extra simulation seeds (101, 202) of the weak (Figure 2) domain-mix data for
# cell / bin16um / spot, so every Figure 3 bar can show mean ± SD over seeds
# 2025/101/202. Seed 2025 already exists and is not touched.
# Flags are the seed-2025 flags exactly (checked against each file's
# uns['sim_params']); only --seed differs. bin8um is left at one seed.
#
#   usual batch_sigma (cell 1.5 / bin16um 0.7 / spot 0.3), next to Figure 2's
#     data/figure_2/smaller_sphere/data/<Figure 2 tag>_seed<n>/simulation_<mod>_z{,_qc}.h5ad
#   batch_sigma 0, next to generate_weakmix_bs0.sh's seed-2025 files
#     data/figure_3/weak_domain_mix/data/<modality>_seed<n>/<stem>_seed<n>{,_qc}.h5ad
#
# Task = modality (3) x batch (usual, 0) x seed (101, 202). Existing outputs are
# skipped. DRY_RUN=1 prints the commands instead of running them.
set -euo pipefail

source /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/_env.sh
cd /dcs04/hicks/data/Jan/sim_project
RUN=(); [[ -n "${DRY_RUN:-}" ]] && RUN=(echo)

TASK="${SLURM_ARRAY_TASK_ID:?run as an array task}"
MODS=(cell bin16um spot); BATCHES=(usual 0); SEEDS=(101 202)
NAME="${MODS[$((TASK / 4))]}"; BATCH="${BATCHES[$(((TASK / 2) % 2))]}"; SEED="${SEEDS[$((TASK % 2))]}"

case "${NAME}" in
    cell) MOD=cell; SIGMA=1.5
          TAG="log_mu_-2.5_theta_0.40_jitter0.15_bsigma15"
          STEM0="log_mu_-2.5_theta_0.40_jitter0.15_bsigma0"
          FLAGS=(--sphere-r-um 2050 --n-cells 24207 --base-gene-lognormal -2.5 0.7
                 --theta 0.40 --theta-jitter 0.15) ;;
    bin16um) MOD=bin; SIGMA=0.7
          TAG="packing_pf0p04_bin16um_log_mu_-2.5_jitter0.6_bsigma07"
          STEM0="packing_pf0p04_bin16um_log_mu_-2.5_jitter0.6_bsigma0"
          FLAGS=(--sphere-r-um 2050 --bin-size-um 16 --base-gene-lognormal -2.5 0.7
                 --theta 2.0 --theta-jitter 0.6) ;;
    spot) MOD=spot; SIGMA=0.3
          TAG="packing_pf0p04_log_mu_-2.25_sigma1.0_theta_0.25_jitter0.10_bsigma03"
          STEM0="packing_pf0p04_log_mu_-2.25_sigma1.0_theta_0.25_jitter0.10_bsigma0"
          FLAGS=(--sphere-r-um 2050 --base-gene-lognormal -2.25 1.0
                 --theta 0.25 --theta-jitter 0.10) ;;
esac

if [[ "${BATCH}" == usual ]]; then
    OUT="sim_paper/data/figure_2/smaller_sphere/data/${TAG}_seed${SEED}"
    STEM="simulation_${MOD}_z"
    BS="${SIGMA}"
else
    OUT="sim_paper/data/figure_3/weak_domain_mix/data/${NAME}_seed${SEED}"
    STEM="${STEM0}_seed${SEED}"
    BS=0
fi
SIM_RAW="${OUT}/${STEM}.h5ad"
SIM_QC="${OUT}/${STEM}_qc.h5ad"
echo "[task ${TASK}] ${NAME} weak mix, batch_sigma=${BS}, seed=${SEED} -> ${OUT}"
"${RUN[@]}" mkdir -p "${OUT}"

if [[ ! -f "${SIM_RAW}" ]]; then
    "${RUN[@]}" "${PYTHON_BIN}" sim_paper/code/data/generate_simulation_noisy.py \
        --modality "${MOD}" "${FLAGS[@]}" --batch-sigma "${BS}" --sync-unaligned-seed \
        --seed "${SEED}" --output-dir "${OUT}" --output-stem "${STEM}"
fi
if [[ ! -f "${SIM_QC}" ]]; then
    "${RUN[@]}" "${PYTHON_BIN}" sim_paper/code/clustering/step00_qc_filter.py \
        --modality "${MOD}" --input "${SIM_RAW}" --output "${SIM_QC}"
fi
echo "[done] ${OUT}"
