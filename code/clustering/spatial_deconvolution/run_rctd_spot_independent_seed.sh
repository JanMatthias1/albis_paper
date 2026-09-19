#!/usr/bin/env bash
# RCTD spot deconvolution, rerun against the independently-seeded cell
# reference from generate_rctd_reference_independent_seed.sh, to verify the
# canonical run's result (overall Pearson r=0.730) isn't an artifact of
# reference/query sharing a template via the nominal (but, per that script's
# header, likely already-decorrelated) shared seed=2025.
#
# Writes to a SEPARATE output dir (RCTD/spot_independent_seed/) so the
# canonical RCTD/spot/ results are untouched -- compare metrics_summary.json
# between the two directly.
#
# Array over the same 3 seeds as generate_rctd_reference_independent_seed.sh
# -- submit with --dependency=aftercorr:<that array job's ID> so each RCTD
# task waits only on its own matching reference-generation task, not the
# whole array.
#SBATCH --job-name=rctd_spot_indep
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/spatial_deconvolution/logs/rctd_spot_indep_%A_%a.out
#SBATCH --array=0-2
#SBATCH --time=04:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --partition=shared
set -euo pipefail

mkdir -p /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/spatial_deconvolution/logs

SEEDS=(999999 314159 271828)
SEED="${SEEDS[${SLURM_ARRAY_TASK_ID:-0}]}"
export RCTD_CELL_H5AD="/dcs04/hicks/data/Jan/sim_project/sim_paper/data/figure_2/smaller_sphere/data/log_mu_-2.5_theta_0.40_jitter0.15_bsigma15_rctdref_seed${SEED}/simulation_cell_z_qc.h5ad"
export RCTD_OUT_DIR="/dcs04/hicks/data/Jan/sim_project/sim_paper/data/figure_3/spatial_deconvolution/RCTD/spot_seed${SEED}"

RCTD_ENV="/dcs04/hicks/data/Jan/sim_project/sim_paper/env/rctd"
export LD_LIBRARY_PATH="${RCTD_ENV}/lib:${LD_LIBRARY_PATH:-}"
export R_HOME="${RCTD_ENV}/lib/R"
export RETICULATE_PYTHON="${RCTD_ENV}/bin/python"

"${RCTD_ENV}/bin/Rscript" \
    /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/spatial_deconvolution/run_rctd_spot.R

echo "[done] RCTD (independent-seed reference) -> ${RCTD_OUT_DIR}"
