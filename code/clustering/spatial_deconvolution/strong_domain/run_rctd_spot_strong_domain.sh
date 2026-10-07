#!/usr/bin/env bash
# RCTD spot deconvolution on the STRONG-domain-mix data -- the counterpart of
# ../weak_domain/run_rctd_spot{,_independent_seed}.sh (Figure 2's weak mix).
# Same fitting + scoring code (../weak_domain/run_rctd_spot.R, full mode),
# only the inputs change:
#   query (spot):     strong-mix spot at its canonical batch_sigma 0.3, seed 2025
#   reference (cell): strong-mix cell at its canonical batch_sigma 1.5, one
#                     array task per reference seed:
#     task 0  seed 2025  cell/bs1.5/          -> RCTD/strong_mix/
#     task 1  seed 101   cell/bs1.5_seed101/  -> RCTD/strong_mix_seed101/
#     task 2  seed 202   cell/bs1.5_seed202/  -> RCTD/strong_mix_seed202/
# then average_rctd_seeds.py --config strong_mix averages all 3 into
# RCTD/strong_mix_seed_avg/ (seed 2025 shares the query's seed; the weak-mix
# replication showed this does not matter).
# All inputs are Figure 2 configs + --strong-domain-mix and nothing else
# (strong_domain_mix/generate/check_matches_figure2.py), written into
# data/figure_3/strong_domain_mix/batch_sigma_slide/ by strong_domain_mix/generate/generate_strong_mix_
# {cell,spot}.sh (seed 2025) and strong_domain_mix/batch_sigma_slide/
# batch_sigma_slide_seeds.sh (seeds 101/202) -- so weak vs strong RCTD differ
# ONLY in how strongly each domain is dominated by 1-2 cell types.
#
# run_figure3_global.sh submits this after those jobs; run standalone only
# once they are done.
#SBATCH --job-name=rctd_spot_strong
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/albis_paper/code/clustering/spatial_deconvolution/logs/rctd_spot_strong_%A_%a.out
#SBATCH --array=0-2
#SBATCH --time=04:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --partition=shared
set -euo pipefail

mkdir -p /dcs04/hicks/data/Jan/sim_project/albis_paper/code/clustering/spatial_deconvolution/logs

SLIDE=/dcs04/hicks/data/Jan/sim_project/albis_paper/data/figure_3/strong_domain_mix/batch_sigma_slide
RCTD_ROOT=/dcs04/hicks/data/Jan/sim_project/albis_paper/data/figure_3/spatial_deconvolution/RCTD
CELL_BODY="log_mu_-2.5_theta_0.40_jitter0.15_strongmix_bsigma15"

REF_SEEDS=("" 101 202)   # "" = default seed 2025 (not written in names)
SEED="${REF_SEEDS[${SLURM_ARRAY_TASK_ID:-0}]}"
export RCTD_CELL_H5AD="${SLIDE}/cell/bs1.5${SEED:+_seed${SEED}}/${CELL_BODY}${SEED:+_seed${SEED}}_qc.h5ad"
export RCTD_SPOT_H5AD="${SLIDE}/spot/bs0.3/packing_pf0p04_log_mu_-2.25_sigma1.0_theta_0.25_jitter0.10_strongmix_bsigma03_qc.h5ad"
export RCTD_OUT_DIR="${RCTD_ROOT}/strong_mix${SEED:+_seed${SEED}}"

for f in "${RCTD_CELL_H5AD}" "${RCTD_SPOT_H5AD}"; do
    [[ -f "${f}" ]] || { echo "[error] missing input ${f} -- run the strong-mix generators first" >&2; exit 1; }
done

RCTD_ENV="/dcs04/hicks/data/Jan/sim_project/albis_paper/env/rctd"
export LD_LIBRARY_PATH="${RCTD_ENV}/lib:${LD_LIBRARY_PATH:-}"
export R_HOME="${RCTD_ENV}/lib/R"
export RETICULATE_PYTHON="${RCTD_ENV}/bin/python"

"${RCTD_ENV}/bin/Rscript" \
    /dcs04/hicks/data/Jan/sim_project/albis_paper/code/clustering/spatial_deconvolution/weak_domain/run_rctd_spot.R

echo "[done] RCTD (strong domain mix, reference seed ${SEED:-2025}) -> ${RCTD_OUT_DIR}"

# Produce the manuscript panels after fitting; the plotter reads which spot
# truth to use from this run's metrics_summary.json.
/dcs04/hicks/data/Jan/sim_project/albis_paper/env/albis-tutorial/bin/python \
    /dcs04/hicks/data/Jan/sim_project/albis_paper/code/clustering/spatial_deconvolution/plot_rctd_results.py \
    "${RCTD_OUT_DIR}"
