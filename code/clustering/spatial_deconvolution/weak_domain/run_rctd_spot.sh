#!/usr/bin/env bash
# RCTD spatial deconvolution of Figure 2's spot data (see run_rctd_spot.R
# for the full rationale: cell -> reference, spot -> query, scored against
# spot's obsm['cell_type_frac_true']).
#SBATCH --job-name=rctd_spot
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/albis_paper/code/clustering/spatial_deconvolution/logs/rctd_spot_%j.out
#SBATCH --time=04:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --partition=shared
set -euo pipefail

mkdir -p /dcs04/hicks/data/Jan/sim_project/albis_paper/code/clustering/spatial_deconvolution/logs

RCTD_ENV="/dcs04/hicks/data/Jan/sim_project/albis_paper/env/rctd"
export LD_LIBRARY_PATH="${RCTD_ENV}/lib:${LD_LIBRARY_PATH:-}"
export R_HOME="${RCTD_ENV}/lib/R"
export RETICULATE_PYTHON="${RCTD_ENV}/bin/python"

"${RCTD_ENV}/bin/Rscript" \
    /dcs04/hicks/data/Jan/sim_project/albis_paper/code/clustering/spatial_deconvolution/weak_domain/run_rctd_spot.R

# Produce the manuscript panels after fitting; no manual style-refresh step.
/dcs04/hicks/data/Jan/sim_project/albis_paper/env/albis-tutorial/bin/python \
    /dcs04/hicks/data/Jan/sim_project/albis_paper/code/clustering/spatial_deconvolution/plot_rctd_results.py \
    "${RCTD_OUT_DIR:-/dcs04/hicks/data/Jan/sim_project/albis_paper/data/figure_3/spatial_deconvolution/RCTD/weak_mix}"
