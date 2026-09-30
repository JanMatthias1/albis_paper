#!/usr/bin/env bash
#SBATCH --job-name=rctd_pilot_plots
#SBATCH --partition=shared
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
set -euo pipefail
project=/dcs04/hicks/data/Jan/sim_project
export MPLBACKEND=Agg
"$project/comparison_methods/env/analysis/bin/python" "$project/sim_paper/code/comparison_methods/clustering/spatial_deconvolution/plot_pilot.py" "${1:?Supply pilot run directory}"
