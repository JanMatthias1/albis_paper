#!/usr/bin/env bash
#SBATCH --job-name=fig5c_final_summary
#SBATCH --partition=shared
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:15:00
set -euo pipefail
export MPLBACKEND=Agg MPLCONFIGDIR=/tmp/fig5c_final_summary
/dcs04/hicks/data/Jan/sim_project/sim_paper/env/albis-tutorial/bin/python /dcs04/hicks/data/Jan/sim_project/sim_paper/code/comparison_methods/clustering/cell_type/aggregate_final_600k.py
