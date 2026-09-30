#!/usr/bin/env bash
#SBATCH --job-name=spider_pool_summary
#SBATCH --partition=shared
#SBATCH --mem=1G
#SBATCH --time=00:05:00
set -euo pipefail
/dcs04/hicks/data/Jan/sim_project/sim_paper/env/albis-tutorial/bin/python /dcs04/hicks/data/Jan/sim_project/sim_paper/code/comparison_methods/clustering/cell_type/summarize_spider_pool.py
