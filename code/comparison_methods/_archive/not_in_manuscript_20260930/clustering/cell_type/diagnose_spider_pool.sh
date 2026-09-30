#!/usr/bin/env bash
#SBATCH --job-name=spider_pool_check
#SBATCH --partition=shared
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=02:00:00
set -euo pipefail
export OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4 NUMBA_NUM_THREADS=4 PYTHONDONTWRITEBYTECODE=1
/dcs04/hicks/data/Jan/sim_project/sim_paper/env/albis-tutorial/bin/python -u /dcs04/hicks/data/Jan/sim_project/sim_paper/code/comparison_methods/clustering/cell_type/diagnose_spider_pool.py --task "$SLURM_ARRAY_TASK_ID"
