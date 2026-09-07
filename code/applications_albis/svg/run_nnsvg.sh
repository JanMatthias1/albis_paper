#!/bin/bash
#SBATCH --job-name=nnsvg_svg
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/applications_albis/svg/logs_nnsvg/nnsvg_%A_%a.out
#SBATCH --array=0-2
#SBATCH --time=03:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --partition=shared
#
# Figure 4A -- nnSVG SVG identification, one array task per modality, single
# representative slice (nnSVG has no 3D mode -- see 3D_nnsvg.py docstring).
# Cost is much higher than scBSP/SPARK-X; nnSVG parallelises across genes, so
# wall time is ~ (n_genes / n_threads) * per-gene BRISC fit (~10s/gene at 50k
# obs). Empirical single-slice, n_threads=4, 2026-09-05: spot ~130s, cell
# ~28min. bin16um's first attempt ran ~4x slow on a contended shared node
# (killed at 1h48m); bumped to 8 cpus + one BiocParallel task per gene with a
# progress bar (see run_nnsvg.R) so a slow run is at least visible.
#
# Extra args are forwarded to 3D_nnsvg.py, so:
#   sbatch run_nnsvg.sh                          # canonical weak-mix, post-batch
#   sbatch run_nnsvg.sh --strongmix              # strong-domain-mix, post-batch
#   sbatch run_nnsvg.sh --strongmix --use-pre-batch   # strongmix oracle
# Output dirs are data/figure_4/svg/nnsvg/<dataset>[_strongmix][_prebatch]/.

set -euo pipefail

mkdir -p /dcs04/hicks/data/Jan/sim_project/sim_paper/code/applications_albis/svg/logs_nnsvg

ENV_PREFIX="/dcs04/hicks/data/Jan/sim_project/sim_paper/env/stagate-pyg"
PYTHON_BIN="${ENV_PREFIX}/bin/python"
export PYTHONUNBUFFERED=1

cd /dcs04/hicks/data/Jan/sim_project

echo "Job started: $(date)"
echo "Host: $(hostname)"
echo "Python: ${PYTHON_BIN}"

DATASETS=(bin16um spot cell)
DATASET="${DATASETS[${SLURM_ARRAY_TASK_ID:-0}]}"

echo "Dataset: ${DATASET}"
echo "Command: ${PYTHON_BIN} sim_paper/code/applications_albis/svg/3D_nnsvg.py --dataset ${DATASET} --n-threads ${SLURM_CPUS_PER_TASK:-8} $*"
"${PYTHON_BIN}" sim_paper/code/applications_albis/svg/3D_nnsvg.py --dataset "${DATASET}" --n-threads "${SLURM_CPUS_PER_TASK:-8}" "$@"

echo "Job finished: $(date)"
