#!/usr/bin/env bash
#SBATCH --job-name=celltype_comparator
#SBATCH --time=04:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --partition=shared
set -euo pipefail
METHOD="${1:?Specify spider or sccube}"
case "$METHOD" in spider|sccube) ;; *) exit 2 ;; esac
PAPER=/dcs04/hicks/data/Jan/sim_project/sim_paper
source "$PAPER/code/clustering/_env.sh"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"
export OPENBLAS_NUM_THREADS="$OMP_NUM_THREADS"
export MKL_NUM_THREADS="$OMP_NUM_THREADS"
export NUMBA_NUM_THREADS="$OMP_NUM_THREADS"
SOURCE="${CELL_CLUSTER_SOURCE:-$PAPER/data/comparison_methods/overview_full_35831924}"
RUN_PREFIX="${CELL_CLUSTER_PREFIX:-pilot}"
EXTRA_ARGS=()
if [[ -n "${CELL_CLUSTER_N_GENES:-}" ]]; then
    EXTRA_ARGS=(--expected-n-genes "$CELL_CLUSTER_N_GENES")
fi
"$PYTHON_BIN" "$PAPER/code/comparison_methods/clustering/cell_type/run_pilot.py" \
  --method "$METHOD" \
  --source "$SOURCE" \
  --out "$PAPER/data/comparison_methods/clustering/cell_type/${RUN_PREFIX}_${METHOD}_${SLURM_JOB_ID:-local_$(date +%Y%m%d_%H%M%S)}" "${EXTRA_ARGS[@]}"
