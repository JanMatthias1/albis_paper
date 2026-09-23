#!/bin/bash
#SBATCH --job-name=albis_compute
#SBATCH --array=0-14
#SBATCH --time=2-00:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --partition=shared
#SBATCH --output=compute_%A_%a.out
set -euo pipefail
# SLURM copies this file; use the project path unless explicitly overridden.
PROJECT="${SIM_PROJECT_ROOT:-/dcs04/hicks/data/Jan/sim_project}"
SCRIPT_DIR="$PROJECT/sim_paper/code/comparison_methods/compute"
PYTHON="$PROJECT/comparison_methods/env/analysis/bin/python"
if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    COUNTS=(10000 20000 40000 50000 100000)
    TECHNOLOGIES=(cell bin spot)
    IDX="$SLURM_ARRAY_TASK_ID"
    if (( IDX < 0 || IDX > 14 )); then
        echo 'Array index must be 0–14' >&2
        exit 2
    fi
    N="${COUNTS[$((IDX / 3))]}"
    TECH="${TECHNOLOGIES[$((IDX % 3))]}"
    OUTPUT_ROOT="${COMPUTE_OUTPUT_ROOT:-$PROJECT/sim_paper/data/comparison_methods/compute/raw}"
    exec "$PYTHON" "$SCRIPT_DIR/run_compute.py" --n-cells "$N" --technology "$TECH" \
        --out-dir "$OUTPUT_ROOT/n${N}_${TECH}" "$@"
else
    exec "$PYTHON" "$SCRIPT_DIR/run_compute.py" "$@"
fi
