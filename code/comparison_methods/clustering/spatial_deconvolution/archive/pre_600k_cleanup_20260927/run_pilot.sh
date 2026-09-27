#!/usr/bin/env bash
#SBATCH --job-name=rctd_comparator
#SBATCH --time=04:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --partition=shared
set -euo pipefail
METHOD="${1:?Specify spider or sccube}"
case "$METHOD" in spider|sccube) ;; *) exit 2 ;; esac
ROOT=/dcs04/hicks/data/Jan/sim_project
PAPER="$ROOT/sim_paper"
CODE="$PAPER/code/comparison_methods/clustering/spatial_deconvolution"
SOURCE="${RCTD_SOURCE:-$PAPER/data/comparison_methods/overview_full_35831924}"
GENE_ARGS=()
RUN_PREFIX="${RCTD_RUN_PREFIX:-pilot}"
if [[ -n "${RCTD_N_GENES:-}" ]]; then
    GENE_ARGS=(--n-genes "$RCTD_N_GENES")
    RUN_PREFIX="genes${RCTD_N_GENES}"
fi
RUN="$PAPER/data/comparison_methods/spatial_deconvolution/${RUN_PREFIX}_${METHOD}_${SLURM_JOB_ID:-local_$(date +%Y%m%d_%H%M%S)}"
mkdir -p "$RUN"
"$ROOT/comparison_methods/env/analysis/bin/python" "$CODE/prepare_pilot.py" \
  --source "$SOURCE" --out "$RUN/inputs" --method "$METHOD" "${GENE_ARGS[@]}"
RCTD_ENV="$PAPER/env/rctd"
export LD_LIBRARY_PATH="$RCTD_ENV/lib:${LD_LIBRARY_PATH:-}"
export R_HOME="$RCTD_ENV/lib/R"
export RETICULATE_PYTHON="$RCTD_ENV/bin/python"
export RCTD_CELL_H5AD="$RUN/inputs/cell.h5ad"
export RCTD_SPOT_H5AD="$RUN/inputs/spot.h5ad"
export RCTD_OUT_DIR="$RUN/results"
SCRIPT="$PAPER/code/clustering/spatial_deconvolution/weak_domain/run_rctd_spot.R"
sha256sum "$SCRIPT" > "$RUN/rctd_script.sha256"
"$RCTD_ENV/bin/Rscript" -e 'set.seed(20260926); source(commandArgs(trailingOnly=TRUE)[1])' "$SCRIPT"
echo "Completed $METHOD: $RUN"

"$ROOT/comparison_methods/env/analysis/bin/python" "$CODE/plot_pilot.py" "$RUN"
