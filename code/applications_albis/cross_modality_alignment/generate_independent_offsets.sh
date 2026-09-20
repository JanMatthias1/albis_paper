#!/bin/bash
#SBATCH --job-name=crossmod_generate
#SBATCH --array=0-2
#SBATCH --partition=shared
#SBATCH --mem=250G
#SBATCH --cpus-per-task=4
#SBATCH --time=12:00:00
set -euo pipefail
ROOT=/dcs04/hicks/data/Jan/sim_project
cd "$ROOT"
source "$ROOT/sim_paper/code/count_distribution/_env.sh"
export MPLCONFIGDIR="/tmp/crossmod-mpl-${SLURM_JOB_ID:-local}"
export XDG_CACHE_HOME="/tmp/crossmod-cache-${SLURM_JOB_ID:-local}"
TECHS=(bin16um spot cell)
MODS=(bin spot cell)
THETAS=(2.0 0.25 0.40)
JITTERS=(1.0 0.10 0.15)
LOGMUS=(-2.5 -2.0 -2.5)
BATCH=(0.7 0.3 1.5)
NCELLS=(600000 600000 24207)
BINS=(16 8 8)
i="${SLURM_ARRAY_TASK_ID:?Array task required}"
OUT="$ROOT/sim_paper/data/figure_4/cross_modality_alignment/independent_offsets/data/${TECHS[$i]}"
mkdir -p "$OUT"
# Intentionally OMIT --sync-unaligned-seed: native modality-specific RNG offsets.
# Keep seed 2025 for tissue/expression, changing only perturbation synchronization.
if [[ ! -f "$OUT/simulation_${MODS[$i]}_z.h5ad" ]]; then
 "$PYTHON_BIN" sim_paper/code/data/generate_simulation_noisy.py \
  --modality "${MODS[$i]}" --slice-axis Z --seed 2025 \
  --sphere-r-um 2050 --n-cells "${NCELLS[$i]}" --cell-r-mean 7.5 \
  --theta "${THETAS[$i]}" --theta-jitter "${JITTERS[$i]}" \
  --noise-scale 1.3 --marker-foldchange 3.5 --shared-marker-foldchange 2.5 \
  --base-gene-lognormal "${LOGMUS[$i]}" 0.7 --batch-sigma "${BATCH[$i]}" \
  --bin-size-um "${BINS[$i]}" --core-fuzz-width-um 102.5 --max-shift 1025 \
  --output-dir "$OUT"
fi
"$PYTHON_BIN" sim_paper/code/clustering/00_qc_filter.py --modality "${MODS[$i]}" \
 --input "$OUT/simulation_${MODS[$i]}_z.h5ad" --output "$OUT/simulation_${MODS[$i]}_z_qc.h5ad"
