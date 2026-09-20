#!/bin/bash
# Native ALBIS generation -> STAIR -> approved plots + both metric workflows.
# Usage: bash run_native_offsets.sh submit [OUTDIR] [MAX_SHIFT] [OFFSET_SEED] [TISSUE_SEED] [REPLACE_TARGET]
# Change OFFSET_SEED for new rotations/directions; changing MAX_SHIFT alone scales
# the same native draws. Defaults reproduce the historical shift3x condition.
set -euo pipefail
# SLURM runs a spool copy of this file, so BASH_SOURCE is not the project path.
ROOT=/dcs04/hicks/data/Jan/sim_project
CODE="$ROOT/sim_paper/code/applications_albis/cross_modality_alignment"
PYTHON="$ROOT/sim_paper/env/albis-tutorial/bin/python"
MODE=${1:-submit}
OUT=${2:-"$ROOT/sim_paper/data/figure_4/cross_modality_alignment/native_offsets_3075_seed12345"}
SHIFT=${3:-3075}
OFFSET_SEED=${4:-12345}
TISSUE_SEED=${5:-2025}
REPLACE_TARGET=${6:-}
export PYTHONDONTWRITEBYTECODE=1
export MPLCONFIGDIR="/tmp/native-crossmod-mpl-${SLURM_JOB_ID:-local}"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"
export OPENBLAS_NUM_THREADS="$OMP_NUM_THREADS"
export MKL_NUM_THREADS="$OMP_NUM_THREADS"
case "$MODE" in
  submit)
    mkdir -p "$OUT/logs"
    # Snapshot source because a working-tree file can change while jobs queue.
    test ! -e "$OUT/submitted_jobs.tsv"
    tar --exclude=__pycache__ --exclude=logs_cross_tech_stair -czf "$OUT/source_snapshot.tar.gz" -C "$ROOT" \
      albis/albis \
      sim_paper/code/applications_albis/cross_modality_alignment sim_paper/code/manuscript_style.py
    "$PYTHON" "$CODE/generate_native_offsets.py" --modality bin16um \
      --outdir "$OUT" --max-shift "$SHIFT" --base-seed-unaligned "$OFFSET_SEED" --seed "$TISSUE_SEED" \
      --print-config > "$OUT/native_config_bin16um.json"
    generation=$(sbatch --parsable --job-name=crossmod_native --partition=shared --array=0-2 \
      --mem=250G --cpus-per-task=4 --time=12:00:00 --output="$OUT/logs/generate_%A_%a.out" \
      "$CODE/run_native_offsets.sh" generate "$OUT" "$SHIFT" "$OFFSET_SEED" "$TISSUE_SEED")
    printf 'generation\t%s\n' "$generation" > "$OUT/submitted_jobs.tsv"
    alignment=$(sbatch --parsable --job-name=crossmod_native_stair --partition=gpu --gres=gpu:l40s:1 \
      --mem=100G --cpus-per-task=4 --time=02:00:00 --dependency="afterok:$generation" \
      --output="$OUT/logs/stair_%j.out" "$CODE/run_native_offsets.sh" align "$OUT" "$SHIFT" "$OFFSET_SEED" "$TISSUE_SEED")
    printf 'alignment\t%s\n' "$alignment" >> "$OUT/submitted_jobs.tsv"
    printf 'Native generation: %s\nSTAIR + plots + metrics: %s\n' "$generation" "$alignment"
    if [[ -n "$REPLACE_TARGET" ]]; then
      publish=$(sbatch --parsable --job-name=crossmod_replace --partition=shared --mem=2G \
        --cpus-per-task=1 --time=00:10:00 --dependency="afterok:$alignment" \
        --output="$OUT/logs/replace_%j.out" "$CODE/run_native_offsets.sh" publish \
        "$OUT" "$SHIFT" "$OFFSET_SEED" "$TISSUE_SEED" "$REPLACE_TARGET")
      printf 'replacement\t%s\n' "$publish" >> "$OUT/submitted_jobs.tsv"
      printf 'Replace previous experiment after success: %s\n' "$publish"
    fi
    ;;
  generate)
    techs=(bin16um spot cell)
    "$PYTHON" -u "$CODE/generate_native_offsets.py" \
      --modality "${techs[${SLURM_ARRAY_TASK_ID:?Array task required}]}" --outdir "$OUT" \
      --max-shift "$SHIFT" --max-deg 270 --base-seed-unaligned "$OFFSET_SEED" --seed "$TISSUE_SEED"
    ;;
  align)
    source /jhpce/shared/jhpce/core/anaconda3/2023.03/etc/profile.d/conda.sh
    conda activate /dcs04/hicks/data/multi-sample-alignment-benchmark/envs/STAIR
    python "$CODE/cross_tech_stair.py" --slice 5 --input-root "$OUT/data" --output-base "$OUT/STAIR/cross_tech"
    "$PYTHON" "$CODE/plot_independent_offsets.py" --base "$OUT"
    "$PYTHON" "$CODE/evaluation_gene_metrics.py" --root "$OUT"
    ;;
  publish)
    test -n "$REPLACE_TARGET"
    "$PYTHON" "$CODE/promote_native_offsets.py" --staging "$OUT" --target "$REPLACE_TARGET"
    ;;
  *) echo 'Mode must be submit, generate, align, or publish' >&2; exit 2 ;;
esac
