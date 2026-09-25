#!/bin/bash
# Strong-domain-mix shift3x: shared-tissue generation -> section-5 STAIR ->
# plots + reference RMSE + gene metrics (same downstream as ../run_native_offsets.sh).
# Usage: bash run_strongmix_offsets.sh submit [OUTDIR]
set -euo pipefail
# SLURM runs a spool copy of this file, so BASH_SOURCE is not the project path.
ROOT=/dcs04/hicks/data/Jan/sim_project
CODE="$ROOT/sim_paper/code/applications_albis/cross_modality_alignment"
HERE="$CODE/strong_domain_mix"
PYTHON="$ROOT/sim_paper/env/albis-tutorial/bin/python"
MODE=${1:-submit}
OUT=${2:-"$ROOT/sim_paper/data/figure_4/cross_modality_alignment/strong_domain_mix_shift3x"}
export PYTHONDONTWRITEBYTECODE=1
export MPLCONFIGDIR="/tmp/strongmix-crossmod-mpl-${SLURM_JOB_ID:-local}"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"
export OPENBLAS_NUM_THREADS="$OMP_NUM_THREADS"
export MKL_NUM_THREADS="$OMP_NUM_THREADS"
case "$MODE" in
  submit)
    test ! -e "$OUT/submitted_jobs.tsv"
    mkdir -p "$OUT/logs"
    tar --exclude=__pycache__ --exclude=logs_cross_tech_stair -czf "$OUT/source_snapshot.tar.gz" -C "$ROOT" \
      albis/albis sim_paper/code/applications_albis/cross_modality_alignment sim_paper/code/manuscript_style.py
    "$PYTHON" "$HERE/generate_strongmix_offsets.py" --outdir "$OUT" --print-config > "$OUT/config.json"
    generation=$(sbatch --parsable --job-name=crossmod_strongmix --partition=shared \
      --mem=150G --cpus-per-task=4 --time=06:00:00 --output="$OUT/logs/generate_%j.out" \
      "$HERE/run_strongmix_offsets.sh" generate "$OUT")
    printf 'generation\t%s\n' "$generation" > "$OUT/submitted_jobs.tsv"
    alignment=$(sbatch --parsable --job-name=crossmod_strongmix_stair --partition=gpu --gres=gpu:l40s:1 \
      --mem=100G --cpus-per-task=4 --time=03:00:00 --dependency="afterok:$generation" \
      --output="$OUT/logs/stair_%j.out" "$HERE/run_strongmix_offsets.sh" align "$OUT")
    printf 'alignment\t%s\n' "$alignment" >> "$OUT/submitted_jobs.tsv"
    printf 'Generation: %s\nSTAIR + plots + metrics: %s\n' "$generation" "$alignment"
    ;;
  generate)
    "$PYTHON" -u "$HERE/generate_strongmix_offsets.py" --outdir "$OUT"
    ;;
  align)
    source /jhpce/shared/jhpce/core/anaconda3/2023.03/etc/profile.d/conda.sh
    # sbatch exports the submitting shell's env; an active venv (e.g. splatter)
    # would stay first on PATH and shadow the STAIR python (job 35909602 failed so).
    if [[ -n "${VIRTUAL_ENV:-}" ]]; then PATH="${PATH//$VIRTUAL_ENV\/bin:/}"; unset VIRTUAL_ENV; fi
    STAIR_ENV=/dcs04/hicks/data/multi-sample-alignment-benchmark/envs/STAIR
    conda activate "$STAIR_ENV"
    "$STAIR_ENV/bin/python" "$CODE/cross_tech_stair.py" --slice 5 --input-root "$OUT/data" --output-base "$OUT/STAIR/cross_tech"
    "$PYTHON" "$CODE/plot_independent_offsets.py" --base "$OUT"
    "$PYTHON" "$CODE/evaluation_gene_metrics.py" --root "$OUT"
    ;;
  *) echo 'Mode must be submit, generate, or align' >&2; exit 2 ;;
esac
