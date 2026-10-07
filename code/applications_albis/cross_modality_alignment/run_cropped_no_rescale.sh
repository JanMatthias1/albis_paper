#!/bin/bash
# Cross-modality STAIR with bin and spot cropped to a 2221 um square capture
# area (Visium's 6.5 mm scaled by the sphere shrink 2050/6000), cell uncropped,
# and no rescaling (the cross_tech_stair.py default), so cell keeps its true size
# (~4.1 mm disc) with the bin/spot squares inside it. Same tissue settings as
# Figure 4C (strong_domain_mix/run_strongmix_offsets.sh); only the crop differs.
# Usage: bash run_cropped_no_rescale.sh submit [OUTDIR]
#   submits generation, then STAIR + plots + metrics (afterok)
#   DRY_RUN=1 prints the sbatch commands without submitting.
set -euo pipefail
# SLURM runs a spool copy of this file, so BASH_SOURCE is not the project path.
ROOT=/dcs04/hicks/data/Jan/sim_project
CODE="$ROOT/albis_paper/code/applications_albis/cross_modality_alignment"
PYTHON="$ROOT/albis_paper/env/albis-tutorial/bin/python"
STAIR_ENV=/dcs04/hicks/data/multi-sample-alignment-benchmark/envs/STAIR
MODE=${1:-submit}
OUT=${2:-"$ROOT/albis_paper/data/figure_4/cross_modality_alignment/_internal_cropped_bin_spot_no_rescale"}
GEN_ARGS=(--crop-modalities bin,spot --crop-window-um 2221)
export PYTHONDONTWRITEBYTECODE=1
export MPLCONFIGDIR="/tmp/crossmod-norescale-mpl-${SLURM_JOB_ID:-local}"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"
export OPENBLAS_NUM_THREADS="$OMP_NUM_THREADS"
export MKL_NUM_THREADS="$OMP_NUM_THREADS"
case "$MODE" in
  submit)
    gen=(sbatch --parsable --job-name=crossmod_crop_gen --partition=shared --mem=150G --cpus-per-task=4
         --time=06:00:00 --output="$OUT/logs/generate_%j.out" "$CODE/run_cropped_no_rescale.sh" generate "$OUT")
    align=(sbatch --parsable --job-name=crossmod_norescale --partition=gpu --gres=gpu:l40s:1 --mem=100G
           --cpus-per-task=4 --time=03:00:00 --output="$OUT/logs/stair_%j.out"
           "$CODE/run_cropped_no_rescale.sh" align "$OUT")
    if [[ -n "${DRY_RUN:-}" ]]; then printf '%q ' "${gen[@]}"; echo; printf '%q ' "${align[@]}"; echo; exit 0; fi
    test ! -e "$OUT" || { echo "Output exists (archive it first): $OUT" >&2; exit 1; }
    mkdir -p "$OUT/logs"
    tar --exclude=__pycache__ --exclude=logs_cross_tech_stair -czf "$OUT/source_snapshot.tar.gz" -C "$ROOT" \
      albis_paper/code/applications_albis/cross_modality_alignment albis_paper/code/manuscript_style.py
    "$PYTHON" "$CODE/strong_domain_mix/generate_strongmix_offsets.py" --outdir "$OUT" "${GEN_ARGS[@]}" \
      --print-config > "$OUT/config.json"
    generation=$("${gen[@]}")
    alignment=$("${align[@]:0:2}" --dependency="afterok:$generation" "${align[@]:2}")
    printf 'generation\t%s\nalignment\t%s\n' "$generation" "$alignment" > "$OUT/submitted_jobs.tsv"
    printf 'Generation: %s\nSTAIR + plots + metrics: %s\n' "$generation" "$alignment"
    ;;
  generate)
    "$PYTHON" -u "$CODE/strong_domain_mix/generate_strongmix_offsets.py" --outdir "$OUT" "${GEN_ARGS[@]}"
    ;;
  align)
    source /jhpce/shared/jhpce/core/anaconda3/2023.03/etc/profile.d/conda.sh
    # An active venv from the submitting shell would shadow the STAIR python.
    if [[ -n "${VIRTUAL_ENV:-}" ]]; then PATH="${PATH//$VIRTUAL_ENV\/bin:/}"; unset VIRTUAL_ENV; fi
    conda activate "$STAIR_ENV"
    "$STAIR_ENV/bin/python" "$CODE/cross_tech_stair.py" --slice 5 --input-root "$OUT/data" \
      --output-base "$OUT/STAIR/cross_tech"
    "$PYTHON" "$CODE/plot/plot_independent_offsets.py" --base "$OUT"
    "$PYTHON" "$CODE/evaluation_gene_metrics.py" --root "$OUT"
    ;;
  *) echo 'Mode must be submit, generate or align' >&2; exit 2 ;;
esac
