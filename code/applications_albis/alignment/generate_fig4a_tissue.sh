#!/bin/bash
# Figure 4A input: generate the strong-domain-mix shared tissue (r = 4600 um,
# 1M cells; bin16um, spot and cell from one ALBIS call) that 3D_stair.py reads.
# Same generator as Figure 4C (cross_modality_alignment/strong_domain_mix/
# generate_strongmix_offsets.py), with a larger sphere and more cells.
#
# Usage: bash generate_fig4a_tissue.sh [OUTDIR]
#   prints the generation job id; then run STAIR gated on it:
#   bash run_stair_3D.sh <job id>
#   DRY_RUN=1 prints the sbatch command without submitting.
set -euo pipefail
ROOT=/dcs04/hicks/data/Jan/sim_project
CM="$ROOT/sim_paper/code/applications_albis/cross_modality_alignment"
PYTHON="$ROOT/sim_paper/env/albis-tutorial/bin/python"
OUT=${1:-"$ROOT/sim_paper/data/figure_4/alignment/strong_domain_mix_shift3x_r4600"}
GEN_ARGS=(--sphere-r-um 4600 --n-cells 1000000)

cmd=(sbatch --parsable --job-name=fig4a_r4600_gen --partition=shared --mem=150G --cpus-per-task=4
     --time=06:00:00 --output="$OUT/logs/generate_%j.out"
     "$CM/strong_domain_mix/run_strongmix_offsets.sh" generate "$OUT" "${GEN_ARGS[@]}")
if [[ -n "${DRY_RUN:-}" ]]; then printf '%q ' "${cmd[@]}"; echo; exit 0; fi

test ! -e "$OUT" || { echo "Output exists (archive it first): $OUT" >&2; exit 1; }
mkdir -p "$OUT/logs"
# record the exact code and settings used
tar --exclude=__pycache__ --exclude=logs_cross_tech_stair -czf "$OUT/source_snapshot.tar.gz" -C "$ROOT" \
    sim_paper/code/applications_albis/cross_modality_alignment sim_paper/code/applications_albis/alignment
"$PYTHON" "$CM/strong_domain_mix/generate_strongmix_offsets.py" --outdir "$OUT" "${GEN_ARGS[@]}" \
    --print-config > "$OUT/config.json"
job=$("${cmd[@]}")
printf 'generation\t%s\n' "$job" > "$OUT/submitted_jobs.tsv"
echo "$job"
