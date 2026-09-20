#!/bin/bash
# Dedicated Figure 4 cross-modality experiment; no Figure 2 outputs are modified.
set -euo pipefail
ROOT=/dcs04/hicks/data/Jan/sim_project
CODE="$ROOT/sim_paper/code/applications_albis/cross_modality_alignment"
OUT="$ROOT/sim_paper/data/figure_4/cross_modality_alignment/independent_offsets"
mkdir -p "$OUT/logs"
generation=$(sbatch --parsable --output="$OUT/logs/generate_%A_%a.out" "$CODE/generate_independent_offsets.sh")
printf 'generation\t%s\n' "$generation" >> "$OUT/submitted_jobs.tsv"
alignment=$(sbatch --parsable --dependency="afterok:${generation}" --output="$OUT/logs/stair_%j.out" "$CODE/run_independent_offsets_stair.sh")
printf 'alignment\t%s\n' "$alignment" >> "$OUT/submitted_jobs.tsv"
printf 'Generation array: %s\nDependent STAIR + plots: %s\n' "$generation" "$alignment"
