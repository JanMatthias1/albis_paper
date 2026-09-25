#!/bin/bash
#SBATCH --job-name=crossmod_shift3x
#SBATCH --partition=gpu
#SBATCH --gres=gpu:l40s:1
#SBATCH --mem=100G
#SBATCH --cpus-per-task=4
#SBATCH --time=02:00:00
set -euo pipefail
ROOT=/dcs04/hicks/data/Jan/sim_project
CODE="$ROOT/sim_paper/code/applications_albis/cross_modality_alignment"
OUT="$ROOT/sim_paper/data/figure_4/cross_modality_alignment/independent_offsets_shift3x"
cd "$ROOT"
test -f "$OUT/perturbation_provenance.json"
source /jhpce/shared/jhpce/core/anaconda3/2023.03/etc/profile.d/conda.sh
conda activate /dcs04/hicks/data/multi-sample-alignment-benchmark/envs/STAIR
python "$CODE/cross_tech_stair.py" --slice 5 --input-root "$OUT/data" --output-base "$OUT/STAIR/cross_tech"
export MPLCONFIGDIR="/tmp/crossmod-mpl-${SLURM_JOB_ID}"
export XDG_CACHE_HOME="/tmp/crossmod-cache-${SLURM_JOB_ID}"
"$ROOT/sim_paper/env/albis-tutorial/bin/python" "$CODE/plot_independent_offsets.py" --base "$OUT"
