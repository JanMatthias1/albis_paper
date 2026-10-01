#!/usr/bin/env bash
#SBATCH --job-name=fig3_noharmony
#SBATCH --partition=shared
#SBATCH --cpus-per-task=8
#SBATCH --mem=120G
#SBATCH --time=12:00:00
set -euo pipefail
project=/dcs04/hicks/data/Jan/sim_project
code="$project/sim_paper/code/clustering/no_harmony_domain"
phase="${1:?phase: sweep/final/select/summary}"
export OPENBLAS_NUM_THREADS=8 OMP_NUM_THREADS=8 NUMBA_NUM_THREADS=8 PYTHONDONTWRITEBYTECODE=1
export MPLCONFIGDIR="/tmp/fig3_noharmony_${SLURM_JOB_ID:-$$}" XDG_CACHE_HOME="/tmp/fig3_noharmony_cache_${SLURM_JOB_ID:-$$}"
mkdir -p "$MPLCONFIGDIR" "$XDG_CACHE_HOME"
tutorial="$project/sim_paper/env/albis-tutorial/bin/python"
if [[ "$phase" == select ]]; then exec "$tutorial" -u "$code/select_parameters.py"; fi
if [[ "$phase" == summary ]]; then
    "$tutorial" -u "$code/summarize.py"
    exec bash "$project/sim_paper/code/clustering/ari_recovery_summary/plot_recovery_no_harmony.sh"
fi
[[ "$phase" == sweep || "$phase" == final ]]
task="${SLURM_ARRAY_TASK_ID:?submit as an array}"
# Only expression controls use the tutorial environment for embedding; BANKSY
# feature generation needs its existing isolated package versions.
embed_python="$project/sim_paper/env/albis-banksy/bin/python"
if [[ "$phase" == final ]] && (( task % 2 == 1 )); then embed_python="$tutorial"; fi
"$embed_python" -u "$code/run.py" --phase "$phase" --task "$task" --stage embed
"$tutorial" -u "$code/run.py" --phase "$phase" --task "$task" --stage cluster
