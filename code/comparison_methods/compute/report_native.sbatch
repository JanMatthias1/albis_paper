#!/usr/bin/env bash
#SBATCH --job-name=native_compute_report
#SBATCH --partition=shared
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:10:00
set -euo pipefail
project=/dcs04/hicks/data/Jan/sim_project
jobs=$(awk -F '\t' '$1=="task" {printf "%s%s", sep, $3; sep=","}' "$RUN_ROOT"/submissions_*.tsv)
sacct -j "$jobs" --parsable2 --noheader --format=JobID,State,ExitCode,ElapsedRaw,MaxRSS,NodeList > "$RUN_ROOT/slurm_accounting.psv"
exec "$project/comparison_methods/env/analysis/bin/python" "$project/sim_paper/code/comparison_methods/compute/report_native.py" --root "$RUN_ROOT"
