#!/usr/bin/env bash
#SBATCH --partition=shared
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
# Tables (measurements.csv, summary.csv, RESULTS.md) and figure for a multi-seed run.
# Submitted by submit_seeds.sh after every benchmark job; can also be rerun by hand:
#   sbatch --export=ALL,RUN_ROOT=<run> report_seeds.sbatch
# Env: RUN_ROOT. Runs the live report/plot code (they only read the measurements).
set -euo pipefail
project=/dcs04/hicks/data/Jan/sim_project
code="$project/albis_paper/code/comparison_methods/compute"
export MPLCONFIGDIR="$RUN_ROOT/cache/report/matplotlib" XDG_CACHE_HOME="$RUN_ROOT/cache/report"
exec "$project/comparison_methods/env/analysis/bin/python" "$code/report_seeds.py" --root "$RUN_ROOT"
