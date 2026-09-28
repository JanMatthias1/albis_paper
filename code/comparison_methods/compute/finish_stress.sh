#!/bin/bash
#SBATCH --job-name=compute_stress_plots
#SBATCH --time=00:30:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --partition=shared
set -euo pipefail
project=/dcs04/hicks/data/Jan/sim_project
root="$1"
job="$2"
sacct -j "$job" --parsable2 --noheader --format=JobIDRaw,State,ExitCode,ElapsedRaw,MaxRSS > "$root/slurm_accounting.psv"
export MPLCONFIGDIR="$root/cache/plots" XDG_CACHE_HOME="$root/cache"
"$project/comparison_methods/env/analysis/bin/python" "$project/sim_paper/code/comparison_methods/compute/finalize_stress.py" --root "$root"
"$project/comparison_methods/env/analysis/bin/python" "$project/sim_paper/code/comparison_methods/compute/plot_compute.py" --input-dir "$root/raw" --output-dir "$root/figures"
