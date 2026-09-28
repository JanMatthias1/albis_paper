#!/usr/bin/env bash
set -euo pipefail
project=/dcs04/hicks/data/Jan/sim_project
root="${1:?benchmark root}"
phase="${2:-final}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export MPLCONFIGDIR="$root/cache/plots" XDG_CACHE_HOME="$root/cache"
mkdir -p "$MPLCONFIGDIR"
"$project/comparison_methods/env/analysis/bin/python" "$project/sim_paper/code/comparison_methods/compute/summarize_native.py" --root "$root"
if [[ "$phase" == pilot ]]; then
    "$project/comparison_methods/env/analysis/bin/python" "$project/sim_paper/code/comparison_methods/compute/advance_native.py" --root "$root"
fi
