#!/usr/bin/env bash
set -euo pipefail
project=/dcs04/hicks/data/Jan/sim_project
root="${1:?benchmark root}"
kind="${2:?reference or task}"
key="${3:?seed or task key}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
export PYTHONDONTWRITEBYTECODE=1
if [[ "$kind" == reference ]]; then
    exec "$project/comparison_methods/env/analysis/bin/python" "$project/sim_paper/code/comparison_methods/compute/native_benchmark.py" reference --root "$root" --seed "$key"
elif [[ "$kind" == task ]]; then
    exec "$project/comparison_methods/env/analysis/bin/python" "$project/sim_paper/code/comparison_methods/compute/native_benchmark.py" task --root "$root" --key "$key"
else
    echo "Unknown task kind: $kind" >&2
    exit 2
fi
