#!/usr/bin/env bash
set -euo pipefail
project=/dcs04/hicks/data/Jan/sim_project
exec "$project/comparison_methods/env/analysis/bin/python" "$project/sim_paper/code/comparison_methods/compute/summarize_native.py" "$@"
