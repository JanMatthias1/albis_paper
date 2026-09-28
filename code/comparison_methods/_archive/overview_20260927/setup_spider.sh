#!/usr/bin/env bash
# Fetch unmodified upstream Python sources into an isolated dependency folder.
# Reuses the existing Spider environment; does not install or modify packages.
set -euo pipefail
script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
source_dir="$script_dir/../../../env/spider-overview-src"
commit=6ccd4da77257f2807c430f8f42fbe2dc175991de
mkdir -p "$source_dir/spider"
for module in __init__ api core neighbors solver sim_expr simulate_10X Annealing enhance utils scsim sim_naive random_based_utils data_op joint_simulator opt random_generator run_scsim; do
    curl -L --fail --silent --show-error \
        "https://raw.githubusercontent.com/YANG-ERA/Spider/$commit/spider/$module.py" \
        -o "$source_dir/spider/$module.py"
done
printf '%s\n' "$commit" > "$source_dir/COMMIT"
echo "Pinned native Spider sources: $source_dir"
