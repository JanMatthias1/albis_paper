#!/usr/bin/env bash
#SBATCH --job-name=splatter_native556
#SBATCH --partition=shared
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=02:00:00
set -euo pipefail
project=/dcs04/hicks/data/Jan/sim_project
out="${1:?Supply new overview directory}"
seed="${2:-20260921}"
"$project/comparison_methods/env/splatter/bin/Rscript" "$project/comparison_methods/code/workflows/splatter_expression.R" \
 --contract "$out/splatter_contract.json" --seed "$seed" --out-dir "$out/splatter"
"$project/comparison_methods/env/analysis/bin/python" - "$out" <<'PY'
import json,sys
from pathlib import Path
p=Path(sys.argv[1])/'splatter'
m=json.loads((p/'manifest.json').read_text())
assert m['status']=='ok',m
names=(p/'genes.tsv').read_text().splitlines()
assert len(names)==len(set(names))==556
print('Verified fresh Splatter pool: 556 unique genes')
PY
