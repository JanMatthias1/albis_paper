#!/usr/bin/env bash
# Scheduler submission only. Each job executes its own direct-native script.
set -euo pipefail
root="${1:?prepared run directory}"
phase="${2:-pilot}"
project=/dcs04/hicks/data/Jan/sim_project
code="$project/sim_paper/code/comparison_methods/compute"
if [[ -e "$root/submissions_${phase}.tsv" ]]; then
    echo "Submission receipt already exists; inspect it before resubmitting." >&2
    exit 2
fi
if [[ "$phase" != pilot && "$phase" != smoke ]]; then
    echo "Only validated pilot/smoke submissions supported at this stage." >&2
    exit 2
fi
# Check frozen source copies; this performs no method execution or adaptation.
"$project/comparison_methods/env/analysis/bin/python" - "$root" "$code" <<'PY'
import hashlib,json,sys
from pathlib import Path
root,code=map(Path,sys.argv[1:]);p=json.loads((root/'protocol.json').read_text())
for name,expected in p['source_hashes'].items():
    assert hashlib.sha256((code/name).read_bytes()).hexdigest()==expected,name
PY
printf 'kind\tkey\tjob_id\n' > "$root/submissions_${phase}.tsv"
count=10000
[[ "$phase" == smoke ]] && count=32
export SETTINGS="$root/settings/n${count}_seed2025.json" REFERENCE="$root/references/seed2025" OUTPUT="$root/reference_job_seed2025" PYTHONHASHSEED=2025
reference_dependencies=()
[[ -n "${VALIDATION_JOBS:-}" ]] && reference_dependencies+=(--dependency="afterok:$VALIDATION_JOBS" --kill-on-invalid-dep=yes)
reference_job=$(sbatch --parsable "${reference_dependencies[@]}" --mem=16G --time=00:30:00 --job-name=native_reference --output="$root/logs/reference_%j.out" "$code/reference_native.sbatch")
reference_job="${reference_job%%;*}"
printf 'reference\tseed2025\t%s\n' "$reference_job" >> "$root/submissions_${phase}.tsv"
declare -A previous
while IFS=$'\t' read -r method n seed settings output reference; do
    [[ "$method" == method ]] && continue
    [[ "$seed" != 2025 ]] && continue
    if [[ "$phase" == pilot && "$n" != 10000 && "$n" != 100000 ]]; then continue; fi
    if [[ "$phase" == smoke && "$n" != 32 ]]; then continue; fi
    export SETTINGS="$settings" OUTPUT="$output" REFERENCE="$reference" PYTHONHASHSEED="$seed"
    dependency="${VALIDATION_JOBS:-}"
    if [[ "$method" != albis ]]; then dependency="${dependency:+$dependency:}$reference_job"; fi
    if [[ -n "${previous[$method]:-}" ]]; then dependency="${dependency:+$dependency:}${previous[$method]}"; fi
    args=()
    [[ "$phase" == smoke ]] && args+=(--mem=16G --time=00:30:00)
    [[ -n "$dependency" ]] && args+=(--dependency="afterok:$dependency" --kill-on-invalid-dep=yes)
    job=$(sbatch --parsable "${args[@]}" --job-name="native_${method}_${n}" --output="$root/logs/${method}_n${n}_%j.out" "$code/${method}_native.sbatch")
    job="${job%%;*}";previous[$method]="$job"
    printf 'task\t%s_n%s_seed%s\t%s\n' "$method" "$n" "$seed" "$job" >> "$root/submissions_${phase}.tsv"
    echo "$method $n seed$seed: $job"
done < "$root/tasks.tsv"
echo "Receipt: $root/submissions_${phase}.tsv"
export RUN_ROOT="$root"
report_dependencies=$(awk -F '\t' '$1=="task" {printf "%s%s", sep, $3; sep=":"}' "$root/submissions_${phase}.tsv")
report_job=$(sbatch --parsable --dependency="afterany:$report_dependencies" --output="$root/logs/report_%j.out" "$code/report_native.sbatch")
printf 'report\t%s\t%s\n' "$phase" "${report_job%%;*}" >> "$root/submissions_${phase}.tsv"
echo "Report job: ${report_job%%;*}"
