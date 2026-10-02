#!/usr/bin/env bash
# Multi-seed native compute benchmark (seeds 2025, 101, 202; 10k-5M cells; 1 thread).
#   bash submit_seeds.sh [RUN_ROOT]   (default data/figure_5/compute/native_seeds_<date>; must not exist)
# Per seed: Splatter reference -> SPIDER per size and 3 scCube batches; ALBIS per size
# (no reference needed). Jobs run the scripts in this folder. Last, report_seeds.sbatch
# (after all jobs) writes measurements.csv, summary.csv, RESULTS.md and figures/compute/.
# Job IDs: RUN_ROOT/submissions.tsv.
set -euo pipefail
project=/dcs04/hicks/data/Jan/sim_project
code="$project/sim_paper/code/comparison_methods/compute"
root="${1:-$project/sim_paper/data/figure_5/compute/native_seeds_$(date +%Y%m%d)}"
"$project/comparison_methods/env/analysis/bin/python" "$code/prepare_seeds.py" --root "$root" >/dev/null
printf 'kind\tmethod\tn_cells\tseed\tjob_id\n' > "$root/submissions.tsv"
log() { printf '%s\t%s\t%s\t%s\t%s\n' "$@" >> "$root/submissions.tsv"; }
sb() { sbatch --parsable --export=ALL,RUN_ROOT="$root" "$@"; }

for seed in 2025 101 202; do
    ref_dir="$root/references/seed$seed"
    ref=$(sb --export=ALL,RUN_ROOT="$root",SETTINGS="$root/settings/n600000_seed$seed.json",REFERENCE="$ref_dir",OUTPUT="$ref_dir" \
        --mem=16G --time=00:30:00 --job-name="ref_seed$seed" --output="$root/logs/reference_seed${seed}_%j.out" \
        "$code/reference_native.sbatch")
    log reference splatter 10000 "$seed" "$ref"
    for n in 10000 100000 200000 400000 600000 1000000 2000000 5000000; do
        mem=32G; time=02:00:00
        [[ $n -ge 5000000 ]] && mem=96G && time=04:00:00
        for method in albis spider; do
            dep=(); [[ $method == spider ]] && dep=(--dependency="afterok:$ref")
            job=$(sb "${dep[@]}" --export=ALL,RUN_ROOT="$root",METHOD="$method",SETTINGS="$root/settings/n${n}_seed$seed.json",OUTPUT="$root/raw/${method}_n${n}_seed$seed",REFERENCE="$ref_dir" \
                --mem="$mem" --time="$time" --job-name="${method}_n${n}_s$seed" \
                --output="$root/logs/${method}_n${n}_seed${seed}_%j.out" "$code/native_extension.sbatch")
            log task "$method" "$n" "$seed" "$job"
        done
    done
    for batch in "b10k_100k:10000 100000" "b200k_600k:200000 400000 600000" "b1m_5m:1000000 2000000 5000000"; do
        name=${batch%%:*}; sizes=${batch#*:}
        job=$(sb --dependency="afterok:$ref" --export=ALL,RUN_ROOT="$root",SEED="$seed",BATCH="$name",SIZES="$sizes" \
            --job-name="sccube_${name}_s$seed" --output="$root/logs/sccube_${name}_seed${seed}_%j.out" \
            "$code/sccube_seed_batch.sbatch")
        log task sccube "${sizes// /,}" "$seed" "$job"
    done
done
all_jobs=$(awk 'NR > 1 {print $5}' "$root/submissions.tsv" | paste -sd:)
report=$(sb --dependency="afterany:$all_jobs" --job-name=compute_report \
    --output="$root/logs/report_%j.out" "$code/report_seeds.sbatch")
log report report all all "$report"
echo "Run: $root"; column -t "$root/submissions.tsv" | head -70
