#!/bin/bash
# Tabulate sweep_cell_kgeom_r2050.sh: domain_true ARI + slice_id leakage per
# k_geom, at cell's new r=2050/n=24207 scale (bs=0, no batch effect). Lower
# k_geom = smaller spatial neighborhood; pick the k_geom with the best
# domain ARI AND low slice_id leakage (a high-ARI point with high leakage
# just means clustering == slice_id, not real domain recovery -- same
# validity check composition_recovery.py's slice_id_leakage_ari exists for
# everywhere else in this project). Read-only; no jobs.
set -euo pipefail
cd /dcs04/hicks/data/Jan/sim_project/sim_paper/data/figure_3/cellbin_batch_sigma_slide/misc/cell_kgeom_sweep_r2050/scores

python3 - <<'PY'
import json, os

KGS = [15, 30, 60, 100, 150, 200]

hdr = f"{'k_geom':>7} | {'domain ARI':>10} | {'celltype ARI':>12} | {'slice leakage':>13}"
print(hdr); print("-" * len(hdr))

rows = []
for kg in KGS:
    p = f"composition_recovery_cell_kg{kg}.json"
    if not os.path.exists(p):
        print(f"{kg:>7} | MISSING ({p})")
        continue
    with open(p) as fh:
        j = json.load(fh)
    dom = j["levels"]["domain"]
    ct = j["levels"]["cell_type"]
    print(f"{kg:>7} | {dom['hard_label']['ari']:>10.4f} | {ct['hard_label']['ari']:>12.4f} | {dom['slice_id_leakage_ari']:>13.4f}")
    rows.append((kg, dom['hard_label']['ari'], dom['slice_id_leakage_ari']))

print()
print("ranked by domain ARI (lower slice_id leakage = more trustworthy):")
for kg, ari, leak in sorted(rows, key=lambda r: -r[1]):
    flag = "  <- high leakage, may be a slice_id artifact" if leak > 0.3 else ""
    print(f"  k_geom={kg:>4}  domain_ARI={ari:.4f}  leakage={leak:.4f}{flag}")
PY
