#!/bin/bash
# Tabulate sweep_spot_theta_rebracket.sh: theta in {0.15,0.20,0.25,0.30,0.35}
# at the CORRECTED theta_jitter=0.10, log_mu=-2.0 fixed. Same composite
# convention as summary_spot_logmu_theta_joint.sh (sum of |ln(sim/real
# ratio)| over theta_hat/total_counts_median/matrix_zero_frac, summed across
# both probe refs, lower = better). Read-only; no jobs.
set -euo pipefail
cd /dcs04/hicks/data/Jan/sim_project/sim_paper/data/count_distribution/sweeps

python3 - <<'PY'
import json, math, os

BASE = "spot_theta_rebracket_probe"
REALS = ["tonsil_visium", "lymph_node_visium"]
THETAS = ["0.15", "0.20", "0.25", "0.30", "0.35"]

def load(path):
    with open(path) as fh:
        return json.load(fh)

hdr = f"{'theta':>6} {'ref':>18} | {'thetaHat s/r (ratio)':>26} | {'tcMed s/r (ratio)':>24} | {'gpc s/r':>13} | {'zero s/r (ratio)':>22} | {'compΣ3':>7}"
print(hdr); print("-" * len(hdr))

grid_totals = {}
for th in THETAS:
    tag = f"spot_theta_rebracket_{th}"
    per_ref_comp = []
    for real in REALS:
        p = f"{BASE}/{tag}_hvg_matched_vs_{real}/comparison_summary.json"
        if not os.path.exists(p):
            print(f"{th:>6} {real:>18} | MISSING ({p})")
            per_ref_comp.append(None); continue
        j = load(p)
        def sr(key):
            d = j[key]; s = d["spot"]; r = [v for k, v in d.items() if k != "spot"][0]
            return s, r, (s / r if r else float("nan"))
        ths, thr, thratio = sr("theta_hat")
        tcs, tcr, tcratio = sr("total_counts_median")
        gps, gpr, _ = sr("genes_per_cell_median")
        zs, zr, zratio = sr("matrix_zero_frac")
        comp = sum(abs(math.log(x)) for x in (thratio, tcratio, zratio))
        per_ref_comp.append(comp)
        print(f"{th:>6} {real:>18} | "
              f"{ths:6.3f}/{thr:6.3f} ({thratio:4.2f}) | "
              f"{tcs:8.0f}/{tcr:8.0f} ({tcratio:4.2f}) | "
              f"{gps:5.0f}/{gpr:5.0f} | "
              f"{zs:5.3f}/{zr:5.3f} ({zratio:4.2f}) | "
              f"{comp:7.3f}")
    if all(c is not None for c in per_ref_comp):
        grid_totals[th] = sum(per_ref_comp)

print()
print("composite summed across BOTH refs (lower = better):")
for th, tot in sorted(grid_totals.items(), key=lambda kv: kv[1]):
    print(f"  theta={th:>5}  ->  {tot:.3f}")
PY
