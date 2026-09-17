#!/bin/bash
# Tabulate the spot log_mu x theta joint mini-grid: the 3 new log_mu=-2.0
# points from sweep_spot_logmu_theta_joint_probe.sh + the 3 already-computed
# log_mu=-2.5 points from sweep_spot_theta_probe.sh (byte-identical config).
# Prints one row per (log_mu, theta, real ref) with sim/real values and ratios
# for theta_hat, total_counts_median, genes_per_cell_median, matrix_zero_frac,
# plus a composite = sum of |ln(ratio)| over {theta_hat, total_counts_median,
# matrix_zero_frac}, and that composite summed across BOTH refs (the score the
# probe sweeps optimize). Read-only; no jobs.
set -euo pipefail
cd /dcs04/hicks/data/Jan/sim_project/sim_paper/data/count_distribution/sweeps

python3 - <<'PY'
import json, math, os

NEW = "spot_joint_logmu_theta_probe"          # log_mu=-2.0 row (this sweep)
OLD = "spot_theta_sweep_probe"                # log_mu=-2.5 row (theta sweep)
REALS = ["tonsil_visium", "lymph_node_visium"]

def load(path):
    with open(path) as fh:
        return json.load(fh)

def rowdir(base, tag, real):
    return f"{base}/{tag}_hvg_matched_vs_{real}/comparison_summary.json"

points = []  # (log_mu, theta, base, tagfmt)
for th in ("0.25", "0.35", "0.5"):
    points.append((-2.5, th, OLD, f"spot_theta_sweep_{th}"))
    points.append((-2.0, th, NEW, f"spot_joint_logmu_-2.0_theta_{th}"))

hdr = f"{'log_mu':>7} {'theta':>6} {'ref':>18} | {'thetaHat s/r (ratio)':>26} | {'tcMed s/r (ratio)':>24} | {'gpc s/r':>13} | {'zero s/r (ratio)':>22} | {'compΣ3':>7}"
print(hdr); print("-" * len(hdr))

grid_totals = {}
for log_mu, th, base, tag in sorted(points):
    per_ref_comp = []
    for real in REALS:
        p = rowdir(base, tag, real)
        if not os.path.exists(p):
            print(f"{log_mu:>7} {th:>6} {real:>18} | MISSING ({p})")
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
        print(f"{log_mu:>7} {th:>6} {real:>18} | "
              f"{ths:6.3f}/{thr:6.3f} ({thratio:4.2f}) | "
              f"{tcs:8.0f}/{tcr:8.0f} ({tcratio:4.2f}) | "
              f"{gps:5.0f}/{gpr:5.0f} | "
              f"{zs:5.3f}/{zr:5.3f} ({zratio:4.2f}) | "
              f"{comp:7.3f}")
    if all(c is not None for c in per_ref_comp):
        grid_totals[(log_mu, th)] = sum(per_ref_comp)

print()
print("composite summed across BOTH refs (lower = better):")
for (log_mu, th), tot in sorted(grid_totals.items(), key=lambda kv: kv[1]):
    print(f"  log_mu={log_mu:>5}  theta={th:>4}  ->  {tot:.3f}")
PY
