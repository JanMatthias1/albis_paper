#!/bin/bash
# Tabulate the cell log_mu x theta joint grid (sweep_cell_logmu_theta_joint_xenium.sh)
# against the two Xenium real refs (lung_cancer, non_diseased_lung). Prints
# one row per (log_mu, theta, real ref) with sim/real values and ratios for
# theta_hat, total_counts_median, genes_per_cell_median, matrix_zero_frac,
# plus a composite = sum of |ln(ratio)| over {theta_hat, total_counts_median,
# matrix_zero_frac}, and that composite summed across BOTH refs (the same
# scoring summary_spot_logmu_theta_joint.sh uses for spot). Read-only; no jobs.
set -euo pipefail
cd /dcs04/hicks/data/Jan/sim_project/sim_paper/data/count_distribution/sweeps

python3 - <<'PY'
import json, math, os

BASE = "cell_joint_logmu_theta_xenium"
REALS = ["lung_cancer", "non_diseased_lung"]
LOG_MUS = ["-2.5", "-2.3", "-2.0"]
THETAS = ["0.30", "0.40", "0.50"]

def load(path):
    with open(path) as fh:
        return json.load(fh)

def rowdir(log_mu, theta, real):
    tag = f"cell_joint_logmu_{log_mu}_theta_{theta}"
    return f"{BASE}/{tag}_hvg_matched_vs_{real}/comparison_summary.json"

hdr = f"{'log_mu':>7} {'theta':>6} {'ref':>18} | {'thetaHat s/r (ratio)':>26} | {'tcMed s/r (ratio)':>24} | {'gpc s/r':>13} | {'zero s/r (ratio)':>22} | {'compΣ3':>7}"
print(hdr); print("-" * len(hdr))

grid_totals = {}
for log_mu in LOG_MUS:
    for th in THETAS:
        per_ref_comp = []
        for real in REALS:
            p = rowdir(log_mu, th, real)
            if not os.path.exists(p):
                print(f"{log_mu:>7} {th:>6} {real:>18} | MISSING ({p})")
                per_ref_comp.append(None); continue
            j = load(p)
            def sr(key):
                d = j[key]; s = d["cell"]; r = [v for k, v in d.items() if k != "cell"][0]
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
    marker = "  <- current Figure 2 config" if (log_mu, th) == ("-2.3", "0.40") else ""
    print(f"  log_mu={log_mu:>5}  theta={th:>4}  ->  {tot:.3f}{marker}")
PY
