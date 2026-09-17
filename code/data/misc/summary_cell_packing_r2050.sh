#!/bin/bash
# Tabulate sweep_cell_packing_r2050.sh: --n-cells 75000/200000/600000 at the
# fixed sphere_r_um=2050 disc, cell's dispersion knobs otherwise held fixed
# (log_mu=-2.3, theta=0.40, jitter=0.15, batch_sigma=1.5).
#
# Baseline row (n_cells=24207, the current production config) is read
# straight from data/count_distribution/figure_2/cell_vs_*/ rather than
# regenerated. One row per (n_cells, ref) with sim/real values + ratios for
# theta_hat, total_counts_median, genes_per_cell_median, matrix_zero_frac,
# and a composite = sum of |ln(ratio)| over {theta_hat, total_counts_median,
# matrix_zero_frac} (same scoring convention as the other Figure 2 sweeps).
#
# n_obs is also shown per row -- the whole point of this sweep is whether a
# bigger per-slice sample (currently only 3,544 cells at slice_id=5) changes
# the picture.
#
# Read-only; no jobs.
set -euo pipefail
cd /dcs04/hicks/data/Jan/sim_project/sim_paper/data

python3 - <<'PY'
import json, math, os

SWEEP = "count_distribution/sweeps/cell_packing_sweep_r2050"
BASELINE_DIR = "count_distribution/figure_2"
N_CELLS_GRID = ["24207", "75000", "200000", "600000"]
REALS = ["lung_cancer", "non_diseased_lung"]

def load(path):
    with open(path) as fh:
        return json.load(fh)

def summ(n_cells, mode, real):
    if n_cells == "24207":
        return f"{BASELINE_DIR}/cell_vs_{real}/{mode}/comparison_summary.json"
    return f"{SWEEP}/cell_packing_sweep_ncells{n_cells}_{mode}_vs_{real}/comparison_summary.json"

def sr(j, key):
    d = j[key]; s = d["cell"]; r = [v for k, v in d.items() if k != "cell"][0]
    return s, r, (s / r if r else float("nan"))

hdr = (f"{'n_cells':>9} {'ref':>18} | {'n_obs (slice)':>13} | "
       f"{'thetaHat s/r (rt)':>22} | {'tcMed s/r (rt)':>22} | "
       f"{'gpc s/r':>13} | {'zero s/r (rt)':>20} | {'compΣ3':>7}")
print(hdr); print("-" * len(hdr))

grid_totals = {}
for n_cells in N_CELLS_GRID:
    per_ref_comp = []
    for real in REALS:
        p = summ(n_cells, "hvg_matched", real)
        if not os.path.exists(p):
            print(f"{n_cells:>9} {real:>18} | MISSING ({p})")
            per_ref_comp.append(None); continue
        j = load(p)
        nobs = j["n_obs"]["cell"]
        ths, thr, thratio = sr(j, "theta_hat")
        tcs, tcr, tcratio = sr(j, "total_counts_median")
        gps, gpr, _ = sr(j, "genes_per_cell_median")
        zs, zr, zratio = sr(j, "matrix_zero_frac")
        comp = sum(abs(math.log(x)) for x in (thratio, tcratio, zratio))
        per_ref_comp.append(comp)

        print(f"{n_cells:>9} {real:>18} | {nobs:>13} | "
              f"{ths:6.3f}/{thr:6.3f} ({thratio:4.2f}) | "
              f"{tcs:7.0f}/{tcr:7.0f} ({tcratio:4.2f}) | "
              f"{gps:5.0f}/{gpr:5.0f} | "
              f"{zs:5.3f}/{zr:5.3f} ({zratio:4.2f}) | "
              f"{comp:7.3f}")
    if all(c is not None for c in per_ref_comp):
        grid_totals[n_cells] = sum(per_ref_comp)

print()
print("composite (hvg_matched) summed across BOTH refs (lower = better):")
for n_cells, tot in sorted(grid_totals.items(), key=lambda kv: kv[1]):
    print(f"  n_cells={n_cells:>7}  ->  {tot:.3f}")
print()
print("Also eyeball qc_filtered/mean_variance_compare.png and mean_dropout_compare.png")
print(f"per n_cells under {SWEEP}/ -- the composite doesn't capture visual scatter/noise.")
PY
