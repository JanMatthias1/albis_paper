#!/bin/bash
# Tabulate sweep_bin16um_sphere_r6000.sh: --n-cells 3008662/5641241/15043310
# at sphere_r_um=6000 (packing fractions 0.8%/1.5%/4.0%), bin16um's
# dispersion knobs otherwise held fixed (log_mu=-2.5, batch_sigma=0.7).
#
# Baseline row (current production: n_cells=600000 @ sphere_r_um=2050, same
# 4.0% fraction as the top of this grid but realized at the small disc) is
# read straight from data/figure_2/packing_pf0p04_bin16um_log_mu_-2.5_bsigma07/
# and data/count_distribution/figure_2/bin_vs_*_16um/ rather than regenerated.
#
# empty_frac (n_fully_empty / n_obs_before_qc, from each tag's
# simulation_bin_z_qc_summary.json) is the KEY column: the whole point of
# packing_pf0p04 was fixing bin's ~77% empty-bin problem at low packing
# fraction post-realwindow (see generate_simulation_noisy.py header comment
# + DATA_VERSIONS.md). If growing the disc to sphere_r_um=6000 at matched
# packing fraction reproduces the baseline's empty_frac, the "grow the
# shared disc instead of shrinking cell" approach is validated; if it climbs
# back toward 66-77%, the fix does NOT carry over and this approach should
# be dropped. theta_hat/total_counts_median (qc_filtered vs
# breast_cancer_visium_hd/human_pancreas_visium_hd) are shown as the
# secondary count-fidelity check.
#
# Read-only; no jobs.
set -euo pipefail
cd /dcs04/hicks/data/Jan/sim_project/sim_paper/data

python3 - <<'PY'
import json, math, os

SWEEP_NOISY = "noisy"
SWEEP_CMP = "count_distribution/sweeps/bin16um_sphere_r6000"
BASELINE_QC = "figure_2/packing_pf0p04_bin16um_log_mu_-2.5_bsigma07/simulation_bin_z_qc_summary.json"
BASELINE_CMP_DIR = "count_distribution/figure_2"
GRID = [("600000@2050(base)", None), ("3008662@6000", "3008662"),
        ("5641241@6000", "5641241"), ("15043310@6000", "15043310")]
REALS = ["breast_cancer_visium_hd", "human_pancreas_visium_hd"]
BASELINE_REAL_SUFFIX = "_16um"  # baseline OUT_ROOT naming quirk, see bin16um_vs_*.sh

def load(path):
    with open(path) as fh:
        return json.load(fh)

def qc_path(n_cells):
    if n_cells is None:
        return BASELINE_QC
    return f"{SWEEP_NOISY}/bin16um_sphere_r6000_ncells{n_cells}/simulation_bin_z_qc_summary.json"

def cmp_path(n_cells, real):
    if n_cells is None:
        return f"{BASELINE_CMP_DIR}/bin_vs_{real}{BASELINE_REAL_SUFFIX}/qc_filtered/comparison_summary.json"
    return f"{SWEEP_CMP}/bin16um_sphere_r6000_ncells{n_cells}_qc_filtered_vs_{real}/comparison_summary.json"

def sr(j, key):
    d = j[key]; s = d["bin"]; r = [v for k, v in d.items() if k != "bin"][0]
    return s, r, (s / r if r > 0 else float("nan"))

hdr = f"{'config':>20} | {'empty_frac (n_before->after)':>34}"
print(hdr); print("-" * len(hdr))
empty_fracs = {}
for label, n_cells in GRID:
    p = qc_path(n_cells)
    if not os.path.exists(p):
        print(f"{label:>20} | MISSING ({p})")
        continue
    j = load(p)
    frac = j["n_fully_empty"] / j["n_obs_before_qc"]
    empty_fracs[label] = frac
    print(f"{label:>20} | {frac:6.1%}  ({j['n_obs_before_qc']} -> {j['n_obs_after_qc']})")

print()
hdr2 = (f"{'config':>20} {'ref':>25} | {'thetaHat s/r (rt)':>22} | "
        f"{'tcMed s/r (rt)':>22} | {'compΣ2':>7}")
print(hdr2); print("-" * len(hdr2))
grid_totals = {}
for label, n_cells in GRID:
    per_ref_comp = []
    for real in REALS:
        p = cmp_path(n_cells, real)
        if not os.path.exists(p):
            print(f"{label:>20} {real:>25} | MISSING ({p})")
            per_ref_comp.append(None); continue
        j = load(p)
        ths, thr, thratio = sr(j, "theta_hat")
        tcs, tcr, tcratio = sr(j, "total_counts_median")
        comp = sum(abs(math.log(x)) for x in (thratio, tcratio) if x and x > 0)
        per_ref_comp.append(comp)
        print(f"{label:>20} {real:>25} | "
              f"{ths:6.3f}/{thr:6.3f} ({thratio:4.2f}) | "
              f"{tcs:7.1f}/{tcr:7.1f} ({tcratio:4.2f}) | "
              f"{comp:7.3f}")
    if all(c is not None for c in per_ref_comp):
        grid_totals[label] = sum(per_ref_comp)

print()
print("composite (qc_filtered, theta_hat+total_counts) summed across BOTH refs (lower = better):")
for label, tot in sorted(grid_totals.items(), key=lambda kv: kv[1]):
    print(f"  {label:>20}  ->  {tot:.3f}")
print()
print(f"PASS/FAIL signal: empty_frac at 15043310@6000 should land near the")
print(f"600000@2050 baseline ({empty_fracs.get('600000@2050(base)', float('nan')):.1%}) if the")
print("packing-fraction fix carries over to the bigger shared disc.")
PY
