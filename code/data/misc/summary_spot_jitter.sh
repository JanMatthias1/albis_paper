#!/bin/bash
# Tabulate the spot theta_jitter sweep (sweep_spot_jitter_probe.sh):
# --theta-jitter 0.10 / 0.15 / 0.25 at the shared spot config otherwise
# fixed (log_mu=-2.0, theta=0.25, batch_sigma=0.3). One row per (jitter, ref)
# with sim/real values + ratios for theta_hat, total_counts_median,
# genes_per_cell_median, matrix_zero_frac, a composite = sum of |ln(ratio)|
# over {theta_hat, total_counts_median, matrix_zero_frac}, and that composite
# summed across BOTH refs (same score the other probe sweeps optimize).
#
# theta_hat is shown for BOTH modes: hvg_matched (the apples-to-apples number
# the composite uses) and qc_filtered (the mode the figure actually reads
# spot's mean_variance panel from). The point of the sweep is the two-cloud
# collapse, which is visual -- check the qc_filtered mean_variance_compare.png
# for each jitter by eye; this table is the supporting quantitative read.
#
# Baseline for reference (current shared tag, jitter=1.0):
#   qc_filtered  theta_hat: spot 0.866 | tonsil 4.30 | lymph_node 2.43
#   hvg_matched  theta_hat: spot 0.167 | tonsil 0.344 | lymph_node 0.762
#
# Read-only; no jobs.
set -euo pipefail
cd /dcs04/hicks/data/Jan/sim_project/sim_paper/data/count_distribution/sweeps

python3 - <<'PY'
import json, math, os

BASE = "spot_jitter_sweep_probe"
JITTERS = ["0.10", "0.15", "0.25"]
REALS = ["tonsil_visium", "lymph_node_visium"]

def load(path):
    with open(path) as fh:
        return json.load(fh)

def summ(jit, mode, real):
    return f"{BASE}/spot_jitter_sweep_{jit}_{mode}_vs_{real}/comparison_summary.json"

def sr(j, key):
    d = j[key]; s = d["spot"]; r = [v for k, v in d.items() if k != "spot"][0]
    return s, r, (s / r if r else float("nan"))

hdr = (f"{'jitter':>6} {'ref':>18} | {'thetaHat hvg s/r (rt)':>26} | "
       f"{'thetaHat qcf s/r (rt)':>26} | {'tcMed s/r (rt)':>22} | "
       f"{'gpc s/r':>13} | {'zero s/r (rt)':>20} | {'compΣ3':>7}")
print(hdr); print("-" * len(hdr))

grid_totals = {}
for jit in JITTERS:
    per_ref_comp = []
    for real in REALS:
        ph = summ(jit, "hvg_matched", real)
        pf = summ(jit, "qc_filtered", real)
        if not os.path.exists(ph):
            print(f"{jit:>6} {real:>18} | MISSING ({ph})")
            per_ref_comp.append(None); continue
        jh = load(ph)
        ths, thr, thratio = sr(jh, "theta_hat")
        tcs, tcr, tcratio = sr(jh, "total_counts_median")
        gps, gpr, _ = sr(jh, "genes_per_cell_median")
        zs, zr, zratio = sr(jh, "matrix_zero_frac")
        comp = sum(abs(math.log(x)) for x in (thratio, tcratio, zratio))
        per_ref_comp.append(comp)

        if os.path.exists(pf):
            jf = load(pf)
            fs, fr, fratio = sr(jf, "theta_hat")
            qcf = f"{fs:6.3f}/{fr:6.3f} ({fratio:4.2f})"
        else:
            qcf = "  (missing qc_filtered)  "

        print(f"{jit:>6} {real:>18} | "
              f"{ths:6.3f}/{thr:6.3f} ({thratio:4.2f}) | "
              f"{qcf:>26} | "
              f"{tcs:7.0f}/{tcr:7.0f} ({tcratio:4.2f}) | "
              f"{gps:5.0f}/{gpr:5.0f} | "
              f"{zs:5.3f}/{zr:5.3f} ({zratio:4.2f}) | "
              f"{comp:7.3f}")
    if all(c is not None for c in per_ref_comp):
        grid_totals[jit] = sum(per_ref_comp)

print()
print("composite (hvg_matched) summed across BOTH refs (lower = better):")
for jit, tot in sorted(grid_totals.items(), key=lambda kv: kv[1]):
    print(f"  jitter={jit:>5}  ->  {tot:.3f}")
PY
