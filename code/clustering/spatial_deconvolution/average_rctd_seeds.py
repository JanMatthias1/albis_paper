#!/usr/bin/env python
"""
Figure 3 -- RCTD spot deconvolution, 3-seed average.

Averages the per-spot, per-cell-type estimated fractions across the 3
independent-reference-seed RCTD runs (RCTD/spot_seed{271828,314159,999999}/,
from run_rctd_spot_independent_seed.sh) into a single combined estimate, then
recomputes the same metrics run_rctd_spot.R computes (per-spot/per-cell-type
Pearson r + RMSE, overall summary) against that average -- a single reference
draw is noisy (canonical r=0.73 vs the 3 seeds' individual 0.68/0.67/0.72),
and averaging estimates across independent reference draws is the standard
way to get a less noisy point estimate of RCTD's true performance on this
query, the same way the seeds themselves were used to bound it.

All 3 seed runs scored the exact same 9889 spot_ids (verified directly,
2026-09-18) -- query (spot) and its QC/min-UMI-driven dropouts are identical
across runs, only the reference differs -- so the average is a clean
elementwise mean, no alignment/intersection needed.

Output, under RCTD/spot_seed_avg/ (parallel to the individual spot_seed*/
dirs and the canonical spot/):
  estimated_fractions_wide.csv  mean of the 3 seeds' est. fractions, same
                                 wide format as run_rctd_spot.R's output
  per_spot_metrics.csv, per_celltype_metrics.csv, metrics_summary.json
                                 recomputed from the averaged estimate,
                                 identical formulas to run_rctd_spot.R
  est_vs_true_scatter.png       same diagnostic scatter as the R script's,
                                 rebuilt in matplotlib (no per-run R rerun)

Run this first, then plot_rctd_results.py <this OUT_DIR> for the spatial
zoom / error boxplot figures (it only needs estimated_fractions_wide.csv +
per_celltype_metrics.csv, both written here).
"""

import json
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import anndata as ad

RCTD_DIR = "/dcs04/hicks/data/Jan/sim_project/sim_paper/data/figure_3/spatial_deconvolution/RCTD"
SEEDS = ["271828", "314159", "999999"]
SEED_DIRS = [os.path.join(RCTD_DIR, f"spot_seed{s}") for s in SEEDS]
OUT_DIR = os.path.join(RCTD_DIR, "spot_seed_avg")
SPOT_H5AD = "/dcs04/hicks/data/Jan/sim_project/sim_paper/data/figure_2/smaller_sphere/data/packing_pf0p04_log_mu_-2.25_sigma1.0_theta_0.25_jitter0.10_dsf_bsigma03/simulation_spot_z_qc.h5ad"
CELL_TYPES = [f"type{i}" for i in range(1, 9)]

os.makedirs(OUT_DIR, exist_ok=True)


def safe_cor(a, b):
    """np.corrcoef(a, b) is NaN-on-nothing but returns [[nan,nan],[nan,nan]]
    when either input is constant -- mirror run_rctd_spot.R's safe_cor
    (NA if either vector has zero sd) explicitly rather than relying on
    that propagating correctly."""
    if np.std(a) == 0 or np.std(b) == 0:
        return np.nan
    return float(np.corrcoef(a, b)[0, 1])


def load_average_estimates():
    frames = []
    for d, seed in zip(SEED_DIRS, SEEDS):
        df = pd.read_csv(os.path.join(d, "estimated_fractions_wide.csv"))
        df = df.set_index("spot_id")[CELL_TYPES]
        frames.append(df)
    ids = [set(f.index) for f in frames]
    assert all(s == ids[0] for s in ids), "spot_id sets differ across seeds -- cannot average elementwise"
    print(f"[load] {len(frames[0])} spots, {len(SEEDS)} seeds ({', '.join(SEEDS)}), spot_id sets identical")

    stacked = np.stack([f.loc[frames[0].index].to_numpy() for f in frames], axis=0)
    avg = pd.DataFrame(stacked.mean(axis=0), index=frames[0].index, columns=CELL_TYPES)
    avg.index.name = "spot_id"
    return avg


def load_true_fractions(spot_ids):
    spot = ad.read_h5ad(SPOT_H5AD)
    true_frac = pd.DataFrame(
        np.asarray(spot.obsm["cell_type_frac_true"]), columns=CELL_TYPES, index=spot.obs_names
    )
    # Same "spot_<1-based positional index>" mapping run_rctd_spot.R /
    # plot_rctd_results.py use to join RCTD's spot_ids back to obs_names.
    positional_idx = pd.Index(spot_ids).str.replace("spot_", "", regex=False).astype(int) - 1
    true_frac_ordered = true_frac.iloc[positional_idx.to_numpy()]
    true_frac_ordered.index = spot_ids
    return true_frac_ordered


def compute_and_write_metrics(est, truth):
    est_arr = est[CELL_TYPES].to_numpy()
    truth_arr = truth[CELL_TYPES].to_numpy()

    per_spot_cor = np.array([safe_cor(est_arr[i], truth_arr[i]) for i in range(est_arr.shape[0])])
    per_spot_rmse = np.sqrt(np.mean((est_arr - truth_arr) ** 2, axis=1))
    per_spot_df = pd.DataFrame({
        "spot_id": est.index,
        "pearson_r": per_spot_cor,
        "rmse": per_spot_rmse,
    })
    per_spot_df.to_csv(os.path.join(OUT_DIR, "per_spot_metrics.csv"), index=False, quoting=2)

    per_type_cor = np.array([safe_cor(est_arr[:, j], truth_arr[:, j]) for j in range(len(CELL_TYPES))])
    per_type_rmse = np.sqrt(np.mean((est_arr - truth_arr) ** 2, axis=0))
    per_type_df = pd.DataFrame({
        "cell_type": CELL_TYPES,
        "pearson_r": per_type_cor,
        "rmse": per_type_rmse,
    })
    per_type_df.to_csv(os.path.join(OUT_DIR, "per_celltype_metrics.csv"), index=False, quoting=2)

    overall_cor = safe_cor(est_arr.ravel(), truth_arr.ravel())
    overall_rmse = float(np.sqrt(np.mean((est_arr - truth_arr) ** 2)))
    # Read from the seed runs (same query for all) instead of a hardcoded
    # count -- was 10034, which silently went stale on the 2026-09-23 spot retune.
    seed_totals = {json.load(open(os.path.join(d, "metrics_summary.json")))["n_spots_total"] for d in SEED_DIRS}
    assert len(seed_totals) == 1, f"seed runs disagree on n_spots_total: {seed_totals}"
    summary = {
        "n_spots_scored": int(est_arr.shape[0]),
        "n_spots_total": seed_totals.pop(),
        "doublet_mode": "full",
        "n_seeds_averaged": len(SEEDS),
        "seeds_averaged": SEEDS,
        "overall_pearson_r": round(overall_cor, 4),
        "overall_rmse": round(overall_rmse, 4),
        "mean_per_spot_pearson_r": round(float(np.nanmean(per_spot_cor)), 4),
        "mean_per_spot_rmse": round(float(np.mean(per_spot_rmse)), 4),
    }
    with open(os.path.join(OUT_DIR, "metrics_summary.json"), "w") as f:
        json.dump(summary, f, indent=2)

    print("\n=== 3-SEED AVERAGE RESULTS ===")
    print(json.dumps(summary, indent=2))
    print(per_type_df.to_string(index=False))
    return per_type_df


def plot_scatter(est, truth):
    fig, axes = plt.subplots(2, 4, figsize=(9, 5), sharex=True, sharey=True)
    for ax, ct in zip(axes.flat, CELL_TYPES):
        ax.scatter(truth[ct], est[ct], s=3, alpha=0.15, color="#1f77b4", linewidths=0)
        ax.plot([0, 1], [0, 1], linestyle="--", color="grey", linewidth=1)
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.set_aspect("equal")
        ax.set_title(ct, fontsize=10)
    for ax in axes[-1, :]:
        ax.set_xlabel("true fraction")
    for ax in axes[:, 0]:
        ax.set_ylabel("RCTD estimated\nfraction (3-seed avg)")
    fig.suptitle("RCTD spot deconvolution (3-seed reference average): estimated vs. true cell-type fraction",
                  fontsize=11, fontweight="bold")
    fig.tight_layout()
    out = os.path.join(OUT_DIR, "est_vs_true_scatter.png")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print("wrote", out)


if __name__ == "__main__":
    avg_est = load_average_estimates()
    avg_est.reset_index().to_csv(os.path.join(OUT_DIR, "estimated_fractions_wide.csv"), index=False, quoting=2)
    print("wrote", os.path.join(OUT_DIR, "estimated_fractions_wide.csv"))

    truth = load_true_fractions(avg_est.index)
    compute_and_write_metrics(avg_est, truth)
    plot_scatter(avg_est, truth)
