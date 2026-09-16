#!/usr/bin/env python3
"""
Figure 4A (step 2) -- quantify how well STAIR does relative to SLAT, SPACEL,
and STAligner on the pairwise slice-alignment translation-error task.

Recycles already-computed results from the multi-sample-alignment-benchmark
repo -- nothing is rerun here. For every algorithm x sectioning-axis, reads
    evaluation_result/simulated_pairwise_relative_metrics/<algo>/
        cell_obs_sectioned_<axis>_relative_alignment_metrics_all_pairs.csv
(45 rows = all C(10,2) slice pairs; produced by evaluation_moving_metrics.py
/ code/vani/INSPIRE/run_relative_metrics_all_simulated.py). All four
algorithms here output aligned coordinates in the micron frame (anchor_scale
~= 1), so moving_relative_mse_anchor_corrected_xy is already the correct
cross-method column -- no GPSA-style true-frame rescaling needed.

Produces, under sim_paper/data/figure_4/alignment_translation_error/ :
    figure4a_metric_summary.csv      per algorithm x axis: median/mean/IQR of
                                      relative RMSE (um) and |translation| (um)
    figure4a_metric_comparison.png   boxplots of per-pair relative RMSE,
                                      one panel per axis, algorithms compared
                                      within each panel
"""

import argparse
import os

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

METRIC_ROOT = (
    "/dcs04/hicks/data/multi-sample-alignment-benchmark/evaluation_result/"
    "simulated_pairwise_relative_metrics"
)
OUT = "/dcs04/hicks/data/Jan/sim_project/sim_paper/data/figure_4/alignment_translation_error"

ALGO_DIRNAME = {"STAIR": "stair", "SLAT": "slat", "SPACEL": "spacel", "STAligner": "staligner"}
ALGO_COLOR = {"STAIR": "#4878d0", "SLAT": "#ee854a", "SPACEL": "#6acc64", "STAligner": "#d65f5f"}
ALGOS = ["STAIR", "SLAT", "SPACEL", "STAligner"]
AXES = ["X", "Y", "Z"]
DATASET = "cell_obs_sectioned"


def load(dataset, axis):
    rows = []
    for algo in ALGOS:
        p = os.path.join(
            METRIC_ROOT, ALGO_DIRNAME[algo],
            f"{dataset}_{axis}_relative_alignment_metrics_all_pairs.csv",
        )
        df = pd.read_csv(p)
        df = df[df["status"] == "success"].copy()
        df["algorithm_label"] = algo
        df["axis"] = axis
        rows.append(df)
    return pd.concat(rows, ignore_index=True)


def summarize(df):
    df = df.copy()
    df["relative_rmse_um"] = np.sqrt(df["moving_relative_mse_anchor_corrected_xy"])
    g = df.groupby(["axis", "algorithm_label"])
    summary = g.agg(
        n_pairs=("relative_rmse_um", "size"),
        median_relative_rmse_um=("relative_rmse_um", "median"),
        mean_relative_rmse_um=("relative_rmse_um", "mean"),
        q25_relative_rmse_um=("relative_rmse_um", lambda x: x.quantile(0.25)),
        q75_relative_rmse_um=("relative_rmse_um", lambda x: x.quantile(0.75)),
        median_translation_norm_um=("moving_relative_translation_norm", "median"),
        median_rotation_deg=("moving_relative_rotation_degrees", lambda x: x.abs().median()),
    ).reset_index()
    return summary.sort_values(["axis", "median_relative_rmse_um"])


def plot_boxplots(df, outpath):
    df = df.copy()
    df["relative_rmse_um"] = np.sqrt(df["moving_relative_mse_anchor_corrected_xy"])

    fig, axes = plt.subplots(1, len(AXES), figsize=(4.2 * len(AXES), 4.6), sharey=True)
    for ax_i, axis in enumerate(AXES):
        ax = axes[ax_i]
        sub = df[df["axis"] == axis]
        data = [sub.loc[sub["algorithm_label"] == a, "relative_rmse_um"].values for a in ALGOS]
        bp = ax.boxplot(data, tick_labels=ALGOS, showfliers=True, patch_artist=True, widths=0.6)
        for patch, algo in zip(bp["boxes"], ALGOS):
            patch.set_facecolor(ALGO_COLOR[algo])
            patch.set_alpha(0.75)
        for median in bp["medians"]:
            median.set_color("black")
        ax.set_title(f"axis {axis}", fontsize=12)
        ax.set_xlabel("")
        ax.tick_params(axis="x", rotation=30)
        if ax_i == 0:
            ax.set_ylabel("relative alignment RMSE (µm)\n(anchor-corrected, moving slice vs. truth)")
        ax.set_yscale("log")
    fig.suptitle("Figure 4A -- pairwise slice-alignment translation error, cell resolution "
                  "(45 slice pairs / algorithm / axis)", fontsize=12)
    fig.patch.set_facecolor("white")
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(outpath, dpi=200, facecolor="white")
    plt.close(fig)
    print("wrote", outpath)


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", default=DATASET)
    ap.add_argument("--axes", nargs="+", default=AXES, choices=AXES)
    ap.add_argument("--outdir", default=OUT)
    return ap.parse_args()


if __name__ == "__main__":
    args = parse_args()
    AXES = args.axes
    os.makedirs(args.outdir, exist_ok=True)

    all_df = pd.concat([load(args.dataset, axis) for axis in AXES], ignore_index=True)

    summary = summarize(all_df)
    summary_path = os.path.join(args.outdir, "figure4a_metric_summary.csv")
    summary.to_csv(summary_path, index=False)
    print("wrote", summary_path)
    print(summary.to_string(index=False))

    plot_boxplots(all_df, os.path.join(args.outdir, "figure4a_metric_comparison.png"))
