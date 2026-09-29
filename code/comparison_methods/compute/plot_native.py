"""Plot recorded direct-native measurements; never execute simulation methods."""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--root", type=Path, required=True)
args = parser.parse_args()
root = args.root
data = pd.read_csv(root / "measurements.csv")
ok = data.loc[data.status.eq("ok")].copy()
if ok.empty:
    raise SystemExit("No successful measurements to plot")
ok["measured_total_seconds"] = ok[
    ["generation_seconds", "input_seconds", "reference_generation_seconds"]
].sum(axis=1, min_count=3)
out = root / "figures"
out.mkdir(exist_ok=True)
ok.to_csv(out / "plotted_measurements.csv", index=False)
methods = {"albis": ("ALBIS", "#4C8FD5"),
           "sccube": ("scCube", "#B9C0C7"),
           "spider": ("SPIDER", "#8E63C7")}
plt.rcParams.update({"font.size": 9, "axes.titlesize": 10,
                     "axes.labelsize": 9, "legend.fontsize": 8,
                     "axes.spines.top": False, "axes.spines.right": False,
                     "savefig.dpi": 300, "svg.fonttype": "none",
                     "pdf.fonttype": 42})


def save(fig, name):
    for ext in ("png", "pdf", "svg"):
        fig.savefig(out / f"{name}.{ext}", bbox_inches="tight")
    plt.close(fig)


fig, axes = plt.subplots(2, 2, figsize=(10, 7))
panels = [("generation_seconds", "Native generation", "Seconds (log scale)"),
          ("measured_total_seconds", "Generation + input + Splatter", "Seconds (log scale)"),
          ("process_peak_rss_gib", "Whole method process peak memory", "Peak RSS (GiB)"),
          ("n_molecules", "ALBIS molecular output", "Molecules (millions)")]
for ax, (column, title, ylabel) in zip(axes.flat, panels):
    for method, (label, color) in methods.items():
        values = ok.loc[ok.method.eq(method)].dropna(subset=[column])
        if values.empty:
            continue
        scale = 1e6 if column == "n_molecules" else 1
        summary = values.groupby("n_cells")[column].agg(["mean", "std", "count"])
        ax.plot(summary.index, summary["mean"] / scale, color=color, lw=1.5, label=label)
        ax.scatter(values.n_cells, values[column] / scale, color=color, s=28, zorder=3)
        replicated = summary.loc[summary["count"].gt(1)]
        if not replicated.empty:
            ax.errorbar(replicated.index, replicated["mean"] / scale,
                        yerr=replicated["std"] / scale, fmt="none", color=color, capsize=3)
    ax.set_xscale("log")
    sizes = sorted(ok.n_cells.unique())
    ax.set_xticks(sizes, [f"{n / 1000:g}k" for n in sizes])
    ax.minorticks_off()
    if column.endswith("seconds"):
        ax.set_yscale("log")
    else:
        ax.set_ylim(bottom=0)
    ax.set(title=title, xlabel="Requested tissue cells", ylabel=ylabel)
    ax.grid(axis="y", color="#DDDDDD", lw=0.5)
    ax.legend(frameon=False)
fig.suptitle("Preliminary native compute benchmark · 556 genes · one CPU thread", y=0.99)
fig.text(0.06, 0.015,
         "Points: successful runs; lines: means (descriptive only). Single-seed pilot has no error bars.\n"
         "scCube includes training. Total excludes imports/startup; memory excludes separate reference process.\n"
         "Native outputs differ: ALBIS molecules, scCube cell expression, SPIDER native expression view. No aggregation/export.",
         fontsize=8)
fig.tight_layout(rect=(0, 0.11, 1, 0.96))
save(fig, "native_scaling")

# Use the largest size completed by all methods; 600k is unavailable in the pilot.
shared = set.intersection(*(set(ok.loc[ok.method.eq(m), "n_cells"]) for m in methods))
if shared:
    size = max(shared)
    stages = ["reference_generation_seconds", "input_seconds", "generation_seconds"]
    fig, axes = plt.subplots(1, 3, figsize=(10, 3.9))
    for ax, (method, (label, color)) in zip(axes, methods.items()):
        values = ok.loc[ok.method.eq(method) & ok.n_cells.eq(size), stages].mean()
        bars = ax.barh(["Splatter reference", "Input / preprocessing", "Native generation"],
                       values, color=["#B9C0C7", "#7393B3", color])
        ax.bar_label(bars, labels=[f"{v:.2f} s" for v in values], padding=4, fontsize=8)
        ax.set_xlim(0, max(values) * 1.4)
        ax.invert_yaxis()
        ax.set(title=f"{label} · total {values.sum():.2f} s", xlabel="Seconds (panel-specific scale)")
        ax.grid(axis="x", color="#DDDDDD", lw=0.5)
    fig.suptitle(f"Preliminary runtime breakdown · {size:,} cells", y=0.99)
    fig.text(0.03, 0.025,
             "Measured stages only; imports/startup excluded. scCube training and synthesis are combined.\n"
             "Splatter cost is charged once to each competitor. 600k breakdown awaits completed runs.", fontsize=8)
    fig.tight_layout(rect=(0, 0.13, 1, 0.94))
    save(fig, "native_stage_breakdown")

counts = data.groupby(["method", "status"]).size().unstack(fill_value=0)
counts.to_csv(out / "run_status_counts.csv")
(out / "README.md").write_text(
    "# Preliminary direct-native compute plots\n\n"
    "Source: ../measurements.csv. Only status=ok records are plotted.\n"
    "Run status counts, including unsubmitted points, are in run_status_counts.csv.\n"
    "PNG, PDF and SVG versions are provided. No extrapolation is performed.\n\n"
    "Runtime total = native generation + input/preprocessing + Splatter generation. "
    "It excludes imports and process startup. scCube generation includes training. "
    "Memory is the whole method process peak RSS, excluding the separate reference process. "
    "Outputs are different native representations, not equivalent outputs. "
    "SPIDER retains its native expression view; no forced materialization is timed.\n\n"
    "The pilot has one seed per size, so uncertainty cannot be estimated. "
    "The breakdown uses the largest size completed by all three methods; "
    "the planned 600k panel must wait for those measurements.\n")
print(out)
