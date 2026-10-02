"""Semantic manuscript colors and Figure 3 typography.

No style is applied on import, so Figure 1/2 remain independent.
"""
import re

import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment

DOMAIN_COLORS = dict(zip([f"D{i}" for i in range(6)],
    ["#8E63C7", "#4C8FD5", "#43B7A5", "#72C69C", "#F2A65A", "#E65F5C"]))
CELLTYPE_COLORS = dict(zip([f"Cell Type {i}" for i in range(1, 9)],
    ["#6E86D6", "#D98B5F", "#7CBD7A", "#B085D6", "#D46F92", "#8FA06A", "#D1B85A", "#68B9C0"]))
MODALITY_COLORS = {"Cell": "#8E63C7", "Bin (16 µm)": "#43B7A5", "Spot": "#F2A65A"}
# Both bin sizes share the bin hue; their resolution remains explicit in labels.
MODALITY_LOOKUP = {"cell": MODALITY_COLORS["Cell"], "bin": MODALITY_COLORS["Bin (16 µm)"],
                   "bin16um": MODALITY_COLORS["Bin (16 µm)"], "spot": MODALITY_COLORS["Spot"]}
SLICE_COLORS = ["#7393B3", "#D79A83", "#8FB9A8", "#B29BC8", "#D49BAD",
                "#A7A7A7", "#C6B879", "#82B6BA", "#A9A0CB", "#B5C98A"]
BOX_EDGE, BOX_MEDIAN, GRID_COLOR = "#4A4A4A", "#2F2F2F", "#DDDDDD"
UNMATCHED_COLOR = "#B9C0C7"
TITLE_SIZE, LABEL_SIZE, TICK_SIZE, LEGEND_SIZE = 17, 16, 12, 12
LEGEND_MARKERSIZE, LEGEND_TITLE_SIZE = 12, 13
PANEL_FIGSIZE = (12, 6.8)
# Preserve the full Figure 3 canvas, including four-row mapped legends.
PANEL_EXPORT_BOTTOM = 0.0
PANEL_MARGINS = dict(left=0.08, right=0.98, top=0.90, bottom=0.34, wspace=0.28)


def apply_style():
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": LABEL_SIZE,
        "axes.titlesize": TITLE_SIZE, "axes.titleweight": "bold", "axes.labelsize": LABEL_SIZE,
        "xtick.labelsize": TICK_SIZE, "ytick.labelsize": TICK_SIZE,
        "legend.fontsize": LEGEND_SIZE, "legend.title_fontsize": LEGEND_TITLE_SIZE,
        "figure.titlesize": TITLE_SIZE, "figure.titleweight": "bold", "savefig.dpi": 300})


def pretty_label(key):
    if key.startswith("leiden_") or key.startswith("louvain_"):
        return "Cluster label"
    return {"slice_id": "Slice ID", "domain_true": "Domain", "cell_type_true": "Cell type",
            "cluster_label": "Cluster label"}.get(key, key.replace("_", " ").title())


def category_label(value, key):
    value = str(value)
    if key == "cell_type_true":
        match = re.fullmatch(r"(?:type|Cell Type\s*)([1-8])", value, re.I)
        if match:
            return f"Cell Type {match[1]}"
    return value


def category_color(value, key):
    value = str(value)
    if value.lower() in {"unassigned", "nan", "none", "-1"}:
        return UNMATCHED_COLOR
    if key == "domain_true":
        return DOMAIN_COLORS[value]
    if key == "cell_type_true":
        return CELLTYPE_COLORS[category_label(value, key)]
    if key == "slice_id":
        idx = int(float(value))
        if not 0 <= idx < len(SLICE_COLORS):
            raise ValueError(f"No manuscript slice color defined for {value}")
        return SLICE_COLORS[idx]
    return UNMATCHED_COLOR


def category_order(values):
    def order(s):
        match = re.search(r"\d+(?:\.\d+)?$", s)
        return (0, float(match[0]), s) if match else (1, 0, s)
    return sorted(set(map(str, values)), key=order)


def matched_labels(obs, pred_key, true_key):
    """Maximum-overlap one-to-one matching, used for display only.

    Extra/unassigned clusters remain gray. Never modifies obs or the ARI.
    """
    table = pd.crosstab(obs[pred_key].astype(str), obs[true_key].astype(str))
    truth = [v for v in category_order(table.columns)
             if v.lower() not in {"unassigned", "nan", "none", "-1"}]
    pred = [v for v in category_order(table.index)
            if v.lower() not in {"unassigned", "nan", "none", "-1"}]
    table = table.reindex(index=pred, columns=truth)
    rows, cols = linear_sum_assignment(-table.to_numpy())
    return {table.index[i]: table.columns[j] for i, j in zip(rows, cols) if table.iloc[i, j] > 0}


def scatter_colors(obs, key, true_key=None):
    """Return matplotlib scatter kwargs, legend labels and exact color swatches."""
    values = obs[key]
    if key not in {"slice_id", "domain_true", "cell_type_true", "cluster_label"} and true_key is None \
            and pd.api.types.is_numeric_dtype(values) and values.nunique() > 20:
        return dict(c=values.to_numpy(), cmap="viridis"), None, None
    categories = category_order(values.astype(str))
    labels = [category_label(v, key) for v in categories]
    if true_key is not None:
        mapping = matched_labels(obs, key, true_key)
        # Show the correspondence explicitly, as in the saved Figure 3 panels.
        labels = [f"{v} → {category_label(mapping[v], true_key)}" if v in mapping
                  else f"{v} (unmatched)" for v in categories]
        colors = [category_color(mapping[v], true_key) if v in mapping else UNMATCHED_COLOR
                  for v in categories]
    else:
        colors = [category_color(v, key) for v in categories]
    cmap = ListedColormap(colors)
    codes = pd.Categorical(values.astype(str), categories=categories).codes
    return dict(c=codes, cmap=cmap, norm=BoundaryNorm(np.arange(len(colors) + 1) - 0.5,
                                                     len(colors))), labels, colors


def legend_handles(labels, colors):
    return [plt.Line2D([], [], linestyle="none", marker="o", markerfacecolor=color,
                      markeredgecolor="white", markersize=LEGEND_MARKERSIZE, label=label)
            for label, color in zip(labels, colors)]
