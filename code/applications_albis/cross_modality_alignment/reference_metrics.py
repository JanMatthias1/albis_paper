"""Figure 4C: reference-only rigid correction and held-out modality errors.

Coordinates use row vectors: corrected = coordinates @ rotation + translation.
Truth is paired within each observation; no cross-platform point matching is used.
"""
import argparse
import csv
import json
from pathlib import Path

import numpy as np

STAGES = {"unaligned": "spatial", "stair_init": "transform_init",
          "stair_fine": "transform_fine"}
TECHS = ("bin16um", "spot", "cell")


def rigid_fit(source, target):
    source, target = np.asarray(source, float), np.asarray(target, float)
    if source.shape != target.shape or source.ndim != 2 or source.shape[1] != 2:
        raise ValueError("Expected matching N x 2 coordinates")
    if not np.isfinite(source).all() or not np.isfinite(target).all():
        raise ValueError("Coordinates must be finite")
    x, y = source - source.mean(0), target - target.mean(0)
    if len(x) < 3 or min(np.linalg.matrix_rank(x), np.linalg.matrix_rank(y)) < 2:
        raise ValueError("Reference needs at least three non-collinear points")
    u, _, vt = np.linalg.svd(x.T @ y)
    correction = np.eye(2)
    correction[-1, -1] = np.linalg.det(u @ vt)
    rotation = u @ correction @ vt
    return rotation, target.mean(0) - source.mean(0) @ rotation


def evaluate(truth, coordinates, technology, reference):
    truth = np.asarray(truth, float)
    technology = np.asarray(technology, str)
    if truth.shape != (len(technology), 2) or not np.isfinite(truth).all():
        raise ValueError("Truth must contain finite N x 2 coordinates")
    if set(technology) != set(TECHS) or reference not in TECHS:
        raise ValueError("Expected bin16um, spot and cell, and a valid reference")
    ref = technology == reference
    targets = [t for t in TECHS if t != reference]
    result = {"reference_modality": reference, "scored_modalities": targets,
              "units": "um", "allow_scale": False, "allow_reflection": False,
              "combined_definition": "sqrt(mean(per-target-modality MSE)); reference excluded",
              "stages": {}}
    corrected = {}
    for stage, values in coordinates.items():
        values = np.asarray(values, float)
        if values.shape != truth.shape or not np.isfinite(values).all():
            raise ValueError(f"Invalid coordinates for {stage}")
        rotation, translation = rigid_fit(values[ref], truth[ref])
        xy = values @ rotation + translation
        corrected[stage] = xy
        distances = np.linalg.norm(xy - truth, axis=1)
        per_tech = {}
        for tech in TECHS:
            d = distances[technology == tech]
            per_tech[tech] = {"n_obs": len(d), "rmse_um": float(np.sqrt(np.mean(d**2))),
                              "median_um": float(np.median(d)),
                              "p95_um": float(np.quantile(d, 0.95))}
        result["stages"][stage] = {
            "rotation": rotation.tolist(), "translation_um": translation.tolist(),
            "per_modality": per_tech,
            "balanced_target_rmse_um": float(np.sqrt(np.mean([
                per_tech[t]["rmse_um"]**2 for t in targets]))),
            "pooled_target_rmse_um": float(np.sqrt(np.mean(distances[~ref]**2)))}
    return result, corrected


def plot_results(result, corrected, truth, technology, domains, outdir, slice_id):
    import sys
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
    from manuscript_style import DOMAIN_COLORS

    labels = {"bin16um": "Bin (16 µm)", "spot": "Spot", "cell": "Cell"}
    plt.rcParams.update({"font.size": 9, "axes.titlesize": 10, "axes.labelsize": 9,
                         "pdf.fonttype": 42, "ps.fonttype": 42})
    ref = result["reference_modality"]
    stem = f"figure4c_reference_{ref}_slice_{slice_id}"
    fig, ax = plt.subplots(figsize=(6.3, 3.7))
    groups = result["scored_modalities"] + ["balanced"]
    x = np.arange(len(groups))
    for offset, stage, label, color in [(-.19, "unaligned", "Unaligned", "#B9C0C7"),
                                       (.19, "stair_fine", "STAIR aligned", "#4C8FD5")]:
        row = result["stages"][stage]
        values = [row["per_modality"][t]["rmse_um"] for t in groups[:-1]]
        values.append(row["balanced_target_rmse_um"])
        bars = ax.bar(x + offset, values, .36, color=color, label=label)
        ax.bar_label(bars, fmt="%.1f", padding=3, fontsize=8)
    ax.set_xticks(x, [labels[t] for t in groups[:-1]] + ["Equal-weight\ncombined"])
    ax.set_ylabel("Reference-corrected RMSE (µm)")
    ax.set_title(f"Reference: {labels[ref]} · section {slice_id}")
    ax.spines[["top", "right"]].set_visible(False)
    ax.margins(y=.22)
    ax.legend(frameon=False)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(outdir / f"{stem}_rmse.{ext}", dpi=300)
    plt.close(fig)

    # Per-modality rows prevent dense bin points from hiding cell/spot errors.
    panels = [corrected["unaligned"], corrected["stair_fine"], truth]
    all_xy = np.concatenate(panels)
    lo, hi = all_xy.min(0), all_xy.max(0)
    center = (lo + hi) / 2
    half = (hi - lo).max() * .53
    fig, axes = plt.subplots(3, 3, figsize=(9, 9))
    colors = np.array([DOMAIN_COLORS[str(d)] for d in domains])
    for r, tech in enumerate(TECHS):
        sel = np.flatnonzero(technology == tech)
        for c, xy in enumerate(panels):
            ax = axes[r, c]
            ax.scatter(xy[sel, 0], xy[sel, 1], c=colors[sel],
                       s=1 if tech == "bin16um" else 4, linewidths=0, rasterized=True)
            ax.set(xlim=(center[0]-half, center[0]+half),
                   ylim=(center[1]-half, center[1]+half), aspect="equal")
            ax.set_xticks([]); ax.set_yticks([])
            for spine in ax.spines.values():
                spine.set_visible(False)
            if r == 0:
                ax.set_title(["Unaligned", "STAIR aligned", "Ground truth"][c])
            if c == 0:
                ax.set_ylabel(labels[tech] + (" (reference)" if tech == ref else ""))
            if c < 2:
                stage = ["unaligned", "stair_fine"][c]
                error = result["stages"][stage]["per_modality"][tech]["rmse_um"]
                ax.text(.5, .01, f"RMSE {error:.1f} µm", ha="center", transform=ax.transAxes)
            if r == 2 and c == 2:
                left, bottom = center - half * .85
                ax.plot([left, left + 1000], [bottom, bottom], color="#4A4A4A", lw=2)
                ax.text(left + 500, bottom + half*.05, "1 mm", ha="center", fontsize=8)
    handles = [Line2D([], [], marker="o", ls="", color=color, label=key)
               for key, color in DOMAIN_COLORS.items()]
    fig.legend(handles=handles, loc="lower center", ncol=6, frameon=False)
    fig.suptitle(f"Section {slice_id} · shared rigid correction from {labels[ref]}")
    fig.tight_layout(rect=(0, .045, 1, .965))
    for ext in ("png", "pdf"):
        fig.savefig(outdir / f"{stem}_domains.{ext}", dpi=300)
    plt.close(fig)


def run(input_path, outdir, reference, slice_id=5):
    import anndata as ad
    input_path, outdir = Path(input_path).resolve(), Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    a = ad.read_h5ad(input_path, backed="r")
    try:
        truth = np.asarray(a.obsm["spatial_true"], float)
        coordinates = {stage: np.asarray(a.obsm[key], float) for stage, key in STAGES.items()}
        technology = a.obs["technology"].astype(str).to_numpy()
        domains = a.obs["domain_true"].astype(str).to_numpy()
    finally:
        a.file.close()
    result, corrected = evaluate(truth, coordinates, technology, reference)
    result.update({"input_h5ad": str(input_path), "input_size_bytes": input_path.stat().st_size,
                   "input_mtime_ns": input_path.stat().st_mtime_ns,
                   "script": str(Path(__file__).resolve()), "slice_id": slice_id,
                   "coordinate_keys": STAGES,
                   "coordinate_frame": "Saved STAIR spatial_true; includes preprocessing rescaling"})
    stem = f"reference_metrics_{reference}_slice_{slice_id}"
    (outdir / f"{stem}.json").write_text(json.dumps(result, indent=2) + "\n")
    with (outdir / f"{stem}.csv").open("w") as fh:
        writer = csv.writer(fh)
        writer.writerow(["stage", "modality", "is_reference", "n_obs", "rmse_um", "median_um", "p95_um"])
        for stage, row in result["stages"].items():
            for tech, stats in row["per_modality"].items():
                writer.writerow([stage, tech, tech == reference, *[stats[k] for k in
                                  ("n_obs", "rmse_um", "median_um", "p95_um")]])
    plot_results(result, corrected, truth, technology, domains, outdir, slice_id)
    print(json.dumps(result, indent=2))
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--reference", choices=TECHS, required=True)
    parser.add_argument("--slice", type=int, default=5)
    args = parser.parse_args()
    run(args.input, args.outdir, args.reference, args.slice)
