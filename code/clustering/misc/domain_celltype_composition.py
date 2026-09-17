#!/usr/bin/env python
"""
Why bin/spot cell_type_true ARI is capped: aggregation makes each observation
a MIXTURE of cell types, not one type, so no method can assign it a single
discrete label well.

For each ground-truth spatial domain (D0..D5), and each modality
(cell / bin 8um / bin 16um / spot), this looks at the observations assigned
to that domain (obs['domain_true'] == D) and reports:

  * composition -- the mean cell-type fraction vector within the domain
    (obsm['cell_type_frac_true']; for `cell` every observation is one type so
    this is just the category counts). Shows what the domain is made of.
  * purity      -- per observation, max_t cell_type_frac_true[t]: the share of
    the single dominant cell type. cell == 1.0 by construction; bin/spot < 1.
    Mean purity within a domain is roughly the ceiling on discrete
    cell-type-label accuracy there, hence on cell_type_true ARI.
  * entropy     -- per observation, Shannon entropy of the composition vector
    (0 = pure, log2(8)=3 = maximally mixed). Complementary mixing measure.

Outputs (data/figure_3/domain_celltype_composition/ by default):
  domain_composition_<D>.png    -- 1x2: stacked composition bar + purity violin, per modality
  mean_purity_by_domain.png     -- grouped bars, x=domain, one bar per modality (the headline)
  summary.csv                   -- modality x domain: n_obs, mean_purity, mean_entropy

Environment: conda activate /dcs04/hicks/data/Jan/sim_project/albis/env/albis-tutorial
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
SIM_PAPER_DIR = SCRIPT_DIR.parents[1]
FIG2 = SIM_PAPER_DIR / "data" / "figure_2"

# modality label -> (h5ad path, modality name for the sim file)
DATASETS = {
    "cell": FIG2 / "log_mu_-2.3_theta_0.40_jitter0.15_bsigma15" / "simulation_cell_z_qc.h5ad",
    "bin 8µm": FIG2 / "packing_pf0p04_log_mu_0.0_bsigma08" / "simulation_bin_z_qc.h5ad",
    "bin 16µm": FIG2 / "packing_pf0p04_bin16um_log_mu_-2.5_bsigma07" / "simulation_bin_z_qc.h5ad",
    "spot": FIG2 / "packing_pf0p04_log_mu_-2.5_bsigma03" / "simulation_spot_z_qc.h5ad",
}
CELL_TYPES = [f"type{i}" for i in range(1, 9)]
DOMAINS = [f"D{i}" for i in range(6)]


def frac_matrix(a: ad.AnnData) -> np.ndarray:
    """(n_obs, 8) cell-type fraction matrix; synthesised one-hot for single-cell data."""
    if "cell_type_frac_true" in a.obsm:
        return np.asarray(a.obsm["cell_type_frac_true"], dtype=float)
    codes = a.obs["cell_type_true"].cat.codes.to_numpy()
    m = np.zeros((a.n_obs, len(CELL_TYPES)))
    m[np.arange(a.n_obs), codes] = 1.0
    return m


def entropy(frac: np.ndarray) -> np.ndarray:
    p = np.clip(frac, 1e-12, None)
    p = p / p.sum(axis=1, keepdims=True)
    return -(p * np.log2(p)).sum(axis=1)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--domains", nargs="+", default=DOMAINS, choices=DOMAINS)
    ap.add_argument("--output-dir", type=Path, default=SIM_PAPER_DIR / "data" / "figure_3" / "domain_celltype_composition")
    ap.add_argument("--violin-sample", type=int, default=20_000, help="Max obs per modality/domain for the purity violin.")
    args = ap.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(0)

    # load once: per modality -> dict(domain_labels, frac, purity, entropy)
    data = {}
    for label, path in DATASETS.items():
        if not path.is_file():
            print(f"[skip] {label}: {path} not found")
            continue
        print(f"[load] {label}: {path.relative_to(SIM_PAPER_DIR)}")
        a = ad.read_h5ad(path)
        frac = frac_matrix(a)
        data[label] = {
            "domain": a.obs["domain_true"].astype(str).to_numpy(),
            "frac": frac,
            "purity": frac.max(axis=1),
            "entropy": entropy(frac),
        }
    labels = list(data)
    cmap = plt.get_cmap("tab10")

    rows = []
    purity_by_domain = {lab: [] for lab in labels}
    for D in args.domains:
        comp, purity_violin = {}, {}
        for lab in labels:
            sel = data[lab]["domain"] == D
            n = int(sel.sum())
            pur = data[lab]["purity"][sel]
            ent = data[lab]["entropy"][sel]
            comp[lab] = data[lab]["frac"][sel].mean(axis=0) if n else np.zeros(len(CELL_TYPES))
            v = pur
            if v.size > args.violin_sample:
                v = rng.choice(v, args.violin_sample, replace=False)
            purity_violin[lab] = v
            mp = float(pur.mean()) if n else float("nan")
            purity_by_domain[lab].append(mp)
            rows.append({"modality": lab, "domain": D, "n_obs": n,
                         "mean_purity": round(mp, 4),
                         "mean_entropy": round(float(ent.mean()) if n else float("nan"), 4)})

        fig, (axc, axp) = plt.subplots(1, 2, figsize=(12, 4.6))
        x = np.arange(len(labels))
        bottom = np.zeros(len(labels))
        for t, ct in enumerate(CELL_TYPES):
            vals = np.array([comp[lab][t] for lab in labels])
            axc.bar(x, vals, 0.6, bottom=bottom, label=ct, color=cmap(t % 10))
            bottom += vals
        axc.set_xticks(x); axc.set_xticklabels(labels)
        axc.set_ylabel("mean cell-type fraction"); axc.set_ylim(0, 1)
        axc.set_title(f"Cell-type composition within domain {D}")
        axc.legend(fontsize=7, ncol=2, loc="upper right")

        parts = axp.violinplot([purity_violin[lab] for lab in labels], showmeans=True, showextrema=False)
        for pc in parts["bodies"]:
            pc.set_alpha(0.6)
        axp.set_xticks(x + 1); axp.set_xticklabels(labels)
        axp.set_ylabel("per-observation purity  (max cell-type fraction)")
        axp.set_ylim(0, 1.02)
        axp.set_title(f"Observation purity within domain {D}\n(≈ ceiling on discrete cell-type accuracy)")
        fig.tight_layout()
        out = args.output_dir / f"domain_composition_{D}.png"
        fig.savefig(out, dpi=200, bbox_inches="tight")
        plt.close(fig)
        print(f"[save] {out}")

    # headline: mean purity, x=domain, grouped by modality
    fig, ax = plt.subplots(figsize=(8, 4.5))
    x = np.arange(len(args.domains))
    w = 0.8 / max(len(labels), 1)
    for i, lab in enumerate(labels):
        bars = ax.bar(x + (i - (len(labels) - 1) / 2) * w, purity_by_domain[lab], w, label=lab)
        ax.bar_label(bars, fmt="%.2f", fontsize=7, padding=2)
    ax.set_xticks(x); ax.set_xticklabels(args.domains)
    ax.set_ylabel("mean per-observation purity"); ax.set_ylim(0, 1.05)
    ax.set_title("Mean cell-type purity per ground-truth domain\n"
                 "cell ≈ 1.0 (one type/obs); bin/spot < 1 → discrete cell_type ARI is capped")
    ax.legend(loc="lower right", fontsize=8)
    fig.tight_layout()
    out = args.output_dir / "mean_purity_by_domain.png"
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"[save] {out}")

    csv_path = args.output_dir / "summary.csv"
    with csv_path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["modality", "domain", "n_obs", "mean_purity", "mean_entropy"])
        w.writeheader(); w.writerows(rows)
    print(f"[save] {csv_path}")


if __name__ == "__main__":
    main()
