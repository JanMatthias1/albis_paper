#!/usr/bin/env python
"""
Figure 4A follow-up -- domain-celltype coupling strength sweep.

Both scBSP and SPARK-X land at near-chance AUROC on the canonical (weak-mix)
`spot` dataset for SVG recovery, and SPARK-X performs no better than scBSP on
the `strongmix` variant either (see SVG_PLAN.md / 3D_scbsp.py history). The
suspected root cause: `domain_true` only weakly determines `cell_type_true`
in the canonical `domain_type_mix` (see generate_simulation_noisy.py's
MANUSCRIPT_DOMAIN_TYPE_MIX docstring), so marker-gene expression -- which IS a
strong, non-spatial cell-type signal -- only weakly tracks *spatial* domain
structure. This script tests that hypothesis directly: sweep the
domain-celltype coupling strength and measure whether SVG-recovery AUROC
climbs toward 1.0 as coupling strengthens (confirms a coupling-strength
effect, not a broken detector) or stays flat (points to something else).

Uses the exact `spot` sim_params already locked for the canonical Figure 2 tag
`packing_pf0p04_log_mu_-2.5_bsigma03` (read from that h5ad's `uns['sim_params']`
and hardcoded below), varying ONLY `domain_type_mix`. Sweep matrices share
STRONG_DOMAIN_TYPE_MIX's topology (each of 6 domains has 2 "signature" cell
types) but vary the signature probability `p_sig`: p_sig=0.125 is exactly
uniform (8 types x 0.125 = no coupling at all, a true chance floor);
p_sig=0.30 matches the already-tested STRONG_DOMAIN_TYPE_MIX; higher values
push further than anything tested so far, up to near-deterministic.

For each p_sig, generates spot data, QC-filters it (mirrors 00_qc_filter.py's
total_counts>0 & n_genes>=3 filter on post-batch counts, inline rather than
shelling out), and writes SPARK-X-ready flat files (counts.mtx genes x cells,
locations.csv, gene/cell name lists, gene_meta.csv) plus a JSON of the
achieved domain/celltype crosstab enrichment (max deviation from the 12.5%
uniform rate) as a model-free x-axis check that p_sig actually drove
coupling strength as intended.

Expected environment:
    conda activate /dcs04/hicks/data/Jan/sim_project/sim_paper/env/albis-tutorial

Example:
    python domain_mix_sweep_generate.py --p-sig 0.125 --out-tag psig0p125
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.io as sio
from scipy import sparse


def find_repo_root(start: Path) -> Path:
    for candidate in [start, *start.parents]:
        nested = candidate / "albis"
        if (nested / "pyproject.toml").is_file() and (nested / "albis" / "__init__.py").is_file():
            return nested
    raise RuntimeError("Could not find the albis repository root.")


REPO_ROOT = find_repo_root(Path.cwd().resolve())
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
sys.modules.pop("albis", None)
import albis as ab  # noqa: E402

SIM_PAPER_DIR = Path("/dcs04/hicks/data/Jan/sim_project/sim_paper")
OUT_ROOT = SIM_PAPER_DIR / "data" / "figure_4" / "svg" / "domain_mix_sweep"

# Exact sim_params locked for the canonical spot Fig2 tag
# (packing_pf0p04_log_mu_-2.5_bsigma03), read from that h5ad's uns['sim_params'].
FIXED = dict(
    sphere_R_um=2050.0,
    capture_window_um=(6500.0, 6500.0),
    n_domains=6,
    core_frac=0.55,
    core_bump_amp=0.25,
    wedge_angle_amp_deg=25.0,
    noise_terms=16,
    noise_freq_range=(3.0, 6.0),
    boundary_fuzz_width_deg=6.0,
    boundary_fuzz_flip_prob=0.15,
    core_fuzz_width_um=102.5,
    core_fuzz_flip_prob=0.25,
    n_cells=600_000,
    cell_radius_kwargs=dict(radius_dist="lognormal", r_mean=7.5, r_sigma=0.28, r_min=4.0, r_max=14.0),
    allow_cell_overlap=False,
    n_cell_types=8,
    n_slices=10,
    bin_size_um=8.0,
    theta=2.0,
    theta_jitter=1.0,
    noise_scale=1.3,
    marker_foldchange=3.5,
    shared_marker_foldchange=2.5,
    base_gene_lognormal=(-2.5, 0.7),
    batch_sigma=0.3,
    domain_size_factors=None,
    max_deg=270.0,
    max_shift=1025.0,
    output_modalities=("spot",),
    slice_axes=("Z",),
)

# STRONG_DOMAIN_TYPE_MIX's domain -> 2 signature cell-type indices (0-based),
# from generate_simulation_noisy.py.
SIGNATURE_TYPES = {
    0: (0, 1),
    1: (2, 3),
    2: (4, 5),
    3: (6, 7),
    4: (0, 4),
    5: (2, 7),
}


def build_domain_type_mix(p_sig: float, n_domains: int = 6, n_cell_types: int = 8) -> np.ndarray:
    """p_sig=0.125 (=1/8) is exactly uniform (no coupling); higher p_sig
    concentrates more probability mass on each domain's 2 signature types."""
    if not (1.0 / n_cell_types <= p_sig < 0.5):
        raise ValueError(f"p_sig must be in [{1/n_cell_types:.4f}, 0.5), got {p_sig}")
    other_p = (1.0 - 2 * p_sig) / (n_cell_types - 2)
    mix = np.full((n_domains, n_cell_types), other_p, dtype=float)
    for d, (t1, t2) in SIGNATURE_TYPES.items():
        mix[d, t1] = p_sig
        mix[d, t2] = p_sig
    assert np.allclose(mix.sum(axis=1), 1.0)
    return mix


def crosstab_enrichment(adata) -> dict:
    """Max |mean per-domain cell_type_frac_true - uniform 1/n_cell_types|, a
    model-free check that p_sig actually moved realized coupling strength as
    intended. Uses the continuous per-spot composition (obsm), not the
    per-spot majority-vote `cell_type_true` label -- spot's per-observation
    cell count is small enough that a majority-vote crosstab is dominated by
    argmax quantization noise rather than the underlying coupling strength."""
    frac = np.asarray(adata.obsm["cell_type_frac_true"])
    domain = adata.obs["domain_true"].astype(str).to_numpy()
    uniform = 1.0 / frac.shape[1]
    max_dev = 0.0
    for d in np.unique(domain):
        mean_frac = frac[domain == d].mean(axis=0)
        max_dev = max(max_dev, float(np.abs(mean_frac - uniform).max()))
    return {"max_abs_dev_from_uniform": max_dev, "uniform_rate": uniform}


def qc_filter(adata):
    X = adata.layers["counts"] if "counts" in adata.layers else adata.X
    n_genes = np.asarray(X.getnnz(axis=1) if sparse.issparse(X) else np.count_nonzero(X, axis=1)).ravel()
    total_counts = np.asarray(X.sum(axis=1)).ravel()
    keep = (total_counts > 0) & (n_genes >= 3)
    print(f"[qc] keeping {keep.sum()}/{len(keep)} ({keep.mean()*100:.1f}%)")
    return adata[keep].copy()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--p-sig", type=float, required=True)
    parser.add_argument("--out-tag", required=True)
    parser.add_argument("--boundary-fuzz-flip-prob", type=float, default=None,
                         help="Override FIXED's boundary_fuzz_flip_prob (default 0.15) -- e.g. 0.0 to "
                              "test whether clean (non-fuzzy) domain edges close the AUROC gap left "
                              "after the p_sig sweep.")
    parser.add_argument("--core-fuzz-flip-prob", type=float, default=None,
                         help="Override FIXED's core_fuzz_flip_prob (default 0.25), same rationale.")
    args = parser.parse_args()

    domain_type_mix = build_domain_type_mix(args.p_sig)
    print(f"[sweep] p_sig={args.p_sig}\n{domain_type_mix}")

    fixed = dict(FIXED)
    if args.boundary_fuzz_flip_prob is not None:
        fixed["boundary_fuzz_flip_prob"] = args.boundary_fuzz_flip_prob
    if args.core_fuzz_flip_prob is not None:
        fixed["core_fuzz_flip_prob"] = args.core_fuzz_flip_prob
    print(f"[sweep] boundary_fuzz_flip_prob={fixed['boundary_fuzz_flip_prob']} "
          f"core_fuzz_flip_prob={fixed['core_fuzz_flip_prob']}")

    simulation = ab.simulate_3d_molecule_sphere_multires(domain_type_mix=domain_type_mix, **fixed)
    adata = simulation["spot_adatas"]["Z"]
    print(f"[gen] raw shape: {adata.n_obs} x {adata.n_vars}")

    enrich = crosstab_enrichment(adata)
    print(f"[gen] crosstab max_abs_dev_from_uniform={enrich['max_abs_dev_from_uniform']:.4f}")

    adata = qc_filter(adata)

    outdir = OUT_ROOT / args.out_tag
    outdir.mkdir(parents=True, exist_ok=True)

    sio.mmwrite(str(outdir / "counts.mtx"), adata.X.T.tocoo())
    (outdir / "cell_names.txt").write_text("\n".join(adata.obs_names))
    (outdir / "gene_names.txt").write_text("\n".join(adata.var_names))
    adata.var[["is_noise"]].to_csv(outdir / "gene_meta.csv")
    coords = adata.obsm["spatial_3d"]
    pd.DataFrame(coords, columns=["x", "y", "z"], index=adata.obs_names).to_csv(outdir / "locations.csv")

    meta = {"p_sig": args.p_sig, "out_tag": args.out_tag, "n_obs_qc": int(adata.n_obs),
            "n_genes": int(adata.n_vars),
            "boundary_fuzz_flip_prob": fixed["boundary_fuzz_flip_prob"],
            "core_fuzz_flip_prob": fixed["core_fuzz_flip_prob"], **enrich}
    (outdir / "sweep_meta.json").write_text(json.dumps(meta, indent=2))
    print(f"[save] {outdir}")


if __name__ == "__main__":
    main()
