#!/usr/bin/env python3
"""Strong-domain-mix version of the independent-offset shift3x experiment.

ONE ALBIS call builds one tissue (cells, domains, cell types, molecules) from
the bin16um strong-mix parameters and aggregates it into bin16um, spot and
cell outputs, so all three modalities share the same underlying tissue
composition -- unlike ../generate_native_offsets.py, which ran one call per
modality with modality-specific n_cells/dispersion (different tissues).

Tissue/expression = Figure 3 strong-mix bin16um config
(code/clustering/strong_mix/generate_strong_mix_bin16um.sh, canonical point):
r=2050 um, 600k cells, log_mu=-2.5, theta=2.0, jitter=0.6, batch_sigma=0.7,
STRONG_DOMAIN_TYPE_MIX. Offsets = shift3x: max_shift=3075, max_deg=270,
base_seed_unaligned=12345, sync_unaligned_seed=False (each modality draws its
own rotation/translation per slice).

    python generate_strongmix_offsets.py --outdir OUT            # writes OUT/data/<tech>/
    python generate_strongmix_offsets.py --outdir OUT --print-config
"""
import argparse
import hashlib
import inspect
import json
import platform
import subprocess
import sys
from importlib.metadata import version
from pathlib import Path

import numpy as np
from scipy import sparse

ROOT = Path(__file__).resolve().parents[5]
REPO = ROOT / "albis"
sys.dont_write_bytecode = True
sys.path.insert(0, str(REPO))
import albis as ab
from albis.simulation_sphere import _to_serializable

# Copied from code/data/generate_simulation_noisy.py (that module parses argv
# at import time, so it can't be imported).
STRONG_DOMAIN_TYPE_MIX = [
    [0.30, 0.30, 0.0667, 0.0667, 0.0667, 0.0667, 0.0667, 0.0667],
    [0.0667, 0.0667, 0.30, 0.30, 0.0667, 0.0667, 0.0667, 0.0667],
    [0.0667, 0.0667, 0.0667, 0.0667, 0.30, 0.30, 0.0667, 0.0667],
    [0.0667, 0.0667, 0.0667, 0.0667, 0.0667, 0.0667, 0.30, 0.30],
    [0.30, 0.0667, 0.0667, 0.0667, 0.30, 0.0667, 0.0667, 0.0667],
    [0.0667, 0.0667, 0.30, 0.0667, 0.0667, 0.0667, 0.0667, 0.30],
]
# output folder -> (ALBIS modality, result key)
TECHS = {"bin16um": ("bin", "bin_adatas"), "spot": ("spot", "spot_adatas"),
         "cell": ("cell", "adata_cell_sectioned")}


def simulator_parameters(max_shift=3075., max_deg=270., perturbation_seed=12345,
                         tissue_seed=2025, slice_axis="Z"):
    parameters = {name: item.default for name, item in
                  inspect.signature(ab.simulate_3d_molecule_sphere_multires).parameters.items()
                  if item.default is not inspect.Parameter.empty}
    parameters.update(
        sphere_R_um=2050., capture_window_um=(6500., 6500.), xenium_capture_window_um=(12000., 24000.),
        n_domains=6, core_frac=.55, core_bump_amp=.25, wedge_angle_amp_deg=25., noise_terms=16,
        noise_freq_range=(3., 6.), boundary_fuzz_width_deg=6., boundary_fuzz_flip_prob=.15,
        core_fuzz_width_um=102.5, core_fuzz_flip_prob=.25,
        n_cells=600000, allow_cell_overlap=False,
        cell_radius_kwargs=dict(radius_dist="lognormal", r_mean=7.5, r_sigma=.28, r_min=4., r_max=14.),
        n_cell_types=8, domain_type_mix=STRONG_DOMAIN_TYPE_MIX,
        base_gene_lognormal=(-2.5, .7), theta=2., theta_jitter=.6, noise_scale=1.3,
        marker_foldchange=3.5, shared_marker_foldchange=2.5,
        n_slices=10, batch_sigma=.7, bin_size_um=16.,
        max_shift=max_shift, max_deg=max_deg, base_seed_unaligned=perturbation_seed,
        sync_unaligned_seed=False, seed=tissue_seed,
        output_modalities=tuple(m for m, _ in TECHS.values()), slice_axes=(slice_axis,))
    return parameters


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--max-shift", type=float, default=3075.)
    parser.add_argument("--max-deg", type=float, default=270.)
    parser.add_argument("--base-seed-unaligned", type=int, default=12345)
    parser.add_argument("--seed", type=int, default=2025)
    parser.add_argument("--n-cells", type=int, default=None, help="Override (smoke tests only)")
    parser.add_argument("--print-config", action="store_true")
    args = parser.parse_args()
    params = simulator_parameters(args.max_shift, args.max_deg, args.base_seed_unaligned, args.seed)
    if args.n_cells:
        params["n_cells"] = args.n_cells
    resolved = _to_serializable(params)
    if args.print_config:
        print(json.dumps(resolved, indent=2))
        return

    destinations = {t: args.outdir / "data" / t for t in TECHS}
    for d in destinations.values():
        d.mkdir(parents=True, exist_ok=False)
    source_files = [Path(__file__).resolve(), REPO / "albis/simulation_sphere.py"]
    revision = subprocess.run(["git", "-C", str(REPO), "rev-parse", "HEAD"], capture_output=True, text=True)
    provenance = {"status": "started", "shared_tissue": "single ALBIS call, all three modalities",
        "simulator_parameters": resolved,
        "albis_git_revision": revision.stdout.strip() if revision.returncode == 0 else None,
        "source_sha256": {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in source_files},
        "versions": {name: version(name) for name in ("numpy", "scipy", "anndata")},
        "python": platform.python_version(), "argv": sys.argv,
        "qc": "post-batch X total counts > 0 and detected genes >= 3"}
    for d in destinations.values():
        (d / "generation_manifest.json").write_text(json.dumps(provenance, indent=2) + "\n")
    print(json.dumps(provenance, indent=2), flush=True)

    sim = ab.simulate_3d_molecule_sphere_multires(**params)
    for tech, (modality, key) in TECHS.items():
        a = sim[key]["Z"]
        a.uns["native_generation_parameters"] = resolved
        stem = f"simulation_{modality}_z"
        a.write_h5ad(destinations[tech] / f"{stem}.h5ad")
        counts = np.asarray(a.X.sum(axis=1)).ravel()
        genes = np.asarray(a.X.getnnz(axis=1) if sparse.issparse(a.X) else np.count_nonzero(a.X, axis=1)).ravel()
        qc = a[(counts > 0) & (genes >= 3)].copy()
        qc.write_h5ad(destinations[tech] / f"{stem}_qc.h5ad")
        record = dict(provenance, status="completed", n_obs_raw=a.n_obs, n_obs_qc=qc.n_obs, n_genes=a.n_vars,
                      slices=sorted(qc.obs.slice_id.astype(str).unique().tolist(), key=int))
        (destinations[tech] / "generation_manifest.json").write_text(json.dumps(record, indent=2) + "\n")
        print(f"Completed {tech}: {a.n_obs} raw -> {qc.n_obs} QC observations, {a.n_vars} genes", flush=True)


if __name__ == "__main__":
    main()
