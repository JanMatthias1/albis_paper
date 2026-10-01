#!/usr/bin/env python3
"""Strong-domain-mix version of the independent-offset shift3x experiment.

ONE ALBIS call builds one tissue (cells, domains, cell types, molecules) from
the bin16um strong-mix parameters and aggregates it into bin16um, spot and
cell outputs, so all three modalities share the same underlying tissue
composition -- unlike the retired code/misc/applications_albis/cross_modality_alignment/weak_domain_mix/generate_native_offsets.py, which ran one call per
modality with modality-specific n_cells/dispersion (different tissues).

Tissue/expression = Figure 3 strong-mix bin16um config
(code/clustering/strong_domain_mix/generate/generate_strong_mix_bin16um.sh, canonical point):
r=2050 um, 600k cells, log_mu=-2.5, theta=2.0, jitter=0.6, batch_sigma=0.7,
STRONG_DOMAIN_TYPE_MIX. Offsets = shift3x: max_shift=3075, max_deg=270,
base_seed_unaligned=12345, sync_unaligned_seed=False (each modality draws its
own rotation/translation per slice).

    python generate_strongmix_offsets.py --outdir OUT            # writes OUT/data/<tech>/
    python generate_strongmix_offsets.py --outdir OUT --print-config
    python generate_strongmix_offsets.py --outdir OUT --crop-modalities bin,spot --crop-window-um 2221
    python generate_strongmix_offsets.py --outdir OUT --sphere-r-um 4600 --n-cells 1000000   # Figure 4A

--sphere-r-um sets the sphere size; the core fuzz width (0.05 R) and the default
max shift (1.5 R) scale with it, so the defaults (2050 um, 600k cells) give exactly
the configuration above. Figure 4A uses r=4600 um with 1M cells so the central
sections are cut square by the 6.5 mm bin/spot window (lower cell density than
Figure 2/3; illustrative). Figure 4C and Figure 5A use the defaults, uncropped.

--crop-modalities/--crop-window-um crop only the listed bin/spot outputs to a
W x W capture square (cell always keeps its window). The intact sphere is built
once and sectioned twice -- listed modalities with the small window, the rest
unchanged -- so the tissue stays shared and uncropped modalities match a run
without the flag. 2221 um = Visium's 6.5 mm scaled by the sphere shrink
(2050/6000), which crops the central sections to a full square.
"""
import argparse
import hashlib
import inspect
import json
import platform
import sys
from importlib.metadata import version
from pathlib import Path

import numpy as np
from scipy import sparse

ROOT = Path(__file__).resolve().parents[5]
sys.dont_write_bytecode = True
import albis as ab
import albis.simulation_sphere
from albis.simulation_sphere import _to_serializable

REQUIRED_ALBIS_VERSION = "0.1.2"
if ab.__version__ != REQUIRED_ALBIS_VERSION:
    raise RuntimeError(
        f"albis {ab.__version__} found at {ab.__file__}; this script requires albis "
        f"{REQUIRED_ALBIS_VERSION} (pip install albis=={REQUIRED_ALBIS_VERSION})."
    )

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
# section_3d_molecule_sphere arguments; the rest go to the base (intact sphere)
# call, mirroring how simulate_3d_molecule_sphere_multires splits them.
SECTION_ONLY = ("n_slices", "batch_sigma", "bin_size_um", "spot_spacing_um", "spot_radius_um", "max_deg",
                "max_shift", "base_seed_unaligned", "sync_unaligned_seed", "slice_axes")
CAPTURE = ("capture_window_um", "capture_window_center_um", "xenium_capture_window_um")


def simulate_cropped(params, crop_modalities, window_um):
    """Build the sphere once; section crop_modalities with a window_um square."""
    base = ab.simulate_3d_molecule_sphere_base(
        **{k: v for k, v in params.items() if k not in SECTION_ONLY + ("_section",)})
    section = {k: params[k] for k in SECTION_ONLY + CAPTURE}
    keep = [m for m, _ in TECHS.values() if m not in crop_modalities]
    sim = ab.section_3d_molecule_sphere(base, **section, output_modalities=keep)
    cropped = ab.section_3d_molecule_sphere(base, **dict(section, capture_window_um=(window_um, window_um)),
                                            output_modalities=crop_modalities)
    for modality, key in TECHS.values():
        if modality in crop_modalities:
            sim[key] = cropped[key]
    return sim


def simulator_parameters(max_shift=None, max_deg=270., perturbation_seed=12345,
                         tissue_seed=2025, slice_axis="Z", sphere_r_um=2050., n_cells=600000):
    """Size-dependent settings scale with the sphere radius as at r=2050:
    core fuzz width 0.05 R (102.5 um) and, unless given, max shift 1.5 R (3075 um)."""
    if max_shift is None:
        max_shift = 1.5 * sphere_r_um
    parameters = {name: item.default for name, item in
                  inspect.signature(ab.simulate_3d_molecule_sphere_multires).parameters.items()
                  if item.default is not inspect.Parameter.empty}
    parameters.update(
        sphere_R_um=sphere_r_um, capture_window_um=(6500., 6500.), xenium_capture_window_um=(12000., 24000.),
        n_domains=6, core_frac=.55, core_bump_amp=.25, wedge_angle_amp_deg=25., noise_terms=16,
        noise_freq_range=(3., 6.), boundary_fuzz_width_deg=6., boundary_fuzz_flip_prob=.15,
        core_fuzz_width_um=.05 * sphere_r_um, core_fuzz_flip_prob=.25,
        n_cells=n_cells, allow_cell_overlap=False,
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
    parser.add_argument("--sphere-r-um", type=float, default=2050.,
                        help="Sphere radius; core fuzz width (0.05 R) and default max shift (1.5 R) scale with it")
    parser.add_argument("--max-shift", type=float, default=None, help="Default 1.5 x --sphere-r-um (3075 at 2050)")
    parser.add_argument("--max-deg", type=float, default=270.)
    parser.add_argument("--base-seed-unaligned", type=int, default=12345)
    parser.add_argument("--seed", type=int, default=2025)
    parser.add_argument("--n-cells", type=int, default=600000)
    parser.add_argument("--crop-modalities", default="", help="Comma list from {bin,spot} to crop")
    parser.add_argument("--crop-window-um", type=float, default=None, help="Side of the square crop (um)")
    parser.add_argument("--print-config", action="store_true")
    args = parser.parse_args()
    crop = [m for m in args.crop_modalities.split(",") if m]
    if not set(crop) <= {"bin", "spot"} or bool(crop) != (args.crop_window_um is not None):
        parser.error("--crop-modalities (bin and/or spot) and --crop-window-um go together")
    params = simulator_parameters(args.max_shift, args.max_deg, args.base_seed_unaligned, args.seed,
                                  sphere_r_um=args.sphere_r_um, n_cells=args.n_cells)
    resolved = _to_serializable(params)
    if crop:
        resolved["crop"] = {"modalities": crop, "window_um": [args.crop_window_um] * 2}
    if args.print_config:
        print(json.dumps(resolved, indent=2))
        return

    destinations = {t: args.outdir / "data" / t for t in TECHS}
    for d in destinations.values():
        d.mkdir(parents=True, exist_ok=False)
    script = Path(__file__).resolve()
    simulator_source = Path(albis.simulation_sphere.__file__)
    provenance = {"status": "started", "shared_tissue": "single ALBIS call, all three modalities",
        "simulator_parameters": resolved,
        "albis_version": ab.__version__,
        "source_sha256": {str(script.relative_to(ROOT)): hashlib.sha256(script.read_bytes()).hexdigest(),
                          "albis/simulation_sphere.py": hashlib.sha256(simulator_source.read_bytes()).hexdigest()},
        "versions": {name: version(name) for name in ("numpy", "scipy", "anndata")},
        "python": platform.python_version(), "argv": sys.argv,
        "qc": "post-batch X total counts > 0 and detected genes >= 3"}
    for d in destinations.values():
        (d / "generation_manifest.json").write_text(json.dumps(provenance, indent=2) + "\n")
    print(json.dumps(provenance, indent=2), flush=True)

    if crop:
        sim = simulate_cropped(params, crop, args.crop_window_um)
    else:
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
