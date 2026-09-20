#!/usr/bin/env python3
"""Reproduce the manuscript's strong-offset inputs using native ALBIS parameters.

Run one modality per process (bin16um and spot need paper-scale memory):
  python generate_native_offsets.py --modality bin16um --outdir OUTPUT

All ten sections are generated from scratch. No fitted/post-hoc transformation
or existing simulated data is used. Outputs include raw and count-QC AnnData,
the full resolved simulator arguments, source hashes, and package versions.
The tissue seed and perturbation seed are independent. Changing only max_shift
with a fixed perturbation seed preserves rotation/direction and scales translation.
"""
import argparse
import hashlib
import inspect
import json
from pathlib import Path
import platform
import subprocess
import sys
from importlib.metadata import version

import numpy as np
from scipy import sparse

ROOT = Path(__file__).resolve().parents[4]
REPO = ROOT / "albis"
sys.dont_write_bytecode = True
sys.path.insert(0, str(REPO))
import albis as ab
from albis.simulation_sphere import _to_serializable

DOMAIN_TYPE_MIX = [
    [.18, .18, .13, .12, .11, .10, .09, .09],
    [.11, .12, .18, .18, .13, .10, .09, .09],
    [.10, .11, .12, .13, .18, .18, .09, .09],
    [.10, .10, .11, .12, .13, .13, .16, .15],
    [.14, .13, .12, .11, .12, .13, .13, .12],
    [.125] * 8,
]
PRESETS = {
    'bin16um': dict(modality='bin', n_cells=600000, theta=2., theta_jitter=1.,
                   base_gene_lognormal=(-2.5, .7), batch_sigma=.7, bin_size_um=16.),
    'spot': dict(modality='spot', n_cells=600000, theta=.25, theta_jitter=.10,
                 base_gene_lognormal=(-2., .7), batch_sigma=.3, bin_size_um=8.),
    'cell': dict(modality='cell', n_cells=24207, theta=.40, theta_jitter=.15,
                 base_gene_lognormal=(-2.5, .7), batch_sigma=1.5, bin_size_um=8.),
}


def simulator_parameters(technology, max_shift=3075., max_deg=270., perturbation_seed=12345,
                         tissue_seed=2025, slice_axis='Z'):
    parameters = {name: item.default for name, item in
                  inspect.signature(ab.simulate_3d_molecule_sphere_multires).parameters.items()
                  if item.default is not inspect.Parameter.empty}
    preset = PRESETS[technology].copy()
    modality = preset.pop('modality')
    parameters.update(preset)
    parameters.update(sphere_R_um=2050., capture_window_um=(6500., 6500.),
        xenium_capture_window_um=(12000., 24000.), n_domains=6, core_frac=.55,
        core_bump_amp=.25, wedge_angle_amp_deg=25., noise_terms=16,
        noise_freq_range=(3., 6.), boundary_fuzz_width_deg=6., boundary_fuzz_flip_prob=.15,
        core_fuzz_width_um=102.5, core_fuzz_flip_prob=.25, n_cell_types=8,
        cell_radius_kwargs=dict(radius_dist='lognormal', r_mean=7.5, r_sigma=.28, r_min=4., r_max=14.),
        allow_cell_overlap=False, domain_type_mix=DOMAIN_TYPE_MIX, n_slices=10,
        noise_scale=1.3, marker_foldchange=3.5, shared_marker_foldchange=2.5,
        max_shift=max_shift, max_deg=max_deg, base_seed_unaligned=perturbation_seed,
        sync_unaligned_seed=False, seed=tissue_seed, output_modalities=(modality,), slice_axes=(slice_axis,))
    return parameters


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--modality', choices=PRESETS, required=True)
    parser.add_argument('--outdir', type=Path, required=True)
    parser.add_argument('--max-shift', type=float, default=3075.)
    parser.add_argument('--max-deg', type=float, default=270.)
    parser.add_argument('--base-seed-unaligned', type=int, default=12345)
    parser.add_argument('--seed', type=int, default=2025)
    parser.add_argument('--slice-axis', choices=['X', 'Y', 'Z'], default='Z')
    parser.add_argument('--print-config', action='store_true', help='Print resolved arguments without generating or writing')
    args = parser.parse_args()
    if any(not np.isfinite(v) or v < 0 for v in (args.max_shift, args.max_deg)) or min(args.seed, args.base_seed_unaligned) < 0:
        parser.error('Shift/rotation bounds must be finite and nonnegative; seeds must be nonnegative')
    params = simulator_parameters(args.modality, args.max_shift, args.max_deg,
                                  args.base_seed_unaligned, args.seed, args.slice_axis)
    resolved = _to_serializable(params)
    if args.print_config:
        print(json.dumps(resolved, indent=2))
        return
    destination = args.outdir / 'data' / args.modality
    destination.mkdir(parents=True, exist_ok=False)
    source_files = [Path(__file__).resolve(), REPO / 'albis/simulation_sphere.py']
    revision = subprocess.run(['git', '-C', str(REPO), 'rev-parse', 'HEAD'], capture_output=True, text=True)
    provenance = {'status': 'started', 'simulator_parameters': resolved,
        'albis_git_revision': revision.stdout.strip() if revision.returncode == 0 else None,
        'source_sha256': {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in source_files},
        'versions': {name: version(name) for name in ('numpy', 'scipy', 'anndata')},
        'python': platform.python_version(), 'argv': sys.argv,
        'qc': 'post-batch X total counts > 0 and detected genes >= 3'}
    manifest = destination / 'generation_manifest.json'
    manifest.write_text(json.dumps(provenance, indent=2) + '\n')
    print(json.dumps(provenance, indent=2), flush=True)
    sim = ab.simulate_3d_molecule_sphere_multires(**params)
    modality = PRESETS[args.modality]['modality']
    key = {'bin': 'bin_adatas', 'spot': 'spot_adatas', 'cell': 'adata_cell_sectioned'}[modality]
    a = sim[key][args.slice_axis]
    a.uns['native_generation_parameters'] = resolved
    stem = f'simulation_{modality}_{args.slice_axis.lower()}'
    a.write_h5ad(destination / f'{stem}.h5ad')
    counts = np.asarray(a.X.sum(axis=1)).ravel()
    genes = np.asarray(a.X.getnnz(axis=1) if sparse.issparse(a.X) else np.count_nonzero(a.X, axis=1)).ravel()
    keep = (counts > 0) & (genes >= 3)
    qc = a[keep].copy()
    qc.write_h5ad(destination / f'{stem}_qc.h5ad')
    provenance.update(status='completed', n_obs_raw=a.n_obs, n_obs_qc=qc.n_obs, n_genes=a.n_vars,
                      slices=sorted(qc.obs.slice_id.astype(str).unique().tolist(), key=int))
    manifest.write_text(json.dumps(provenance, indent=2) + '\n')
    print(f'Completed {args.modality}: {a.n_obs} raw → {qc.n_obs} QC observations, {a.n_vars} genes', flush=True)


if __name__ == '__main__':
    main()
