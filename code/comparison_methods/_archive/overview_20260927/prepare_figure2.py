"""Prepare a fresh overview using exact Figure 2 ALBIS files and synthetic reference."""
import argparse
import json
from pathlib import Path
import anndata as ad
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
SOURCES = {
    'cell': 'log_mu_-2.5_theta_0.40_jitter0.15_bsigma15',
    'bin': 'packing_pf0p04_bin16um_log_mu_-2.5_bsigma07',
    'spot': 'packing_pf0p04_log_mu_-2.0_theta_0.25_jitter0.10_bsigma03',
}
p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--out', type=Path, required=True)
p.add_argument('--settings', type=Path, default=HERE/'settings_figure2.json')
a = p.parse_args()
a.out.mkdir(parents=True, exist_ok=True)
if any(a.out.iterdir()):
    raise FileExistsError(f'Refusing nonempty output: {a.out}')
cfg = json.loads(a.settings.read_text())
(a.out/'settings.json').write_text(json.dumps(cfg, indent=2)+'\n')
(a.out/'albis').mkdir()
report = dict(method='albis', source='Exact Figure 2 raw outputs, linked without modification',
              separate_modality_simulations=True, display_offset_um=[2050,2050,2050], outputs={})
for modality, tag in SOURCES.items():
    source = ROOT/'sim_paper/data/figure_2/smaller_sphere/data'/tag/f'simulation_{modality}_z.h5ad'
    x = ad.read_h5ad(source, backed='r')
    if set(x.obs.slice_id.astype(int)) != set(range(cfg['n_slices'])):
        raise ValueError(f'Unexpected sections: {source}')
    if x.uns['sim_params']['sphere_r_um']*2 != cfg['extent_um']:
        raise ValueError(f'Unexpected physical extent: {source}')
    report['outputs'][modality] = dict(source=str(source), bytes=source.stat().st_size,
        n_observations=x.n_obs, n_genes=x.n_vars,
        simulated_cells=int(x.uns['sim_params']['n_cells']),
        slices=x.obs.slice_id.value_counts().sort_index().to_dict())
    x.file.close()
    (a.out/'albis'/f'{modality}.h5ad').symlink_to(source)
(a.out/'albis/manifest.json').write_text(json.dumps(report, indent=2)+'\n')
reference = ROOT/cfg['splatter_source']
if json.loads((reference/'manifest.json').read_text())['status'] != 'ok':
    raise ValueError('Synthetic reference failed')
(a.out/'splatter').symlink_to(reference)
(a.out/'input_provenance.json').write_text(json.dumps(dict(albis=report,
    splatter_source=str(reference), splatter_reference_cells=10000,
    sccube_requested_cells=cfg['n_cells'], selected_slice_id=5), indent=2)+'\n')
print(a.out, flush=True)
