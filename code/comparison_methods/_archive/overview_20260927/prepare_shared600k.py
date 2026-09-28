"""Link the shared Figure 4 ALBIS tissue and completed 600k scCube realization."""
import argparse
import json
from pathlib import Path

import anndata as ad
import numpy as np

ROOT = Path(__file__).resolve().parents[4]
SOURCE = ROOT / 'sim_paper/data/figure_4/cross_modality_alignment/strong_domain_mix_shift3x_cropped_bin/data'
COMPARATORS = ROOT / 'sim_paper/data/comparison_methods/overview_figure2_600k_35881874'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists():
        raise FileExistsError(args.out)
    cfg = json.loads((COMPARATORS / 'settings.json').read_text())
    cfg.update(albis_source=str(SOURCE),
               display_limits_um=[0, 4100],
               overview_note='600,000 cells per tissue; 10 sections. ALBIS bins: 16 µm, cropped field. '
                             'scCube targets: 3 cells/bin, 10 cells/spot.',
               unavailable_methods={'spider': 'No 600,000-cell output\nPinned native 3D implementation fails'})
    report = {'shared_tissue': True, 'coordinate_key': 'spatial_3d',
              'display_offset_um': [2050, 2050, 2050],
              'bin_crop_window_um': [2221, 2221], 'outputs': {}}
    links = {}
    for modality, folder in [('cell', 'cell'), ('bin', 'bin16um'), ('spot', 'spot')]:
        source = SOURCE / folder / f'simulation_{modality}_z.h5ad'
        manifest = json.loads((source.parent / 'generation_manifest.json').read_text())
        assert manifest['shared_tissue'] == 'single ALBIS call, all three modalities'
        assert manifest['simulator_parameters']['n_cells'] == cfg['n_cells'] == 600000
        assert manifest['simulator_parameters']['bin_size_um'] == 16
        with_data = ad.read_h5ad(source, backed='r')
        assert set(with_data.obs.slice_id.astype(int)) == set(range(10))
        assert np.isfinite(with_data.obsm['spatial_3d']).all()
        report['outputs'][modality] = dict(source=str(source), n_observations=with_data.n_obs,
            n_genes=with_data.n_vars, slices=with_data.obs.slice_id.value_counts().sort_index().to_dict())
        if modality == 'cell':
            assert with_data.n_obs == 600000
        with_data.file.close()
        links[modality] = source
    sc_manifest = json.loads((COMPARATORS / 'sccube/manifest.json').read_text())
    assert sc_manifest['outputs']['cell']['n_observations'] == 600000
    for key in ['n_cells', 'n_slices', 'extent_um', 'sccube_cells_per_bin', 'sccube_cells_per_spot']:
        assert sc_manifest['settings'][key] == cfg[key]
    for modality in links:
        source = COMPARATORS / 'sccube' / f'{modality}.h5ad'
        data = ad.read_h5ad(source, backed='r')
        assert data.n_obs == sc_manifest['outputs'][modality]['n_observations']
        data.file.close()
    args.out.mkdir(parents=True)
    (args.out / 'albis').mkdir()
    for modality, source in links.items():
        (args.out / 'albis' / f'{modality}.h5ad').symlink_to(source)
    (args.out / 'sccube').symlink_to(COMPARATORS / 'sccube', target_is_directory=True)
    (args.out / 'settings.json').write_text(json.dumps(cfg, indent=2) + '\n')
    (args.out / 'albis/manifest.json').write_text(json.dumps(report, indent=2) + '\n')
    (args.out / 'input_provenance.json').write_text(json.dumps(dict(albis=report,
        sccube_source=str(COMPARATORS / 'sccube'), selected_slice_id=4,
        caveats=['ALBIS: 556 genes; scCube: 2000 genes.',
                 'ALBIS realized type proportions differ from scCube equal requested proportions.',
                 'ALBIS bins retain the requested source crop; cell and spot fields are not cropped.',
                 'scCube capture units are occupancy-driven, not fixed 16 um bins or circular spots.',
                 'SPIDER has no successful native 600k 3D output in the pinned implementation.']), indent=2) + '\n')
    print(args.out)


if __name__ == '__main__':
    main()
