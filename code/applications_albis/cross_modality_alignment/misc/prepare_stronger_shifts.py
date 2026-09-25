"""Create a section-5 experiment by increasing existing rigid translations.

Keep the original expression, truth, rotations and QC membership. This changes
the actual STAIR input, not just the display. Never overwrite an experiment.
"""
import argparse
import json
from pathlib import Path

import anndata as ad
import numpy as np
from reference_metrics import rigid_fit, TECHS


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--outdir', type=Path, required=True)
    parser.add_argument('--factor', type=float, default=3.)
    args = parser.parse_args()
    if not np.isfinite(args.factor) or args.factor <= 1:
        raise ValueError('Translation factor must exceed 1')
    if args.outdir.exists():
        raise FileExistsError(f'Refusing to overwrite {args.outdir}')
    args.outdir.mkdir(parents=True)
    provenance = {'translation_factor': args.factor, 'slice_id': 5,
                  'source': str(args.source.resolve()),
                  'definition': 'new_xy = old_unaligned_xy + (factor - 1) * fitted_original_translation',
                  'modalities': {}}
    for tech in TECHS:
        modality = 'bin' if tech == 'bin16um' else tech
        name = f'simulation_{modality}_z_qc.h5ad'
        source = args.source / 'data' / tech / name
        backed = ad.read_h5ad(source, backed='r')
        try:
            a = backed[backed.obs['slice_id'].astype(str).to_numpy() == '5'].to_memory()
        finally:
            backed.file.close()
        truth = np.asarray(a.obsm['spatial'], float)[:, :2]
        original = np.asarray(a.obsm['spatial_unaligned'], float).copy()
        rotation, translation = rigid_fit(truth, original[:, :2])
        updated = original.copy()
        updated[:, :2] += (args.factor - 1) * translation
        fitted = truth @ rotation + args.factor * translation
        residual = float(np.sqrt(np.mean(np.sum((updated[:, :2] - fitted)**2, axis=1))))
        if residual > .01:
            raise ValueError(f'{tech}: unexpected non-rigid residual {residual}')
        a.obsm['spatial_unaligned'] = updated
        # Retain original provenance explicitly; it no longer describes new offsets.
        if 'rigid_perturb_inplane' in a.uns:
            a.uns['source_rigid_perturb_inplane'] = a.uns.pop('rigid_perturb_inplane')
        record = {'source_h5ad': str(source.resolve()), 'n_obs': a.n_obs,
                  'rotation': rotation.tolist(), 'original_translation_um': translation.tolist(),
                  'new_translation_um': (args.factor * translation).tolist(),
                  'rigid_residual_um': residual}
        a.uns['stronger_shifts'] = {'translation_factor': args.factor,
                                   'source_h5ad': str(source.resolve()),
                                   'rotation': rotation, 'translation_um': args.factor * translation}
        destination = args.outdir / 'data' / tech / name
        destination.parent.mkdir(parents=True)
        a.write_h5ad(destination)
        provenance['modalities'][tech] = record
        print(tech, json.dumps(record), flush=True)
    (args.outdir / 'perturbation_provenance.json').write_text(json.dumps(provenance, indent=2) + '\n')


if __name__ == '__main__':
    main()
