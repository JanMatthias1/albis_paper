"""Replace a previous experiment only after native generation and evaluation succeed."""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path


def promote(staging, target):
    if staging.is_symlink() or target.is_symlink():
        raise ValueError('Staging and target must be real directories')
    staging, target = staging.resolve(), target.resolve()
    if staging == target or staging.parent != target.parent:
        raise ValueError('Staging and target must be distinct sibling experiment directories')
    if not target.is_dir():
        raise FileNotFoundError(target)
    for tech, modality in [('bin16um', 'bin'), ('spot', 'spot'), ('cell', 'cell')]:
        manifest = json.loads((staging / 'data' / tech / 'generation_manifest.json').read_text())
        if manifest['status'] != 'completed' or manifest['slices'] != [str(i) for i in range(10)]:
            raise ValueError(f'{tech}: incomplete native generation')
        for suffix in ('', '_qc'):
            path = staging / 'data' / tech / f'simulation_{modality}_z{suffix}.h5ad'
            if not path.is_file() or path.stat().st_size == 0:
                raise FileNotFoundError(path)
    required = [
        'STAIR/cross_tech/slice_5/adata_results/Sim_CrossTech_STAIR_slice_5.h5ad',
        'plots/slice5_sphere_3d_domain_true.png',
        'plots/reference_metrics_bin16um_slice_5.json',
        'gene_similarity/gene_similarity_ridges_average_all.png',
        'gene_similarity/gene_similarity_ridges_average_all.pdf',
        'gene_similarity/gene_similarity_ridges_average_nz.png',
        'gene_similarity/gene_similarity_average_per_gene.csv',
        'gene_similarity/provenance.json',
    ]
    for relative in required:
        path = staging / relative
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError(path)
    stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    archive = target.with_name(target.name + '_before_native_' + stamp)
    record = {'replaced_at_utc': stamp, 'target': str(target), 'archived_previous': str(archive),
              'staging_alias': str(staging), 'note': 'Staging path remains an alias for recorded provenance paths.'}
    (staging / 'replacement_manifest.json').write_text(json.dumps(record, indent=2) + '\n')
    target.rename(archive)
    try:
        staging.rename(target)
    except Exception:
        archive.rename(target)
        raise
    # Saved input/log paths remain resolvable after promotion.
    staging.symlink_to(target.name, target_is_directory=True)
    print(json.dumps(record, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--staging', type=Path, required=True)
    parser.add_argument('--target', type=Path, required=True)
    args = parser.parse_args()
    promote(args.staging, args.target)
