"""Adapt completed native overview outputs to the unchanged Albis RCTD interface."""
import argparse
import hashlib
import json
from pathlib import Path
import anndata as ad
import numpy as np
from scipy.sparse import csr_matrix

TYPES = [f'type{i}' for i in range(1, 9)]

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--method', choices=['spider', 'sccube'], required=True)
    parser.add_argument('--n-genes', type=int, default=None)
    parser.add_argument('--gene-seed', type=int, default=20260926)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    report = {'method': args.method, 'source': str(args.source.resolve()),
              'truth': 'contributing-cell fractions, not molecule fractions',
              'reference': 'same simulated tissue as query; optimistic pilot',
              'gene_selection': 'all native genes; identical ordered reference/query panel',
              'rounding': 'nearest integer (numpy.rint) on separate copies; no depth normalization'}
    datasets = {}
    for modality in ['cell', 'spot']:
        source = args.source / args.method / f'{modality}.h5ad'
        a = ad.read_h5ad(source)
        assert a.obs_names.is_unique and a.var_names.is_unique
        if args.n_genes is not None:
            if not 0 < args.n_genes <= a.n_vars:
                raise ValueError('n-genes must be between 1 and the native gene count')
            # Same sorted candidate pool and seed give the same panel across methods.
            genes = np.sort(np.random.default_rng(args.gene_seed).choice(
                np.sort(a.var_names.to_numpy()), args.n_genes, replace=False))
            a = a[:, genes].copy()
            report['gene_selection'] = 'uniform random native-gene subset; identical reference/query panel'
            report['gene_seed'] = args.gene_seed
            report['selected_genes'] = genes.tolist()
        a.X = csr_matrix(a.X, dtype=np.float64)
        assert np.isfinite(a.X.data).all() and (a.X.data >= 0).all()
        before = np.asarray(a.X.sum(axis=1)).ravel()
        fractional = float(np.mean(a.X.data != np.rint(a.X.data)))
        if fractional and args.method != 'sccube':
            raise ValueError('Unexpected fractional SPIDER expression')
        delta = np.abs(a.X.data - np.rint(a.X.data)).sum()
        a.X.data = np.rint(a.X.data)
        a.X.eliminate_zeros()
        a.X = a.X.astype(np.int64)
        totals = np.asarray(a.X.sum(axis=1)).ravel()
        detected = np.asarray((a.X > 0).sum(axis=1)).ravel()
        report[modality] = dict(source=str(source.resolve()),
            source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
            n_obs=a.n_obs, n_genes=a.n_vars, fractional_nonzero_share=fractional,
            rounding_absolute_change=float(delta),
            total_before_quantiles=np.percentile(before, [0,25,50,75,100]).tolist(),
            total_after_quantiles=np.percentile(totals,[0,25,50,75,100]).tolist(),
            detected_genes_quantiles=np.percentile(detected,[0,25,50,75,100]).tolist())
        if modality == 'cell':
            assert set(a.obs.cell_type_true.astype(str)) == set(TYPES)
            assert all((a.obs.loc[totals >= 100, 'cell_type_true'].astype(str) == t).sum() >= 25 for t in TYPES)
        else:
            # capture_data() constructs cell_type_counts columns in TYPES order.
            counts = np.asarray(a.obsm['cell_type_counts'], dtype=float)
            assert counts.shape == (a.n_obs, len(TYPES))
            assert np.isfinite(counts).all() and (counts >= 0).all()
            occupied = counts.sum(axis=1) > 0
            assert np.array_equal(np.asarray(TYPES)[counts[occupied].argmax(axis=1)],
                                  a.obs.cell_type_true.astype(str).to_numpy()[occupied])
            keep = occupied & (totals > 0)
            report[modality]['excluded_empty_or_zero_expression'] = int((~keep).sum())
            a = a[keep].copy()
            counts = counts[keep]
            a.obsm['cell_type_frac_true'] = counts / counts.sum(axis=1, keepdims=True)
            assert np.allclose(a.obsm['cell_type_frac_true'].sum(axis=1), 1)
            assert np.isfinite(a.obsm['spatial']).all()
            a.uns['cell_type_frac_true_columns'] = TYPES
        datasets[modality] = a
    assert datasets['cell'].var_names.equals(datasets['spot'].var_names)
    for modality, a in datasets.items():
        a.write_h5ad(args.out / f'{modality}.h5ad', compression='gzip')
        # The unchanged R script uses sequential IDs: preserve their original identity.
        a.obs.assign(rctd_id=[f'{"ref" if modality == "cell" else "spot"}_{i+1}'
                             for i in range(a.n_obs)]).to_csv(args.out / f'{modality}_id_map.csv')
    (args.out / 'input_audit.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2), flush=True)

if __name__ == '__main__':
    main()
