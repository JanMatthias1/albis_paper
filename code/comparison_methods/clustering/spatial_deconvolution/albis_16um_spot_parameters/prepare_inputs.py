"""Validate native counts and retain nonempty cells/spots for unchanged RCTD."""
import argparse
import json
from pathlib import Path
import anndata as ad
import numpy as np
from scipy import sparse
from generate import OUT, SEEDS


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--seed', type=int, choices=SEEDS, required=True)
    seed = p.parse_args().seed
    run = OUT / f'seed{seed}'
    dest = run / 'inputs'
    dest.mkdir(exist_ok=False)
    records = {}
    genes = None
    for modality in ['cell', 'spot']:
        source = run / f'data/{modality}/simulation_{modality}_z.h5ad'
        a = ad.read_h5ad(source)
        assert a.n_vars == 556 and a.obs.slice_id.nunique() == 10
        if modality == 'cell':
            assert a.n_obs == 600000
            genes = a.var_names.copy()
        else:
            assert a.var_names.equals(genes)
        vals = a.X.data if sparse.issparse(a.X) else a.X.ravel()
        assert np.isfinite(vals).all() and (vals >= 0).all() and np.equal(vals, np.rint(vals)).all()
        n_raw = a.n_obs
        keep = np.asarray(a.X.sum(axis=1)).ravel() > 0
        if 'is_empty' in a.obs:
            keep &= ~a.obs.is_empty.to_numpy().astype(bool)
        keep &= a.obs.cell_type_true.astype(str).to_numpy() != 'unassigned'
        a = a[keep].copy()
        assert a.obs_names.is_unique
        if modality == 'spot':
            truth = np.asarray(a.obsm['cell_type_frac_true'])
            assert truth.shape == (a.n_obs, 8) and np.isfinite(truth).all() and (truth >= 0).all()
            assert np.allclose(truth.sum(axis=1), 1, atol=1e-5)
        else:
            assert set(a.obs.cell_type_true.astype(str)) == {f'type{i}' for i in range(1, 9)}
        a.X = sparse.csr_matrix(a.X, dtype=np.float64)
        a.obs.to_csv(dest / f'{modality}_obs.csv.gz')
        a.write_h5ad(dest / f'{modality}.h5ad', compression='gzip')
        records[modality] = {'source': str(source.resolve()), 'raw_observations': n_raw, 'retained_observations': a.n_obs}
        print(modality, records[modality], flush=True)
    (dest / 'input_audit.json').write_text(json.dumps({'simulation_seed': seed, 'sources': records, 'reference': 'Matching cell output from the same tissue and simulation seed; not independent.', 'truth': 'Native captured-molecule fractions, ordered type1..type8.', 'qc': 'Remove zero-total, explicitly empty or unassigned observations. RCTD internal QC remains default.', 'expression': 'Native integer counts; dtype conversion only, no normalization or rounding.', 'rctd_seed': 0}, indent=2) + '\n')


if __name__ == '__main__':
    main()
