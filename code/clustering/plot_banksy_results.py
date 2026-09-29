#!/usr/bin/env python
"""Render BANKSY PCA/UMAP and matched-cluster panels from saved results.

Uses the same plotting functions as the plain PCA workflow. Reads only obs and
embeddings, preserves scientific h5ads, and caches display UMAPs alongside plots.
No dependence on misc/ or the one-off Figure 3 style-refresh job directory.
"""
from __future__ import annotations
import argparse
import importlib.util
from pathlib import Path
import numpy as np
import h5py
import anndata as ad
from anndata._io.specs import read_elem
import step03_cluster_and_plot as cluster
from pc_pairs import sampled_indices

SPEC = importlib.util.spec_from_file_location('pca_plotting', Path(__file__).with_name('step01_pca_harmony.py'))
pca = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(pca)


def read_plot_data(path):
    with h5py.File(path) as handle:
        obs = read_elem(handle['obs'])
        embeddings = {key: read_elem(handle['obsm'][key]) for key in handle['obsm']
                      if key.startswith(('X_pca', 'X_umap'))}
    return ad.AnnData(obs=obs, obsm=embeddings)


def render(input_path, cluster_path=None, output_dir=None, max_obs=50000, seed=0, neighbors=15):
    a = read_plot_data(input_path)
    a.obsm['X_pca_post_harmony'] = a.obsm['X_pca_harmony'].copy()
    out = output_dir / 'banksy_matrix/plots' if output_dir else input_path.parent / 'plots'
    out.mkdir(parents=True, exist_ok=True)
    for key in ('slice_id', 'domain_true', 'cell_type_true'):
        if key in a.obs:
            for suffix in ('pre', 'post'):
                pca.plot_two_dims(a, f'X_pca_{suffix}_harmony', key, out / f'pca_{suffix}harmony_by_{key}.png')
    indices = sampled_indices(a.n_obs, max_obs, False, seed)
    u = a[indices].copy()
    cache_dir = out / 'umap_cache'
    cache_dir.mkdir(exist_ok=True)
    for suffix in ('pre', 'post'):
        rep = f'X_pca_{suffix}_harmony'
        key = f'X_umap_pca_{suffix}_harmony'
        if key in u.obsm:
            continue
        cache = cache_dir / f'{suffix}_seed{seed}_n{neighbors}.npz'
        reused = False
        if cache.exists():
            with np.load(cache) as saved:
                if np.array_equal(saved['obs'], np.asarray(u.obs_names, str)) and np.array_equal(saved['pca'], u.obsm[rep]):
                    u.obsm[key] = saved['coords'].copy()
                    reused = True
        if not reused:
            pca.compute_umap(u, rep, key, neighbors, seed)
            np.savez_compressed(cache, obs=np.asarray(u.obs_names, str), pca=u.obsm[rep], coords=u.obsm[key])
    pca.save_umap_plots(u, out, ['slice_id', 'domain_true', 'cell_type_true'], 4, .85)
    if cluster_path is not None:
        labels = read_plot_data(cluster_path)[u.obs_names].copy()
        if not np.allclose(labels.obsm['X_pca_harmony'], u.obsm['X_pca_harmony']):
            raise ValueError('Clustering and BANKSY input embeddings differ')
        labels.obsm['X_umap'] = u.obsm['X_umap_pca_post_harmony'].copy()
        out = output_dir / 'ari/plots' if output_dir else cluster_path.parent / 'plots'
        out.mkdir(parents=True, exist_ok=True)
        for truth in ('domain_true', 'cell_type_true'):
            pred = f'leiden_{truth}'
            if truth not in labels.obs or pred not in labels.obs:
                raise ValueError(f'Missing {truth} or {pred} in {cluster_path}')
            cluster.plot_umap_true_vs_predicted(labels, truth, pred, out / f'umap_true_vs_predicted_{truth}.png')
            cluster.plot_contingency_heatmap(labels, truth, pred, out / f'contingency_{truth}.png')
    print(f'[save] BANKSY plots: {out}')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--cluster-input', type=Path, help='ARI h5ad with leiden_domain_true / leiden_cell_type_true')
    parser.add_argument('--output-dir', type=Path, help='Alternate root for banksy_matrix/plots and ari/plots')
    parser.add_argument('--max-obs', type=int, default=50000)
    parser.add_argument('--random-state', type=int, default=0)
    parser.add_argument('--n-neighbors', type=int, default=15)
    args = parser.parse_args()
    if args.max_obs <= 0 or args.n_neighbors < 2:
        parser.error('--max-obs must be positive and --n-neighbors at least 2')
    render(args.input, args.cluster_input, args.output_dir, args.max_obs, args.random_state, args.n_neighbors)


if __name__ == '__main__':
    main()
