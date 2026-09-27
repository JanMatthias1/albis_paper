"""Cell-type recovery on the same 600k ALBIS simulations; no Harmony."""
import argparse
import gc
import json
import sys
import hashlib
from pathlib import Path
import numpy as np
import anndata as ad
import scanpy as sc
from scipy.sparse.csgraph import connected_components
from sklearn.metrics import adjusted_rand_score
from generate import OUT, ROOT, SEEDS
sys.path.insert(0, str(ROOT / 'sim_paper/code/comparison_methods/clustering/cell_type'))
import cluster


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--seed', type=int, choices=SEEDS, required=True)
    seed = p.parse_args().seed
    run = OUT / f'seed{seed}'
    dest = run / 'clustering'
    dest.mkdir(exist_ok=False)
    source = run / 'data/cell/simulation_cell_z.h5ad'
    x = ad.read_h5ad(source)
    assert x.shape == (600000, 556) and x.obs.slice_id.nunique() == 10
    keep = np.asarray(x.X.sum(axis=1)).ravel() > 0
    x = x[keep].copy() if not keep.all() else x
    np.random.seed(0)
    cluster.load('prep', '01.2_pca_harmony.py').preprocess_for_pca(x, 30, 0, 10000)
    a = ad.AnnData(obs=x.obs.copy())
    a.obsm['X_pca_pre_harmony'] = x.obsm['X_pca_pre_harmony'].copy()
    for key in ['spatial', 'spatial_3d']:
        if key in x.obsm:
            a.obsm[key] = x.obsm[key].copy()
    del x
    gc.collect()
    print('Building Gaussian K256 graph', flush=True)
    sc.pp.neighbors(a, n_neighbors=256, use_rep='X_pca_pre_harmony', method='gauss', random_state=0)
    assert np.isfinite(a.obsp['connectivities'].data).all()
    nc, labels = connected_components(a.obsp['connectivities'], directed=False)
    report = dict(method='albis', seed=seed, genes=556, n_cells_generated=600000, n_cells_clustered=a.n_obs, n_slices=10, batch_sigma=.7, harmony_applied=False, n_neighbors=256, neighbor_weighting='gauss', n_pcs=30, clustering_seed=0, connected_components=int(nc), largest_components=np.sort(np.bincount(labels))[-10:][::-1].tolist(), input=str(source.resolve()), shared_code_sha256={n: hashlib.sha256((cluster.SHARED/n).read_bytes()).hexdigest() for n in ['01.2_pca_harmony.py', '02_leiden_resolution_sweep.py']})
    (dest / 'graph_diagnostic.json').write_text(json.dumps(report, indent=2) + '\n')
    resolution, k, trials = cluster.load('sweep', '02_leiden_resolution_sweep.py').find_resolution_for_k(a, 8, .01, 3, 15, 0, 'leiden_cell_type')
    report.update(cell_type_ari=float(adjusted_rand_score(a.obs.cell_type_true.astype(str), a.obs.leiden_cell_type.astype(str))), achieved_clusters=int(k), resolution=float(resolution), trials=trials)
    a.obs.to_csv(dest / 'labels.csv.gz')
    a.write_h5ad(dest / 'embeddings_labels.h5ad')
    (dest / 'metrics.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report), flush=True)


if __name__ == '__main__':
    main()
