"""Reuse the Albis PCA/Harmony and count-matched Leiden pipeline on native cells."""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import anndata as ad
import numpy as np
from scipy import sparse

CODE = Path(__file__).resolve().parents[3] / 'clustering'

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--method', required=True, choices=['spider', 'sccube'])
    parser.add_argument('--source', required=True, type=Path)
    parser.add_argument('--out', required=True, type=Path)
    parser.add_argument("--expected-n-genes", type=int)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    source = args.source / args.method / 'cell.h5ad'
    a = ad.read_h5ad(source)
    assert a.obs_names.is_unique and a.var_names.is_unique
    if args.expected_n_genes is not None:
        assert a.n_vars == args.expected_n_genes, (a.n_vars, args.expected_n_genes)
    values = a.X.data if sparse.issparse(a.X) else np.asarray(a.X)
    assert np.isfinite(values).all() and (values >= 0).all()
    assert (np.asarray(a.X.sum(axis=1)).ravel() > 0).all()
    assert set(a.obs.cell_type_true.astype(str)) == {f'type{i}' for i in range(1,9)}
    assert a.obs.slice_id.notna().all() and a.obs.slice_id.nunique() > 1
    scripts = ['01.2_pca_harmony.py', '02_leiden_resolution_sweep.py', 'step03_cluster_and_plot.py']
    audit = dict(method=args.method, source=str(source.resolve()),
                 source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                 n_cells=a.n_obs, n_genes=a.n_vars, n_slices=int(a.obs.slice_id.nunique()),
                 cell_types=a.obs.cell_type_true.value_counts().to_dict(),
                 preprocessing='all genes; normalize_total 10000; log1p; scale max_value 10',
                 n_pcs=30, n_neighbors=15, batch_key='slice_id', random_state=0,
                 resolution_search=dict(lo=0.01, hi=3.0, max_iter=15, target_k=8,
                     selection='closest cluster count, not maximum ARI'),
                 script_sha256={s:hashlib.sha256((CODE/s).read_bytes()).hexdigest() for s in scripts})
    (args.out/'input_audit.json').write_text(json.dumps(audit,indent=2)+'\n')
    print(json.dumps(audit,indent=2),flush=True)
    del a
    pca = args.out/'pca_harmony'/'cell_pca_harmony.h5ad'
    subprocess.run([sys.executable,str(CODE/scripts[0]),'--modality','cell',
        '--input',str(source),'--output',str(pca),'--no-umap-sample',
        '--plot-colors','slice_id','cell_type_true'],check=True)
    # Competitors have no domain_true. Reuse the exact search implementation,
    # selecting only the real cell-type labels rather than fabricating domains.
    spec = importlib.util.spec_from_file_location('albis_resolution_sweep',CODE/scripts[1])
    sweep = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(sweep)
    sweep.GROUND_TRUTH_COLS = ('cell_type_true',)
    sys.argv = [str(CODE/scripts[1]),'--modality','cell','--input',str(pca),
                '--output-dir',str(args.out/'ari_recovery')]
    sweep.main()
    result = json.loads((args.out/'ari_recovery'/'ari_summary_cell.json').read_text())[0]
    subprocess.run([sys.executable,str(CODE/scripts[2]),'--modality','cell',
        '--input',str(pca),'--output-dir',str(args.out/'leiden_celltype_matched'),
        '--pipeline','pca_harmony','--resolution',str(result['resolution']),
        '--plot-colors','slice_id','cell_type_true','cluster_label'],check=True)
    print('[done]',json.dumps(result),flush=True)

if __name__ == '__main__':
    main()
