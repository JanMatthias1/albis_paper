"""Figure 3 cell-type recovery at batch zero, reusing the domain experiment."""
import argparse
import hashlib
import json
import shutil
import sys
from pathlib import Path

CODE = Path(__file__).resolve().parent
sys.path.insert(0, str(CODE.parent / 'no_harmony_domain'))
from common import ROOT, SHARED, SEEDS, MODALITIES, K_GEOM, LAMBDAS, save, load_shared, task_for

BASE = ROOT / 'albis_paper/data/figure_3/no_harmony_cell_type'
DOMAIN = BASE.parent / 'no_harmony_domain'


def prepare(root):
    root.mkdir(exist_ok=False)
    for name in ['logs', 'summary', 'source_snapshot']:
        (root/name).mkdir()
    for phase in ['sweep', 'final']:
        rows = json.loads((DOMAIN/f'{phase}_tasks.json').read_text())
        for task in rows:
            task['directory'] = str(root / Path(task['directory']).relative_to(DOMAIN))
            task['target'] = 'cell_type_true'
        save(root/f'{phase}_tasks.json', rows)
    for name in ['input_audit.json', 'shared_source_sha256.json']:
        shutil.copy2(DOMAIN/name, root/name)
    protocol = json.loads((DOMAIN/'protocol.json').read_text())
    protocol.update(target='cell_type_true', selection='Maximum cell-type ARI among exact-eight-cluster runs on strong seed2025; exact ties choose smaller lambda.',
                    reference_comparison='Previous batch-zero Figure3 cell-type scores with Harmony; BANKSY lambda is independently retuned.',
                    aggregate_labels='Bins and spots use their stored dominant cell-type annotation; this is not mixture recovery.',
                    reuse='Reuse uncorrected PCA embeddings and K15 graphs from completed domain runs when input and parameters match. Recluster to eight types. Generate missing selected-lambda embeddings using the unchanged domain embedding function.')
    protocol['leiden']['target_clusters'] = 8
    protocol['selection_caveat'] = 'Seed2025 participates in selection and the displayed three-seed summary; not held out.'
    save(root/'protocol.json', protocol)
    for directory, prefix in [(CODE, ''), (CODE.parent/'no_harmony_domain', 'domain_')]:
        for f in directory.iterdir():
            if f.is_file():
                shutil.copy2(f, root/'source_snapshot'/f'{prefix}{f.name}')
    save(root/'dependency_sha256.json', {str(p): hashlib.sha256(p.read_bytes()).hexdigest()
         for p in [CODE.parent/'no_harmony_domain'/n for n in ['common.py', 'run.py']]})
    shutil.copy2(CODE/'README.md', root/'README.md')
    print('Prepared 21 cell-type lambda sweep tasks and 36 final comparisons.', flush=True)


def check(root):
    for path, digest in json.loads((root/'dependency_sha256.json').read_text()).items():
        assert hashlib.sha256(Path(path).read_bytes()).hexdigest() == digest, path
    import run as domain_run
    domain_run.check_sources(root)


def reuse_source(task):
    if task['phase'] == 'sweep':
        return DOMAIN/f"sweep/strong/{task['modality']}/seed2025/lam{task['lambda']:g}"
    source = DOMAIN/f"final/{task['mix']}/{task['modality']}/seed{task['seed']}/{task['pipeline']}"
    report = json.loads((source/'embedding_provenance.json').read_text())
    old = report['task']
    if all(old[k] == task[k] for k in ['input', 'pipeline', 'k_geom', 'lambda']):
        return source
    return None


def reuse_final(root, task):
    # With fixed parameters (no_harmony_domain/use_fixed_parameters.py) there is no sweep to reuse.
    if task['phase'] == 'final' and task['pipeline'] == 'banksy' and task['mix'] == 'strong' and task['seed'] == 2025:
        source = root/f"sweep/strong/{task['modality']}/seed2025/lam{task['lambda']:g}"
        if not (source/'metrics.json').exists():
            return False
        dest = Path(task['directory'])
        dest.parent.mkdir(parents=True, exist_ok=True)
        if not dest.exists():
            dest.symlink_to(source, target_is_directory=True)
        assert dest.resolve() == source.resolve()
        return True
    return False


def embed(root, task):
    if reuse_final(root, task):
        return
    dest = Path(task['directory'])
    if (dest/'embedding.h5ad').exists():
        assert json.loads((dest/'embedding_provenance.json').read_text())['task'] == task
        return
    source = reuse_source(task)
    if source is None:
        import run as domain_run
        domain_run.embed(root, task)
        return
    report = json.loads((source/'embedding_provenance.json').read_text())
    assert not report['harmony_applied']
    assert all(report['task'][k] == task[k] for k in ['input', 'pipeline', 'lambda', 'k_geom'])
    dest.mkdir(parents=True, exist_ok=True)
    (dest/'embedding.h5ad').symlink_to((source/'embedding.h5ad').resolve())
    save(dest/'embedding_provenance.json', dict(task=task, harmony_applied=False,
         reused_embedding=str((source/'embedding.h5ad').resolve()),
         graph_source=str((source/'embeddings_labels.h5ad').resolve()),
         original_provenance=str((source/'embedding_provenance.json').resolve())))
    print('Reused embedding and graph:', source, flush=True)


def cluster(root, task):
    if reuse_final(root, task):
        return
    import anndata as ad
    import numpy as np
    import scanpy as sc
    from scipy.sparse.csgraph import connected_components
    from sklearn.metrics import adjusted_rand_score
    dest = Path(task['directory'])
    if (dest/'metrics.json').exists():
        assert json.loads((dest/'metrics.json').read_text())['task'] == task
        return
    report = json.loads((dest/'embedding_provenance.json').read_text())
    a = ad.read_h5ad(report.get('graph_source', dest/'embedding.h5ad'))
    assert not a.uns['harmony_applied'] and 'X_pca_harmony' not in a.obsm
    assert a.obs.cell_type_true.nunique() == 8 and a.obs.slice_id.nunique() == 10
    save(dest/'status.json', dict(stage='clustering', status='running', task=task))
    if 'graph_source' not in report:
        sc.pp.neighbors(a, n_neighbors=15, use_rep='X_pca_pre_harmony', random_state=0, method='umap')
    else:
        params = a.uns['neighbors']['params']
        assert params['n_neighbors'] == 15 and params['use_rep'] == 'X_pca_pre_harmony'
    assert np.isfinite(a.obsp['connectivities'].data).all()
    nc, labels = connected_components(a.obsp['connectivities'], directed=False)
    diagnostic = dict(connected_components=int(nc), largest_components=np.sort(np.bincount(labels))[-10:][::-1].tolist(),
                      expression_graph_k=15, weighting='umap', harmony_applied=False, reused='graph_source' in report)
    sweep = load_shared('leiden_sweep', 'step02_leiden_resolution_sweep.py')
    res, k, trials = sweep.find_resolution_for_k(a, 8, .01, 3., 15, 0, 'leiden_cell_type_true')
    # The saved result contains only this experiment's inferred labels.
    if 'leiden_domain_true' in a.obs:
        del a.obs['leiden_domain_true']
    result = dict(task=task, n_obs=a.n_obs, n_slices=int(a.obs.slice_id.nunique()),
                  cell_type_ari=float(adjusted_rand_score(a.obs.cell_type_true.astype(str), a.obs.leiden_cell_type_true.astype(str))),
                  slice_ari=float(adjusted_rand_score(a.obs.slice_id.astype(str), a.obs.leiden_cell_type_true.astype(str))),
                  target_clusters=8, achieved_clusters=int(k), resolution=float(res), trials=trials,
                  harmony_applied=False, graph=diagnostic)
    a.obs.to_csv(dest/'labels.csv.gz')
    a.write_h5ad(dest/'embeddings_labels.h5ad')
    save(dest/'graph_diagnostic.json', diagnostic)
    save(dest/'metrics.json', result)
    save(dest/'status.json', dict(stage='clustering', status='complete' if k == 8 else 'cluster_count_mismatch'))
    print(json.dumps(result), flush=True)


def select(root):
    import pandas as pd
    rows = []
    for task in json.loads((root/'sweep_tasks.json').read_text()):
        path = Path(task['directory'])/'metrics.json'
        m = json.loads(path.read_text())
        assert m['task'] == task and not m['harmony_applied']
        rows.append(dict(modality=task['modality'], **{'lambda': task['lambda']},
                         cell_type_ari=m['cell_type_ari'], achieved_clusters=m['achieved_clusters'], path=str(path)))
    pd.DataFrame(rows).to_csv(root/'summary/sweep_metrics.csv', index=False)
    choices = {}
    for modality in MODALITIES:
        eligible = [r for r in rows if r['modality'] == modality and r['achieved_clusters'] == 8]
        if not eligible:
            raise RuntimeError(f'No exact-eight-cluster result for {modality}')
        choices[modality] = dict(sorted(eligible, key=lambda r: (-r['cell_type_ari'], r['lambda']))[0], k_geom=K_GEOM[modality])
    result = dict(selection_seed=2025, selection_mix='strong', target='cell_type_true', modalities=choices)
    path = root/'selected_parameters.json'
    if path.exists():
        assert json.loads(path.read_text()) == result, 'Frozen parameters differ'
    else:
        save(path, result)
    print(json.dumps(result, indent=2), flush=True)


def summarize(root):
    import pandas as pd
    rows = []
    for task in json.loads((root/'final_tasks.json').read_text()):
        resolved = task_for(root, 'final', task['task'])
        path = Path(task['directory'])/'metrics.json'
        m = json.loads(path.read_text())
        assert not m['harmony_applied']
        assert all(m['task'][k] == resolved[k] for k in ['input', 'pipeline', 'lambda', 'k_geom', 'target'])
        rows.append(dict(mix=task['mix'], modality=task['modality'], seed=task['seed'], pipeline=task['pipeline'],
                         **{'lambda': resolved['lambda']}, n_obs=m['n_obs'], cell_type_ari=m['cell_type_ari'],
                         achieved_clusters=m['achieved_clusters'], status='complete' if m['achieved_clusters'] == 8 else 'cluster_count_mismatch', path=str(path)))
    df = pd.DataFrame(rows)
    df.to_csv(root/'summary/metrics_by_seed.csv', index=False)
    df.groupby(['mix', 'modality', 'pipeline']).agg(mean_ari=('cell_type_ari','mean'), sample_sd=('cell_type_ari','std'), n_seeds=('seed','size')).to_csv(root/'summary/metrics_mean_sd.csv')
    sys.path.insert(0, str(SHARED/'ari_recovery_summary'))
    import plot_banksy_vs_pca_recovery as previous
    old, lambdas = previous.load()
    comparison = []
    for r in rows:
        # previous Harmony runs may be archived; then the comparison is left empty (as in no_harmony_domain/summarize.py)
        value = old.get((r['mix'], 'banksy' if r['pipeline'] == 'banksy' else 'genes_bs0', r['modality'], 'cell_type_true'), {}).get(r['seed'])
        comparison.append(dict(r, previous_harmony_ari=value, ari_difference=None if value is None else r['cell_type_ari']-value,
                               previous_lambda=lambdas[(r['modality'],'cell_type_true')] if r['pipeline'] == 'banksy' else None,
                               interpretation='BANKSY lambda also retuned; not an isolated Harmony effect.' if r['pipeline'] == 'banksy' else 'Same expression pipeline with Harmony removed.'))
    pd.DataFrame(comparison).to_csv(root/'summary/comparison_to_previous_harmony.csv', index=False)
    assert len(df) == 36 and (df.status == 'complete').all(), 'Inspect cluster-count mismatches before plotting'
    from plot import plot
    plot(root)
    print('Validated all 36 cell-type comparisons and generated summary figures.', flush=True)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('action', choices=['prepare','embed','cluster','select','summary'])
    p.add_argument('--phase', choices=['sweep','final'])
    p.add_argument('--task', type=int)
    p.add_argument('--root', type=Path, default=BASE)
    args = p.parse_args()
    if args.action == 'prepare':
        prepare(args.root)
        return
    check(args.root)
    if args.action in ['embed', 'cluster']:
        assert args.phase is not None and args.task is not None
        globals()[args.action](args.root, task_for(args.root, args.phase, args.task))
    elif args.action == 'select':
        select(args.root)
    else:
        summarize(args.root)


if __name__ == '__main__':
    main()
