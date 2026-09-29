"""Validate all 18 existing batch-zero datasets and prepare the sweep/final manifests."""
import csv
import hashlib
import json
import shutil
import anndata as ad
from common import *

def main():
    BASE.mkdir(exist_ok=False)
    (BASE / 'logs').mkdir()
    (BASE / 'summary').mkdir()
    inputs = {}
    audit = []
    for mix in ['strong', 'weak']:
        table = SHARED / f'{mix}_domain_mix/banksy_harmony_batch_zero/domain/final_tasks.tsv'
        for row in csv.DictReader(table.open(), delimiter='\t'):
            modality, seed = row['modality'], int(row['seed'])
            path = (ROOT / row['qc_h5ad']).resolve(strict=True)
            a = ad.read_h5ad(path, backed='r')
            assert a.n_vars == 556 and a.obs.slice_id.nunique() == 10
            assert a.obs.domain_true.nunique() == 6 and a.obs.cell_type_true.nunique() == 8
            assert a.obs_names.is_unique and a.var_names.is_unique
            assert 'bsigma0' in path.name
            audit.append(dict(mix=mix, modality=modality, seed=seed, source=str(path), n_obs=a.n_obs, n_vars=a.n_vars, bytes=path.stat().st_size, mtime_ns=path.stat().st_mtime_ns))
            inputs[(mix, modality, seed)] = str(path)
            a.file.close()
    assert len(inputs) == 18
    sweep, final = [], []
    for m in MODALITIES:
        for lam in LAMBDAS:
            sweep.append(dict(task=len(sweep), phase='sweep', mix='strong', modality=m, seed=2025, pipeline='banksy', k_geom=K_GEOM[m], **{'lambda':lam}, input=inputs[('strong',m,2025)], directory=str(BASE/f'sweep/strong/{m}/seed2025/lam{lam:g}')))
    for mix in ['strong','weak']:
        for m in MODALITIES:
            for seed in SEEDS:
                for pipeline in ['banksy','expression']:
                    final.append(dict(task=len(final), phase='final', mix=mix, modality=m, seed=seed, pipeline=pipeline, k_geom=K_GEOM[m], **{'lambda':None}, input=inputs[(mix,m,seed)], directory=str(BASE/f'final/{mix}/{m}/seed{seed}/{pipeline}')))
    save(BASE/'sweep_tasks.json',sweep)
    save(BASE/'final_tasks.json',final)
    save(BASE/'input_audit.json',audit)
    with (BASE/'runs.tsv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(sweep[0]),delimiter='\t');w.writeheader();w.writerows(sweep+final)
    save(BASE/'protocol.json',dict(batch_sigma=0,harmony=False,seeds=SEEDS,sweep_seed=2025,lambda_grid=LAMBDAS,k_geom=K_GEOM,max_m=1,nbr_weight_decay='scaled_gaussian',stagger_scale=5,n_pcs=30,expression_graph=dict(k=15,method='umap',random_state=0),leiden=dict(target_clusters=6,resolution_bounds=[.01,3.],max_trials=15,random_state=0,flavor='igraph',n_iterations=2),selection='Maximum domain ARI among exact-six-cluster runs on strong seed2025; exact ARI ties choose smaller lambda. No eligible run stops the dependency chain.',expression_preprocessing='normalize_total10000,log1p,scale max10,arpack PCA30; original Figure3 helper.',banksy_preprocessing='normalize_total10000 without log; all556 genes; original BANKSY feature construction and seeded PCA.',reference_comparison='Tuned-batch cell-type panel remains existing PCA/Harmony/Leiden results.',selection_caveat='Seed2025 participates in selection; not held out.',seed2025_reuse='Final strong BANKSY seed2025 links to the selected completed sweep result, avoiding a redundant fit.',spatial_graph_check='Verify no BANKSY edges connect distinct slices after original staggering.'))
    snap=BASE/'source_snapshot';snap.mkdir()
    for f in CODE.iterdir():
        if f.is_file():shutil.copy2(f,snap/f.name)
    for name in ['step01_build_banksy_matrix.py','step01_pca_harmony.py','step02_leiden_resolution_sweep.py']:
        shutil.copy2(SHARED/name,snap/name)
    save(BASE/'shared_source_sha256.json',{n:hashlib.sha256((SHARED/n).read_bytes()).hexdigest() for n in ['step01_build_banksy_matrix.py','step01_pca_harmony.py','step02_leiden_resolution_sweep.py']})
    shutil.copy2(CODE/'README.md',BASE/'README.md')
    print(f'Validated {len(inputs)} inputs; prepared {len(sweep)} sweep tasks and {len(final)} final tasks.',flush=True)

if __name__=='__main__':main()
