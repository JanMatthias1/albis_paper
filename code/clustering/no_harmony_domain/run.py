"""Build uncorrected embeddings, then cluster in the existing Leiden environment."""
import argparse
import gc
import hashlib
import importlib.metadata
import json
import time
from pathlib import Path
import numpy as np
import anndata as ad
import scanpy as sc
from common import *

def check_sources(root):
    for name,digest in json.loads((root/'shared_source_sha256.json').read_text()).items():
        assert hashlib.sha256((SHARED/name).read_bytes()).hexdigest()==digest, f'Shared source changed: {name}'

def embed(root,task):
    dest=Path(task['directory'])
    if (dest/'embedding.h5ad').exists():
        assert json.loads((dest/'embedding_provenance.json').read_text())['task']==task
        return
    dest.mkdir(parents=True,exist_ok=True)
    save(dest/'status.json',dict(stage='embedding',status='running',task=task))
    a=ad.read_h5ad(task['input'])
    assert a.n_vars==556 and a.obs.slice_id.nunique()==10
    obs=a.obs.copy();spatial=np.asarray(a.obsm['spatial']).copy()
    report=dict(task=task,harmony_applied=False,n_obs=a.n_obs,n_genes=a.n_vars,versions={n:importlib.metadata.version(n) for n in ['scanpy','anndata','numpy','scipy','scikit-learn']})
    started=time.monotonic()
    if task['pipeline']=='expression':
        load_shared('plain_pca','01.2_pca_harmony.py').preprocess_for_pca(a,30,0,10000)
        pcs=a.obsm['X_pca_pre_harmony'].copy()
    else:
        from banksy.initialize_banksy import initialize_banksy
        from banksy.embed_banksy import generate_banksy_matrix
        from banksy_utils.umap_pca import pca_umap
        xy=spatial[:,:2].astype(float).copy()
        slices=sorted(a.obs.slice_id.unique())
        masks=[(a.obs.slice_id==s).to_numpy() for s in slices]
        step=max(float(np.ptp(xy[m,0])) for m in masks)*5
        for i,mask in enumerate(masks):xy[mask,0]=xy[mask,0]-xy[mask,0].min()+i*step
        a.obsm['spatial']=xy;a.obs['x_stagger']=xy[:,0];a.obs['y_stagger']=xy[:,1]
        sc.pp.normalize_total(a,target_sum=10000)
        b=initialize_banksy(a,('x_stagger','y_stagger','spatial'),task['k_geom'],nbr_weight_decay='scaled_gaussian',max_m=1,plt_edge_hist=False,plt_nbr_weights=False,plt_agf_angles=False,plt_theta=False)
        codes=a.obs.slice_id.astype(str).to_numpy()
        audit={}
        for m,w in b['scaled_gaussian']['weights'].items():
            w=w.tocsr();cross=0
            for start in range(0,a.n_obs,10000):
                chunk=w[start:start+10000].tocoo()
                cross+=int(np.count_nonzero(codes[chunk.row+start]!=codes[chunk.col]))
            audit[str(m)]=dict(edges=int(w.nnz),cross_slice_edges=cross)
            assert cross==0, f'Cross-slice spatial edges at m={m}: {cross}'
        report['spatial_graph']=audit
        b,_=generate_banksy_matrix(a,b,[task['lambda']],1)
        matrix=b['scaled_gaussian'][task['lambda']]['adata']
        b={'scaled_gaussian':{task['lambda']:{'adata':matrix}}}
        np.random.seed(0)
        pca_umap(b,pca_dims=[30],add_umap=False,plt_remaining_var=False)
        pcs=matrix.obsm['reduced_pc_30'].copy()
        del matrix,b
    assert pcs.shape==(len(obs),30) and np.isfinite(pcs).all()
    embedding=ad.AnnData(obs=obs)
    embedding.obsm['X_pca_pre_harmony']=pcs
    embedding.obsm['spatial']=spatial
    embedding.uns['harmony_applied']=False
    del a;gc.collect()
    embedding.write_h5ad(dest/'embedding.h5ad')
    report['seconds']=time.monotonic()-started
    save(dest/'embedding_provenance.json',report)
    save(dest/'status.json',dict(stage='embedding',status='complete'))
    print('Embedding saved',dest,flush=True)

def cluster(root,task):
    from scipy.sparse.csgraph import connected_components
    from sklearn.metrics import adjusted_rand_score
    dest=Path(task['directory'])
    if (dest/'metrics.json').exists():
        assert json.loads((dest/'metrics.json').read_text())['task']==task
        return
    a=ad.read_h5ad(dest/'embedding.h5ad')
    assert 'X_pca_harmony' not in a.obsm and not a.uns['harmony_applied']
    save(dest/'status.json',dict(stage='clustering',status='running'))
    sc.pp.neighbors(a,n_neighbors=15,use_rep='X_pca_pre_harmony',random_state=0,method='umap')
    assert np.isfinite(a.obsp['connectivities'].data).all()
    nc,lab=connected_components(a.obsp['connectivities'],directed=False)
    diagnostic=dict(connected_components=int(nc),largest_components=np.sort(np.bincount(lab))[-10:][::-1].tolist(),expression_graph_k=15,weighting='umap',harmony_applied=False)
    save(dest/'graph_diagnostic.json',diagnostic)
    sweep=load_shared('leiden_sweep','02_leiden_resolution_sweep.py')
    target=int(a.obs.domain_true.nunique());assert target==6
    res,k,trials=sweep.find_resolution_for_k(a,target,.01,3.,15,0,'leiden_domain_true')
    result=dict(task=task,n_obs=a.n_obs,n_slices=int(a.obs.slice_id.nunique()),domain_ari=float(adjusted_rand_score(a.obs.domain_true.astype(str),a.obs.leiden_domain_true.astype(str))),slice_ari=float(adjusted_rand_score(a.obs.slice_id.astype(str),a.obs.leiden_domain_true.astype(str))),target_clusters=target,achieved_clusters=int(k),resolution=float(res),trials=trials,harmony_applied=False,graph=diagnostic)
    a.obs.to_csv(dest/'labels.csv.gz')
    a.write_h5ad(dest/'embeddings_labels.h5ad')
    save(dest/'metrics.json',result)
    save(dest/'status.json',dict(stage='clustering',status='complete' if k==target else 'cluster_count_mismatch'))
    print(json.dumps(result),flush=True)

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--phase',choices=['sweep','final'],required=True)
    p.add_argument('--task',type=int,required=True)
    p.add_argument('--stage',choices=['embed','cluster'],required=True)
    p.add_argument('--root',type=Path,default=BASE)
    args=p.parse_args();check_sources(args.root)
    task=task_for(args.root,args.phase,args.task)
    # Retain the exact selection-seed result in the final summary, without re-fitting.
    if args.phase=='final' and task['pipeline']=='banksy' and task['mix']=='strong' and task['seed']==2025:
        dest=Path(task['directory']);selected=args.root/f"sweep/strong/{task['modality']}/seed2025/lam{task['lambda']:g}"
        assert (selected/'metrics.json').exists()
        dest.parent.mkdir(parents=True,exist_ok=True)
        if not dest.exists():dest.symlink_to(selected,target_is_directory=True)
        assert dest.resolve()==selected.resolve()
        print('Reusing selected seed2025 sweep result:',selected,flush=True);return
    if args.stage=='embed':embed(args.root,task)
    else:cluster(args.root,task)

if __name__=='__main__':main()
