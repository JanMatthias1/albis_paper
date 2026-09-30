"""Common all-cell PCA/Harmony/Leiden analysis; reuse Figure 3 implementations."""
import argparse,json,importlib.util,hashlib,gc,sys
from pathlib import Path
import numpy as np
import anndata as ad
import scanpy as sc
from scipy import sparse
from sklearn.metrics import adjusted_rand_score
ROOT=Path(__file__).resolve().parents[5]
SHARED=ROOT/'sim_paper/code/clustering'
sys.path.insert(0,str(SHARED))

def load(name,filename):
    spec=importlib.util.spec_from_file_location(name,SHARED/filename)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module

def main():
    p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True)
    p.add_argument('--method',choices=['albis','spider','sccube'],required=True)
    p.add_argument('--genes',type=int,required=True);a=p.parse_args()
    dest=a.run/'clustering';dest.mkdir(exist_ok=False)
    source=a.run/a.method/'cell.h5ad'
    print('Loading',source,flush=True);x=ad.read_h5ad(source)
    assert x.n_obs==600000 and x.n_vars==a.genes
    assert x.obs.cell_type_true.nunique()==8 and x.obs.slice_id.nunique()==10
    assert x.obs_names.is_unique and x.var_names.is_unique
    vals=x.X.data if sparse.issparse(x.X) else np.asarray(x.X)
    assert np.isfinite(vals).all() and (vals>=0).all()
    del vals
    raw_n=x.n_obs;keep=np.asarray(x.X.sum(axis=1)).ravel()>0
    if not keep.all():x=x[keep].copy()
    x.obs['slice_id']=x.obs.slice_id.astype(str).astype('category')
    prep=load('shared_pca','01.2_pca_harmony.py')
    sweep=load('shared_leiden','02_leiden_resolution_sweep.py')
    np.random.seed(0)
    print('PCA',flush=True);prep.preprocess_for_pca(x,30,0,10000)
    print('Harmony',flush=True);prep.run_harmony(x,'slice_id','X_pca_pre_harmony','X_pca_harmony')
    # Save embeddings without another copy of the large expression matrix.
    embedding=ad.AnnData(obs=x.obs.copy())
    for key in ['X_pca_pre_harmony','X_pca_harmony']:embedding.obsm[key]=x.obsm[key].copy()
    for key in ['spatial','spatial_3d']:
        if key in x.obsm:embedding.obsm[key]=x.obsm[key].copy()
    del x;gc.collect()
    print('Neighbor graph and count-matched Leiden',flush=True)
    sc.pp.neighbors(embedding,n_neighbors=15,use_rep='X_pca_harmony',random_state=0)
    resolution,k,trials=sweep.find_resolution_for_k(embedding,8,.01,3.,15,0,'leiden_cell_type')
    ari=adjusted_rand_score(embedding.obs.cell_type_true.astype(str),embedding.obs.leiden_cell_type.astype(str))
    embedding.obs.to_csv(dest/'labels.csv.gz')
    embedding.write_h5ad(dest/'embeddings_labels.h5ad')
    cfg=json.loads((a.run/'settings.json').read_text())
    result={'method':a.method,'genes':a.genes,'simulation_seed':cfg['seed'],'clustering_seed':0,
        'n_cells_generated':raw_n,'n_cells_clustered':embedding.n_obs,'n_slices':10,
        'batch_sigma':cfg.get('batch_sigma',0),'cell_type_ari':float(ari),
        'target_clusters':8,'achieved_clusters':int(k),'resolution':float(resolution),'trials':trials,
        'input':str(source.resolve()),'normalization':10000,'log1p':True,'scale_clip':10,
        'n_pcs':30,'n_neighbors':15,'harmony_batch':'slice_id',
        'selection':'Closest known cluster count, not maximum ARI',
        'shared_code_sha256':{name:hashlib.sha256((SHARED/name).read_bytes()).hexdigest() for name in ['01.2_pca_harmony.py','02_leiden_resolution_sweep.py']}}
    (dest/'metrics.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result),flush=True)
if __name__=='__main__':main()
