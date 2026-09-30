"""Cluster original pools and audit exact profile reuse; no production results modified."""
import argparse,json,hashlib,gc
from pathlib import Path
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
from scipy.io import mmread
from scipy.sparse.csgraph import connected_components
from sklearn.metrics import adjusted_rand_score
import cluster
BASE=cluster.ROOT/'sim_paper/data/figure_5/clustering/cell_type'
OUT=BASE/'spider_pool_diagnostic_20260928'
def hashes(x):
 x=x.tocsr();x.sum_duplicates();x.eliminate_zeros();x.sort_indices()
 return [hashlib.sha256(x.indices[x.indptr[i]:x.indptr[i+1]].astype('<i8').tobytes()+x.data[x.indptr[i]:x.indptr[i+1]].astype('<f8').tobytes()).hexdigest() for i in range(x.shape[0])]
def main():
 p=argparse.ArgumentParser();p.add_argument('--task',type=int,required=True);task=p.parse_args().task
 genes,seed=[(g,s) for g in [556,2000] for s in [2025,101,202]][task]
 dest=OUT/f'genes{genes}_seed{seed}';dest.mkdir(exist_ok=False)
 pool=BASE/f'experiment_600k/pools/genes{genes}_seed{seed}/splatter'
 meta=pd.read_csv(pool/'cells.tsv',sep='\t').set_index('cell_id')
 a=ad.AnnData(mmread(pool/'counts.mtx').T.tocsr(),obs=meta,var=pd.DataFrame(index=(pool/'genes.tsv').read_text().splitlines()))
 keys=hashes(a.X);lookup={}
 for h,t in zip(keys,a.obs.cell_type_true):lookup.setdefault(h,set()).add(str(t))
 assert a.shape==(10000,genes)
 report=dict(genes=genes,seed=seed,pool=str(pool),pool_cells=a.n_obs,pool_unique_profiles=len(lookup),pool_profiles_with_multiple_labels=sum(len(v)>1 for v in lookup.values()))
 np.random.seed(0);cluster.load('prep','step01_pca_harmony.py').preprocess_for_pca(a,30,0,10000)
 results=[]
 for k,method in [(256,'gauss'),(15,'umap')]:
  sc.pp.neighbors(a,n_neighbors=k,use_rep='X_pca_pre_harmony',method=method,random_state=0)
  nc,comp=connected_components(a.obsp['connectivities'],directed=False)
  res,n,trials=cluster.load('sweep','step02_leiden_resolution_sweep.py').find_resolution_for_k(a,8,.01,3,15,0,'leiden_cell_type')
  result=dict(k=k,weighting=method,components=int(nc),component_ari=float(adjusted_rand_score(a.obs.cell_type_true,comp)),clusters=int(n),ari=float(adjusted_rand_score(a.obs.cell_type_true,a.obs.leiden_cell_type)),resolution=float(res),trials=trials)
  results.append(result);print(result,flush=True)
  a.obs.to_csv(dest/f'pool_labels_{method}{k}.csv.gz')
 report['pool_clustering']=results
 (dest/'pool_metrics.json').write_text(json.dumps(report,indent=2)+'\n')
 del a;gc.collect()
 source=BASE/f'experiment_600k/runs/spider_genes{genes}_seed{seed}/spider/cell.h5ad'
 a=ad.read_h5ad(source);assert a.shape==(600000,genes)
 expanded=hashes(a.X);counts={};unmatched=0;wrong=0
 for h,t in zip(expanded,a.obs.cell_type_true.astype(str)):
  counts[h]=counts.get(h,0)+1
  if h not in lookup:unmatched+=1
  elif t not in lookup[h]:wrong+=1
 freq=np.array(list(counts.values()))
 report.update(spider_input=str(source),spider_cells=a.n_obs,spider_unique_profiles=len(counts),spider_profiles_not_in_pool=len(set(counts)-set(lookup)),spider_cells_not_in_pool=unmatched,spider_cells_with_label_mismatch=wrong,pool_profiles_used=len(set(counts)&set(lookup)),profile_multiplicity_quantiles=np.quantile(freq,[0,.25,.5,.75,1]).tolist())
 final=BASE/f'final_600k_batchzero_noharmony_20260928/spider_genes{genes}_seed{seed}/metrics.json'
 m=json.loads(final.read_text());report['expanded_clustering']={k:m[k] for k in ['cell_type_ari','achieved_clusters','connected_components','n_neighbors','neighbor_weighting']}
 (dest/'diagnostic.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report),flush=True)
if __name__=='__main__':main()
