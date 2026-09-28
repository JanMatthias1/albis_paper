"""Full-cell KNN graph sensitivity and common-pipeline rerun; no profile collapsing."""
import argparse,json,gc
from pathlib import Path
import numpy as np
import anndata as ad
import scanpy as sc
from scipy.sparse.csgraph import connected_components
from sklearn.metrics import adjusted_rand_score
import cluster
ROOT=cluster.ROOT
OUT=ROOT/'sim_paper/data/comparison_methods/clustering/cell_type/experiment_600k'
CANDIDATES=[(128,'umap'),(256,'umap'),(256,'gauss')]
def build(run,k,method,pilot):
 old=json.loads((run/'clustering/metrics.json').read_text())
 basis='X_pca_harmony' if old['method']=='albis' else 'X_pca_pre_harmony'
 dest=run/(f'knn_pre_harmony_pilot_{method}{k}' if pilot else 'clustering_full_knn_method_harmony');dest.mkdir(exist_ok=False)
 a=ad.read_h5ad(run/'clustering/embeddings_labels.h5ad')
 assert a.n_obs==old['n_cells_clustered']
 a.obsp.clear();a.uns.pop('neighbors',None);gc.collect()
 print('Building full graph',a.n_obs,k,method,flush=True)
 sc.pp.neighbors(a,n_neighbors=k,use_rep=basis,method=method,random_state=0)
 weights=a.obsp['connectivities'];assert np.isfinite(weights.data).all()
 n,components=connected_components(weights,directed=False);sizes=np.sort(np.bincount(components))[::-1]
 report=dict(old,graph_basis=basis,harmony_applied=old['method']=='albis',harmony_batch='slice_id' if old['method']=='albis' else None,n_neighbors=k,neighbor_weighting=method,connected_components=int(n),largest_components=sizes[:12].tolist(),n_graph_nodes=a.n_obs,collapsed_profiles=False)
 (dest/'graph_diagnostic.json').write_text(json.dumps(report,indent=2)+'\n');print('Components',n,'largest',sizes[:12],flush=True)
 if pilot:return
 res,achieved,trials=cluster.load('sweep','02_leiden_resolution_sweep.py').find_resolution_for_k(a,8,.01,3,15,0,'leiden_cell_type')
 report.update(resolution=float(res),achieved_clusters=int(achieved),trials=trials,cell_type_ari=float(adjusted_rand_score(a.obs.cell_type_true.astype(str),a.obs.leiden_cell_type.astype(str))))
 a.obs.to_csv(dest/'labels.csv.gz');a.write_h5ad(dest/'embeddings_labels.h5ad')
 (dest/'metrics.json').write_text(json.dumps(report,indent=2)+'\n')
 print(json.dumps(report),flush=True)
def main():
 p=argparse.ArgumentParser();p.add_argument('--pilot',type=int);p.add_argument('--task',type=int);p.add_argument('--select',action='store_true');args=p.parse_args()
 pilotrun=OUT/'runs/spider_genes556_seed2025'
 if args.pilot is not None:build(pilotrun,*CANDIDATES[args.pilot],True);return
 if args.select:
  diagnostics=[]
  for k,m in CANDIDATES:
   path=pilotrun/f'knn_pre_harmony_pilot_{m}{k}/graph_diagnostic.json'
   if not path.exists():raise RuntimeError(f'Missing pilot: {path}')
   d=json.loads(path.read_text());diagnostics.append(d)
  acceptable=[d for d in diagnostics if d['connected_components']<=8]
  if not acceptable:raise RuntimeError('No pilot has <=8 components. Inspect diagnostics; do not launch final graph.')
  chosen=acceptable[0]
  report={'n_neighbors':chosen['n_neighbors'],'method':chosen['neighbor_weighting'],'selection':'First candidate with <=8 connected components; no labels/ARI used. Preference: UMAP128, UMAP256, Gaussian256.','candidates':[{'k':d['n_neighbors'],'method':d['neighbor_weighting'],'components':d['connected_components']} for d in diagnostics],'all_cells_retained':True,'profile_collapse':False}
  (OUT/'knn_protocol_method_harmony.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report));return
 if args.task is not None:
  r=json.loads((OUT/'runs.json').read_text())[args.task];protocol=json.loads((OUT/'knn_protocol_method_harmony.json').read_text());build(Path(r['directory']),protocol['n_neighbors'],protocol['method'],False);return
 p.error('Specify --pilot, --select, or --task')
if __name__=='__main__':main()
