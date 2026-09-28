"""Paired pre/post-Harmony K=15 graphs on identical saved full-cell embeddings."""
import argparse,json,gc,importlib.metadata
from pathlib import Path
import numpy as np
import anndata as ad
import scanpy as sc
from scipy.sparse.csgraph import connected_components
ROOT=Path(__file__).resolve().parents[5]
OUT=ROOT/'sim_paper/data/comparison_methods/clustering/cell_type/experiment_600k'

def stats(g,slices):
 g=g.copy();g.eliminate_zeros()
 n,labels=connected_components(g,directed=False);sizes=np.bincount(labels)
 degrees=np.diff(g.indptr)
 # Unique slice/component pairs avoid a large dense component-by-cell matrix.
 pairs=np.unique(np.column_stack([labels,slices]),axis=0)
 ns=np.bincount(pairs[:,0],minlength=n)
 return {'components':int(n),'largest_components':np.sort(sizes)[-10:][::-1].tolist(),'median_component_cells':float(np.median(sizes)),'singletons':int((sizes==1).sum()),'cells_in_components_le100':int(sizes[sizes<=100].sum()),'components_spanning_multiple_slices':int((ns>1).sum()),'fraction_cells_in_multislice_components':float(sizes[ns>1].sum()/len(labels)),'degree_quantiles':np.quantile(degrees,[0,.25,.5,.75,1]).tolist(),'nonzero_edges_directed':int(g.nnz)}

def main():
 p=argparse.ArgumentParser();p.add_argument('--genes',type=int,required=True);a=p.parse_args()
 run=OUT/f'runs/spider_genes{a.genes}_seed2025';dest=OUT/'harmony_graph_diagnostic';dest.mkdir(exist_ok=True)
 x=ad.read_h5ad(run/'clustering/embeddings_labels.h5ad');assert x.n_obs==600000
 slices=x.obs.slice_id.astype(int).to_numpy()
 report={'run':str(run),'n_cells':x.n_obs,'n_genes':a.genes,'n_neighbors':15,'method':'umap','metric':'euclidean','random_state':0,'scanpy_version':importlib.metadata.version('scanpy'),'design':'Identical cells/order, graph settings and software; change only pre/post Harmony embedding. No clustering or parameter selection.'}
 report['saved_post_harmony']=stats(x.obsp['connectivities'],slices)
 for basis in ['X_pca_pre_harmony','X_pca_harmony']:
  graph=ad.AnnData(obs=x.obs[[]].copy());graph.obsm['embedding']=x.obsm[basis].copy()
  sc.pp.neighbors(graph,n_neighbors=15,use_rep='embedding',method='umap',metric='euclidean',random_state=0)
  report[basis]=stats(graph.obsp['connectivities'],slices)
  print(basis,json.dumps(report[basis]),flush=True)
  (dest/f'spider_genes{a.genes}_seed2025.json').write_text(json.dumps(report,indent=2)+'\n')
  del graph;gc.collect()
 print('DONE',flush=True)
if __name__=='__main__':main()
