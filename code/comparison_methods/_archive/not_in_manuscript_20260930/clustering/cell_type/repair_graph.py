"""Duplicate-aware neighbor graph: exact raw profiles within slice, no truth in grouping."""
import argparse,json,hashlib,gc
from pathlib import Path
import numpy as np
import anndata as ad
import scanpy as sc
from scipy import sparse
from scipy.sparse.csgraph import connected_components
from sklearn.metrics import adjusted_rand_score
import cluster

def main():
 p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);p.add_argument('--group-across-slices',action='store_true');a=p.parse_args()
 run=a.run;old=json.loads((run/'clustering/metrics.json').read_text());method=old['method']
 dest=run/('clustering_unique_profiles' if a.group_across_slices else 'clustering_duplicate_aware');dest.mkdir(exist_ok=False)
 raw=ad.read_h5ad(run/method/'cell.h5ad');x=raw.X.tocsr();x.sort_indices();x.eliminate_zeros()
 codes={};inverse=np.empty(raw.n_obs,dtype=np.int32);first=[]
 for i,sid in enumerate(raw.obs.slice_id.astype(str)):
  lo,hi=x.indptr[i:i+2];h=hashlib.sha256(x.indices[lo:hi].tobytes()+x.data[lo:hi].tobytes()).digest();key=h if a.group_across_slices else (sid,h)
  if key not in codes:codes[key]=len(first);first.append(i)
  inverse[i]=codes[key]
 ids=raw.obs_names.copy();del raw,x,codes;gc.collect()
 emb=ad.read_h5ad(run/'clustering/embeddings_labels.h5ad');assert ids.equals(emb.obs_names)
 idx=np.array(first);n=len(idx);counts=np.bincount(inverse)
 if n==emb.n_obs:
  report=dict(old,n_unique_profile_slice=n,duplicate_cells=0,graph_policy='No exact raw duplicates: original 15NN/PCA/Harmony/Leiden results retained')
  (dest/'metrics.json').write_text(json.dumps(report,indent=2)+'\n');print('No duplicates; retained original analysis',flush=True);return
 # Mean numerical PCA/Harmony representations for identical raw-profile/slice groups.
 coords=np.zeros((n,emb.obsm['X_pca_harmony'].shape[1]));np.add.at(coords,inverse,emb.obsm['X_pca_harmony']);coords/=counts[:,None]
 graph=ad.AnnData(obs=emb.obs.iloc[idx].copy());graph.obsm['X_pca_harmony']=coords
 sc.pp.neighbors(graph,n_neighbors=15,use_rep='X_pca_harmony',random_state=0)
 nc,components=connected_components(graph.obsp['connectivities'],directed=False)
 report=dict(old,n_unique_profile_slice=n,duplicate_cells=emb.n_obs-n,connected_components=nc,largest_component=int(np.bincount(components).max()),graph_policy='Exact raw-expression SHA256'+(' across slices' if a.group_across_slices else ' + slice_id')+'; mean corrected embedding per group; 15NN; unweighted unique-profile graph; broadcast labels to original cells')
 (dest/'diagnostic.json').write_text(json.dumps(report,indent=2)+'\n');print('Unique',n,'components',nc,flush=True)
 res,k,trials=cluster.load('sweep','step02_leiden_resolution_sweep.py').find_resolution_for_k(graph,8,.01,3,15,0,'leiden_cell_type')
 labels=graph.obs.leiden_cell_type.astype(str).to_numpy()[inverse]
 report.update(cell_type_ari=float(adjusted_rand_score(emb.obs.cell_type_true.astype(str),labels)),achieved_clusters=int(k),resolution=float(res),trials=trials)
 emb.obs['leiden_cell_type']=labels;emb.obs.to_csv(dest/'labels.csv.gz');graph.write_h5ad(dest/'unique_profile_graph.h5ad');np.save(dest/'cell_to_profile.npy',inverse)
 (dest/'metrics.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report),flush=True)
if __name__=='__main__':main()
