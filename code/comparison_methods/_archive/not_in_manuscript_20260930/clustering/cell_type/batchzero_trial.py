"""All-method no-Harmony Gaussian K256 trial, with requested ALBIS batch-zero data."""
import argparse,json,gc
from pathlib import Path
import numpy as np,pandas as pd,anndata as ad,scanpy as sc
from scipy.sparse.csgraph import connected_components
from sklearn.metrics import adjusted_rand_score
import cluster
ROOT=cluster.ROOT
BASE=ROOT/'sim_paper/data/figure_5/clustering/cell_type'
OUT=BASE/'batchzero_noharmony_gauss256'
ALBIS=ROOT/'sim_paper/data/figure_3/strong_domain_mix/batch_sigma_slide/cell'

def prepare():
 OUT.mkdir(exist_ok=False);runs=[]
 old=json.loads((BASE/'experiment_600k/runs.json').read_text())
 for row in old:
  r=dict(row);r['old_directory']=r.pop('directory');r['directory']=str(OUT/f"{r['method']}_genes{r['genes']}_seed{r['seed']}")
  if r['method']=='albis':
   folder=ALBIS/('bs0' if r['seed']==2025 else f"bs0_seed{r['seed']}")
   paths=list(folder.glob('*_qc.h5ad'));assert len(paths)==1,paths;r['input']=str(paths[0])
  else:r['input']=str(Path(r['old_directory'])/r['method']/'cell.h5ad')
  runs.append(r)
 (OUT/'runs.json').write_text(json.dumps(runs,indent=2)+'\n')
 (OUT/'protocol.json').write_text(json.dumps({'harmony':False,'albis_batch_sigma':0,'albis_source':str(ALBIS),'n_neighbors':256,'edge_weighting':'gauss','n_pcs':30,'clustering_seed':0,'target_clusters':8,'profile_collapse':False,'caveat':'ALBIS requested source has about 24k cells per seed versus 600k competitor cells. All observations retained except zero-total cells. Same K/weighting across methods; not an equal-cell-number benchmark.'},indent=2)+'\n')

def run(task):
 r=json.loads((OUT/'runs.json').read_text())[task];dest=Path(r['directory']);dest.mkdir(exist_ok=False)
 saved=Path(r['old_directory'])/'clustering/embeddings_labels.h5ad'
 if r['method']!='albis' and saved.exists():
  a=ad.read_h5ad(saved);a.obsp.clear();a.uns.pop('neighbors',None)
  a.obsm.pop('X_pca_harmony',None);origin='Saved pre-Harmony PCA from the original all-cell preprocessing'
 else:
  x=ad.read_h5ad(r['input']);keep=np.asarray(x.X.sum(1)).ravel()>0
  if not keep.all():x=x[keep].copy()
  cluster.load('prep','step01_pca_harmony.py').preprocess_for_pca(x,30,0,10000)
  a=ad.AnnData(obs=x.obs.copy());a.obsm['X_pca_pre_harmony']=x.obsm['X_pca_pre_harmony'].copy();del x;gc.collect();origin='Fresh PCA using shared preprocessing; no Harmony'
 assert a.n_vars==0 and a.obs.cell_type_true.nunique()==8
 print('Building',r['method'],a.n_obs,'K256 Gaussian',flush=True)
 sc.pp.neighbors(a,n_neighbors=256,use_rep='X_pca_pre_harmony',method='gauss',random_state=0)
 assert np.isfinite(a.obsp['connectivities'].data).all()
 nc,lab=connected_components(a.obsp['connectivities'],directed=False)
 report=dict(r,n_cells_clustered=a.n_obs,n_slices=int(a.obs.slice_id.nunique()),harmony_applied=False,batch_sigma=0,n_neighbors=256,neighbor_weighting='gauss',graph_basis='X_pca_pre_harmony',connected_components=int(nc),largest_components=np.sort(np.bincount(lab))[-10:][::-1].tolist(),pca_origin=origin)
 (dest/'graph_diagnostic.json').write_text(json.dumps(report,indent=2)+'\n')
 resolution,k,trials=cluster.load('sweep','step02_leiden_resolution_sweep.py').find_resolution_for_k(a,8,.01,3,15,0,'leiden_cell_type')
 report.update(cell_type_ari=float(adjusted_rand_score(a.obs.cell_type_true.astype(str),a.obs.leiden_cell_type.astype(str))),achieved_clusters=int(k),resolution=float(resolution),trials=trials)
 a.obs.to_csv(dest/'labels.csv.gz');a.write_h5ad(dest/'embeddings_labels.h5ad')
 (dest/'metrics.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report),flush=True)

def aggregate():
 import matplotlib
 matplotlib.use('Agg')
 import matplotlib.pyplot as plt
 rows=[];status=[]
 for r in json.loads((OUT/'runs.json').read_text()):
  p=Path(r['directory'])/'metrics.json';d=json.loads(p.read_text()) if p.exists() else None
  status.append(dict(r,status='missing' if d is None else ('complete' if d['achieved_clusters']==8 else 'cluster_count_mismatch')))
  if d:rows.append(d)
 pd.DataFrame(status).to_csv(OUT/'status.csv',index=False)
 if rows:pd.DataFrame(rows).drop(columns=['trials']).to_csv(OUT/'metrics.csv',index=False)
 if len(rows)!=15 or any(d['achieved_clusters']!=8 for d in rows):raise RuntimeError('Incomplete or count-mismatched trial; inspect status.csv')
 frame=pd.DataFrame(rows);conditions=[('albis',556),('sccube',556),('sccube',2000),('spider',556),('spider',2000)]
 fig,ax=plt.subplots(figsize=(9,5))
 ax.boxplot([frame[(frame.method==m)&(frame.genes==g)].cell_type_ari for m,g in conditions],showfliers=False)
 for j,seed in enumerate([2025,101,202]):
  vals=[frame[(frame.method==m)&(frame.genes==g)&(frame.seed==seed)].cell_type_ari.item() for m,g in conditions]
  ax.scatter(np.arange(1,6)+(j-1)*.08,vals,label=f'Seed {seed}',zorder=3)
 ax.set(xticks=range(1,6),xticklabels=['ALBIS\n556','scCube\n556','scCube\n2000','SPIDER\n556','SPIDER\n2000'],ylabel='Cell-type ARI',ylim=(-.05,1.05),title='No Harmony · Gaussian KNN (K=256)');ax.legend(frameon=False)
 fig.text(.5,.015,'Gene counts shown; ALBIS ~24k cells, competitors 600k cells per seed',ha='center');fig.tight_layout(rect=(0,.04,1,1))
 for ext in ['png','pdf','svg']:fig.savefig(OUT/f'cell_type_ari_boxplot.{ext}',dpi=220)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--prepare',action='store_true');p.add_argument('--task',type=int);p.add_argument('--aggregate',action='store_true');a=p.parse_args()
 if a.prepare:prepare()
 elif a.aggregate:aggregate()
 elif a.task is not None:run(a.task)
 else:p.error('Choose prepare, task, or aggregate')
