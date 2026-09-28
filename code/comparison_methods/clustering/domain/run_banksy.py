"""Single-slice BANKSY domain pilot; select Leiden resolution by count, never ARI."""
import argparse,json
from pathlib import Path
import numpy as np
import scanpy as sc
import anndata as ad
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn.metrics import adjusted_rand_score
from banksy.initialize_banksy import initialize_banksy
from banksy.embed_banksy import generate_banksy_matrix
from banksy_utils.umap_pca import pca_umap


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--slice-id',type=int,default=4)
    p.add_argument('--seed',type=int,default=0)
    p.add_argument('--lambda-value',type=float,default=.8)
    p.add_argument('--k-geom',type=int,default=60)
    a=p.parse_args()
    if a.out.exists():raise FileExistsError(a.out)
    a.out.mkdir(parents=True)
    source=ad.read_h5ad(a.input,backed='r')
    keep=source.obs.slice_id.astype(int).to_numpy()==a.slice_id
    x=source[keep].to_memory();source.file.close()
    assert {'domain_true','cell_type_true'}<=set(x.obs)
    original_n=x.n_obs
    sc.pp.filter_cells(x,min_counts=1);sc.pp.filter_genes(x,min_cells=1)
    assert x.obs.domain_true.notna().all() and x.obs.domain_true.nunique()>1
    k=x.obs.domain_true.nunique()
    x.obs['x']=x.obsm['spatial'][:,0];x.obs['y']=x.obsm['spatial'][:,1]
    sc.pp.normalize_total(x,target_sum=10000)
    decay='scaled_gaussian'
    b=initialize_banksy(x,('x','y','spatial'),a.k_geom,nbr_weight_decay=decay,max_m=1,
        plt_edge_hist=False,plt_nbr_weights=False,plt_agf_angles=False,plt_theta=False)
    b,_=generate_banksy_matrix(x,b,[a.lambda_value],1)
    y=b[decay][a.lambda_value]['adata'];y.obs=x.obs.copy();y.obsm['spatial']=x.obsm['spatial'].copy()
    np.random.seed(a.seed)
    pcs=min(30,y.n_obs-1,y.n_vars-1)
    pca_umap(b,pca_dims=[pcs],add_umap=False,plt_remaining_var=False)
    sc.pp.neighbors(y,n_neighbors=15,use_rep=f'reduced_pc_{pcs}',random_state=a.seed)
    # Fixed bracket and budget, including endpoints. Ground truth only supplies k.
    lo,hi=.01,3.;history=[];best=None
    for r in [lo,hi]:
        sc.tl.leiden(y,resolution=r,random_state=a.seed,key_added='cluster')
        labels=y.obs.cluster.astype(str).copy();n=labels.nunique();history.append({'resolution':r,'clusters':n})
        if best is None or abs(n-k)<best[0]:best=(abs(n-k),r,labels)
    for _ in range(15):
        r=(lo+hi)/2
        sc.tl.leiden(y,resolution=r,random_state=a.seed,key_added='cluster')
        labels=y.obs.cluster.astype(str).copy();n=labels.nunique();history.append({'resolution':r,'clusters':n})
        if abs(n-k)<best[0]:best=(abs(n-k),r,labels)
        if n==k:break
        if n<k:lo=r
        else:hi=r
    y.obs['banksy_domain']=best[2]
    ari=adjusted_rand_score(y.obs.domain_true.astype(str),best[2])
    labels=y.obs[['slice_id','domain_true','cell_type_true','banksy_domain']].copy()
    labels['x']=x.obs.x;labels['y']=x.obs.y;labels.to_csv(a.out/'labels.csv.gz')
    result={'source':str(a.input.resolve()),'slice_id':a.slice_id,'n_before_qc':original_n,
        'n_cells':y.n_obs,'n_genes':x.n_vars,'lambda':a.lambda_value,'k_geom':a.k_geom,'max_m':1,
        'preprocessing':'positive-count cells/genes; normalize_total 10000, no log; BANKSY PCA; no Harmony',
        'seed':a.seed,'n_pcs':pcs,'target_domains':int(k),'recovered_clusters':int(best[2].nunique()),
        'resolution':best[1],'domain_ari':ari,'resolution_search':history,
        'selection':'closest known domain count, first candidate on ties; no ARI tuning'}
    (a.out/'metrics.json').write_text(json.dumps(result,indent=2)+'\n')
    fig,axs=plt.subplots(1,2,figsize=(9,4))
    for ax,key in zip(axs,['domain_true','banksy_domain']):
        values=labels[key].astype('category');ax.scatter(labels.x,labels.y,c=values.cat.codes,cmap='tab10',s=1,rasterized=True)
        ax.set_title('True domains' if key=='domain_true' else f'BANKSY clusters (ARI {ari:.3f})')
        ax.set_aspect('equal');ax.axis('off')
    fig.tight_layout()
    for ext in ['png','pdf']:fig.savefig(a.out/f'domain_recovery.{ext}',dpi=250)
    plt.close(fig)
    print(json.dumps(result),flush=True)

if __name__=='__main__':main()
