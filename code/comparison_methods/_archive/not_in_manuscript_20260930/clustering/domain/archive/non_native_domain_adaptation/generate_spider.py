"""Compose native SPIDER gyrus and zone-specific cell-type solvers."""
import argparse, json, importlib.metadata
from dataclasses import asdict
from pathlib import Path
import numpy as np
import pandas as pd
import anndata as ad
from scipy.io import mmread
from scipy.sparse.csgraph import connected_components
import spider
from spider.api import simulate_cells, make_transition_matrix
from spider.solver import AnnealingConfig, solve_cell_types_by_zone
from spider.sim_expr import get_sim_cell_level_expr


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--pool',type=Path,required=True)
    p.add_argument('--n-cells',type=int,default=10000)
    p.add_argument('--seed',type=int,default=20260921)
    p.add_argument('--iterations',type=int,default=50000)
    p.add_argument('--n-slices',type=int,default=10)
    a=p.parse_args()
    assert importlib.metadata.version('st-spider')=='1.2.0'
    if a.out.exists():raise FileExistsError(a.out)
    a.out.mkdir(parents=True)
    np.random.seed(a.seed)
    config=AnnealingConfig(max_iter=a.iterations,random_state=a.seed,verbose=True)
    print(f'Generating {a.n_cells} cells: four 3D gyrus groups',flush=True)
    tissue=simulate_cells(n_cells=a.n_cells,n_celltypes=4,prior=[.25]*4,
        pattern='gyrus',dimensions=3,plate_shape=[4100]*3,config=config,random_state=a.seed)
    domains=tissue.labels.copy()
    pd.DataFrame(tissue.history).to_csv(a.out/'domain_optimization_history.csv',index=False)
    print(f'Gyrus completed; transition loss {tissue.loss}. Assigning eight cell types within zones.',flush=True)
    priors={i:np.array([.30 if j//2==i else .4/6 for j in range(8)]) for i in range(4)}
    transitions={i:np.tile(priors[i],(8,1)) for i in range(4)}
    zoned=solve_cell_types_by_zone(tissue.adjacency,domains,transitions,8,priors=priors,config=config)
    labels=zoned['labels']
    print('Cell-type assignment completed; sampling expression.',flush=True)
    assert set(domains)==set(range(4)) and set(labels)==set(range(8))
    meta=pd.read_csv(a.pool/'cells.tsv',sep='\t')
    types=[f'type{i}' for i in range(1,9)]
    reference=ad.AnnData(mmread(a.pool/'counts.mtx').T.tocsr(),
        obs=meta.set_index('cell_id'),var=pd.DataFrame(index=(a.pool/'genes.tsv').read_text().splitlines()))
    expression=get_sim_cell_level_expr(celltype_assignment=labels,adata=reference,
        Num_celltype=8,Num_ct_sample=np.bincount(labels,minlength=8),match_list=types,ct_key='cell_type_true')
    assert np.array_equal(expression.obs.cell_type_true.astype(str).to_numpy(),np.array(types)[labels])
    obs=pd.DataFrame({'domain_true':pd.Categorical([f'D{i}' for i in domains]),
        'cell_type_true':pd.Categorical(np.array(types)[labels]),
        'slice_id':np.minimum((tissue.coordinates[:,2]/4100*a.n_slices).astype(int),a.n_slices-1)},
        index=[f'cell{i}' for i in range(a.n_cells)])
    out=ad.AnnData(expression.X.copy(),obs=obs,var=expression.var.copy())
    out.obsm['spatial_3d']=tissue.coordinates;out.obsm['spatial']=tissue.coordinates[:,:2]
    assert np.isfinite(out.obsm['spatial_3d']).all() and out.obs.slice_id.nunique()==a.n_slices
    assert out.n_obs==a.n_cells and out.obs_names.is_unique and out.var_names.is_unique
    out.obs.to_csv(a.out/'annotations.csv.gz')
    components={}
    for i in range(4):
        idx=np.flatnonzero(domains==i)
        n,lab=connected_components(tissue.adjacency[idx][:,idx],directed=False)
        components[f'D{i}']={'components':int(n),'largest_fraction':float(np.bincount(lab).max()/len(idx))}
    composition=pd.crosstab(obs.domain_true,obs.cell_type_true)
    composition.to_csv(a.out/'domain_celltype_counts.csv')
    audit={'version':importlib.metadata.version('st-spider'),'n_cells':a.n_cells,'n_genes':out.n_vars,
        'seed':a.seed,'n_domains':4,'n_celltypes':8,'n_slices':a.n_slices,'extent_um':4100,
        'slice_counts':out.obs.slice_id.value_counts().sort_index().to_dict(),
        'pattern':'gyrus','pattern_strength':.8,'max_iterations':a.iterations,
        'solver_config':asdict(config),
        'domain_tolerance_reached':bool(tissue.loss <= config.tol),
        'domain_prior':[.25]*4,'celltype_priors':{str(k):v.tolist() for k,v in priors.items()},
        'domain_target':tissue.target_transition.tolist(),'domain_observed':tissue.observed_transition.tolist(),
        'domain_loss':tissue.loss,'components':components,'pool':str(a.pool.resolve()),
        'workflow':'Two native operations: gyrus group generation, then zone-specific cell-type assignment.',
        'limitations':['Domains are generated groups, not guaranteed connected layers.',
            'Zone solvers omit cross-zone edges.'],
        'zones':{str(k):{'converged':v.converged,'loss':v.loss,'iterations':v.n_iter} for k,v in zoned['zone_results'].items()}}
    out.uns['generation_manifest']=json.dumps(audit)
    print('Writing cell.h5ad with domain_true, cell_type_true and slice_id.',flush=True)
    out.write_h5ad(a.out/'cell.h5ad')
    (a.out/'manifest.json').write_text(json.dumps(audit,indent=2)+'\n')
    print(json.dumps(audit),flush=True)

if __name__=='__main__':main()
