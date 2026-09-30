"""Native circular SPIDER aggregation of saved 600k cells; reference is source Splatter pool."""
import argparse,json,sys
from pathlib import Path
import anndata as ad,numpy as np,pandas as pd
from scipy.io import mmread
from scipy.sparse import csr_matrix
from spider.sim_expr import get_sim_spot_level_expr
ROOT=Path(__file__).resolve().parents[5]
sys.path.insert(0,str(ROOT/'sim_paper/code/comparison_methods/overview'))
from generate import capture_data,TYPES

def main():
 p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
 a.out.mkdir(parents=True,exist_ok=False);cfg=json.loads((a.run/'settings.json').read_text())
 cells=ad.read_h5ad(a.run/'spider/cell.h5ad');assert cells.n_obs==600000
 sections=[]
 for sid in range(10):
  sub=cells[cells.obs.slice_id.astype(int)==sid].copy();codes=pd.Categorical(sub.obs.cell_type_true,categories=TYPES).codes
  expr,xy,membership=get_sim_spot_level_expr(Num_sample=sub.n_obs,spot_diameter=2*cfg['spot_radius_um'],image_width=cfg['extent_um'],image_height=cfg['extent_um'],celltype_assignment=codes,cell_spatial=sub.obsm['spatial'],sim_cell_expr=sub,gap=cfg['spot_spacing_um'],coord_type='generic',spot_generate_type='circle',cell_coord_type='generic')[:3]
  counts=np.asarray(membership@np.eye(8,dtype=int)[codes]);sec=capture_data(expr,cells.var_names,xy,counts,sid,cfg)
  keep=(counts.sum(1)>0)&(np.asarray(sec.X.sum(1)).ravel()>0);sec=sec[keep].copy();counts=counts[keep]
  sec.obsm['cell_type_frac_true']=counts/counts.sum(1,keepdims=True);sec.uns['cell_type_frac_true_columns']=TYPES
  assert np.allclose(sec.obsm['cell_type_frac_true'].sum(1),1)
  sections.append(sec);print('Aggregated slice',sid,sec.n_obs,flush=True)
 spots=ad.concat(sections,index_unique='-slice',merge='same');spots.uns['cell_type_frac_true_columns']=TYPES
 pool=a.run/'splatter';meta=pd.read_csv(pool/'cells.tsv',sep='\t').set_index('cell_id');genes=(pool/'genes.tsv').read_text().splitlines()
 ref=ad.AnnData(mmread(pool/'counts.mtx').T.tocsr(),obs=meta,var=pd.DataFrame(index=genes))
 assert ref.var_names.equals(spots.var_names) and set(ref.obs.cell_type_true)==set(TYPES)
 for name,x in [('cell',ref),('spot',spots)]:
  x.X=csr_matrix(x.X,dtype=np.float64);assert np.isfinite(x.X.data).all() and (x.X.data>=0).all() and np.equal(x.X.data,np.rint(x.X.data)).all()
  x.write_h5ad(a.out/f'{name}.h5ad',compression='gzip')
 report={'source':str(a.run),'seed':cfg['seed'],'genes':ref.n_vars,'underlying_cells':600000,'reference_cells':ref.n_obs,'query_spots':spots.n_obs,'reference':'Original synthetic Splatter source pool; profiles contributed to query via sampling. Not independent/held-out.','truth':'Contributing-cell fractions','expression':'Native integer counts, no rounding, normalization or depth adjustment','capture':'Native SPIDER circular spots, diameter55um, spacing100um; all ten slices','limitations':'Matched gene count does not match marker strength or depth; ALBIS molecule-fraction truth differs.'}
 (a.out/'input_audit.json').write_text(json.dumps(report,indent=2)+'\n')
if __name__=='__main__':main()
