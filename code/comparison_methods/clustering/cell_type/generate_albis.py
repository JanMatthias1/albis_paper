"""Native ALBIS cell output with Figure 5A parameters and batch_sigma=1.5."""
import argparse,json,sys,hashlib
from pathlib import Path
ROOT=Path(__file__).resolve().parents[5]
sys.path.insert(0,str(ROOT/'albis'))
import albis
import numpy as np

def main():
    p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);a=p.parse_args()
    cfg=json.loads((a.run/'settings.json').read_text())
    dest=a.run/'albis';dest.mkdir(exist_ok=False)
    assert cfg['n_cells']==600000 and cfg['batch_sigma']==1.5
    assert cfg['theta']==2 and cfg['theta_jitter']==.6
    x=albis.simulate_3d_molecule_sphere_multires(**cfg)['adata_cell_sectioned']['Z']
    assert x.n_obs==600000 and x.n_vars==556 and x.obs.cell_type_true.nunique()==8
    assert x.obs.slice_id.nunique()==10
    x.write_h5ad(dest/'cell.h5ad')
    report={'method':'albis','parameters':cfg,'n_cells':x.n_obs,'n_genes':x.n_vars,
        'counts_positive':int((np.asarray(x.X.sum(axis=1)).ravel()>0).sum()),
        'albis_source_sha256':hashlib.sha256((ROOT/'albis/albis/simulation_sphere.py').read_bytes()).hexdigest()}
    (dest/'manifest.json').write_text(json.dumps(report,indent=2)+'\n')
    print('ALBIS generation complete',flush=True)
if __name__=='__main__':main()
