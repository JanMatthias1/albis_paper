"""Fail before alignment if independently generated modality shifts coincide."""
from pathlib import Path
import json
import h5py
import numpy as np
from anndata._io.specs import read_elem
ROOT=Path(__file__).resolve().parents[3]/'data/figure_4/cross_modality_alignment/independent_offsets'
rows=[]
for tech,mod in [('bin16um','bin'),('spot','spot'),('cell','cell')]:
    path=ROOT/f'data/{tech}/simulation_{mod}_z_qc.h5ad'
    with h5py.File(path) as f:
        obs=read_elem(f['obs'])
        true=np.asarray(read_elem(f['obsm/spatial']))[:,:2]
        shifted=np.asarray(read_elem(f['obsm/spatial_unaligned']))[:,:2]
        params=read_elem(f['uns/sim_params'])
    assert not params['sync_unaligned_seed'], path
    for sid in sorted(obs.slice_id.astype(str).unique(),key=float):
        mask=(obs.slice_id.astype(str)==sid).to_numpy()
        x,y=true[mask],shifted[mask]
        xc,yc=x-x.mean(0),y-y.mean(0)
        u,_,vt=np.linalg.svd(xc.T@yc);r=u@vt;t=y.mean(0)-x.mean(0)@r
        residual=float(np.sqrt(np.mean(np.sum((x@r+t-y)**2,axis=1))))
        assert np.linalg.det(r)>0 and residual<.01, (tech,sid,residual)
        rows.append(dict(technology=tech,slice_id=sid,n_obs=int(mask.sum()),
                         angle_deg=float(np.degrees(np.arctan2(r[0,1],r[0,0]))),
                         translation_um=t.tolist(),rigid_fit_rmse_um=residual))
for sid in sorted({r['slice_id'] for r in rows},key=float):
    group=[r for r in rows if r['slice_id']==sid]
    assert len(group)==3
    for i in range(3):
        for j in range(i):
            assert not np.allclose(group[i]['translation_um'],group[j]['translation_um'],atol=.01)
            assert not np.isclose(group[i]['angle_deg'],group[j]['angle_deg'],atol=.01)
(ROOT/'perturbation_validation.json').write_text(json.dumps(rows,indent=2))
print('Validated distinct rotations and translations across modalities in all ten sections.')
