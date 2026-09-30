"""Render the Albis RCTD plot layouts for saved comparator pilot estimates."""
import argparse
import importlib.util
import json
from pathlib import Path
import sys
import anndata as ad
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('run', type=Path)
parser.add_argument('--slice-id', default='5')
args = parser.parse_args()
results = args.run/'results'
summary = json.loads((results/'metrics_summary.json').read_text())
audit = json.loads((args.run/'inputs/input_audit.json').read_text())
source = Path(__file__).resolve().parents[3]/'clustering/spatial_deconvolution/plot_rctd_results.py'
sys.argv = [str(source), str(results)]
spec = importlib.util.spec_from_file_location('albis_rctd_plots', source)
plots = importlib.util.module_from_spec(spec)
spec.loader.exec_module(plots)
spot = ad.read_h5ad(summary['spot_h5ad'])
wide = pd.read_csv(results/'estimated_fractions_wide.csv')
assert wide.spot_id.is_unique
indices = wide.spot_id.str.removeprefix('spot_').astype(int).to_numpy()-1
assert ((indices>=0)&(indices<spot.n_obs)).all()
assert len(indices)==summary['n_spots_scored']
truth = pd.DataFrame(np.asarray(spot.obsm['cell_type_frac_true'])[indices],columns=plots.CELL_TYPES)
est = wide[plots.CELL_TYPES].copy()
assert np.isfinite(est.to_numpy()).all() and np.isfinite(truth.to_numpy()).all()
assert np.isclose(np.sqrt(np.mean((est.to_numpy()-truth.to_numpy())**2)),summary['overall_rmse'],atol=0.0001)
spatial = pd.DataFrame(np.asarray(spot.obsm['spatial'])[indices,:2],columns=['x','y'])
meta = spot.obs.iloc[indices][['slice_id']].reset_index(drop=True)
selected = meta.slice_id.astype(str).to_numpy()==args.slice_id
if selected.sum()<2: raise ValueError('Requested slice has fewer than two scored spots')
xy = spatial.loc[selected].to_numpy()
pitch = float(np.median(cKDTree(xy).query(xy,k=2)[0][:,1]))
assert pitch>0
extent = json.loads(spot.uns['overview_settings_json'])['extent_um']
plots.SLICE_ID=args.slice_id
plots.ZOOM_XMIN=plots.ZOOM_YMIN=-pitch/2
plots.ZOOM_XMAX=plots.ZOOM_YMAX=extent+pitch/2
# Existing plot multiplies supplied radius by 1.9. Here radius is only a glyph size.
radius = .32*pitch/1.9
method = {'spider':'SPIDER','sccube':'scCube'}[audit['method']]
original = plots.plt.Figure.savefig
def save_with_caption(fig, path, *a, **kw):
    fig.text(.5,.002,f'{method} native-data pilot | truth: contributing-cell fractions | circles: display markers',
             ha='center',va='bottom',fontsize=8)
    original(fig,path,*a,**kw)
    original(fig,Path(path).with_suffix('.pdf'),**{k:v for k,v in kw.items() if k!='dpi'})
plots.plt.Figure.savefig=save_with_caption
for types,pairs,name in [
    (plots.CELL_TYPES,1,'rctd_spatial_zoom.png'),
    (plots.CELL_TYPES,2,'rctd_spatial_zoom_4x4.png'),
    (['type1','type4'],1,'rctd_spatial_zoom_main.png'),
    (plots.CELL_TYPES[:4],2,'rctd_spatial_zoom_4types.png')]:
    plots.plot_spatial_zoom(truth,est,spatial,meta,radius,cell_types=types,pairs_per_row=pairs,out_name=name)
plots.plot_error_boxplot(truth,est)
(results/'plot_provenance.json').write_text(json.dumps(dict(
    method=audit['method'],spot_h5ad=summary['spot_h5ad'],slice_id=args.slice_id,
    n_scored=len(indices),n_spatial_panel=int(selected.sum()),
    glyph_radius=.32*pitch,spatial_markers='display circles, not physical capture footprints',
    truth='contributing-cell fractions',
    error_plot='Albis layout; violin values capped at pooled 99th percentile, y-axis cropped; metrics use full errors'
),indent=2)+'\n')
