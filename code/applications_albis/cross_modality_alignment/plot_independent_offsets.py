"""Render only the dedicated cross-modality independent-offset experiment."""
from pathlib import Path
import argparse
import plot_cross_tech_stair as overlay
import plot_slice5_3d_domain as spatial
import plot_slice5_tech_breakdown as breakdown
from reference_metrics import run as evaluate_reference
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--base', type=Path, default=Path(__file__).resolve().parents[3]/'data/figure_4/cross_modality_alignment/independent_offsets')
base = parser.parse_args().base
for module in (overlay,spatial,breakdown):
    module.BASE=str(base)
out=base/'plots';out.mkdir(exist_ok=True)
a,metrics=overlay.load(5)
for key in ('technology','domain_true'):
    overlay.plot_overlay_2d(a,key,True,str(out),5)
# Legacy pooled/per-modality fits are retained in metrics.json for provenance.
# Figure 4C evaluation uses one reference-only transform for every modality.
evaluate_reference(base/'STAIR/cross_tech/slice_5/adata_results/Sim_CrossTech_STAIR_slice_5.h5ad',
                   out, reference='bin16um', slice_id=5)
spatial.plot(a,5,str(out),reference='bin16um',common_limits=True)
breakdown.plot(a,str(out),5)
