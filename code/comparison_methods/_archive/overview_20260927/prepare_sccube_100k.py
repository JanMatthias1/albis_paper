"""Prepare an isolated 100k-cell scCube rerun using the original expression pool."""
import argparse
import json
from pathlib import Path

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--out', type=Path, required=True)
args = parser.parse_args()
source = Path(__file__).resolve().parents[3] / 'data/comparison_methods/overview_full_35831924'
cfg = json.loads((source/'settings.json').read_text())
cfg['n_cells'] = 100000
cfg['unavailable_methods'] = {m: 'Not rerun at 100,000 cells' for m in ['albis','spider']}
cfg['plot_marker_sizes'] = {'cell': .3, 'bin': .6, 'spot': 3}
args.out.mkdir(parents=True, exist_ok=False)
(args.out/'settings.json').write_text(json.dumps(cfg,indent=2)+'\n')
(args.out/'splatter').symlink_to(source/'splatter', target_is_directory=True)
(args.out/'input_provenance.json').write_text(json.dumps(dict(
    source_overview=str(source), requested_cells=100000,
    expression_training_pool='Original 10000-cell, 2000-gene Splatter pool',
    changed_generation_settings={'n_cells':100000},
    note='Same 600 um extent, 10 sections, 3 cells/bin and 10 cells/spot targets; native geometry.'
),indent=2)+'\n')
