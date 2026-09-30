"""Figure 5D: scCube's native spot memberships for Figure 5A slice 5 (stored slice_id 4).

The saved Figure 5A spot.h5ad has scCube's spot centres, compositions and cell
counts, but not which cells went into each spot. This calls scCube's own
generate_spot_data_random exactly as generate.py did (Visium, n_cell=10, scCube
grid coordinates; function shared with Supp. Fig. 6), checks that the result
reproduces the saved Figure 5A spots, and writes the memberships.
Run with comparison_methods/env/sccube/bin/python.
"""
import argparse
import sys
from pathlib import Path
import json

import anndata as ad
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "overview/figure_supp_6"))
from plot_figure_supp_6 import BASE, SLICE_ID, native_spots  # noqa: E402

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--out", type=Path, required=True)
args = parser.parse_args()
args.out.mkdir(parents=True, exist_ok=True)

from scCube.sccube import scCube
from scCube.utils import calculate_spot_prop

cfg = json.loads((BASE / "settings.json").read_text())
scale = cfg["extent_um"] / cfg["sccube_grid_size"]
n_cell = cfg["sccube_cells_per_spot"]
a = ad.read_h5ad(BASE / "sccube/cell.h5ad")
a = a[a.obs.slice_id.to_numpy().astype(int) == SLICE_ID].copy()
ids = a.obs_names.to_numpy()
native_xy = np.asarray(a.obsm["spatial_3d_native"])[:, :2]
meta = pd.DataFrame({"Cell": ids, "Cell_type": a.obs.cell_type_true.astype(str).to_numpy(),
                     "point_x": native_xy[:, 0], "point_y": native_xy[:, 1]})
X = a.X.toarray() if hasattr(a.X, "toarray") else np.asarray(a.X)
source = pd.DataFrame(X.T, index=a.var_names, columns=ids)
names, xy, prop, n, membership = native_spots(scCube(), calculate_spot_prop, source, meta, n_cell, scale)

b = ad.read_h5ad(BASE / "sccube/spot.h5ad", backed="r")
sel = b.obs.slice_id.to_numpy().astype(int) == SLICE_ID
saved_xy, saved_prop = np.asarray(b.obsm["spatial"])[sel], np.asarray(b.obsm["sccube_spot_prop"])[sel]
saved_n = b.obs.n_source_cells.to_numpy()[sel]
order = lambda p: np.lexsort((p[:, 1], p[:, 0]))
i, j = order(xy), order(saved_xy)
assert len(xy) == len(saved_xy) and np.allclose(xy[i], saved_xy[j]) and np.allclose(prop[i], saved_prop[j])
assert np.array_equal(n[i], saved_n[j])
print(f"n_cell={n_cell}: {len(xy)} spots reproduce the saved Figure 5A slice {SLICE_ID} spots", flush=True)

pd.DataFrame(dict(spot=names, x_um=xy[:, 0], y_um=xy[:, 1], n_cells=n)).to_csv(args.out / "sccube_spots.csv", index=False)
membership.to_csv(args.out / "sccube_membership.csv.gz", index=False)
