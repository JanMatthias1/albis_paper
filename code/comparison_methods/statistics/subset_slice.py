"""Write one slice of a Figure 5A output to its own h5ad (no other change).

count_distribution.py restricts only --input to a slice; --compare-input is
used as given, so both sides are pre-subset here to the same slice.
"""
import argparse
from pathlib import Path

import anndata as ad

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--input", type=Path, required=True)
parser.add_argument("--output", type=Path, required=True)
parser.add_argument("--slice-id", type=int, default=4, help="4 = 5th of 10 slices from the bottom")
args = parser.parse_args()
a = ad.read_h5ad(args.input, backed="r")
keep = a.obs.slice_id.astype(int).to_numpy() == args.slice_id
if not keep.any():
    raise SystemExit(f"slice_id {args.slice_id} not found in {args.input}")
sub = a[keep].to_memory()
a.file.close()
args.output.parent.mkdir(parents=True, exist_ok=True)
sub.write_h5ad(args.output)
print(f"{args.input} -> {args.output}: slice_id {args.slice_id}, {sub.n_obs} x {sub.n_vars}")
