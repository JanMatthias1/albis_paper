#!/usr/bin/env python
"""
Extract raw Visium HD Human Breast Cancer downloads into per-sample folders
ready for unzip_into_adata() / sc.read_visium().

Adapted from the commented-out extraction block in Vani/JM's
Visium_HD_human_colon_{16,8}_um scripts
(multi-sample-alignment-benchmark/code/01_preprocessing/), split out into its
own script per JM's request rather than left as a commented block.

Expects 10x downloads sitting flat in base_dir, one archive per sample:
    {prefix}_binned_outputs.tar.gz   (required; contains binned_outputs/square_XXXum/...)
    {prefix}_spatial.tar.gz          (optional; older/split Visium layout)
    {prefix}_cloupe_*.cloupe         (optional)
    {prefix}_metrics_summary.csv     (optional)
    {prefix}_molecule_info.h5        (optional)

Unlike the colon downloads (single bin size, spatial/ shipped separately),
current 10x Visium HD bundles ship every square_XXXum folder already
self-contained with its own spatial/ subfolder inside binned_outputs.tar.gz.
The {prefix}_spatial.tar.gz handling is kept only as a defensive fallback in
case a future/older-style download splits it out again.

Run from anywhere, e.g.:
    conda activate /dcs04/hicks/data/multi-sample-alignment-benchmark/envs/preprocessing_JM
    python extract_breast_visium_hd.py
"""

import argparse
import gzip
import shutil
import tarfile
from pathlib import Path

DEFAULT_BASE_DIR = Path("/dcs04/hicks/data/Jan/sim_project/albis_paper/data/real_data/breast_cancer")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Extract raw Visium HD breast cancer downloads.")
    parser.add_argument("--base-dir", type=Path, default=DEFAULT_BASE_DIR)
    parser.add_argument(
        "--pattern",
        default="*_binned_outputs.tar.gz",
        help="Glob (within --base-dir) matching each sample's binned_outputs archive.",
    )
    return parser.parse_args()


def extract_sample(base_dir: Path, binned_tar: Path) -> None:
    prefix = binned_tar.name.replace("_binned_outputs.tar.gz", "")
    sample_dir = base_dir / prefix
    spatial_tar = base_dir / f"{prefix}_spatial.tar.gz"
    cloupe_files = list(base_dir.glob(f"{prefix}_cloupe_*.cloupe"))
    metrics_file = base_dir / f"{prefix}_metrics_summary.csv"
    molecule_file = base_dir / f"{prefix}_molecule_info.h5"

    print(f"Processing {prefix} ...")
    sample_dir.mkdir(parents=True, exist_ok=True)

    # 1. Extract binned_outputs.tar.gz -> sample_dir/binned_outputs/square_XXXum/...
    print(f"  Extracting {binned_tar.name} -> {sample_dir}")
    with tarfile.open(binned_tar, "r:gz") as tar:
        tar.extractall(sample_dir)

    # 2. Extract a separate spatial.tar.gz, if this download shipped one
    if spatial_tar.exists():
        print(f"  Extracting {spatial_tar.name} -> {sample_dir / 'spatial'}")
        with tarfile.open(spatial_tar, "r:gz") as tar:
            tar.extractall(sample_dir / "spatial")

    # 3. Move auxiliary files alongside the sample, if present
    for f in [*cloupe_files, metrics_file, molecule_file]:
        if f.exists():
            shutil.move(str(f), sample_dir / f.name)

    # 4. Decompress any .gz files inside every spatial/ folder under sample_dir
    #    (covers both the per-bin-size spatial/ dirs bundled in binned_outputs,
    #    and the top-level spatial/ dir from step 2 if it was used)
    for spatial_dir in sample_dir.glob("**/spatial"):
        for gz_file in spatial_dir.glob("*.gz"):
            out_file = gz_file.with_suffix("")
            with gzip.open(gz_file, "rb") as f_in, open(out_file, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)
            gz_file.unlink()

    print(f"Finished {prefix}\n")


def main() -> None:
    args = parse_args()
    binned_tars = sorted(args.base_dir.glob(args.pattern))
    if not binned_tars:
        raise SystemExit(f"No archives matching {args.pattern!r} found under {args.base_dir}")

    for binned_tar in binned_tars:
        extract_sample(args.base_dir, binned_tar)


if __name__ == "__main__":
    main()
