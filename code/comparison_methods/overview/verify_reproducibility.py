"""Compare a fresh Figure 5A run with earlier outputs made from the same settings.

Checks only; nothing is modified. Writes OUT/reproducibility_check.json:
- Splatter 556-gene pool vs overview_native556_20260926/splatter (same seed/settings)
- ALBIS cell/bin/spot vs Figure 4 strong_domain_mix_shift3x (same generator/seed)
scCube/SPIDER are not compared: they have no earlier 556-gene 600k run, and
SPIDER's annealer is unseedable (see README).
"""
import argparse
import hashlib
import json
from pathlib import Path

import anndata as ad
import numpy as np
from scipy import sparse

DATA = Path(__file__).resolve().parents[3] / "data"
POOL_REF = DATA / "comparison_methods/overview_native556_20260926/splatter"
ALBIS_REF = DATA / "figure_4/cross_modality_alignment/strong_domain_mix_shift3x/data"
ALBIS_FILES = {"cell": "cell/simulation_cell_z.h5ad", "bin": "bin16um/simulation_bin_z.h5ad",
               "spot": "spot/simulation_spot_z.h5ad"}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def same_matrix(a, b):
    a = a.tocsr() if sparse.issparse(a) else sparse.csr_matrix(a)
    b = b.tocsr() if sparse.issparse(b) else sparse.csr_matrix(b)
    return a.shape == b.shape and (a != b).nnz == 0


def compare_adata(new, ref):
    a, b = ad.read_h5ad(new), ad.read_h5ad(ref)
    return dict(
        n_obs=[a.n_obs, b.n_obs], n_vars=[a.n_vars, b.n_vars],
        X_identical=same_matrix(a.X, b.X),
        obs_names_identical=bool(a.obs_names.equals(b.obs_names)),
        var_names_identical=bool(a.var_names.equals(b.var_names)),
        obs_identical=bool(a.obs.equals(b.obs)),
        obsm_identical={k: bool(k in b.obsm and np.array_equal(np.asarray(a.obsm[k]), np.asarray(b.obsm[k])))
                        for k in a.obsm})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", type=Path, required=True)
    run = parser.parse_args().run.resolve()
    report = {"splatter_pool": {f: sha(run / "splatter" / f) == sha(POOL_REF / f)
                                for f in ["counts.mtx", "genes.tsv", "cells.tsv"]},
              "splatter_reference": str(POOL_REF), "albis_reference": str(ALBIS_REF), "albis": {}}
    for modality, rel in ALBIS_FILES.items():
        report["albis"][modality] = compare_adata(run / "albis" / f"{modality}.h5ad", ALBIS_REF / rel)
        print(modality, report["albis"][modality], flush=True)
    (run / "reproducibility_check.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report["splatter_pool"]), flush=True)


if __name__ == "__main__":
    main()
