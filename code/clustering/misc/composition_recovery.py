#!/usr/bin/env python
"""
Figure 3B, fairer scoring: how well do Leiden clusters recover the ground-truth
structure of an *aggregated* observation, when that observation is a mixture of
cell types rather than a single labelled cell?

Motivation (PI feedback 2026-09-06): the existing metric is hard-label ARI of
resolution-matched Leiden clusters vs `cell_type_true` / `domain_true`. For a
single `cell` that is well posed. For a `bin` or (especially) a `spot`, the unit
is a mixture -- `cell_type_true` there is only the argmax of a ~5-6-way
composition, so hard-label ARI penalises the method for something it structurally
cannot do. The generator already stores the soft truth:
`obsm['cell_type_frac_true']` (n x 8) and `obsm['domain_frac_true']` (n x 6) for
bin/spot; `cell` gets a one-hot synthesised from its hard labels so every
modality is on one table.

Reads each modality's resolution-matched clustering written by
02_leiden_resolution_sweep.py:
    data/figure_3/pca_harmony_single_cell/<mod>/ari_recovery_qc/
        simulation_<mod>_z_ari_recovery.h5ad
        -> obs['leiden_cell_type_true'], obs['leiden_domain_true']
           (Leiden at the resolution that hits k=8 / k=6)
        -> obsm['X_pca_harmony'], obsm['{cell_type,domain}_frac_true']

For each modality x {cell_type, domain} it reports:
  hard-label   : ARI, V-measure, homogeneity, completeness  (vs argmax truth)
  leak         : slice_id_leakage_ari -- ARI(predicted clusters, slice_id).
                 High leak means a cluster's apparent "recovery" is actually
                 just rediscovering slice boundaries, not real structure.

Scores one clustered h5ad per invocation (--h5ad), writing
<out-dir>/composition_recovery_<tag>.json -- every sweep script (batch_sigma
slides, banksy_batch_compare, etc.) calls this per config/point.
slice_id_leakage_ari is read by summary_banksy_lambda_kgeom.py's table; no
current script visualizes it (plot_banksy_batch_compare_3mod.py's dagger
annotation was removed 2026-09-16, see below).

2026-09-16: dropped the other, no-args mode this file used to also support
(score all 4 canonical modalities at once, write a combined
composition_recovery_summary.csv + composition_recovery.png to
data/figure_3/composition_recovery/) -- nothing has called it that way since
the "Proposed Fig 3B" reframing (2026-09-06) it was built for was abandoned;
that output is archived at data/figure_3/_archive_20260914/composition_recovery/,
and plot_composition_recovery.py (the only reader) moved to misc/ as
superseded.

2026-09-16: also dropped the composition/local_knn/oracle/context metrics
(mean JSD/L1/Pearson r, composition R2, kNN concordance vs chance, KMeans-
on-truth oracle ceiling, effective-categories-per-obs) -- computed for the
same abandoned "Proposed Fig 3B" reframing above, and unread by anything
active (only misc/summary_banksy_lambda_kgeom.py, a dormant summary for an
already-decided BANKSY lambda/k_geom sweep, read a few of them -- it'll break
loudly if rerun, which is fine, that sweep's already been superseded).
X_pca_harmony is no longer read either, since it only fed the kNN check.

Expected environment:
    source code/count_distribution/_env.sh   (albis-tutorial: scanpy/sklearn/scipy)

Example:
    python sim_paper/code/clustering/composition_recovery.py \\
        --h5ad data/figure_3/cellbin_batch_sigma_slide/spot/bs0.3/ari/simulation_spot_z_ari_recovery.h5ad \\
        --tag spot_bs0.3 --out-dir data/figure_3/cellbin_batch_sigma_slide/scores
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import h5py
import pandas as pd
from sklearn.metrics import (
    adjusted_rand_score,
    homogeneity_completeness_v_measure,
)

SCRIPT_DIR = Path(__file__).resolve().parent
SIM_PAPER_DIR = SCRIPT_DIR.parents[1]
LEVELS = {  # level -> (frac obsm key, hard-label obs col, k)
    "cell_type": ("cell_type_frac_true", "cell_type_true", 8),
    "domain": ("domain_frac_true", "domain_true", 6),
}


# --------------------------------------------------------------------------- io
def _read_obs_col(f: h5py.File, col: str):
    """Return an obs column as an object/np array, handling anndata categoricals."""
    node = f["obs"][col]
    if isinstance(node, h5py.Group):  # categorical: categories[] + codes[]
        cats = node["categories"][:]
        cats = np.array([c.decode() if isinstance(c, bytes) else c for c in cats], dtype=object)
        codes = node["codes"][:]
        out = np.empty(codes.shape, dtype=object)
        good = codes >= 0
        out[good] = cats[codes[good]]
        out[~good] = None
        return out
    arr = node[:]
    if arr.dtype.kind == "S":
        arr = np.array([x.decode() for x in arr], dtype=object)
    return arr


def load_parts(path: Path, need_obsm: list[str]):
    """Pull only the obs columns / obsm arrays we need (files are up to 12 GB)."""
    obs, obsm = {}, {}
    with h5py.File(path, "r") as f:
        for c in ("leiden_cell_type_true", "leiden_domain_true",
                  "cell_type_true", "domain_true", "slice_id"):
            if c in f["obs"]:
                obs[c] = _read_obs_col(f, c)
        for k in need_obsm:
            if k in f.get("obsm", {}):
                obsm[k] = f["obsm"][k][:]
    return obs, obsm


# ---------------------------------------------------------------------- metrics
def one_hot(labels: np.ndarray, categories: list[str]) -> np.ndarray:
    idx = {c: i for i, c in enumerate(categories)}
    out = np.zeros((labels.shape[0], len(categories)), dtype=np.float64)
    for i, lab in enumerate(labels):
        j = idx.get(lab)
        if j is not None:
            out[i, j] = 1.0
    return out


def codes_from(labels: np.ndarray) -> np.ndarray:
    _, inv = np.unique(labels.astype(str), return_inverse=True)
    return inv.astype(np.int64)


# ------------------------------------------------------------------------- main
def score_path(path: Path, name: str) -> dict:
    """Score one clustered h5ad (needs leiden_{cell_type,domain}_true, the
    *_frac_true obsm or the hard labels, and slice_id for the leakage check)."""
    if not Path(path).is_file():
        raise SystemExit(f"missing: {path}")
    need = [LEVELS[l][0] for l in LEVELS]
    obs, obsm = load_parts(Path(path), need)
    slice_codes = codes_from(obs["slice_id"]) if "slice_id" in obs else None
    n_obs = len(next(iter(obs.values())))
    out = {"modality": name, "n_obs": n_obs, "levels": {}}

    for level, (frac_key, hard_col, k) in LEVELS.items():
        pred = codes_from(obs[f"leiden_{level}_true"])
        if frac_key in obsm:
            frac = np.asarray(obsm[frac_key], dtype=np.float64)
        else:  # cell: synthesise one-hot from the hard label
            cats = sorted(set(obs[hard_col]) - {None})
            frac = one_hot(obs[hard_col], cats)

        row_sum = frac.sum(1)
        keep = row_sum > 1e-9
        n_drop = int((~keep).sum())
        frac = frac[keep] / row_sum[keep, None]
        pred_l = pred[keep]
        hard_truth = frac.argmax(1)

        ari = adjusted_rand_score(hard_truth, pred_l)
        hom, com, vme = homogeneity_completeness_v_measure(hard_truth, pred_l)

        leak = (float(adjusted_rand_score(slice_codes[keep], pred_l))
                if slice_codes is not None else None)

        out["levels"][level] = {
            "k": k,
            "n_dropped_zero_frac": n_drop,
            "hard_label": {
                "ari": float(ari),
                "v_measure": float(vme),
                "homogeneity": float(hom),
                "completeness": float(com),
            },
            "slice_id_leakage_ari": leak,
        }
        lk = "  (leak %.2f)" % leak if leak is not None else ""
        print(f"  [{name}/{level}] ARI={ari:.3f} V={vme:.3f}{lk}")
    return out


def to_rows(results: list[dict]) -> pd.DataFrame:
    rows = []
    for r in results:
        for level, d in r["levels"].items():
            rows.append({
                "modality": r["modality"], "n_obs": r["n_obs"], "level": level,
                "ari": d["hard_label"]["ari"],
                "v_measure": d["hard_label"]["v_measure"],
                "slice_id_leakage_ari": d.get("slice_id_leakage_ari"),
            })
    return pd.DataFrame(rows)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--h5ad", type=Path, required=True,
                    help="score this clustered h5ad; writes composition_recovery_<tag>.json")
    ap.add_argument("--tag", default=None, help="name for the run (default: file stem)")
    ap.add_argument("--out-dir", type=Path, required=True)
    args = ap.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    tag = args.tag or Path(args.h5ad).stem
    print(f"[score] {tag}  <-  {args.h5ad}")
    r = score_path(args.h5ad, tag)
    with open(args.out_dir / f"composition_recovery_{tag}.json", "w") as f:
        json.dump(r, f, indent=2)
    print(to_rows([r]).to_string(index=False))


if __name__ == "__main__":
    main()
