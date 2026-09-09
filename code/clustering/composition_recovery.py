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
ari_vs_ground_truth.py:
    data/figure_3/pca_harmony_single_cell/<mod>/ari_recovery_qc/
        simulation_<mod>_z_ari_recovery.h5ad
        -> obs['leiden_cell_type_true'], obs['leiden_domain_true']
           (Leiden at the resolution that hits k=8 / k=6)
        -> obsm['X_pca_harmony'], obsm['{cell_type,domain}_frac_true']

For each modality x {cell_type, domain} it reports:
  hard-label   : ARI, V-measure, homogeneity, completeness  (vs argmax truth)
  composition  : mean Jensen-Shannon divergence and mean L1 between each obs's
                 true composition and its cluster's mean composition;
                 mean per-category Pearson r (cluster-mean fraction vs true
                 fraction, across obs); composition variance explained
                 R2 = 1 - SS_within / SS_total on the fraction vectors.
  local        : kNN label concordance in X_pca_harmony (fraction of each obs's
                 15 nearest neighbours sharing its argmax label), vs the
                 chance baseline sum(p_k^2). Subsampled to 50k for speed.
  oracle       : cluster the TRUE fraction vectors with KMeans(k) -> the best any
                 expression method could do given the aggregation. Reports
                 ARI_method / ARI_oracle and R2_method / R2_oracle.
  context      : effective number of categories per obs, exp(Shannon entropy of
                 the fraction vector); fraction of obs that are >=80% one
                 category.

Outputs to data/figure_3/composition_recovery/:
    composition_recovery_<mod>.json     full per-modality metrics
    composition_recovery_summary.csv    flat table, all modalities x levels
    composition_recovery.png            grouped-bar panel (replaces the stale
                                        ari_recovery_summary/ for this framing)

Expected environment:
    source code/count_distribution/_env.sh   (albis-tutorial: scanpy/sklearn/scipy)

Example:
    python sim_paper/code/clustering/composition_recovery.py                 # all 4
    python sim_paper/code/clustering/composition_recovery.py --modalities spot
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
from sklearn.cluster import KMeans
from sklearn.neighbors import NearestNeighbors

SCRIPT_DIR = Path(__file__).resolve().parent
SIM_PAPER_DIR = SCRIPT_DIR.parents[1]
FIG3 = SIM_PAPER_DIR / "data" / "figure_3" / "pca_harmony_single_cell"
OUT_DIR = SIM_PAPER_DIR / "data" / "figure_3" / "composition_recovery"

# modality -> h5ad basename stem (bin8 and bin16 share the "bin" stem)
MODS = {
    "cell": "simulation_cell_z_ari_recovery.h5ad",
    "bin": "simulation_bin_z_ari_recovery.h5ad",
    "bin16um": "simulation_bin_z_ari_recovery.h5ad",
    "spot": "simulation_spot_z_ari_recovery.h5ad",
}
LEVELS = {  # level -> (frac obsm key, hard-label obs col, k)
    "cell_type": ("cell_type_frac_true", "cell_type_true", 8),
    "domain": ("domain_frac_true", "domain_true", 6),
}
KNN_SUBSAMPLE = 50_000
KNN_K = 15
SEED = 0


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


def cluster_mean_composition(frac: np.ndarray, pred: np.ndarray) -> np.ndarray:
    K = int(pred.max()) + 1
    D = frac.shape[1]
    counts = np.bincount(pred, minlength=K).astype(np.float64)
    means = np.zeros((K, D))
    for d in range(D):
        means[:, d] = np.bincount(pred, weights=frac[:, d], minlength=K)
    means /= counts[:, None].clip(min=1.0)
    return means[pred]  # (n, D) predicted composition per obs


def jsd_rows(P: np.ndarray, Q: np.ndarray) -> np.ndarray:
    """Jensen-Shannon divergence per row, base 2, in [0, 1]."""
    M = 0.5 * (P + Q)
    with np.errstate(divide="ignore", invalid="ignore"):
        kl_pm = np.where(P > 0, P * np.log2(P / M), 0.0).sum(1)
        kl_qm = np.where(Q > 0, Q * np.log2(Q / M), 0.0).sum(1)
    return 0.5 * kl_pm + 0.5 * kl_qm


def comp_r2(frac: np.ndarray, pred_comp: np.ndarray) -> float:
    ss_within = float(((frac - pred_comp) ** 2).sum())
    ss_total = float(((frac - frac.mean(0, keepdims=True)) ** 2).sum())
    return 1.0 - ss_within / ss_total if ss_total > 0 else float("nan")


def per_category_pearson(frac: np.ndarray, pred_comp: np.ndarray) -> list[float]:
    rs = []
    for d in range(frac.shape[1]):
        a, b = frac[:, d], pred_comp[:, d]
        if a.std() < 1e-12 or b.std() < 1e-12:
            rs.append(float("nan"))
        else:
            rs.append(float(np.corrcoef(a, b)[0, 1]))
    return rs


def knn_concordance(emb: np.ndarray, hard_codes: np.ndarray, rng: np.random.Generator):
    n = emb.shape[0]
    if n > KNN_SUBSAMPLE:
        sel = rng.choice(n, KNN_SUBSAMPLE, replace=False)
        emb, hard_codes = emb[sel], hard_codes[sel]
    nn = NearestNeighbors(n_neighbors=KNN_K + 1).fit(emb)
    _, idx = nn.kneighbors(emb)
    neigh = hard_codes[idx[:, 1:]]  # drop self
    same = (neigh == hard_codes[:, None]).mean()
    _, cnt = np.unique(hard_codes, return_counts=True)
    p = cnt / cnt.sum()
    return float(same), float((p ** 2).sum())  # observed, chance baseline


def eff_categories(frac: np.ndarray):
    with np.errstate(divide="ignore", invalid="ignore"):
        H = -(np.where(frac > 0, frac * np.log(frac), 0.0)).sum(1)
    eff = np.exp(H)
    return float(np.median(eff)), float(eff.mean()), float((frac.max(1) >= 0.8).mean())


# ------------------------------------------------------------------------- main
def score_path(path: Path, name: str) -> dict:
    """Score one clustered h5ad (needs X_pca_harmony, leiden_{cell_type,domain}_true,
    the *_frac_true obsm or the hard labels, and slice_id for the leakage check)."""
    if not Path(path).is_file():
        raise SystemExit(f"missing: {path}")
    need = ["X_pca_harmony"] + [LEVELS[l][0] for l in LEVELS]
    obs, obsm = load_parts(Path(path), need)
    emb = np.asarray(obsm["X_pca_harmony"], dtype=np.float64)[:, :30]
    slice_codes = codes_from(obs["slice_id"]) if "slice_id" in obs else None
    rng = np.random.default_rng(SEED)
    out = {"modality": name, "n_obs": int(emb.shape[0]), "levels": {}}

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
        emb_l = emb[keep]
        hard_truth = frac.argmax(1)

        ari = adjusted_rand_score(hard_truth, pred_l)
        hom, com, vme = homogeneity_completeness_v_measure(hard_truth, pred_l)

        pred_comp = cluster_mean_composition(frac, pred_l)
        jsd = float(jsd_rows(frac, pred_comp).mean())
        l1 = float(np.abs(frac - pred_comp).sum(1).mean())
        rs = per_category_pearson(frac, pred_comp)
        r2 = comp_r2(frac, pred_comp)

        knn_obs, knn_chance = knn_concordance(emb_l, hard_truth.astype(np.int64), rng)

        oracle = KMeans(n_clusters=k, random_state=SEED, n_init=10).fit_predict(frac)
        ari_oracle = adjusted_rand_score(hard_truth, oracle)
        r2_oracle = comp_r2(frac, cluster_mean_composition(frac, oracle))

        leak = (float(adjusted_rand_score(slice_codes[keep], pred_l))
                if slice_codes is not None else None)

        med_eff, mean_eff, pure80 = eff_categories(frac)

        out["levels"][level] = {
            "k": k,
            "n_dropped_zero_frac": n_drop,
            "hard_label": {
                "ari": float(ari),
                "v_measure": float(vme),
                "homogeneity": float(hom),
                "completeness": float(com),
            },
            "composition": {
                "mean_jsd": jsd,
                "mean_l1": l1,
                "per_category_pearson": [None if np.isnan(x) else round(x, 4) for x in rs],
                "mean_pearson": float(np.nanmean(rs)),
                "r2_variance_explained": float(r2),
            },
            "local_knn": {
                "k": KNN_K,
                "subsampled_to": min(int(keep.sum()), KNN_SUBSAMPLE),
                "concordance": knn_obs,
                "chance_baseline": knn_chance,
                "lift_over_chance": knn_obs - knn_chance,
            },
            "oracle": {
                "ari_oracle": float(ari_oracle),
                "ari_method_over_oracle": (float(ari / ari_oracle) if ari_oracle > 0.01 else None),
                "r2_oracle": float(r2_oracle),
                "r2_method_over_oracle": (float(r2 / r2_oracle) if r2_oracle > 0.01 else None),
            },
            "context": {
                "eff_categories_per_obs_median": med_eff,
                "eff_categories_per_obs_mean": mean_eff,
                "frac_obs_ge80pct_one_category": pure80,
            },
            "slice_id_leakage_ari": leak,
        }
        lk = "  (leak %.2f)" % leak if leak is not None else ""
        print(f"  [{name}/{level}] ARI={ari:.3f} V={vme:.3f} | JSD={jsd:.3f} "
              f"meanR={np.nanmean(rs):.3f} R2={r2:.3f} (oracle {r2_oracle:.3f}) | "
              f"kNN={knn_obs:.3f} (chance {knn_chance:.3f}) | eff#cats={med_eff:.2f}{lk}")
    return out


def score_modality(mod: str) -> dict:
    return score_path(FIG3 / mod / "ari_recovery_qc" / MODS[mod], mod)


def to_rows(results: list[dict]) -> pd.DataFrame:
    rows = []
    for r in results:
        for level, d in r["levels"].items():
            rows.append({
                "modality": r["modality"], "n_obs": r["n_obs"], "level": level,
                "ari": d["hard_label"]["ari"],
                "v_measure": d["hard_label"]["v_measure"],
                "comp_mean_pearson": d["composition"]["mean_pearson"],
                "comp_mean_jsd": d["composition"]["mean_jsd"],
                "comp_r2": d["composition"]["r2_variance_explained"],
                "comp_r2_oracle": d["oracle"]["r2_oracle"],
                "comp_r2_frac_of_oracle": d["oracle"]["r2_method_over_oracle"],
                "ari_oracle": d["oracle"]["ari_oracle"],
                "ari_frac_of_oracle": d["oracle"]["ari_method_over_oracle"],
                "knn_concordance": d["local_knn"]["concordance"],
                "knn_chance": d["local_knn"]["chance_baseline"],
                "eff_cats_per_obs_median": d["context"]["eff_categories_per_obs_median"],
                "frac_obs_ge80pct": d["context"]["frac_obs_ge80pct_one_category"],
                "slice_id_leakage_ari": d.get("slice_id_leakage_ari"),
            })
    return pd.DataFrame(rows)


def make_plot(df: pd.DataFrame, out_png: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    mods = ["cell", "bin", "bin16um", "spot"]
    df = df[df.modality.isin(mods)].copy()
    ct = df[df.level == "cell_type"].set_index("modality").reindex(mods)
    dm = df[df.level == "domain"].set_index("modality").reindex(mods)
    x = np.arange(len(mods))
    w = 0.38
    fig, ax = plt.subplots(2, 3, figsize=(15, 8.5))

    def grouped(a, series, labels, title, ylab, ymax=None):
        for i, (s, lab) in enumerate(zip(series, labels)):
            a.bar(x + (i - (len(series) - 1) / 2) * w, s, w, label=lab)
        a.set_xticks(x); a.set_xticklabels(mods)
        a.set_title(title); a.set_ylabel(ylab)
        if ymax is not None:
            a.set_ylim(0, ymax)
        a.legend(fontsize=8); a.grid(axis="y", alpha=0.3)

    grouped(ax[0, 0], [ct.ari, ct.v_measure], ["ARI", "V-measure"],
            "Hard-label recovery (cell_type)", "score", 1.0)
    grouped(ax[0, 1], [ct.comp_mean_pearson, ct.comp_mean_jsd],
            ["mean per-type r", "mean JSD"],
            "Composition recovery (cell_type)", "score", 1.0)
    grouped(ax[0, 2], [ct.comp_r2, ct.comp_r2_oracle], ["method", "oracle (KMeans on truth)"],
            "Composition variance explained (cell_type)", "R²", 1.0)
    grouped(ax[1, 0], [ct.knn_concordance, ct.knn_chance], ["kNN concordance", "chance"],
            "Local kNN label concordance (cell_type)", "fraction", 1.0)
    grouped(ax[1, 1], [dm.ari, dm.v_measure], ["ARI", "V-measure"],
            "Hard-label recovery (domain)", "score", 1.0)
    ax[1, 2].bar(x, ct.eff_cats_per_obs_median, w * 1.6, color="#8c613c")
    ax[1, 2].set_xticks(x); ax[1, 2].set_xticklabels(mods)
    ax[1, 2].set_title("Aggregation confound:\neffective #cell-types per obs (median)")
    ax[1, 2].set_ylabel("exp(Shannon entropy)"); ax[1, 2].grid(axis="y", alpha=0.3)
    ax[1, 2].axhline(1.0, color="k", lw=0.8, ls="--")

    fig.suptitle("Figure 3B (reframed): clustering recovery of aggregated structure", y=1.0)
    fig.tight_layout()
    fig.savefig(out_png, dpi=160, bbox_inches="tight")
    print(f"[plot] {out_png}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--modalities", nargs="+", default=list(MODS),
                    choices=list(MODS))
    ap.add_argument("--out-dir", type=Path, default=OUT_DIR)
    ap.add_argument("--no-plot", action="store_true")
    # single-file mode: score one arbitrary clustered h5ad (e.g. a BANKSY sweep
    # build) instead of the canonical pca_harmony_single_cell/<mod>/ outputs.
    ap.add_argument("--h5ad", type=Path, default=None,
                    help="score just this clustered h5ad; writes composition_recovery_<tag>.json")
    ap.add_argument("--tag", default=None, help="name for the --h5ad run (default: file stem)")
    args = ap.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    if args.h5ad is not None:
        tag = args.tag or Path(args.h5ad).stem
        print(f"[score] {tag}  <-  {args.h5ad}")
        r = score_path(args.h5ad, tag)
        with open(args.out_dir / f"composition_recovery_{tag}.json", "w") as f:
            json.dump(r, f, indent=2)
        print(to_rows([r]).to_string(index=False))
        return

    results = []
    for mod in args.modalities:
        print(f"[score] {mod}")
        r = score_modality(mod)
        results.append(r)
        with open(args.out_dir / f"composition_recovery_{mod}.json", "w") as f:
            json.dump(r, f, indent=2)

    df = to_rows(results)
    csv = args.out_dir / "composition_recovery_summary.csv"
    df.to_csv(csv, index=False)
    print(f"[csv] {csv}")
    print(df.to_string(index=False))

    if not args.no_plot and set(args.modalities) >= {"cell", "bin", "bin16um", "spot"}:
        make_plot(df, args.out_dir / "composition_recovery.png")


if __name__ == "__main__":
    main()
