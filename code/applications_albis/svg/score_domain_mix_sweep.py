#!/usr/bin/env python
"""
Figure 4A follow-up -- score + plot the domain-mix coupling-strength sweep.

Reads each domain_mix_sweep_generate.py output dir under
data/figure_4/svg/domain_mix_sweep/<tag>/ (sweep_meta.json, gene_meta.csv,
sparkx_pvalues.csv written by run_sparkx.R), computes AUROC/AUPRC of
-log10(combinedPval) vs. is_noise, and plots both against the realized
domain/celltype coupling strength (`max_abs_dev_from_uniform`, NOT just the
nominal p_sig -- the model-free x-axis is more defensible).

Expected environment: either sim_paper/env/albis-tutorial or
sim_paper/env/stagate-pyg (both have scanpy/sklearn/matplotlib).

Example:
    python score_domain_mix_sweep.py
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from sklearn.metrics import average_precision_score, roc_auc_score  # noqa: E402

SWEEP_ROOT = Path("/dcs04/hicks/data/Jan/sim_project/sim_paper/data/figure_4/svg/domain_mix_sweep")


def main() -> None:
    rows = []
    for outdir in sorted(SWEEP_ROOT.iterdir()):
        meta_path = outdir / "sweep_meta.json"
        pv_path = outdir / "sparkx_pvalues.csv"
        if not meta_path.is_file():
            continue
        if not pv_path.is_file():
            print(f"[skip] no sparkx_pvalues.csv yet: {outdir}")
            continue
        meta = json.loads(meta_path.read_text())
        pv = pd.read_csv(pv_path)
        gmeta = pd.read_csv(outdir / "gene_meta.csv")
        gmeta.columns = ["gene_names", "is_noise"]
        df = pv.merge(gmeta, on="gene_names", how="left", validate="one_to_one")

        y_true = (~df["is_noise"]).astype(int)
        p = df["combinedPval"].to_numpy(dtype=float)
        score = -np.log10(np.clip(p, np.finfo(float).tiny, None))
        auroc = roc_auc_score(y_true, score)
        auprc = average_precision_score(y_true, score)

        rows.append({
            "tag": meta["out_tag"],
            "p_sig": meta["p_sig"],
            "coupling_strength": meta["max_abs_dev_from_uniform"],
            "n_obs_qc": meta["n_obs_qc"],
            "n_signal": int(y_true.sum()),
            "n_noise": int((1 - y_true).sum()),
            "auroc": auroc,
            "auprc": auprc,
        })
        print(f"[load] {meta['out_tag']}: p_sig={meta['p_sig']:.3f} "
              f"coupling={meta['max_abs_dev_from_uniform']:.4f} AUROC={auroc:.4f} AUPRC={auprc:.4f}")

    if not rows:
        raise SystemExit(f"No scored sweep points found under {SWEEP_ROOT}")

    summary = pd.DataFrame(rows).sort_values("p_sig")
    csv_path = SWEEP_ROOT / "domain_mix_sweep_summary.csv"
    summary.to_csv(csv_path, index=False)
    print(f"[save] {csv_path}")

    # x-axis = p_sig, the cleanly-controlled sweep parameter, not the
    # coupling_strength diagnostic -- that crosstab is computed on the raw
    # pre-QC adata (includes the ~77% off-tissue "unassigned" realwindow
    # spots), which appears to add a p_sig-independent floor (0.125 -> 0.125
    # unchanged at the first two sweep points) rather than cleanly tracking
    # only the intended on-tissue coupling change. Kept in the CSV for
    # reference but not plotted.
    fig, ax = plt.subplots(figsize=(6, 4.5))
    ax.plot(summary["p_sig"], summary["auroc"], "o-", color="#4393c3", label="AUROC")
    ax.plot(summary["p_sig"], summary["auprc"], "s-", color="#f4a582", label="AUPRC")
    ax.axhline(0.5, color="black", linewidth=0.8, linestyle=":", label="chance (AUROC)")
    ax.axvline(1 / 8, color="gray", linewidth=0.8, linestyle="--", label="p_sig=1/8 (no coupling)")
    ax.set_xlabel("Domain signature-type probability p_sig\n(1/8 = uniform/no coupling, 0.5 = near-deterministic)")
    ax.set_ylabel("Score")
    ax.set_title("Figure 4A follow-up: SVG recovery (SPARK-X) vs.\ndomain-celltype coupling strength (spot)")
    ax.set_ylim(0, 1.05)
    ax.legend(loc="lower right", fontsize=8, frameon=False)
    ax.grid(ls=":", alpha=0.4)
    ax.set_axisbelow(True)
    fig.tight_layout()

    out_path = SWEEP_ROOT / "domain_mix_sweep.png"
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"[save] {out_path}")


if __name__ == "__main__":
    main()
