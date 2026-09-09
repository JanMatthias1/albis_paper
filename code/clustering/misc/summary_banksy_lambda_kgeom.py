#!/usr/bin/env python
"""
Tabulate the BANKSY lambda x k_geom sweep (banksy_lambda_kgeom_sweep.sh).

Reads data/figure_3/banksy_lambda_kgeom_sweep/scores/composition_recovery_*.json
(one per {cell,bin} x lambda x k_geom build, written by
composition_recovery.py --h5ad) plus the plain-pipeline and BANKSY-default
baselines already on disk, and prints a grid per modality:

  domain ARI | domain comp-R2 | R2 vs oracle ceiling | domain kNN-lift |
  slice_id leakage ARI | cell_type ARI

A run only counts as a genuine domain-recovery hit if domain ARI is up AND the
slice_id leakage ARI stays low (< ~0.3) -- otherwise the "clusters" are just
tracking the staggered-slice batch.

Usage:  python sim_paper/code/clustering/misc/summary_banksy_lambda_kgeom.py
"""
from __future__ import annotations
import glob
import json
import re
from pathlib import Path

import pandas as pd

ROOT = Path("/dcs04/hicks/data/Jan/sim_project/sim_paper/data/figure_3")
SWEEP = ROOT / "banksy_lambda_kgeom_sweep"

BASELINES = {  # (modality, pipeline) -> ari_summary path
    ("cell", "plain (weak mix)"): ROOT / "pca_harmony_single_cell/cell/ari_recovery_qc/ari_summary_cell.json",
    ("cell", "plain (strong mix)"): ROOT / "../clustering_log_mu_-2.5_theta_0.25_strongdomainmix/cell/ari_recovery_qc/ari_summary_cell.json",
    ("cell", "BANKSY dflt (l0.8 kg15 stag1.5)"): ROOT / "../clustering_log_mu_-2.5_theta_0.25_strongdomainmix/cell/banksy_ari_recovery/ari_summary_cell.json",
    ("bin", "plain (weak mix)"): ROOT / "pca_harmony_single_cell/bin/ari_recovery_qc/ari_summary_bin.json",
    ("bin", "plain (strong mix)"): ROOT / "../clustering_packing_pf0p04_bsigma05_strongdomainmix/bin/ari_recovery_qc/ari_summary_bin.json",
    ("bin", "BANKSY dflt (l0.8 kg15 stag1.5)"): ROOT / "../clustering_packing_pf0p04_bsigma05_strongdomainmix/bin/banksy_ari_recovery/ari_summary_bin.json",
}


def ari_of(path: Path) -> dict:
    try:
        d = json.load(open(path))
    except Exception:
        return {}
    return {r["ground_truth"]: round(r["ari"], 4) for r in d}


def main() -> None:
    # baselines
    print("=" * 78)
    print("BASELINES (hard-label ARI)")
    print(f"{'modality':<6} {'pipeline':<34} {'domain':>8} {'cell_type':>10}")
    print("-" * 62)
    for (mod, pipe), p in BASELINES.items():
        a = ari_of(p)
        print(f"{mod:<6} {pipe:<34} {str(a.get('domain_true','-')):>8} {str(a.get('cell_type_true','-')):>10}")

    # sweep
    rows = []
    for f in sorted(glob.glob(str(SWEEP / "scores" / "composition_recovery_*.json"))):
        m = re.search(r"composition_recovery_(cell|bin)_lam([\d.]+)_kg(\d+)\.json$", f)
        if not m:
            continue
        mod, lam, kg = m.group(1), float(m.group(2)), int(m.group(3))
        d = json.load(open(f))
        lv = d["levels"]
        dom, ct = lv["domain"], lv["cell_type"]
        rows.append(dict(
            modality=mod, lam=lam, k_geom=kg,
            domain_ari=round(dom["hard_label"]["ari"], 4),
            domain_v=round(dom["hard_label"]["v_measure"], 4),
            domain_r2=round(dom["composition"]["r2_variance_explained"], 4),
            domain_r2_oracle=round(dom["oracle"]["r2_oracle"], 4),
            domain_r2_frac=(round(dom["oracle"]["r2_method_over_oracle"], 3)
                            if dom["oracle"]["r2_method_over_oracle"] is not None else None),
            domain_knn=round(dom["local_knn"]["concordance"], 4),
            domain_knn_chance=round(dom["local_knn"]["chance_baseline"], 4),
            leak_ari=(round(dom["slice_id_leakage_ari"], 4)
                      if dom.get("slice_id_leakage_ari") is not None else None),
            celltype_ari=round(ct["hard_label"]["ari"], 4),
            celltype_r2=round(ct["composition"]["r2_variance_explained"], 4),
        ))
    if not rows:
        print("\n(no sweep score files yet under", SWEEP / "scores", ")")
        return
    df = pd.DataFrame(rows).sort_values(["modality", "lam", "k_geom"])
    out_csv = SWEEP / "summary_banksy_lambda_kgeom.csv"
    df.to_csv(out_csv, index=False)

    print("\n" + "=" * 78)
    print("SWEEP  (--stagger-scale 5)")
    for mod in ("cell", "bin"):
        sub = df[df.modality == mod]
        if sub.empty:
            continue
        print(f"\n--- {mod} ---")
        print(sub.drop(columns="modality").to_string(index=False))
        hits = sub[(sub.domain_ari > 0.10) & (sub.leak_ari.fillna(1) < 0.30)]
        if len(hits):
            print(f"  >> genuine domain-recovery hit(s): "
                  + ", ".join(f"l{r.lam}/kg{r.k_geom} (ARI {r.domain_ari}, leak {r.leak_ari})"
                              for r in hits.itertuples()))
        else:
            print("  >> no genuine hit (domain ARI stays <=0.10, or gains are slice_id leakage)")
    print(f"\n[csv] {out_csv}")


if __name__ == "__main__":
    main()
