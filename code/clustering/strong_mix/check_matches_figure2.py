"""Check every strong-mix / seed dataset is the Figure 2 config + strong mix only.

Compares uns['sim_params'] of each
data/figure_3/cellbin_batch_sigma_slide/{cell,bin16um,spot}/bs*/ raw h5ad
(one tree per point since 2026-09-23; previously data/noisy/) against the
canonical Figure 2 dataset of the same modality. The only allowed differences are
strong_domain_mix (must be True here, False in Figure 2), batch_sigma (the
swept variable) and seed (seed replicates). Anything else is reported as a
mismatch and the script exits non-zero.

Usage: python check_matches_figure2.py   (reads only; safe while jobs run)
"""
import sys
from pathlib import Path

import h5py
import numpy as np

DATA = Path(__file__).resolve().parents[3] / "data"
FIG2 = DATA / "figure_2" / "smaller_sphere" / "data"
CANON = {
    "cell": FIG2 / "log_mu_-2.5_theta_0.40_jitter0.15_bsigma15" / "simulation_cell_z.h5ad",
    "bin16um": FIG2 / "packing_pf0p04_bin16um_log_mu_-2.5_jitter0.6_bsigma07" / "simulation_bin_z.h5ad",
    "spot": FIG2 / "packing_pf0p04_log_mu_-2.25_sigma1.0_theta_0.25_jitter0.10_dsf_bsigma03" / "simulation_spot_z.h5ad",
}
SLIDE = DATA / "figure_3" / "cellbin_batch_sigma_slide"
STEM = {"cell": "cell", "bin16um": "bin", "spot": "spot"}
ALLOWED = {"strong_domain_mix", "batch_sigma", "seed"}


def params(path):
    def walk(g, prefix=""):
        out = {}
        for k, v in g.items():
            if isinstance(v, h5py.Group):
                out.update(walk(v, f"{prefix}{k}."))
            else:
                x = v[()]
                out[prefix + k] = x.decode() if isinstance(x, bytes) else np.asarray(x).tolist()
        return out
    with h5py.File(path, "r") as f:
        return walk(f["uns/sim_params"])


def main():
    canon = {m: params(p) for m, p in CANON.items()}
    n_ok = n_bad = n_pending = 0
    for d in sorted(SLIDE.glob("*/bs*")):
        mod = d.parent.name
        if mod not in STEM or not d.is_dir():
            continue
        # Raw file is named by its generation parameters (see generate_strong_mix_*.sh):
        # the one *_strongmix_bsigma*.h5ad that isn't the _qc file.
        raws = [p for p in d.glob("*_strongmix_bsigma*.h5ad") if not p.name.endswith("_qc.h5ad")]
        if len(raws) > 1:
            n_bad += 1
            print(f"MISMATCH {mod}/{d.name}: {len(raws)} raw files: {[p.name for p in raws]}")
            continue
        raw = raws[0] if raws else d / "missing.h5ad"
        if not raw.is_file():
            n_pending += 1
            continue
        try:
            got = params(raw)
        except OSError:  # still being written
            n_pending += 1
            continue
        ref = canon[mod]
        diffs = [(k, ref.get(k), got.get(k)) for k in sorted(set(ref) | set(got))
                 if k not in ALLOWED and ref.get(k) != got.get(k)]
        if got.get("strong_domain_mix") is not True:
            diffs.append(("strong_domain_mix", "True expected", got.get("strong_domain_mix")))
        if diffs:
            n_bad += 1
            print(f"MISMATCH {mod}/{d.name}: " + "; ".join(f"{k}: fig2={a} here={b}" for k, a, b in diffs))
        else:
            n_ok += 1
    print(f"\n{n_ok} match Figure 2 + strong mix, {n_bad} mismatch, {n_pending} not written yet")
    sys.exit(1 if n_bad else 0)


if __name__ == "__main__":
    main()
