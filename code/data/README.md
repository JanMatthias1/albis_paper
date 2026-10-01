# Simulation generation

Generates the ALBIS simulations used across the figures, one modality (cell,
bin or spot) per call. Each run writes an h5ad with the 10-slice z-stack, a
`_summary.json` with the settings used, and diagnostic plots.

| File | Purpose |
|---|---|
| `generate_simulation_noisy.py` | one ALBIS simulation via `albis.simulate_3d_molecule_sphere_multires()` |
| `run_generate_simulation_noisy.sh` | SLURM array over spot, bin and cell; extra flags are passed through (`--check-env` only checks the environment) |
| `plot_qc_slices.py` | per-slice spatial plot of a QC'd Figure 2 simulation |

The figure scripts call `generate_simulation_noisy.py` with their own settings
(for example `code/count_distribution/smaller_sphere/*.sh`). The main options are:

| Option | Controls |
|---|---|
| `--modality`, `--bin-size-um` | cell, bin (8 / 16 µm) or spot output |
| `--base-gene-lognormal`, `--theta`, `--theta-jitter`, `--noise-scale` | gene means and count dispersion |
| `--marker-foldchange`, `--shared-marker-foldchange` | marker-gene strength |
| `--strong-domain-mix`, `--domain-size-factors` | how distinct the spatial domains are |
| `--batch-sigma` | per-slice batch effect |
| `--sphere-r-um`, `--n-cells`, `--cell-r-mean`, `--capture-window-um` | tissue size, cell density and capture area |
| `--core-fuzz-width-um`, `--max-shift` | domain boundaries and slice misalignment |
| `--seed`, `--out-tag`, `--output-dir`, `--output-stem` | reproducibility and output location |

Output goes to `data/noisy/[<out-tag>/]` unless `--output-dir` is set.

Requires `albis` 0.1.2 (`pip install albis==0.1.2`); the script checks the
version at start-up.

```bash
python generate_simulation_noisy.py --modality cell --out-tag my_run
sbatch run_generate_simulation_noisy.sh [flags]
```
