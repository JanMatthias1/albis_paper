# Simulated tissue and capture (Figure 1)

Illustrates how ALBIS builds and captures a tissue: one 3D sphere, a capture
window, sectioning into slices and binning, and the per-slice misalignment.
Both scripts use the same simulation settings as `code/data/generate_simulation_noisy.py`.

| Panel | Shows | Script |
|---|---|---|
| 1A | intact 3D sphere before slicing, cells coloured by spatial domain | `01_plot_intact_sphere.py` |
| 1B | capture window: captured cells coloured by domain, the rest grey | `02_plot_capture_and_bins.py` |
| 1C | the sphere sectioned into slices and aggregated into Visium HD-like bins | `02_plot_capture_and_bins.py` |
| 1D | the same bins at each slice's random rotation and shift, before any alignment | `02_plot_capture_and_bins.py` |

```bash
sbatch run_intact_sphere.sh          # 1A
sbatch run_capture_and_bins.sh       # 1B-D (--max-shift: maximum per-slice shift for 1D, in µm)
```

Options for both: `--max-cells` (cells drawn, for plotting speed), `--dpi`, `--outdir`.
Figures (PNG and PDF) are written to `outputs/`; `outputs/shift1500/` holds a
1B-D run with `--max-shift 1500`.

Environment: `albis-tutorial` (`env/create_tutorial_env.sh`).
