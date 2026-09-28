# Figure5A supplement: scCube spot aggregation

From the project root:

```bash
MPLCONFIGDIR=/tmp/fig5a_spot_aggregation comparison_methods/env/analysis/bin/python sim_paper/code/comparison_methods/overview/spot_aggregation/plot_sccube_spot_aggregation.py
```

Reuses the Figure5A scCube slice5 cells. Displays targets3,100,1000,5000,10000 cells/spot with the same spatial scale; writes PNG/PDF/SVG, per-spot counts, memberships and provenance to data/comparison_methods/figure_5A_sccube_spot_aggregation/.

Targets are requested mean occupancies, not exact per-spot cell counts. Native grid rounding and spatial density determine realized occupancy. The geometry-only implementation is validated against saved native Figure5A spots at target10; no expression is regenerated.
