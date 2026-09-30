# Domain recovery: active code

- `run_banksy.py`: BANKSY analysis of supplied domain-labelled data; does not generate domains.
- `run_pilot.sh albis`: existing ALBIS-only pilot launcher.

The custom two-stage SPIDER domain-generation workflow, launchers, plots and diagnostics have been removed from active code and archived under `archive/non_native_domain_adaptation/`. Related data are archived under `data/figure_5/clustering/domain/archive/non_native_domain_adaptation/` and are not benchmark inputs.

No active SPIDER domain-generation launcher remains. Future SPIDER domain simulations must call its native geometric-pattern functions directly (`layer_cell_level_sim` or `naive_cell_level_sim`) and retain their native region and cell-type outputs. Those functions are 2D; no custom domain construction, solver composition or artificial 3D stacking is permitted. The installed SPIDER package was not modified.
