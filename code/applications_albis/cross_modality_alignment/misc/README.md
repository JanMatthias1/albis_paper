# Retired cross-modality code (moved here 2026-09-25)

Kept for provenance only. None of these scripts are called by the current
pipelines. Their internal paths still assume the old flat layout.

- `generate_independent_offsets.sh`, `submit_independent_offsets.sh`,
  `run_independent_offsets_stair.sh`, `validate_independent_offsets.py`: the
  original 1x (max_shift 1025 um) independent-offset run -> `independent_offsets/`.
- `prepare_stronger_shifts.py`, `run_stronger_shifts_stair.sh`: the post-hoc
  3x translation rescale. Superseded by the native regeneration in
  `../weak_domain_mix/` (archived data: `independent_offsets_shift3x_before_native_*`).
- `run_cross_tech_stair.sh`, `logs_cross_tech_stair/`: the original Fig4E
  cross-tech STAIR runner and its logs.
- `plot_slice5_independent_alignment.py`: diagnostic that aligns each technology
  on its own (per-technology 3D_stair.py runs).
- `plot_slice5_zorder_test.py`: diagnostic for overlay draw order.
