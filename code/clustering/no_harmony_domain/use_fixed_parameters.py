"""Write selected_parameters.json from an earlier selection instead of re-sweeping.

Use after `prepare` (in place of `sweep` + `select`) to rerun the final runs on
regenerated data with the manuscript's BANKSY lambdas kept fixed:

  python use_fixed_parameters.py --experiment domain --source OLD/selected_parameters.json
  python use_fixed_parameters.py --experiment cell_type --source OLD/selected_parameters.json

Only lambda and k_geom are carried over; the old ARIs and paths are not, since
they describe the earlier data.
"""
import argparse
import json
from pathlib import Path
from common import BASE, K_GEOM, MODALITIES, save

ROOTS = dict(domain=BASE, cell_type=BASE.parent / 'no_harmony_cell_type')


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--experiment', choices=ROOTS, required=True)
    p.add_argument('--source', type=Path, required=True)
    args = p.parse_args()
    root = ROOTS[args.experiment]
    old = json.loads(args.source.read_text())
    choices = {}
    for m in MODALITIES:
        lam = old['modalities'][m]['lambda']
        assert old['modalities'][m].get('k_geom', K_GEOM[m]) == K_GEOM[m], f'{m}: k_geom changed'
        choices[m] = dict(modality=m, **{'lambda': lam}, k_geom=K_GEOM[m])
    result = {k: v for k, v in old.items() if k != 'modalities'}
    result.update(modalities=choices, fixed_from=str(args.source.resolve()),
                  note='Lambda fixed from the earlier selection; no sweep or re-selection on this data.')
    path = root / 'selected_parameters.json'
    assert root.is_dir() and (root / 'final_tasks.json').exists(), f'run prepare for {args.experiment} first'
    assert not path.exists(), f'{path} already exists'
    save(path, result)
    print(json.dumps(result, indent=2), flush=True)


if __name__ == '__main__':
    main()
