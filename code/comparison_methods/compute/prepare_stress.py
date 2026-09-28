"""Create a fresh full sweep with 16um bins and 55um circular spots."""
import argparse
import json
from pathlib import Path

PROJECT = Path(__file__).resolve().parents[4]
p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--root', type=Path, required=True)
a = p.parse_args()
a.root.mkdir(parents=True, exist_ok=False)
(a.root/'raw').mkdir()
(a.root/'logs').mkdir()
# Geometry changed: rerun the baseline instead of mixing old square captures.
tasks = []
for n in [10000, 20000, 40000, 50000, 100000, 200000, 500000, 1000000]:
    for tech in ['cell', 'bin', 'spot']:
        for method, step in [('albis','our_method'),('spider','splatter_spider'),('sccube','splatter_sccube')]:
            tasks.append(dict(index=len(tasks), n_cells=n, technology=tech, method=method,
                              step=step, directory=f'n{n}_{tech}_{method}'))
(a.root/'stress_plan.json').write_text(json.dumps(dict(tasks=tasks, memory='64G',
    wall_limit='2-00:00:00', cpus=4, max_concurrent_tasks=3, n_genes=556,
    scenario='random_null', independent_method_jobs=True,
    geometry_version='bin16_circle55_spacing100',
    geometry='ALBIS/SPIDER: 16um square bins, 55um circular spots at 100um spacing; scCube native occupancy grid',
    splatter='Fresh full-size pool for each competitor task; included in metrics'), indent=2)+'\n')
print(a.root)
