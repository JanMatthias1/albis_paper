"""Prepare a multi-seed native compute run: settings for every size x seed, the task
table and the SHA-256 of every run script (recorded in protocol.json; the jobs run the
scripts from this folder, nothing is copied into the run). Settings scale the 600k template exactly as every
earlier run did: extent, ALBIS sphere radius and core fuzz width by (n / 600000)^(1/3),
target cells per type = n / 8; the seed is set in `seed` and `albis.seed`."""
import argparse
import hashlib
import json
from pathlib import Path

SIZES = [10000, 100000, 200000, 400000, 600000, 1000000, 2000000, 5000000]
SEEDS = [2025, 101, 202]
# scCube trains one VAE per batch, as in the single-seed run (10k-100k, 200k-600k, 1M-5M)
SCCUBE_BATCHES = {'b10k_100k': [10000, 100000], 'b200k_600k': [200000, 400000, 600000],
                  'b1m_5m': [1000000, 2000000, 5000000]}
CODE_FILES = ['albis_native.py', 'spider_native.py', 'sccube_split.py', 'sample_rss.py', 'reference_native.R',
              'reference_native.sbatch', 'native_extension.sbatch', 'sccube_seed_batch.sbatch',
              'prepare_seeds.py', 'submit_seeds.sh', 'settings_template_n600000_seed2025.json']

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--root', type=Path, required=True)
args = parser.parse_args()
root = args.root.resolve()
code = Path(__file__).resolve().parent
root.mkdir(parents=True, exist_ok=False)
for folder in ['settings', 'raw', 'logs', 'references', 'sccube']:
    (root / folder).mkdir()

template = json.loads((code / 'settings_template_n600000_seed2025.json').read_text())
rows = ['method\tn_cells\tseed\tsettings\toutput']
for seed in SEEDS:
    for n in SIZES:
        s = json.loads(json.dumps(template))
        factor = (n / 600000) ** (1 / 3)
        s.update(n_cells=n, seed=seed, target_cells_per_type={f'Group{i}': n // 8 for i in range(1, 9)},
                 extent_um=template['extent_um'] * factor)
        s['albis'].update(n_cells=n, seed=seed, sphere_R_um=template['albis']['sphere_R_um'] * factor,
                          core_fuzz_width_um=template['albis']['core_fuzz_width_um'] * factor)
        path = root / 'settings' / f'n{n}_seed{seed}.json'
        path.write_text(json.dumps(s, indent=2) + '\n')
        for method in ['albis', 'spider']:
            rows.append('\t'.join(map(str, [method, n, seed, path, root / 'raw' / f'{method}_n{n}_seed{seed}'])))
    for batch, sizes in SCCUBE_BATCHES.items():
        rows.append('\t'.join(map(str, ['sccube', ','.join(map(str, sizes)), seed, root / 'settings',
                                        root / 'sccube' / f'seed{seed}_{batch}'])))
(root / 'tasks.tsv').write_text('\n'.join(rows) + '\n')
(root / 'protocol.json').write_text(json.dumps(dict(
    seeds=SEEDS, sizes=SIZES, sccube_batches=SCCUBE_BATCHES, n_genes=556, threads=1,
    hardware_constraint='sapphirerapids', reference='one Splatter pool per seed (reference_native.R)',
    template='settings_template_n600000_seed2025.json (copied from the single-seed pilot run)',
    source_hashes={n: hashlib.sha256((code / n).read_bytes()).hexdigest() for n in CODE_FILES}),
    indent=2) + '\n')
print(root)
