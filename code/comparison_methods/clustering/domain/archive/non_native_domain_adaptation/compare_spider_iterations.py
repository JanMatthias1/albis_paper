"""Compare native optimizer diagnostics, without clustering-score tuning."""
import argparse
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baseline', type=Path, required=True)
    parser.add_argument('--current', type=Path, required=True)
    args = parser.parse_args()
    old, new = [json.loads((p/'manifest.json').read_text()) for p in (args.baseline, args.current)]
    for key in ('seed', 'n_cells', 'n_slices', 'domain_prior', 'celltype_priors', 'domain_target'):
        assert old[key] == new[key], f'Unmatched parameter: {key}'
    comparison = {name: {key: data[key] for key in
        ('max_iterations', 'domain_loss', 'components', 'zones')}
        for name, data in [('baseline', old), ('current', new)]}
    comparison['domain_loss_reduction'] = old['domain_loss'] - new['domain_loss']
    comparison['note'] = 'Connectivity is diagnostic; convergence does not guarantee connected layers.'
    (args.current/'iteration_comparison.json').write_text(json.dumps(comparison, indent=2)+'\n')
    history = pd.read_csv(args.current/'domain_optimization_history.csv')
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    axes[0].plot(history.iteration, history.best_loss)
    axes[0].axhline(old['domain_loss'], color='gray', linestyle='--', label='50k run')
    axes[0].axhline(new['solver_config']['tol'], color='black', linestyle=':', label='Tolerance')
    axes[0].set(xlabel='Iteration', ylabel='Best transition loss'); axes[0].legend()
    groups = sorted(new['components'])
    x = list(range(len(groups)))
    for offset, data, label in [(-.2, old, '50k'), (.2, new, '3 million budget')]:
        axes[1].bar([v+offset for v in x], [data['components'][g]['largest_fraction'] for g in groups],
                    width=.4, label=label)
    axes[1].set(xticks=x, xticklabels=groups, ylim=(0,1), ylabel='Largest connected fraction')
    axes[1].legend(); fig.tight_layout()
    fig.savefig(args.current/'iteration_comparison.png', dpi=200)
    plt.close(fig)
    print(json.dumps(comparison), flush=True)


if __name__ == '__main__':
    main()
