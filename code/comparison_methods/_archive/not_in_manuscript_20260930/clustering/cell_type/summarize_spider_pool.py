import json,csv
from pathlib import Path
ROOT=Path('/dcs04/hicks/data/Jan/sim_project/sim_paper/data/figure_5/clustering/cell_type/spider_pool_diagnostic_20260928')
rows=[]
for genes in [556,2000]:
 for seed in [2025,101,202]:
  p=ROOT/f'genes{genes}_seed{seed}/diagnostic.json'
  if not p.exists():continue
  d=json.loads(p.read_text());r={k:d[k] for k in ['genes','seed','pool_unique_profiles','spider_unique_profiles','spider_cells_not_in_pool','spider_cells_with_label_mismatch']}
  for m in d['pool_clustering']:r[f"pool_ari_{m['weighting']}{m['k']}"]=m['ari']
  r['expanded_ari']=d['expanded_clustering']['cell_type_ari'];rows.append(r)
if rows:
 with (ROOT/'summary.csv').open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
lines=['# SPIDER expression-pool diagnostic',f'Completed: {len(rows)}/6 conditions.','', 'Original 10k pools are clustered without Harmony using Gaussian K256 and UMAP K15; Leiden targets eight types. Expanded raw profiles are matched exactly using canonical SHA256 hashes.','']
for r in rows:lines.append(f"- {r['genes']} genes, seed {r['seed']}: pool ARI Gaussian256={r['pool_ari_gauss256']:.4f}, UMAP15={r['pool_ari_umap15']:.4f}; expanded ARI={r['expanded_ari']:.4f}; {r['spider_unique_profiles']} unique expanded profiles; unmatched cells={r['spider_cells_not_in_pool']}; label mismatches={r['spider_cells_with_label_mismatch']}.")
lines+=['','The installed SPIDER get_sim_cell_level_expr samples reference cells with replacement within each assigned type and copies their expression unchanged. Perfect pool ARI demonstrates separability before spatial expansion. Replication can fragment low-K graphs; this does not imply additional biological expression diversity or realism.','No production clustering results have been changed.']
(ROOT/'REPORT.md').write_text('\n'.join(lines)+'\n')
print(f'Summarized {len(rows)}/6 conditions')
