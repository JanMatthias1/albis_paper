# Real data QC (Figure 2 references)

QC of the real 10x datasets that the ALBIS simulations are compared against in
Figure 2. Each script reads the raw 10x output from `data/real_data/<dataset>/`
and writes a QC-filtered AnnData with raw counts, plus QC plots and a
`qc_summary.json`, to `data/real_data_qc/<dataset>/`. No normalization is
applied.

| Technology | Datasets | Script | QC |
|---|---|---|---|
| Xenium (cells) | non-diseased lung, lung cancer | `xenium_qc.py` (`run_real_data_qc.sh`) | 4-MAD low-tail cutoffs on total counts and genes per cell; drops cells with zero count density or control probe/codeword counts |
| Visium HD (8 / 16 µm bins) | breast cancer, human pancreas | `visium_hd_qc.py` (raw archives unpacked first with `extract_breast_visium_hd.py`) | SpotSweeper local outliers on total counts, genes and % mitochondrial |
| Xenium (cells) | human brain: healthy, glioblastoma, Alzheimer's (FFPE, with add-on panel; raw in `data/real_data/brain_samples_xenium/`) | `xenium_qc.py` (`run_brain_xenium_qc.sh`) → `data/real_data_qc/brain_samples_xenium/<slice>/` | same as the lung Xenium QC |
| Visium (spots) | lymph node, tonsil (CytAssist, probe-based) | `visium_qc.py` (`run_visium_lymph_qc.sh`, `run_visium_tonsil_qc.sh`) | SpotSweeper local outliers, as for Visium HD |

Other scripts:

- `plot_qc_supplement.py`: supplementary figure showing kept and removed
  observations for every dataset → `data/real_data_qc/supplementary/`
- `pca_umap.py` (`run_pca_umap.sh`): PCA/UMAP of each QC'd dataset, for
  comparison with the simulated data

Environment: `albis-tutorial`, plus `spotsweeper` for the Visium/Visium HD QC
and `scikit-misc` for `pca_umap.py`.
