#!/usr/bin/env python
# coding: utf-8

# In[1]:


# Preprocessing script for Visium human breast data
# enviroment: preprocessing_JM
# Adapted from Visium_HD_human_colon_8_um.ipynb (multi-sample-alignment-benchmark/code/01_preprocessing/)

from pathlib import Path
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
from statsmodels import robust
import spotsweeper.local_outliers as lo
import spotsweeper.plot_QC as plot_QC
import spotsweeper.plot_QCpdf as pdf
import gzip
import os
import shutil


# In[2]:


# Raw downloads are extracted into per-sample folders by extract_breast_visium_hd.py
# (sim_paper/code/real_data_qc/extract_breast_visium_hd.py) before this script runs:
#   conda activate /dcs04/hicks/data/multi-sample-alignment-benchmark/envs/preprocessing_JM
#   python extract_breast_visium_hd.py


# In[3]:


def unzip_into_adata(data_path, library_id=None):
    """
    Load Visium / Visium HD data into AnnData.
    Expected folder structure:
    sample/
    ├── filtered_feature_bc_matrix.h5
    └── spatial/
        ├── scalefactors_json.json
        ├── tissue_hires_image.png
        ├── tissue_lowres_image.png
        ├── tissue_positions_list.csv (or tissue_positions.parquet)
    """
    data_path = Path(data_path)
    spatial_path = data_path / "spatial"

    # --- 1. Unzip any .gz files in the spatial folder ---
    for gz_file in spatial_path.glob("*.gz"):
        out_file = gz_file.with_suffix('')
        print(f"Unzipping {gz_file.name} -> {out_file.name}")
        with gzip.open(gz_file, 'rb') as f_in, open(out_file, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
        gz_file.unlink()

    # --- 2. Convert tissue_positions.parquet to CSV if needed ---
    parquet_file = spatial_path / "tissue_positions.parquet"
    csv_file = spatial_path / "tissue_positions_list.csv"
    if parquet_file.exists() and not csv_file.exists():
        print(f"Converting {parquet_file.name} → {csv_file.name}")
        df = pd.read_parquet(parquet_file)
        expected_cols = [
            "barcode",
            "in_tissue",
            "array_row",
            "array_col",
            "pxl_row_in_fullres",
            "pxl_col_in_fullres"
        ]
        # enforce column order
        df = df[expected_cols]
        df.to_csv(csv_file, index=False, header=False)

    # --- 3. Load Visium HD data ---
    print(f"Reading Visium HD data from: {data_path}")
    adata = sc.read_visium(
        path=str(data_path),
        count_file="filtered_feature_bc_matrix.h5",
        library_id=library_id,
        load_images=True
    )

    return adata

def mad_based_cutoffs(x, direction="both", log=False, nmads=3):
    """
    From Vani's script

    Returns lower and/or upper cutoffs based on median ± nmads * MAD.
    direction: "both", "lower", or "higher"
    """
    x = np.array(x)
    if log:
        x = np.log1p(x)

    med = np.median(x)
    mad = robust.mad(x, c=1)  # Median Absolute Deviation
    lower = med - nmads * mad
    upper = med + nmads * mad

    if log:
        lower, upper = np.expm1(lower), np.expm1(upper)

    if direction == "lower":
        return lower
    elif direction == "higher":
        return upper
    else:
        return lower, upper

# Helper: remove genes starting with ERCC / MT-
def prefilter_specialgenes(adata, Gene1Pattern="ERCC", Gene2Pattern="MT-"):
    mask1 = ~adata.var_names.str.startswith(Gene1Pattern)
    mask2 = ~adata.var_names.str.startswith(Gene2Pattern)
    keep = mask1 & mask2
    adata._inplace_subset_var(keep)
    return adata


# In[4]:


# Script based on Vani's preprocessing steps, CORRECTED

# Define your sample names and file paths
# CHANGE THE PATH AND FILE NAMES ACCORDINGLY, THE FOLLOWING STEPS WILL PLOT THE QC METRICS,
# 8um bins to match sim_app's bin_size_um=8.0 (see sim_paper/code/data/generate_simulation.py)

sample_files = {
    "Visium_HD_11mm_Human_Breast_Cancer": "/dcs04/hicks/data/Jan/sim_project/sim_paper/data/real_data/breast_cancer/Visium_HD_11mm_Human_Breast_Cancer/binned_outputs/square_008um",
}

# Initialize a dictionary to store the AnnData objects and QC thresholds
adata_dict = {}

for sample_name, data_path in sample_files.items():

    print(f"Loading {sample_name}")
    # Load the data from the above sample_files
    adata = unzip_into_adata(data_path, library_id=None)

    adata.var_names_make_unique()
    print("Raw shape:", adata.shape)

    # Retain only in_tissue spots
    if "in_tissue" not in adata.obs:
        raise ValueError(f"{sample_name}: 'in_tissue' missing in .obs")

    before = adata.n_obs
    adata = adata[adata.obs["in_tissue"] == 1].copy()
    after = adata.n_obs
    print(f"Spots retained after in_tissue filter: {after}/{before}")

    # Make annotation categorical
    #adata.obs["annotation"] = adata.obs["annotation"].astype("category")
    #print("Annotation value counts:")
    #print(adata.obs["annotation"].value_counts())

    # 3. Compute QC metrics and mitochondrial metrics
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)
    # Now have: total_counts, n_genes_by_counts, pct_counts_mt in adata.obs

    # 4. Spotsweeper: detect local outliers & retain QC-passing spots
    adata.obs["region"] = "sample1"
    lo.local_outliers(adata, metric="total_counts", sample_key="region", direction="lower", log=True)
    lo.local_outliers(adata, metric="n_genes_by_counts", sample_key="region", direction="lower", log=True)
    adata.obs.rename(columns={"n_genes_by_counts_outliers": "n_genes_outliers"}, inplace=True)
    lo.local_outliers(adata, metric="pct_counts_mt", sample_key="region", direction="higher", log=False)
    adata.obs.rename(columns={"pct_counts_mt_outliers": "mt_outliers"}, inplace=True)

    # Plot Spotsweeper outliers
    #fig_out = plot_dir / f"{sample_name}_spotsweeper_outliers.png"
    sc.pl.spatial(
        adata,
        img_key="hires",
        color=["total_counts_outliers", "n_genes_outliers", "mt_outliers"],
        size=1.5,
        title=f"{sample_name}: Spotsweeper outliers",
        show=False,
    )
    #plt.savefig(fig_out, dpi=300, bbox_inches="tight")
    plt.show()
    plt.close()

    # Keep only spots that pass all three outlier filters
    before = adata.n_obs
    adata = adata[
        (~adata.obs["total_counts_outliers"]) &
        (~adata.obs["n_genes_outliers"]) &
        (~adata.obs["mt_outliers"]),
        :
    ].copy()
    after = adata.n_obs
    print(f"Remaining spots after Spotsweeper QC: {after}/{before}")

    # 5. Remove lowly expressed genes + ERCC/MT genes
    sc.pp.filter_genes(adata, min_cells=3)
    print("Shape after filter_genes(min_cells=3):", adata.shape)

    prefilter_specialgenes(adata, Gene1Pattern="ERCC", Gene2Pattern="MT-")
    print("Shape after removing ERCC/MT genes:", adata.shape)
    """
    # 6. Normalize, log1p, HVG
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=2000)
    print("Number of HVGs:", adata.var["highly_variable"].sum())
    hvg = adata.var["highly_variable"]
    print("Any ERCC in HVGs?", any(adata.var_names[hvg].str.startswith("ERCC")))
    print("Any MT- in HVGs?", any(adata.var_names[hvg].str.startswith("MT-")))
    """
    # 8. Manual annotation plotting on remaining spots
    #fig_annot = plot_dir / f"{sample_name}_annotation.png"
    sc.pl.spatial(
        adata,
        img_key="hires",
        #color="annotation",
        title=f"{sample_name} annotation (post-QC)",
        show=False,
    )
    plt.show()
    plt.close()
    #plt.savefig(fig_annot, dpi=300, bbox_inches="tight")
    #plt.close()
    print(f"Saved annotation plot")

    # Store in dictionary
    adata_dict[sample_name] = {
        "adata": adata}


# In[5]:


# CHANGE OUTPUT DIRECTORY ACCORDINGLY
output_dir= "/dcs04/hicks/data/Jan/sim_project/sim_paper/data/real_data_qc/breast_cancer"
os.makedirs(output_dir, exist_ok=True)

for sample_name, data in adata_dict.items():
    adata = data["adata"]

    # save file (uncomment to actually write)
    output_path = os.path.join(output_dir, f"{sample_name}_qc_filtered_8_um.h5ad")
    adata.write(output_path)
    print(f"  → Would save to: {output_path}")
