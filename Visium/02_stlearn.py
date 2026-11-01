# ================================
# Load required libraries
# ================================
import os
import random
import warnings
import numpy as np
import pandas as pd
import scanpy as sc
import stlearn as st

# stLearn plotting parameters
st.settings.set_figure_params(dpi=120)

# -------------------------------
# Fix all random seeds for reproducibility
# -------------------------------
seed = 0
np.random.seed(seed)
random.seed(seed)
os.environ["PYTHONHASHSEED"] = str(seed)

# Ignore warnings for cleaner output
warnings.filterwarnings("ignore")

# ================================
# Function: Load Visium data + tissue annotation
# ================================
def load_visium_with_tissue(path, count_file="filtered_feature_bc_matrix.h5",
                            tissue_file=None, load_images=True):
    """
    Load 10X Visium spatial transcriptomics data and attach tissue-region annotation.

    Purpose:
    --------
    1. Read the standard Visium gene expression + spatial coordinates.
    2. If available, load user-provided tissue annotations (normal/border/tumor).
    3. Ensure the 'use_quality' field exists for stLearn image extraction.
    """

    # 1. Load Visium data
    adata = sc.read_visium(
        path=path,
        count_file=count_file,
        load_images=load_images
    )

    # 2. Attach tissue region annotation to adata.obs
    if tissue_file is not None:
        tissue_df = pd.read_csv(tissue_file, sep="\t")
        tissue_df = tissue_df.set_index("barcode")

        # Align annotation order with adata.obs_names
        adata.obs["tissue"] = tissue_df.loc[adata.obs_names, "tissue"].values

    # 3. Get the first library_id (required by Visium format)
    library_id = list(adata.uns["spatial"].keys())[0]

    # 4. Ensure "use_quality" exists (needed for downstream stLearn image operations)
    if "use_quality" not in addata.uns["spatial"][library_id]:
        adata.uns["spatial"][library_id]["use_quality"] = "hires"
        print(f"✅ Added 'use_quality' = 'hires' to library_id: {library_id}")

    return adata


# Load lvmet1 sample
lvmet1 = load_visium_with_tissue(
    path="./Gaoqiang/ST-liver1/",
    tissue_file="./Gaoqiang/ST-liver1.txt"
)


# ================================
# Function: stLearn preprocessing + SME normalization
# ================================
def process_visium_sample(adata, out_dir=None, min_cells=1,
                          pca_comps=50, save_name=None):
    """
    Full preprocessing pipeline for a Visium sample using stLearn.

    Purpose:
    --------
    This function performs the following biological and computational tasks:

    1. Gene filtering (removes rarely detected genes).
    2. Library-size normalization + log1p transformation.
    3. Prepare spatial coordinates (imagerow/imagecol).
    4. Generate image tiles & extract vision features (Deep Learning features).
    5. PCA for dimensionality reduction.
    6. SME normalization (spatial morphological expression normalization),
       which integrates histology with gene expression to improve expression accuracy.
    7. Save the processed AnnData object for downstream analyses.
    """

    # Create output directory if needed
    if out_dir is not None:
        os.makedirs(out_dir, exist_ok=True)
        os.chdir(out_dir)

    # Ensure gene names are unique
    if not adata.var_names.is_unique:
        print(f"⚠️ Duplicate gene names found: {adata.var_names.duplicated().sum()}")
        adata = adata[:, ~adata.var_names.duplicated()].copy()

    # -------------------------------
    # Step 1: Filter genes + normalize + log transform
    # -------------------------------
    st.pp.filter_genes(adata, min_cells=min_cells)
    st.pp.normalize_total(adata)
    st.pp.log1p(adata)

    # -------------------------------
    # Step 2: Add spatial coordinates
    # -------------------------------
    if "spatial" in adata.obsm:
        adata.obs["imagerow"] = adata.obsm["spatial"][:, 1]
        adata.obs["imagecol"] = adata.obsm["spatial"][:, 0]
    else:
        print("⚠️ Spatial coordinates missing in adata.obsm")

    # Ensure use_quality exists
    library_id = list(adata.uns["spatial"].keys())[0]
    if "use_quality" not in adata.uns["spatial"][library_id]:
        adata.uns["spatial"][library_id]["use_quality"] = "hires"

    # -------------------------------
    # Step 3: stLearn image tiling + feature extraction
    # -------------------------------
    st.pp.tiling(adata, out_path="tiling")
    st.pp.extract_feature(adata)

    # -------------------------------
    # Step 4: PCA (before SME)
    # -------------------------------
    st.em.run_pca(adata, n_comps=pca_comps)

    # -------------------------------
    # Step 5: SME normalization
    # -------------------------------
    data_SME = adata.copy()
    st.spatial.SME.SME_normalize(data_SME, use_data="raw")
    data_SME.X = data_SME.obsm["raw_SME_normalized"]

    # PCA again after SME normalization
    st.em.run_pca(data_SME, n_comps=pca_comps)

    # -------------------------------
    # Step 6: Save processed data
    # -------------------------------
    if save_name is None:
        save_name = "sample_processed.h5ad"

    data_SME.write(save_name)
    print(f"✅ Saved SME-normalized data: {save_name}")

    return data_SME


# ================================
# Process all samples
# ================================
samples = {
    "lvmet1": lvmet1,
    "lvmet4": lvmet4,
    "HCC1L": HCC1L,
    "HCC5B": HCC5B,
    "HCC5C": HCC5C,
    "HCC5D": HCC5D
}

for name, adata in samples.items():
    out_dir = f"/hwdata/home/zcshixi7/project/liver/data/ST_processed/{name}"
    process_visium_sample(
        adata,
        out_dir=out_dir,
        save_name=f"{name}_SME.h5ad"
    )

