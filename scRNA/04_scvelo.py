import scvelo as scv
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc
import igraph
import loompy as lmp
import anndata
import os

# === 1. Set scVelo settings ===
# Set verbosity for scVelo. Verbosity 3 will show hints, which can be helpful for debugging.
scv.settings.verbosity = 3  

# Set scVelo to presenter view mode (adjust figure size)
scv.settings.presenter_view = True  

# Set figure parameters to use scVelo's recommended style for visualization
scv.settings.set_figure_params('scvelo')  

# Create a directory to store figures if it does not exist
fig_dir = "./figures1/"
os.makedirs(fig_dir, exist_ok=True)

# === 2. Read data ===
# Load the main Seurat-like AnnData object containing single-cell RNA-seq data
all = sc.read_h5ad("./scvelo/Endo.h5ad")

# Load the list of cell IDs to retain in the dataset
sample_obs = pd.read_csv(
    './cellID_obs.csv',
    header=None
)[0].values

# Filter the original dataset to include only the cells specified in sample_obs
adata = all[np.isin(all.obs.index, sample_obs)]  

# === 3. Read additional metadata files ===
# Read UMAP coordinates for cells (used for visualization)
cell_umap = pd.read_csv(
    './cell_embeddings.csv',
    header=0,
    names=["Cell ID", "UMAP_1", "UMAP_2"]
)

# Read cell type annotations for each cell
cell_celltype = pd.read_csv(
    './cell_celltype.csv',
    header=0,
    names=["Cell ID", "celltype"]
)

# === 4. Process and create index for cell annotations ===
# Create a DataFrame with cell IDs from the filtered data (adata)
sample_one_index = pd.DataFrame(adata.obs.index, columns=["Cell ID"])

# Merge the cell index with UMAP coordinates based on the Cell ID
umap_ordered = sample_one_index.merge(cell_umap, on="Cell ID")

# Retain only the UMAP coordinates (UMAP_1, UMAP_2) for further use
umap_ordered = umap_ordered.iloc[:, 1:]  

# Assign UMAP coordinates to the AnnData object (adata) in the obsm attribute for visualization
adata.obsm["X_umap"] = umap_ordered.values  

# Merge the cell index with cell type information based on Cell ID
celltype_ordered = sample_one_index.merge(cell_celltype, on="Cell ID")

# Assign cell types to the AnnData object (adata)
adata.obs['celltype'] = celltype_ordered["celltype"].values

# === 5. Ensure unique gene names ===
# Make sure gene names are unique across the dataset, which is required for analysis
adata.var_names_make_unique()

# === 6. Preprocessing the data ===
# Filter and normalize the data: keep genes with at least 50 counts and select top 500 genes
scv.pp.filter_and_normalize(adata, min_shared_counts=50, n_top_genes=500)

# Compute moments (first and second order statistics) required for scVelo analysis
scv.pp.moments(adata, n_pcs=50, n_neighbors=50)

# === 7. Fit dynamical model ===
# Recover the dynamics of gene expression for each cell. This is the core of scVelo's dynamical model.
scv.tl.recover_dynamics(adata, n_jobs=24)  # Use 24 threads for faster computation

# === 8. Calculate velocity and velocity graph using dynamical model ===
# Compute RNA velocity based on the recovered dynamics. This represents the directionality of transcriptional change.
scv.tl.velocity(adata, mode="dynamical")

# Compute the velocity graph, which captures cell-to-cell velocity relationships.
scv.tl.velocity_graph(adata)

# === 9. Visualization: Velocity embedding with cell types ===
# Choose a color palette for cell types
ident_colours = ["#9E0142","#D73027","#F46D43","#5E4FA2","#B3DE69","#41AE76"]

# Visualize the RNA velocity stream on UMAP, colored by cell type
scv.pl.velocity_embedding_stream(
    adata, 
    basis='umap',  # Use UMAP for embedding
    color="celltype",  # Color by cell type
    palette=ident_colours,  # Set the color palette for different cell types
    size=20,  # Set the size of points
    alpha=1,  # Set opacity
    fontsize=0.5,  # Adjust font size
    save=os.path.join(fig_dir, 'embedding_grid_dynamical.pdf')  # Save the plot to the specified directory
)

# === 10. Calculate and visualize latent time ===
# Latent time represents the continuous progression of cells along their trajectory
scv.tl.recover_latent_time(adata)  # Compute latent time based on the dynamical model

# Visualize latent time on UMAP, coloring by latent time
scv.pl.velocity_embedding_stream(
    adata,
    basis="umap",  # Use UMAP for embedding
    color="latent_time",  # Color by latent time
    title='Latent Time (dynamical)',  # Title for the plot
    fontsize=20,  # Font size
    legend_fontsize=20,  # Font size for the legend
    min_mass=2,  # Minimum mass threshold for plotting
    color_map="plasma",  # Choose a color map
    save=os.path.join(fig_dir, 'latent_time_dynamical.pdf')  # Save the plot
)
