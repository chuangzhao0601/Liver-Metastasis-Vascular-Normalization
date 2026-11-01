import tempfile
import anndata
import subprocess
import scvi
import sys
import torch
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import scanpy as sc
import seaborn as sns
import os

# Importing cell2location related modules
from cell2location.models import Cell2location, RegressionModel
from cell2location.plt import plot_spatial
from cell2location.utils import select_slide
from cell2location.utils.filtering import filter_genes

# Change working directory
os.chdir('./')
print(os.getcwd())

results_folder = './results/'
ref_run_name = f'{results_folder}/reference_signatures'
run_name = f'{results_folder}/cell2location_map'

# Function to read and preprocess Visium spatial transcriptomics data
def read_and_qc(sample_path, sample_name):
    """
    Read a single Visium spatial transcriptomics dataset, compute QC, and remove duplicate genes.
    
    :param sample_path: Path to the sample directory
    :param sample_name: Name of the sample
    :return: Preprocessed AnnData object
    """
    adata = sc.read_visium(sample_path, count_file='filtered_feature_bc_matrix.h5', load_images=True)
    adata.uns['spatial'] = {sample_name: adata.uns['spatial'].pop(list(adata.uns['spatial'].keys())[0])}
    
    print(f"Processing {sample_name}: {adata.n_vars} genes before removing duplicates")

    # Remove duplicated genes
    if adata.var_names.duplicated().any():
        adata = adata[:, ~adata.var_names.duplicated()].copy()
        print(f"After removing duplicates: {adata.n_vars} genes")
    
    adata.obs['sample'] = sample_name

    # Compute QC metrics and mitochondrial gene fraction
    sc.pp.calculate_qc_metrics(adata, inplace=True)
    adata.var['mt'] = [gene.startswith('MT-') for gene in adata.var_names]
    adata.obs['mt_frac'] = adata[:, adata.var['mt'].tolist()].X.sum(1).A.squeeze() / adata.obs['total_counts']

    adata.obs_names = [f"{sample_name}_{spot}" for spot in adata.obs_names]
    adata.obs.index.name = 'spot_id'

    return adata

# Paths to slices
slice_dirs = [
    "./slice1/",
    "./slice2/",
    "./slice3/",
]
sample_names = ["slice1", "slice2", "slice3"]

# Read and process the slices
slides = [read_and_qc(slice_dir, sample_name) for slice_dir, sample_name in zip(slice_dirs, sample_names)]

# Concatenate all slices into a single AnnData object
adata_vis = slides[0].concatenate(
    slides[1:],
    batch_key="sample",
    uns_merge="unique",
    batch_categories=sample_names,
    index_unique=None  # Spot names are already unique
)

# Ensure gene names are unique
adata_vis.var_names_make_unique()
print(f"Final combined dataset: {adata_vis.n_obs} spots, {adata_vis.n_vars} genes")

# Load the list of excluded genes
genes_df = pd.read_csv("./excludedgenes.csv", sep="\t")
genes_df.columns = genes_df.columns.str.strip()  # Clean up any spaces in column names
excluded_genes = pd.unique(genes_df.values.ravel())
excluded_genes = [g for g in excluded_genes if pd.notna(g)]

print(f"{len(excluded_genes)} genes to exclude")

# Remove duplicate genes from the spatial dataset
adata_vis = adata_vis[:, ~adata_vis.var_names.duplicated()]

# Identify genes that are both in the excluded gene list and the spatial data
excluded_genes_in_spatial = list(set(excluded_genes).intersection(adata_vis.var_names))

# Optional: Save the expression of excluded genes to obsm (for later reference)
adata_vis.obsm['excluded_genes'] = adata_vis[:, excluded_genes_in_spatial].X.toarray()

# Remove excluded genes from the dataset
adata_vis = adata_vis[:, [g for g in adata_vis.var_names if g not in excluded_genes_in_spatial]]

# Read reference data (e.g., for regression model)
adata_ref = sc.read("./adata_ref/all.h5ad")

# Save the original counts to .raw if not already done
if adata_ref.raw is None:
    adata_ref.raw = adata_ref.copy()

adata_ref.X = adata_ref.raw.X.copy()  # Use raw counts as X

# Add gene names to the 'SYMBOL' column
if 'SYMBOL' not in adata_ref.var.columns:
    adata_ref.var['SYMBOL'] = adata_ref.var_names

adata_ref.var_names_make_unique()

# Filter genes based on some criteria (e.g., minimum expression or presence in multiple cells)
selected = filter_genes(adata_ref, cell_count_cutoff=60, cell_percentage_cutoff2=0.00015, nonz_mean_cutoff=1.50)
adata_ref = adata_ref[:, selected].copy()

# Prepare the reference data for regression model setup
cell2location.models.RegressionModel.setup_anndata(adata=adata_ref, 
                        batch_key='batch',  # Specify batch key
                        labels_key='subcelltype1',  # Specify cell type labels
                        categorical_covariate_keys=['tissue']  # Specify tissue covariates
)

# Create the regression model and train it
mod = RegressionModel(adata_ref)
mod.view_anndata_setup()
mod.train(max_epochs=300)

mod.plot_history(20)

# Export the estimated cell abundance from the posterior distribution
adata_ref = mod.export_posterior(adata_ref, sample_kwargs={'num_samples': 1000, 'batch_size': 5000})

# Save the model
mod.save(f"{ref_run_name}", overwrite=True)

# Clean up _index column if present
if "_index" in adata_ref.var.columns:
    adata_ref.var = adata_ref.var.drop(columns=["_index"])

if adata_ref.raw is not None:
    raw_df = adata_ref.raw.to_adata()  # Convert raw to mutable AnnData object
    if "_index" in raw_df.var.columns:
        raw_df.var = raw_df.var.drop(columns=["_index"])
    adata_ref.raw = raw_df  # Replace raw data

# Save AnnData object
adata_file = f"{ref_run_name}/sc.h5ad"
adata_ref.write(adata_file)

# Load the saved model
adata_ref = sc.read_h5ad(adata_file)
mod = cell2location.models.RegressionModel.load(f"{ref_run_name}", adata_ref)

# Export posterior distribution
adata_ref = mod.export_posterior(adata_ref, sample_kwargs={'num_samples': 1000, 'batch_size': 5000})

# Plot QC
mod.plot_QC()

# Export quantiles for confidence intervals
adata_ref = mod.export_posterior(adata_ref, use_quantiles=True, sample_kwargs={'batch_size': 5000})

# Export estimated expression in each cluster
if 'means_per_cluster_mu_fg' in adata_ref.varm.keys():
    inf_aver = adata_ref.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{i}' 
                                    for i in adata_ref.uns['mod']['factor_names']]].copy()
else:
    inf_aver = adata_ref.var[[f'means_per_cluster_mu_fg_{i}' 
                                    for i in adata_ref.uns['mod']['factor_names']]].copy()

# Intersect genes between spatial and reference datasets
intersect = np.intersect1d(adata_vis.var_names, inf_aver.index)
adata_vis = adata_vis[:, intersect].copy()
inf_aver = inf_aver.loc[intersect, :].copy()

# Prepare AnnData for the Cell2location model
cell2location.models.Cell2location.setup_anndata(adata=adata_vis, batch_key="sample")

# Create and train the Cell2location model
mod = cell2location.models.Cell2location(
    adata_vis, 
    cell_state_df=inf_aver,  # Cell state signature (e.g., estimated cell abundance)
    N_cells_per_location=30,  # Expected number of cells per spot
    detection_alpha=200  # Hyperparameter for RNA detection normalization
)

mod.view_anndata_setup()

mod.train(max_epochs=30000, batch_size=None, train_size=1)

# Save the trained model
mod.save(f"{run_name}", overwrite=True)

# Export posterior for estimated cell abundance
adata_vis = mod.export_posterior(adata_vis, sample_kwargs={'num_samples': 1000, 'batch_size': mod.adata.n_obs})

# Save spatial data results
adata_file = f"{run_name}/sp.h5ad"
adata_vis.write(adata_file)

# Load the results
adata_vis = sc.read_h5ad(adata_file)

# Add quantile data for plotting
adata_vis.obs[adata_vis.uns['mod']['factor_names']] = adata_vis.obsm['q05_cell_abundance_w_sf']

# Select a slide and plot in spatial coordinates
slide = select_slide(adata_vis, 'ST-liver1')

# Visualize in spatial coordinates
with mpl.rc_context({'axes.facecolor':  'black', 'figure.figsize': [4.5, 5]}):
    sc.pl.spatial(slide, cmap='magma', color=["Mono_CD14", "Mph_IL1B", "LSEC", "CD8_ISG15", "Fibro_FAP"],
                  ncols=4, size=1.3, img_key='hires', vmin=0, vmax='p99.2')

# Plot spatial results for multiple cell types
clust_labels = ["Fibro_FAP", "CapEC", "cSMC"]
slide = select_slide(adata_vis, 'ST-liver4')

with mpl.rc_context({'figure.figsize': (15, 15)}):
    fig = plot_spatial(
        adata=slide,
        color=clust_labels,
        show_img=True,
        style='fast',
        max_color_quantile=0.992,
        circle_diameter=6,
        colorbar_position='right'
    )
