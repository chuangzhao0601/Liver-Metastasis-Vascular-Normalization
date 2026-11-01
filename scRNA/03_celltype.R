# Set working directory ------------------------------------------------------
setwd("./celltype/")

# Load necessary libraries ---------------------------------------------------
library(Seurat)        # Main package for single-cell RNA-seq analysis
library(harmony)       # Used for batch effect correction using the Harmony algorithm
library(future)        # For parallel processing support (multithreading)
library(dplyr)         # For manipulating metadata and data frames

# Set up parallel processing -------------------------------------------------
plan("multicore", workers = 10)  # Set up parallel processing using 10 CPU cores
options(future.globals.maxSize= 1000 * 1024^3)  # Increase memory limit for parallel processing

# Read data -------------------------------------------------------------------
# Load the preprocessed Seurat object containing the data for analysis
all <- readRDS("./all.rds")

# Filter the dataset to keep only cells from "Endothelial" cell type
all <- subset(all, celltype %in% c("Endothelial"))

# Perform quality control by filtering cells based on feature count, doublet score, and contamination score
all <- subset(all, subset = nFeature_RNA >= 300 & scDblFinder.score < 0.8 & Contamination < 0.8)

# Remove samples with too few cells (less than 10 cells per sample) ---------------
cell_counts <- table(all$batch)  # Count cells per sample
valid_samples <- names(cell_counts[cell_counts >= 10])  # Keep only samples with >= 10 cells

# Split dataset by batch (sample) and retain only valid samples
all_list <- SplitObject(all, split.by = "batch")
all_list <- all_list[valid_samples]  # Keep only samples with sufficient cells

# Calculate highly variable genes (HVGs) for each sample --------------------
hvg_list <- map(all_list, function(x) {
  DefaultAssay(x) <- "RNA"  # Set the RNA assay as the default
  x <- NormalizeData(x, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)  # Normalize data
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000)  # Find top 3000 variable genes
  x@assays[["RNA"]]@var.features  # Extract variable features
}) %>% unlist()  # Combine all variable features from each sample

# Select top 1000 shared highly variable genes
hvg_list <- table(hvg_list) %>% sort(decreasing = TRUE)  # Count frequency of each variable gene
gene_selection <- names(hvg_list)[1:1000]  # Select the top 1000 most variable genes

# Merge the Seurat objects from different samples ---------------------------------
all_merged <- reduce(all_list, merge, merge.data = TRUE)  # Merge samples into a single Seurat object

# Normalize the merged dataset
DefaultAssay(all_merged) <- "RNA"  # Ensure RNA assay is used
all_merged <- NormalizeData(all_merged, normalization.method = "LogNormalize", scale.factor = 10000)

# Cell cycle scoring: Assign S-phase and G2M-phase scores to each cell --------------
s.genes <- cc.genes$s.genes  # Get genes related to the S-phase
g2m.genes <- cc.genes$g2m.genes  # Get genes related to the G2M-phase

# Perform cell cycle scoring and use it for identity assignment (for observation only)
all_merged <- CellCycleScoring(all_merged, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

# Scale data and perform PCA (Principal Component Analysis) ---------------------
all_merged <- ScaleData(all_merged, features = gene_selection, vars.to.regress = c("S.Score", "G2M.Score", "nFeature_RNA"), verbose = FALSE)  # Scale data and regress out unwanted variation (e.g., cell cycle)
all_merged <- RunPCA(all_merged, features = gene_selection, npcs = 20, verbose = FALSE)  # Run PCA using the top 1000 variable genes

# Perform batch effect correction using Harmony ----------------------------------
all_merged <- RunHarmony(
  object = all_merged,
  group.by.vars = c("batch"),  # Correct for batch effects based on the "batch" variable
  plot_convergence = TRUE,  # Plot the convergence of the Harmony algorithm
  assay.use = "RNA",  # Use the RNA assay for correction
  max.iter.harmony = 10  # Limit the number of iterations to 10 for Harmony
)

# Find neighbors and clusters based on Harmony-corrected dimensions ----------------
all_merged <- FindNeighbors(all_merged, reduction = "harmony", dims = 1:20)  # Find neighbors using Harmony-corrected data
resolutions <- seq(0.4, 2.0, by = 0.1)  # Range of resolutions to test for clustering
all_merged <- FindClusters(all_merged, resolution = resolutions, algorithm = 1, verbose = FALSE)  # Perform clustering at different resolutions

# Perform UMAP (Uniform Manifold Approximation and Projection) for visualization ----
all_merged <- RunUMAP(
  object = all_merged,
  reduction = "harmony",  # Use Harmony-corrected data for UMAP
  dims = 1:20,  # Use the first 20 principal components
  reduction.name = "umap_harmony"  # Name the resulting UMAP embedding
)

# Save the final Seurat object with all results ---------------------------------
saveRDS(all_merged, "Endothelial.rds")  # Save the processed Seurat object to an RDS file
