# Set working directory ------------------------------------------------------
setwd("./dataset/")

# Load necessary libraries ---------------------------------------------------
library(Seurat)        # Main package for single-cell RNA-seq analysis
library(harmony)       # Batch effect correction using Harmony
library(future)        # For parallel processing support
library(purrr)         # For working with lists

# Set up parallel processing -------------------------------------------------
plan("multicore", workers = 10)           # Use 10 CPU cores for parallel processing
options(future.globals.maxSize= 1000 * 1024^3)  # Increase memory size for parallel processing

# Read data -------------------------------------------------------------------
# Load the preprocessed Seurat object that contains the data for analysis
all <- readRDS("./combined_seu.rds")

# Set the default assay to RNA (for gene expression analysis)
DefaultAssay(all) <- "RNA"

# Read the list of excluded genes from a file (to remove unwanted genes later)
genes_df <- read.csv("/hwdata/home/zcshixi7/reference/scRNA/excludedgenes.csv", header = TRUE, stringsAsFactors = FALSE,sep = "\t")
excluded_genes <- unique(na.omit(unlist(genes_df)))  # Get a list of unique excluded genes

# Filter out excluded genes from the counts matrix
excluded_genes <- intersect(excluded_genes, rownames(all))  # Ensure the excluded genes are present in the dataset
counts <- all@assays$RNA@counts

# Create a filtered counts matrix excluding the unwanted genes
filtered_counts <- counts[!rownames(counts) %in% excluded_genes, ]

# Update the Seurat object with the filtered counts data
all[["RNA"]] <- CreateAssayObject(counts = filtered_counts)

# Add metadata to the Seurat object
all <- AddMetaData(all, metadata = Matrix::colSums(filtered_counts), col.name = "nCount_RNA")  # Total RNA counts per cell
all <- AddMetaData(all, metadata = Matrix::colSums(filtered_counts > 0), col.name = "nFeature_RNA")  # Number of features (genes) detected per cell

# Filter cells based on feature count and doublet score (cell quality control)
all <- subset(all, subset = nFeature_RNA >= 100 & scDblFinder.score < 0.8)

# Split the dataset by sample identity (for further analysis by sample)
all_list <- SplitObject(all, split.by = "orig.ident")

# Calculate highly variable genes (HVGs) for each sample --------------------
hvg_list <- map(all_list, function(x) {
  DefaultAssay(x) <- "RNA"
  x <- NormalizeData(x, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000)  # Identify the top 3000 HVGs
  x@assays[["RNA"]]@var.features
}) %>% unlist()  # Collect all HVGs from each sample

# Select the top 2000 shared highly variable genes
hvg_list <- table(hvg_list) %>% sort(decreasing = TRUE)  # Count the frequency of each gene being variable
gene_selection <- names(hvg_list)[1:2000]  # Select the top 2000 genes

# Merge the individual Seurat objects (samples) into one integrated object
all_merged <- reduce(all_list, merge, merge.data = TRUE)

# Normalize data for the merged object
DefaultAssay(all_merged) <- "RNA"
all_merged <- NormalizeData(all_merged, normalization.method = "LogNormalize", scale.factor = 10000)

# Perform cell cycle scoring (assign S-phase and G2M-phase scores based on cell cycle genes)
s.genes <- cc.genes$s.genes  # Get genes related to the S-phase
g2m.genes <- cc.genes$g2m.genes  # Get genes related to the G2M-phase

all_merged <- CellCycleScoring(all_merged, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)  # Score each cell for its cell cycle stage

# Scale the data and perform PCA (Principal Component Analysis) -----------------
all_merged <- ScaleData(all_merged, features = gene_selection, vars.to.regress = c("S.Score","G2M.Score","nFeature_RNA"), verbose = FALSE)  # Scale data and regress out unwanted variation
all_merged <- RunPCA(all_merged, features = gene_selection, npcs = 30, verbose = FALSE)  # Perform PCA and retain the first 30 principal components

# Perform batch effect correction using Harmony ----------------------------------
all_merged <- RunHarmony(
  object = all_merged,
  group.by.vars = c("orig.ident"),  # Correct for batch effects based on the "orig.ident" variable (e.g., samples)
  plot_convergence = TRUE,  # Plot convergence of Harmony algorithm
  assay.use = "RNA",  # Use the RNA assay for correction
  max.iter.harmony = 10  # Limit the number of iterations to 10 for Harmony
)

# Find neighbors and clusters based on Harmony-corrected dimensions ----------------
all_merged <- FindNeighbors(all_merged, reduction = "harmony", dims = 1:30)  # Find neighbors using Harmony-corrected data
resolutions <- seq(0.4, 1.6, by = 0.2)  # Range of resolutions for cluster identification
all_merged <- FindClusters(all_merged, resolution = resolutions, algorithm = 1, verbose = FALSE)  # Perform clustering at different resolutions

# Perform UMAP (Uniform Manifold Approximation and Projection) for visualization ----
all_merged <- RunUMAP(
  object = all_merged,
  reduction = "harmony",  # Use Harmony-corrected data for UMAP
  dims = 1:30,  # Use the first 30 principal components
  reduction.name = "umap_harmony"  # Name the resulting UMAP embedding
)

# Save the final Seurat object with all results ---------------------------------
saveRDS(all_merged, "all_merged.rds")  # Save the processed and analyzed Seurat object to a file
