# Load necessary libraries
suppressPackageStartupMessages({
  library(Seurat)          # Seurat for single-cell RNA-seq analysis
  library(SingleCellExperiment)  # For SingleCellExperiment format
  library(scDblFinder)     # For doublet detection
  library(celda)           # For cell contamination correction (decontX)
  library(decontX)         # DecontX for contamination correction
  library(dplyr)           # Data manipulation
  library(ggplot2)         # Plotting
})

# Set base directory for raw data files
base_path <- "./rawdata/"

# List of next batch samples to process
next_samples <- c("HRR573079", "HRR573090", "HRR573091", "HRR573093")

# Initialize an empty list to store processed Seurat objects
seurat_list <- list()

# Iterate over each sample
for (sample_name in next_samples) {
  
  # Define the path to the filtered feature matrix for the current sample
  data_dir <- file.path(base_path, sample_name, sample_name, "outs", "filtered_feature_bc_matrix")
  
  # Check if the directory exists, if not, skip the sample
  if (!dir.exists(data_dir)) {
    message(sprintf("❌ %s: Path does not exist, skipping.", sample_name))
    next
  }
  
  # Read the expression matrix using Read10X function
  message(sprintf("✅ Reading data for %s", sample_name))
  expr_mat <- Read10X(data.dir = data_dir)
  
  # Create a Seurat object with the expression matrix
  seu <- CreateSeuratObject(counts = expr_mat, project = sample_name)
  seu$orig.ident <- sample_name  # Assign the sample name as an identity
  
  # ➤ Detect doublets using scDblFinder (convert Seurat object to SingleCellExperiment)
  message(sprintf("🔍 Running scDblFinder for %s", sample_name))
  sce <- as.SingleCellExperiment(seu)
  sce <- scDblFinder(sce, dbr = 0.05)  # dbr = doublet rate threshold
  
  # ➤ Correct for contamination using decontX (recommended to run on SingleCellExperiment format)
  message(sprintf("💧 Running decontX for %s", sample_name))
  sce <- celda::decontX(sce)  # DecontX contamination correction
  
  # ➤ Extract contamination scores and add them to metadata
  contamination <- sce@metadata$decontX$contamination
  colData(sce)$decontX_contamination <- contamination
  
  # ➤ Convert back to Seurat object and store contamination score in Seurat metadata
  seu <- as.Seurat(sce)
  seu$orig.ident <- sample_name
  seu$Contamination <- seu$decontX_contamination
  
  # ➤ Save the processed Seurat object to the list
  seurat_list[[sample_name]] <- seu
}

# Print message indicating that all samples are processed
message("🎉 All samples processed successfully!")

# Save the list of Seurat objects
saveRDS(seurat_list, "seurat_list.rds")

combined_seu <- merge(seurat_list[[1]], y = seurat_list[-1], add.cell.ids = names(seurat_list))

#Quality control
combined_seu[["percent.mt"]] <- PercentageFeatureSet(combined_seu, pattern = "^MT-|^mt-")
combined_seu[["percent.rb"]] <- PercentageFeatureSet(combined_seu, pattern = "^RP[SL]|^Rp[sl]")

VlnPlot(combined_seu,features = c("nFeature_RNA","nCount_RNA","percent.mt"),group.by = "orig.ident",pt.size = 0,raster=FALSE)

combined_seu <- subset(combined_seu, subset =nFeature_RNA >= 100 & nFeature_RNA < 10000
              & nCount_RNA < 50000)

saveRDS(combined_seu,"combined_seu.rds")
