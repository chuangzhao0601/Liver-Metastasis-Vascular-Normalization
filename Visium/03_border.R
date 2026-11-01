# Load necessary libraries for spatial transcriptomics and single-cell analysis
library(Seurat)         # For single-cell RNA-seq analysis and data processing
library(ggplot2)        # For plotting gene expression trends
library(dplyr)          # For data manipulation and summarization

# Load the spatial transcriptomics data for "lvmet1" sample from a directory
lvmet1 <- Load10X_Spatial(
  data.dir = "./Gaoqiang/ST-liver1/",   # Directory containing the filtered 10X feature matrix
  filename = "filtered_feature_bc_matrix.h5",   # File containing filtered feature-barcode matrix
  assay = "Spatial",   # Specify the assay type as "Spatial" for spatial transcriptomics data
  slice = "ST_liver1"  # Name of the slice in the dataset
)

# Read the barcode metadata file to assign tissue labels for each barcode in lvmet1
barcode_lvmet1 <- read.table("./Gaoqiang/ST-liver1.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
lvmet1$tissue <- barcode_lvmet1$tissue   # Assign the 'tissue' information to the Seurat object

# Visualize the spatial distribution of the samples based on tissue type
SpatialDimPlot(lvmet1, group.by = "tissue", pt.size.factor = 1, alpha = 0.6)

# Add sample labels for each object to distinguish them later
lvmet1$sample <- "lvmet1"
lvmet4$sample <- "lvmet4"
HCC1L$sample <- "HCC1L"
HCC5B$sample <- "HCC5B"
HCC5C$sample <- "HCC5C"
HCC5D$sample <- "HCC5D"

# Combine multiple Seurat objects into one dataset for integrated analysis
combined <- merge(lvmet1, y = list(lvmet4, HCC1L, HCC5B, HCC5C, HCC5D))

# Perform SCTransform normalization on the combined dataset, which standardizes the data
combined <- SCTransform(combined, assay = "Spatial", verbose = FALSE)

# Perform PCA to reduce dimensionality and capture the most important variation in the data
combined <- RunPCA(combined)

# Identify the neighbors of each cell based on PCA dimensions (1 to 30)
combined <- FindNeighbors(combined, dims = 1:30)

# Perform clustering based on PCA and find cell clusters
combined <- FindClusters(combined, resolution = 0.5)

# Run UMAP for dimensionality reduction and visualization in two dimensions
combined <- RunUMAP(combined, dims = 1:30)

# Function to plot gene expression trends across tissue regions and samples
plot_gene_trend_combined <- function(seu, gene, 
                                     tissue_order = c("normal", "border", "tumor"),
                                     slice_col = "slice",
                                     tissue_col = "tissue") {
  # Check if the gene exists in the Seurat object, if not, stop the function
  if (!(gene %in% rownames(seu))) {
    stop(paste("Gene", gene, "not found in Seurat object"))
  }
  
  # Extract expression data and metadata for the specified gene and columns
  expr_df <- FetchData(seu, vars = c(gene, tissue_col, slice_col))
  colnames(expr_df) <- c("expression", "tissue", "sample")  # Rename columns for clarity
  
  # Set the order of tissue types for plotting
  expr_df$tissue <- factor(expr_df$tissue, levels = tissue_order)
  
  # Group the data by sample and tissue, then calculate the mean expression for each group
  mean_expr_df <- expr_df %>%
    group_by(sample, tissue) %>%
    summarise(mean_expr = mean(expression, na.rm = TRUE)) %>%
    ungroup()   # Remove the grouping for further operations
  
  # Plot the gene expression trend using a line plot
  p <- ggplot(mean_expr_df, aes(x = tissue, y = mean_expr, group = sample, color = sample)) +
    geom_line(size = 1) +     # Add lines connecting the points
    geom_point(size = 2) +    # Add points at each data point
    scale_color_manual(values = use_colors) +   # Define custom colors for each sample
    theme_classic() +         # Apply a classic theme for the plot
    labs(
      title = paste0("Expression trend of ", gene),  # Set the title of the plot
      x = "Tissue region",    # Label for the x-axis
      y = paste0("Mean expression of ", gene)   # Label for the y-axis
    ) +
    theme(
      text = element_text(size = 14),  # Set the font size for the text
      legend.title = element_blank(),  # Remove legend title
      legend.position = "right"    # Place the legend on the right side
    )
  
  # Return the plot
  return(p)
}

# Plot the gene expression trends for specific genes across tissue regions
plot_gene_trend_combined(combined, gene = "CCL2", tissue_col = "tissue", slice_col = "sample")
plot_gene_trend_combined(combined, gene = "CCR2", tissue_col = "tissue", slice_col = "sample")
plot_gene_trend_combined(combined, gene = "TNF", tissue_col = "tissue", slice_col = "sample")
plot_gene_trend_combined(combined, gene = "CXCL8", tissue_col = "tissue", slice_col = "sample")
