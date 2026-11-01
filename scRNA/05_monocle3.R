# Set working directory
setwd("./monocle3/")

# Load necessary libraries
library(Seurat)        # Seurat for single-cell RNA-seq analysis
library(dplyr)         # Data manipulation
library(monocle3)      # Monocle3 for trajectory analysis
library(ggplot2)       # Plotting
library(RSpectra)      # For fast singular value decomposition
library(patchwork)     # For combining ggplot objects
library(Matrix)        # For sparse matrix operations
library(irlba)         # For randomized SVD

# Read in the Seurat object
Endo<- readRDS('/all.rds')

# Extract the expression matrix and metadata from the Seurat object
expr_matrix <- as(as.matrix(Endo@assays$RNA@counts), 'sparseMatrix')
cell_metadata <- Endo@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(Endo))

rownames(gene_annotation) <- rownames(Endo)

# Create a Monocle3 CellDataSet object
cds <- new_cell_data_set(expr_matrix,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)

# Preprocess the data: normalize, scale, and perform PCA
cds <- preprocess_cds(cds, num_dim = 20)  # Equivalent to Seurat's NormalizeData + ScaleData + RunPCA
plot_pc_variance_explained(cds)  # Plot the variance explained by each principal component

# Perform dimensionality reduction (PCA)
cds <- reduce_dimension(cds, preprocess_method = "PCA")

# Integrate UMAP embedding from Seurat (harmonized) into Monocle3
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(Endo, reduction = "umap_harmony")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed

# Perform clustering and learn trajectory graph in Monocle3
cds <- cluster_cells(cds, reduction_method = "UMAP")
cds <- learn_graph(cds)

# Save the processed CellDataSet object
saveRDS(cds, "monocle3.rds")

# Read the saved object (optional, just for continuity)
cds <- readRDS("./monocle3.rds")

# Plot cells colored by 'subcelltype1' and show trajectory labels
p <- plot_cells(cds, color_cells_by = "subcelltype1", label_groups_by_cluster = FALSE,
                label_leaves = FALSE, label_branch_points = FALSE)
p

# Plot with additional labels for leaves and branch points of the trajectory
plot_cells(cds,
           color_cells_by = "subcelltype1",
           label_cell_groups = FALSE,
           label_leaves = TRUE,
           label_branch_points = TRUE,
           graph_label_size = 5, group_label_size = 6,
           cell_size = 0.5)
# Black circles represent branch points; light gray circles show different outcomes along the trajectory.

# Function to get the earliest principal node (starting point) for ordering cells by pseudotime
get_earliest_principal_node <- function(cds, celltype = "LSEC") {
  cell_ids <- which(colData(cds)[, "subcelltype1"] == celltype)
  
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  
  root_pr_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                                             (which.max(table(closest_vertex[cell_ids,]))))]
  
  return(root_pr_nodes)
}

# Order cells along the trajectory based on the earliest principal node
cds = order_cells(cds, root_pr_nodes = get_earliest_principal_node(cds))

# Plot cells along pseudotime
p <- plot_cells(cds,
                 color_cells_by = "pseudotime",
                 label_cell_groups = FALSE,
                 label_leaves = FALSE,
                 label_branch_points = FALSE,
                 graph_label_size = 1.5, cell_size = 0.2)
p
ggsave("panec_monocle3.svg", plot = p, device = "svg", height = 8, width = 9, dpi = 36)

# Plot specific genes along pseudotime, colored by celltype
monocle3::plot_genes_in_pseudotime(cds[c("LYVE1", "STAB2", "ITGA6", "LAMB1"),], color_cells_by = "subcelltype1", ncol = 2)

# Custom color scheme for cell types
use_colors <- c(LSEC = "#7FBC41", VenEC = "#08519C", CapEC = "#FEE08B", ArtEC = "#D73027", EndoMT = "#FCBBA1", LEC = "#A6CEE3")

# Select a subset of cells for visualization (30% of the cells)
subset_cells <- sample(colnames(cds), size = round(ncol(cds) * 0.3))
cds_subset <- cds[, subset_cells]

# Plot genes in pseudotime with custom colors
p <- monocle3::plot_genes_in_pseudotime(cds_subset[c("LYVE1", "STAB2", "ITGA6", "LAMB1"),], color_cells_by = "subcelltype1") +
  scale_color_manual(values = use_colors)  # Apply custom colors to the plot
p
ggsave("pseudotime_genes.svg", plot = p, device = "svg", height = 9, width = 7)

# Custom gene list for pseudotime plotting
gene_list <- c("LYVE1", "STAB2")

# Subset of cells for visualization (10% of the cells)
subset_cells <- sample(colnames(cds), size = round(ncol(cds) * 0.1))
cds_subset <- cds[, subset_cells]

# Plot genes in pseudotime with facetting and custom colors
p <- monocle3::plot_genes_in_pseudotime(cds_subset[gene_list, ], color_cells_by = "subcelltype1") +
  facet_wrap(~ factor(gene_short_name, levels = gene_list), scales = "free_y") +
  scale_color_manual(values = use_colors)  # Apply custom colors to the plot
p
ggsave("pseudotime_genes_facet.svg", plot = p, device = "svg", height = 7, width = 8)
