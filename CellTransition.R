setwd("/mnt/storage/qinluo/zhuo")
# Load necessary library
library(limma)

# Load the normalized gene expression matrix (rows are genes, columns are samples)
# Replace 'your_expression_matrix.csv' with the path to your data
expression_matrix <- read.table("TLUAD_TSCLC.txt", 
                                row.names = 1)

s151 = readRDS("/mnt/storage/qinluo/zhuo/s151.rds")


# Load required libraries
library(limma)

# Prepare FPKM data and condition labels
bulk_data <- as.data.frame(expression_matrix)
condition_labels <- gsub("\\d+$", "", colnames(bulk_data))  # Extract condition names (LUAD, SCLC, etc.)

# Design matrix for differential expression
design <- model.matrix(~ 0 + factor(condition_labels))
colnames(design) <- unique(condition_labels)

# Fit linear model
fit <- lmFit(bulk_data, design)
contrast_matrix <- makeContrasts(
  LUADvsOthers = LUAD - (SCLC + T.LUAD + T.SCLC) / 3,
  SCLCvsOthers = SCLC - (LUAD + T.LUAD + T.SCLC) / 3,
  TLUADvsOthers = T.LUAD - (LUAD + SCLC + T.SCLC) / 3,
  TSCLCvsOthers = T.SCLC - (LUAD + SCLC + T.LUAD) / 3,
  levels = design
)
fit <- contrasts.fit(fit, contrast_matrix)
fit <- eBayes(fit)

marker_genes <- list(
  LUAD = rownames(topTable(fit, coef = "LUADvsOthers", number = 200, p.value = 0.05)),
  SCLC = rownames(topTable(fit, coef = "SCLCvsOthers", number = 200, p.value = 0.05)),
  T.LUAD = rownames(topTable(fit, coef = "TLUADvsOthers", number = 200, p.value = 0.05)),
  T.SCLC = rownames(topTable(fit, coef = "TSCLCvsOthers", number = 200, p.value = 0.05)))


all_genes <- unlist(marker_genes)
gene_counts <- table(all_genes)

# Keep only unique genes (appear in one marker set)
marker_genes <- lapply(marker_genes, function(genes) {
  genes[gene_counts[genes] == 1]  # Include only genes that appear once
})



# Step 2: Extract Expression Data
# Extract normalized expression matrix from Seurat object
expr_matrix <- as.matrix(s151@assays$SCT@data)

# Subset to marker genes only
marker_genes_combined <- unique(unlist(marker_gene_sets))
expr_matrix_subset <- expr_matrix[rownames(expr_matrix) %in% marker_genes_combined, ]

# Step 3: Compute Metagene Scores
# Compute average expression for each marker set (metagene) in each spot
metagene_scores <- lapply(marker_gene_sets, function(genes) {
  genes_in_data <- intersect(genes, rownames(expr_matrix_subset))
  if (length(genes_in_data) > 0) {
    colMeans(expr_matrix_subset[genes_in_data, , drop = FALSE])  # Average expression
  } else {
    rep(0, ncol(expr_matrix_subset))  # Assign 0 if no marker genes found
  }
})

# Convert to a data frame
metagene_scores_df <- as.data.frame(t(do.call(rbind, metagene_scores)))
rownames(metagene_scores_df) <- colnames(expr_matrix_subset)  # Spot names
colnames(metagene_scores_df) <- names(marker_gene_sets)       # Metagene names

# Step 4: Normalize Metagene Scores
# Normalize scores (z-score across spots)
metagene_scores_df_norm <- as.data.frame(scale(metagene_scores_df))
colnames(metagene_scores_df_norm) <- paste0("Norm_", colnames(metagene_scores_df))

# Step 5: Add Scores to Seurat Metadata
# Add raw and normalized scores to Seurat object
for (type in names(marker_gene_sets)) {
  s151 <- AddMetaData(s151, metadata = metagene_scores_df[[type]], col.name = paste0("Metagene_", type))
  s151 <- AddMetaData(s151, metadata = metagene_scores_df_norm[[paste0("Norm_", type)]], col.name = paste0("Norm_Metagene_", type))
}

# Step 6: Identify Transition Spots
# Define a co-enrichment threshold
transition_threshold <- 0.2

# Flag transition spots with co-enrichment of T.LUAD and T.SCLC
s151 <- AddMetaData(s151, metadata = metagene_scores_df_norm$Norm_T.LUAD > transition_threshold & 
                      metagene_scores_df_norm$Norm_T.SCLC > transition_threshold, 
                    col.name = "Transition")

# Step 7: Visualize Results

## Spatial Distribution of Metagene Scores
SpatialFeaturePlot(s151, features = paste0("Metagene_", names(marker_gene_sets))) +
  ggtitle("Metagene Scores (Raw)")

SpatialFeaturePlot(s151, features = paste0("Norm_Metagene_", names(marker_gene_sets))) +
  ggtitle("Metagene Scores (Normalized)")

## Visualize Transition Spots
SpatialFeaturePlot(s151, features = "Transition") +
  ggtitle("Transition Spots (T.LUAD and T.SCLC Co-Enrichment)")

## Heatmap of Normalized Scores
library(pheatmap)

# Prepare data for heatmap
heatmap_data <- as.matrix(metagene_scores_df_norm)
rownames(heatmap_data) <- colnames(expr_matrix_subset)  # Spot names

# Plot heatmap
pheatmap(heatmap_data, cluster_rows = TRUE, cluster_cols = TRUE, 
         main = "Normalized Metagene Scores Heatmap")
