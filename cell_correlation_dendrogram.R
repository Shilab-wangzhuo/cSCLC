# Assuming 's151' is your Seurat object containing deconvolved cell type proportions
# Extract proportions for relevant cell types (T cells, NSCLC, SCLC, macrophages)

s151rlist = readRDS("resList_s151.rds")
cell_type_proportions <- s151rlist$mat
colnames(cell_type_proportions) #= c("Plasma_cell","Mast","SCLC-A",
                                    "Neutrophil", "SCLC-N","T_cell",
                                    "Macrophage","B_cell", "SCLC-P",
                                    "Fibroblast","LUAD")

summary(cell_type_proportions)
#keep max >.1 & 3rd >0.05
cell_type_proportions = as.data.frame(cell_type_proportions)
cell_type_proportions = cell_type_proportions[,-c(2,8,9)]
# Calculate the correlation matrix between the cell type proportions


correlation_matrix <- cor(cell_type_proportions, 
                          method = "pearson")
print(correlation_matrix)
# Convert the correlation matrix to a distance matrix
correlation_distance <- as.dist(1 - correlation_matrix)

# Perform hierarchical clustering
hc <- hclust(correlation_distance, method = "ward.D2")
# Using base R to plot the dendrogram
plot(hc, main = "Hierarchical Clustering of Cell-Type Correlations", xlab = "Cell Types", sub = "", ylab = "Correlation Distance")

# Install the dendextend package if not installed
# install.packages("dendextend")
library(dendextend)

# Assuming 'hc' is the hierarchical clustering object obtained earlier
# Convert the hclust object to a dendrogram
dend <- as.dendrogram(hc)
dend <- rev(dend)

# Assign colors to branches
dend <- color_branches(dend, k = 3)  # k specifies the number of clusters/colors; adjust as needed


# You can also manually specify the colors if desired:
dend <- color_branches(dend, k = 3,
        col = c("#87bdd8", "#588c7e",  "#eeac99"))

pdf("s151_dendrogram.pdf",width = 3,height=3)
# Set plot parameters
par(mar = c(5, 2, 2, 10))  # Adjust margins for labels
# Plot the dendrogram with branches facing left and smaller font size
plot( dend, main = "Cell-Type Correlations", 
     ylab = "Correlation Distance", 
     horiz = TRUE,  # Make the dendrogram horizontal
     edgePar = list(lwd = 1),  # Adjust line width for visibility
     cex = 0.6)  # Adjust 'cex' for smaller font size
dev.off()
