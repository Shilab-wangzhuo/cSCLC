#finalref = readRDS("finalref_lcnce.rds")
finalref = readRDS("finalref.rds")
#finalref = readRDS("finalref.rds")
finalref$cell_type [which(finalref$cell_type=="epithelial cell")]="LCNCE"

finalref = subset(finalref,subset=cell_type!="Endothelial"&
                    cell_type!="DC")

finalref$cell_type[which(finalref$cell_type=="LUAD")] = "NSCLC"
finalref$cell_type[which(finalref$cell_type=="LUSC")] = "NSCLC"
finalref$cell_type[which(finalref$cell_type=="LCNCE")] = "NSCLC"

source("./codes/GetuniqueDegs.R")

deg = matrix_combined
counts =  as.data.frame(merged_seurat@assays$SCT@data)
genes = intersect(rownames(counts),deg$gene)
deg = deg[genes,]
table(deg$cluster)


deg = deg[c(deg$gene[which(deg$cluster=="SCLC-A")],
        deg$gene[which(deg$cluster=="SCLC-N")],
         deg$gene[which(deg$cluster=="SCLC-P")],
   deg$gene[which(deg$cluster=="NSCLC")],
  #deg$gene[which(deg$cluster=="LUSC")],
  #deg$gene[which(deg$cluster=="LCNCE")],
  deg$gene[which(deg$cluster=="T cell")],
  deg$gene[which(deg$cluster=="B cell")],
  deg$gene[which(deg$cluster=="Macrophage")],
  deg$gene[which(deg$cluster=="Neutrophil")],
  deg$gene[which(deg$cluster=="Mast")],
  deg$gene[which(deg$cluster=="Plasma cell")],
  deg$gene[which(deg$cluster=="Fibroblast")]),]
merged_seurat = JoinLayers(merged_seurat,assay = "Spatial")
DefaultAssay(merged_seurat) ="Spatial"
counts = merged_seurat@assays$Spatial@layers$counts

#merged_seurat somehow does not work with findmarkers, so we made a new merged_seurat2
rownames(counts) = rownames(merged_seurat)
colnames(counts) =colnames(merged_seurat)
merged_seurat2  = CreateSeuratObject(counts)
merged_seurat2@active.ident = merged_seurat@active.ident
merged_seurat2$seurat_clusters= merged_seurat$seurat_clusters
merged_seurat2 = NormalizeData(merged_seurat2)
merged_seurat2 = ScaleData(merged_seurat2,features=rownames(deg))
deg2 = FindAllMarkers(merged_seurat2,features= rownames(deg))

# refine the deg genes by using the genes show up in the cluster markers
deg = deg[unique(deg2$gene),]


finalref$cell_type = factor(finalref$cell_type,levels=c(
 "SCLC-A","SCLC-N","SCLC-P","NSCLC",#"LUAD","LUSC","LCNCE",
 "T cell","B cell","Macrophage","Neutrophil","Mast",
 "Plasma cell","Fibroblast"))
deg$cluster = factor(deg$cluster ,levels=c(
  "SCLC-A","SCLC-N","SCLC-P","NSCLC",#"LUAD","LUSC","LCNCE",
  "T cell","B cell","Macrophage","Neutrophil","Mast",
  "Plasma cell","Fibroblast"))
deg = deg[order(deg$cluster),]
library(dplyr)
# Calculate average expression per cluster in finalref
ref_clusters <- AggregateExpression(finalref,
                                    group.by = "cell_type",        
                                    return.seurat = TRUE)
ref_clusters$cell_type = factor(ref_clusters$cell_type,levels=c(
  "SCLC-A","SCLC-N","SCLC-P","NSCLC",#"LUAD","LUSC","LCNCE",
  "T cell","B cell","Macrophage","Neutrophil","Mast",
  "Plasma cell","Fibroblast"))
ref_cluster_means <- as.data.frame(ref_clusters@assays$RNA@layers$data)
rownames(ref_cluster_means) = rownames(ref_clusters)
colnames(ref_cluster_means) = colnames(ref_clusters)


DoHeatmap(ref_clusters,features = rownames(deg),
          group.by="cell_type")

# Calculate average expression per cluster in merged_seurat
query_clusters <- AverageExpression(merged_seurat,
                  group.by = "seurat_clusters",        
                          return.seurat = TRUE,assay="SCT")
query_clusters = ScaleData(query_clusters,features=rownames(deg))
DoHeatmap(query_clusters,features = rownames(deg))

query_cluster_means <- as.data.frame(query_clusters@assays$SCT@layers$data)
rownames(query_cluster_means) = rownames(query_clusters)
colnames(query_cluster_means) = colnames(query_clusters)

#cosine similarity
library(lsa)

# Ensure that both data frames have the same genes
ref_cluster_means <- ref_cluster_means[rownames(deg), ]
query_cluster_means <- query_cluster_means[rownames(deg), ]
ref_cluster_means = scale(ref_cluster_means)
query_cluster_means = scale(query_cluster_means)
# Initialize an empty matrix to store cosine similarity values
cosine_similarities <- matrix(0, nrow = ncol(query_cluster_means), ncol = ncol(ref_cluster_means),
                              dimnames = list(colnames(query_cluster_means), colnames(ref_cluster_means)))

# Loop through each cluster (column) in query_cluster_means
for (i in 1:ncol(query_cluster_means)) {
  query_vector <- query_cluster_means[, i]  # Get the expression profile for the i-th query cluster
  
  # Loop through each cluster (column) in ref_cluster_means
  for (j in 1:ncol(ref_cluster_means)) {
    ref_vector <- ref_cluster_means[, j]  # Get the expression profile for the j-th reference cluster
    
    # Calculate cosine similarity between the query and reference vectors
    cosine_similarities[i, j] <- cosine(query_vector, ref_vector)
  }
}

# Display the result
View(cosine_similarities)

#calculate jensen-shannon
library(phangorn)  
ref_cluster_means <- as.data.frame(ref_clusters@assays$RNA@layers$data)
rownames(ref_cluster_means) = rownames(ref_clusters)
colnames(ref_cluster_means) = colnames(ref_clusters)

query_cluster_means <- as.data.frame(query_clusters@assays$SCT@layers$data)
rownames(query_cluster_means) = rownames(query_clusters)
colnames(query_cluster_means) = colnames(query_clusters)

ref_cluster_means <- ref_cluster_means[rownames(deg), ]
query_cluster_means <- query_cluster_means[rownames(deg), ]
# Function to calculate Jensen-Shannon divergence
js_divergence <- function(p, q) {
  # Convert p and q into probability distributions using prop.table
  p <- prop.table(p)
  q <- prop.table(q)
  
  # Replace any 0s with a very small value to avoid log(0) issues
  p[p == 0] <- 1e-10
  q[q == 0] <- 1e-10
  
  # Calculate the average distribution
  m <- 0.5 * (p + q)
  
  # Calculate the JS divergence
  js <- 0.5 * (sum(p * log(p / m)) + sum(q * log(q / m)))
  return(js)
}

# Calculate JSD for each cluster pair
js_results <- matrix(0, nrow = ncol(query_cluster_means), ncol = ncol(ref_cluster_means),
                     dimnames = list(colnames(query_cluster_means), colnames(ref_cluster_means)))
for (i in 1:ncol(query_cluster_means)) {
  for (j in 1:ncol(ref_cluster_means)) {
    js_results[i, j] <- js_divergence(query_cluster_means[, i],
                                      ref_cluster_means[, j])
  }
}
#visualize result
library(ggplot2)
library(reshape2)

# Heatmap for cosine similarity
cosine_df <- melt(cosine_similarities)
ggplot(cosine_df, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "blue", high = "red") +
  labs(title = "Cosine Similarity Heatmap")

# Heatmap for Jensen-Shannon Divergence
js_df <- melt(1-js_results)
ggplot(js_df, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "blue", high = "red") +
  labs(title = "Jensen-Shannon Divergence Heatmap")


####
# Normalize JSD to a 0-1 scale (assuming max JSD = log(2))
max_jsd <- max(js_results)
normalized_js <- 1 - (js_results / max_jsd)  # Inverting to make higher values more similar

#
# Combine metrics with weights
weight_cosine <- 0.7 # Weight for cosine similarity
weight_js <- 0.3     # Weight for Jensen-Shannon divergence

combined_similarity <- weight_cosine * cosine_similarities + 
  weight_js * normalized_js

# Visualize the combined similarity
library(pheatmap)
rownames(combined_similarity) = paste0("cluster_",0:9)
#remove neutrophil and mast cells because none has high similarity to this
#combined_similarity = combined_similarity[,-c(10,11)]

combined_similarity = combined_similarity[c(2,5,3,9,6,10,7,8,1,4),]
pheatmap(t(combined_similarity), color = colorRampPalette(c("royalblue",
            "white", "red"))(50),cluster_rows = F,cluster_cols = F,
         main = "Combined Similarity Heatmap",scale="row")


deg = deg[order(deg$avg_log2FC,decreasing = T),]
deg2 = deg2[which(deg2$avg_log2FC>0),]
deg2 = deg2[order(deg2$p_val_adj,decreasing = F),]
deg2 = deg2[deg2$pct.1>=0.5,]

genesplot5 = deg2$gene[which(deg2$cluster==1)]
#genesplot5 = intersect(genesplot5,deg$gene[which(deg$cluster=="Fibroblast")])
fibroblast_markers <- c("VIM", "ACTA2", "PDGFRA", "PDGFRB", "COL1A1", 
                        "COL1A2", "FAP", "THY1", "SPARC", "S100A4", 
                        "FN1", "DCN", "LOX", "MMP2", "MMP9", 
                        "MMP14", "TNC", "CDH2", "ITGB1", "LUM")
genesplot5 = intersect(fibroblast_markers,genesplot5)
DotPlot(merged_seurat2,features = genesplot5)+coord_flip()

genesplot6 = deg2$gene[which(deg2$cluster==4)]
#genesplot6 = intersect(genesplot6,deg$gene[which(deg$cluster=="T cell")])[1:10]
t_cell_markers <- c("CD3E", "CD3D", "CD3G", "CD4", "CD8A", 
                    "CD8B", "CD2", "CD5", "CD7", "CD28", 
                    "CD27", "CD45 (PTPRC)", "CCR7", "LCK", "FOXP3", 
                    "GZMA", "GZMB", "IFNG", "IL2", "CTLA4")
genesplot6 = intersect(t_cell_markers,genesplot6)
DotPlot(merged_seurat2,features = genesplot6)+coord_flip()

genesplot7 = deg2$gene[which(deg2$cluster==4)]
#genesplot7 = intersect(genesplot7,deg$gene[which(deg$cluster=="Macrophage")])[1:10]
macrophage_markers <- c("CD68", "CD14", "CSF1R", "ITGAM", "ITGAX",
                        "CD163", "CD80", "CD86", "MRC1", "TLR4", 
                        "FCGR1A", "HLA-DRA", "ARG1", "NOS2", 
                        "IL1B", "TNF", "CXCL8", "CCL2", "MARCO", "SIRPA")
genesplot7 = intersect(macrophage_markers,genesplot7)
DotPlot(merged_seurat2,features = genesplot7)+coord_flip()

genesplot8 = deg2$gene[which(deg2$cluster==4)]
neutrophil_markers <- c("CD66b", "CD11b", "CD15", "CSF3R", "LY6G",
                        "MPO", "ELANE", "PRTN3", "LCN2", "S100A8", 
                        "S100A9", "FCGR3B", "CEACAM8", "CXCR1", "CXCR2", 
                        "ITGB2", "RETNLG", "OLFM4", "LTF", "DEFA4")
genesplot8 = intersect(neutrophil_markers,genesplot8)
genesplot8 = intersect(genesplot8,deg$gene[which(deg$cluster=="Neutrophil")])
genesplot8 = c("S100A8","CSF3R","NCF1","TLR4","MNDA")
DotPlot(merged_seurat2,features = genesplot8)+coord_flip()

genesplot10 = deg2$gene[which(deg2$cluster==7)]
genesplot10 = intersect(genesplot10,deg$gene[which(deg$cluster=="Plasma cell")])

genesplot10 = intersect(plasma_cell_markers,genesplot10)
DotPlot(merged_seurat2,features = genesplot10)+coord_flip()

genesplot12 = deg2$gene[which(deg2$cluster==7)]
genesplot12 = intersect(genesplot12,deg$gene[which(deg$cluster=="Mast")])[1:10]
DotPlot(merged_seurat2,features = genesplot12)+coord_flip()

genesplot4 = deg2$gene[which(deg2$cluster==9)]
genesplot4 = intersect(genesplot4,deg$gene[which(deg$cluster=="SCLC-P")])[1:10]
DotPlot(merged_seurat2,features = genesplot4)+coord_flip()

genesplot11 = deg2$gene[which(deg2$cluster==9)]
b_cell_markers <- c("CD19", "CD20", "CD22","HLA-DR","CD40",
                    "CD79A", "CD79B", "CD30", "CD21", "CD27", 
                    "CD95", "CD1d", "CD21", "CD5",
                    "MS4A1", "PAX5", "BLNK", "CD21", "CD24", 
                    "CD40", "CR2", "IGHM", "IGHD", "CXCR5", 
                    "BCL6", "SPIB", "CD72", "FCRL5", "TNFRSF13C")
genesplot11 = intersect(b_cell_markers,genesplot11)
#genesplot11 = intersect(genesplot11,deg$gene[which(deg$cluster=="B cell")])
DotPlot(merged_seurat2,features = genesplot11)+coord_flip()

genesplot9 = deg2$gene[which(deg2$cluster==5)]
genesplot9 = intersect(genesplot9,deg$gene[which(deg$cluster=="LUAD")])[1:10]
DotPlot(merged_seurat2,features = genesplot9)+coord_flip()


genesplot2 = deg2$gene[which(deg2$cluster==6)]
genesplot2 = intersect(genesplot2,deg$gene[which(deg$cluster=="SCLC-N")])[1:10]
DotPlot(merged_seurat2,features = genesplot2)+coord_flip()

genesplot3 = deg2$gene[which(deg2$cluster==8)]
genesplot3 = intersect(genesplot3,deg$gene[which(deg$cluster=="SCLC-N")])[1:10]
DotPlot(merged_seurat2,features = genesplot3)+coord_flip()


genesplot1 = deg2$gene[which(deg2$cluster==0)]
#genesplot1 = intersect(genesplot1,deg$gene[which(deg$cluster=="SCLC-A")])[1:10]
DotPlot(merged_seurat2,features = genesplot1[1:3])+coord_flip()




merged_seurat2$seurat_clusters = factor(merged_seurat2$seurat_clusters,
                              levels = c(1,4,2,7,3,9,5,6,8,0))

DotPlot(merged_seurat2,features = c(genesplot5[1:3], genesplot6[1:3],
                                    genesplot7[1:3],
                                    genesplot8[1:3], genesplot10[1:3],
                                    genesplot12[1:3],
                                    genesplot4[1:3],genesplot11[1:3],
                                    genesplot9[1:3],genesplot2[1:3],
                                   genesplot3[1:3],genesplot1[1:3]
                                    ),
        group.by="seurat_clusters")+
  coord_flip()


library(ggplot2)
library(RColorBrewer)

# Create a custom palette
custom_palette <- colorRampPalette(brewer.pal(9, "Blues"))(100)  # Adjust the palette as needed

pdf("coembedding_dotplot.pdf",width=5,height=8)
# Generate the DotPlot
DotPlot(merged_seurat2, 
                    features = c(genesplot5[1:3], genesplot6[1:3],
                                 genesplot7[1:3], genesplot8[1:3],
                                 genesplot10[1:3], genesplot12[1:3],
                                 genesplot4[1:3], genesplot11[1:3],
                                 genesplot9[1:3], genesplot2[1:3],
                                 genesplot3[1:3], genesplot1[1:3]),
                    group.by = "seurat_clusters") + 
  coord_flip() +  # Flip coordinates
  scale_color_gradientn(colors = custom_palette) +  # Set color limits and custom palette
  theme(
    axis.text.y = element_text(size = 10, face = "italic"),  # Smaller and italic font for gene names
    axis.title.x = element_text(size = 12),  # Adjust x-axis title size
    axis.title.y = element_text(size = 12),  # Adjust y-axis title size
    plot.title = element_text(size = 14, face = "bold")  # Adjust plot title size
  ) + 
  ggtitle("Dot Plot of Gene Expression by Cluster")  # Optional title

dev.off()

