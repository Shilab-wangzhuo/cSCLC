

library(Seurat)
setwd("~/projects/csclc")
matrix_combined = readRDS("matrix_combined.rds")
ga = matrix_combined$gene[which(matrix_combined$cluster=="SCLC-A")] 
gn = matrix_combined$gene[which(matrix_combined$cluster=="SCLC-N")] 
gp = matrix_combined$gene[which(matrix_combined$cluster=="SCLC-P")] 
gb = matrix_combined$gene[c(which(matrix_combined$cluster=="LUAD"),
                            which(matrix_combined$cluster=="LUSC"))] 
merger = readRDS("merger3.rds")
sall = subset(merger,subset=data=="sall")
sall = FindVariableFeatures(sall,nfeatures=5000)
genes = VariableFeatures(sall)

rm(sall)
ga = intersect(ga, genes)
gn = intersect(gn, genes)
gp = intersect(gp, genes)
gb = intersect(gb, genes)
library(AUCell)
expr_matrix <- GetAssayData(merger,assay = "RNA",
                            slot="count")

expr_matrix = expr_matrix# Create gene set list
library(GSEABase)
ga <-  GeneSet(ga, setName="A")
gn <-  GeneSet(gn, setName="N")
gp <-  GeneSet(gp, setName="P")
gb <-  GeneSet(gb, setName="NSCLC")
# Run AUCell to calculate enrichment scores
cells_rankings <- AUCell_buildRankings(expr_matrix, plotStats = F)

#saveRDS(cells_rankings,"cells_rankings.rs")

#cell_rankings = readRDS("cells_rankings.rds")


set.seed(333)
auc_scores1 <- AUCell_calcAUC(ga, cells_rankings)
auc_scores2 <- AUCell_calcAUC(gn, cells_rankings)
auc_scores3 <- AUCell_calcAUC(gp, cells_rankings)
auc_scores4 <- AUCell_calcAUC(gb, cells_rankings)

merger$A =t(auc_scores1@assays@data$AUC)
merger$N =t(auc_scores2@assays@data$AUC)
merger$P =t(auc_scores3@assays@data$AUC)
merger$NSCLC =t(auc_scores4@assays@data$AUC)
#FeaturePlot(merger,features = c("A","N","P"))



merger2 = subset(merger,subset=seurat_clusters!=7&
                   seurat_clusters!=5&
                   seurat_clusters!=2)
merger2$seurat_clusters = as.character(merger2$seurat_clusters)
merger2$seurat_clusters[which(merger2$cell_type_fine=="SCLC-A")]="A"
merger2$seurat_clusters[which(merger2$cell_type_fine=="SCLC-N")]="N"
merger2$seurat_clusters[which(merger2$cell_type_fine=="SCLC-P")]="P"
merger2$seurat_clusters[which(merger2$cell_type_fine=="NSCLC")]="NSCLC"
merger2 = subset(merger2,subset=cell=="SCLC-A"|
                   cell=="SCLC-N"|
                   cell=="SCLC-P"|
                   cell=="LUAD"|
                   cell=="LUSC"|
                   cell=="sall")
                 
#merger2 = subset(merger2,subset=seurat_clusters!=3&
#                   seurat_clusters!=9&
#                   seurat_clusters!=8)
VlnPlot(merger2,features = c("A","N","P","NSCLC"),pt.size = 0,group.by = 
        "seurat_clusters")
saveRDS(merger2,"merger2.rds")
dt =  data.frame(A=merger2$A,N=merger2$N,P=merger2$P,NSCLC=merger2$NSCLC,
                  cell=merger2$seurat_clusters,
                 sample = merger2$batch)
dt$sample[-grep("matrix",dt$sample)]="ANP"
                                               
rownames(dt) = colnames(merger2)
#dt = dt[-which(dt$cell=="Neuroendocrine"),]
dt = dt[-which(dt$cell=="NSCLC"),]

dt$A = (dt$A-min(dt$A))/(max(dt$A) - min(dt$A))
dt$N = (dt$N-min(dt$N))/(max(dt$N)- min(dt$N))
dt$P = (dt$P-min(dt$P))/(max(dt$P)- min(dt$P))
#dt$NSCLC = (dt$NSCLC-min(dt$NSCLC))/(max(dt$NSCLC)- min(dt$NSCLC))
library(dplyr)
library(tibble)
# Convert row names to a column
dt <- dt %>%
  rownames_to_column(var = "rowname")

dt_new 

# Subset the dataframe
dt_sub <- dt %>%
  group_by(cell) %>%  # Group by cell type
  group_modify(~ {
    n_rows <- nrow(.x)  # Number of rows in the group
    if (n_rows <= 1000) {
      return(.x)  # Return all rows if less than or equal to 1000
    } else {
      return(slice_sample(.x, n = 1000))  # Sample 1000 rows if more than 1000
    }
  }) %>%
  ungroup()  # Remove grouping

# Restore row names
dt_sub <- dt_sub %>%
  column_to_rownames(var = "rowname")

# Check the result
head(dt_sub)

dt_sub = as.data.frame(dt_sub)
library(pheatmap)

dt_sub = dt_sub[which(dt_sub$sample=="ANP"),]

# Choose colors for cell types
annotation_colors <- list(
  Cell_Type = c(
    "0" = "orange", "1" = "#117755",
    "3" = "#2288BB", "4" = "#82b74b",
    "6"="#98c1d9",
    "A" = "#db2b39", "N" = "#29335c",
    "P"="#f9d5e5" ))

dt_sub$cell = factor(dt_sub$cell,levels = c("A","N","P","NSCLC",#"Neuroendocrine",
                                            "0","1","3","4","6"))
dt_sub = dt_sub[order(dt_sub$cell),]
#dt_sub = dt_sub[-which(rowSums(dt_sub[,2:5])==0),]
# Create a data frame for row annotations (cell types)
row_annotation <- data.frame(Cell_Type = dt_sub$cell)
rownames(row_annotation) <- rownames(dt_sub)


p=pheatmap(
  dt_sub[,2:5],
  annotation_row =row_annotation,
  color = colorRampPalette(c("#2874A6", "white","#C0392B"))(100)[-c(1:10,50:60)],  # Heatmap color gradient
  cluster_rows = F,  # Cluster rows
  cluster_cols = F,  # Cluster columns
  scale = "column",show_rownames = F,annotation_colors = annotation_colors)
pdf("ANP_NSCLC_plasticity_heatmap.pdf",height = 4,width = 4)
# Generate the heatmap
print(p)
dev.off()



# Install and load required libraries
library(ggtern)

# Normalize A, N, P values so that A + N + P = 1
dt_sub$A_norm <- dt_sub$A / (dt_sub$A + dt_sub$N + dt_sub$P)
dt_sub$N_norm <- dt_sub$N / (dt_sub$A + dt_sub$N + dt_sub$P)
dt_sub$P_norm <- dt_sub$P / (dt_sub$A + dt_sub$N + dt_sub$P)

dt_sub$batch = "CSCLC"
# Create the ternary plot
dt_sub$batch[which(dt_sub$cell=="A"|dt_sub$cell=="N"|dt_sub$cell=="P")] = "SCLC"
#dt_sub$batch[which(dt_sub$cell=="NSCLC")] = "NSCLC"
#dt_sub$batch[which(dt_sub$cell=="Neuroendocrine")] = "NeuroEndocrine"

dt_sub2 = dt_sub
dt_sub2$cell = as.character(dt_sub2$cell)
dt_sub2$cell[which(dt_sub2$batch=="CSCLC")] = "CSCLC"



p=ggtern(data = dt_sub, aes(x = A_norm, y = N_norm, z = P_norm, color = cell)) +
  geom_point(size =0, alpha = 0) +  # Add points with size and transparency
  scale_color_manual(values = c(
    "#db2b39", "#29335c","#fcba04",  "#98c1d9", 
    "#b1cbbb", "#c94c4c", "#ffcccb", "#66cdaa", "#ffa500"  # Customize more colors if needed
  )) +  # Custom colors for cell types
  labs(
    title = "Ternary Plot of A, N, P Values",
    x = "A", y = "N", z = "P", color = "Cell Type"
  ) +
  theme_classic() +  # Minimal theme
  theme(
    legend.position = "right",          # Place legend on the right
    axis.title = element_text(size = 12),  # Customize axis title size
    axis.text = element_text(size = 10)    # Customize axis text size
  )+#facet_wrap(~batch)+  
  stat_density_tern(aes(fill = ..level..), 
                    geom = "polygon", alpha = 0.4,bins=50)
+
  scale_fill_gradientn(
    colors = rev(color_palette))


p=ggtern(data = dt_sub2, aes(x = A_norm, y = N_norm, z = P_norm, color = cell)) +
  geom_point(size =0, alpha = 0) +  # Add points with size and transparency
  scale_color_manual(values = c(
    "#db2b39", "#29335c","#fcba04",  "#98c1d9", 
    "#b1cbbb", "#c94c4c", "#ffcccb", "#66cdaa", "#ffa500"  # Customize more colors if needed
  )) +  # Custom colors for cell types
  labs(
    title = "Ternary Plot of A, N, P Values",
    x = "A", y = "N", z = "P", color = "Cell Type"
  ) +
  theme_classic() +  # Minimal theme
  theme(
    legend.position = "right",          # Place legend on the right
    axis.title = element_text(size = 12),  # Customize axis title size
    axis.text = element_text(size = 10)    # Customize axis text size
  )+#facet_wrap(~batch)+  
  stat_density_tern(aes(fill = ..level..), 
                    geom = "polygon", alpha = 0.4,bins=50)


dt_new = dt_sub2[intersect(rownames(dt_sub2),cells),]


p = ggtern(data = dt_sub2, aes(x = A_norm, y = N_norm, z = P_norm, color = cell)) +
  geom_point(size = 0, alpha = 0) +  # Transparent points to plot density only
  scale_color_manual(values = c(
    "#db2b39", "#29335c", "#fcba04", "#98c1d9", 
    "#b1cbbb", "#c94c4c", "#ffcccb", "#66cdaa", "#ffa500"  # Custom colors
  )) +
  labs(
    title = "Ternary Plot of A, N, P Values",
    x = "A", y = "N", z = "P", color = "Cell Type"
  ) +
  theme_classic() +
  theme(
    legend.position = "right",
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  ) +
  stat_density_tern(aes(fill = after_stat(level)), 
                    geom = "polygon", alpha = 0.4, bins = 50)pdf("ANP_plasticity2.pdf",height = 4,width = 4)
# Generate the heatmap
print(p)
dev.off()


dt_sub3 = dt_sub2[which(dt_sub2$cell=="A"|dt_sub2$cell=="N"|
                          dt_sub2$cell=="P"),]

library(RColorBrewer)
library(viridis)
viridis_colors <- viridis(100)  # Generate 100 colors
subset_colors <- viridis_colors[20:80]  # Middle 60% (index 21 to 80)
pdf("ANP_by_reference.pdf",width=5,height = 5)
ggtern(data = dt_sub3, aes(x = A_norm, y = N_norm, z = P_norm, 
                             color = cell)) +
  geom_point(size =.5, alpha = .7) +  # Add points with size and transparency
  scale_color_manual(values = c(
    "#db2b39",  "#ffa500",  "#98c1d9"  # Customize more colors if needed
  )) +  stat_density_tern(aes(fill = log10(..level..)), 
                         geom = "polygon", alpha = .3,bins=50)+
  scale_fill_gradientn(
    colors = subset_colors)+theme_minimal()
dev.off()




# Subset the dataframe
dt_sub4 <- dt %>%
  group_by(sample) %>%  # Group by cell type
  group_modify(~ {
    n_rows <- nrow(.x)  # Number of rows in the group
    if (n_rows <= 1000) {
      return(.x)  # Return all rows if less than or equal to 1000
    } else {
      return(slice_sample(.x, n = 1000))  # Sample 1000 rows if more than 1000
    }
  }) %>%
  ungroup()  # Remove grouping

# Restore row names
dt_sub4 <- dt_sub4 %>%
  column_to_rownames(var = "rowname")



dt_sub4 = as.data.frame(dt_sub4)
# Normalize A, N, P values so that A + N + P = 1
dt_sub4$A_norm <- dt_sub4$A / (dt_sub4$A + dt_sub4$N + dt_sub4$P)
dt_sub4$N_norm <- dt_sub4$N / (dt_sub4$A + dt_sub4$N + dt_sub4$P)
dt_sub4$P_norm <- dt_sub4$P / (dt_sub4$A + dt_sub4$N + dt_sub4$P)

dt_sub4 = dt_sub4[-which(dt_sub4$sample=="ANP"),]


dt_sub4$cell = "CSCLC"
pdf("ANP_by_sample.pdf",width=15,height = 5)
ggtern(data = dt_sub4, aes(x = A_norm, y = N_norm, z = P_norm)) +
  geom_point(size =.5, alpha = .7,color="black") + 
  stat_density_tern(aes(fill = log10(..level..)), 
                         geom = "polygon", alpha = .1,bins=50)+
  #scale_color_manual(values="grey")+
  scale_fill_gradientn(
    colors = subset_colors)+
    theme_minimal()+facet_wrap(~sample,ncol=4)
dev.off()



subset_colors1 <- viridis(100)[20:80]  # Middle 60% (index 21 to 80)
A=ggtern(data = dt_sub3, aes(x = A_norm, y = N_norm, z = P_norm)) +
 stat_density_tern(data = subset(dt_sub3, cell == "A"), 
                    aes(x = A_norm, y = N_norm, z = P_norm, fill = log10(..level..)), 
                    geom = "polygon", alpha = 0.5, bins = 50) +
  scale_fill_viridis(option = "magma", name = "Density_A") +
  scale_fill_gradientn(colors = subset_colors1)+theme_minimal()

subset_colors2 <- magma(100)[20:80]
N=ggtern(data = dt_sub3, aes(x = A_norm, y = N_norm, z = P_norm)) +
  stat_density_tern(data = subset(dt_sub3, cell == "N"), 
                    aes(x = A_norm, y = N_norm, z = P_norm, fill = log10(..level..)), 
                    geom = "polygon", alpha = 0.5, bins = 50) +
  scale_fill_viridis(option = "viridis", name = "Density_N") +
  scale_fill_gradientn(colors = subset_colors2)+theme_minimal()

subset_colors3 <- cividis(100)[80:100]
P=ggtern(data = dt_sub3, aes(x = A_norm, y = N_norm, z = P_norm)) +
  stat_density_tern(data = subset(dt_sub3, cell == "P"), 
                    aes(x = A_norm, y = N_norm, z = P_norm, fill = log10(..level..)), 
                    geom = "polygon", alpha = 0.5, bins = 30) +
  scale_fill_gradientn(colors = subset_colors3)+theme_minimal()

z= A+N+P
pdf("ANP_by_ref.pdf",width=25,height = 9)
print(z)
dev.off()



#####################


dt =  data.frame(A=merger2$A,N=merger2$N,P=merger2$P,NSCLC=merger2$NSCLC,
                 cell=merger2$seurat_clusters,
                 sample = merger2$batch)
dt$sample[-grep("matrix",dt$sample)]="ANP"

rownames(dt) = colnames(merger2)
dt = dt[-which(dt$cell=="Neuroendocrine"),]
#dt = dt[-which(dt$cell=="NSCLC"),]

dt$SCLC = rowMax(as.matrix(dt[,1:3]))
dt$SCLC = (dt$SCLC-min(dt$SCLC))/(max(dt$SCLC) - min(dt$SCLC))
dt$NSCLC = (dt$NSCLC-min(dt$NSCLC))/(max(dt$NSCLC)- min(dt$NSCLC))
dt$sample[which(dt$sample=="ANP")] = dt$cell[which(dt$sample=="ANP")]

library(dplyr)
library(tibble)
# Convert row names to a column
dt <- dt %>%
  rownames_to_column(var = "rowname")

# Subset the dataframe
dt_sub <- dt %>%
  group_by(sample) %>%  # Group by cell type
  group_modify(~ {
    n_rows <- nrow(.x)  # Number of rows in the group
    if (n_rows <= 1000) {
      return(.x)  # Return all rows if less than or equal to 1000
    } else {
      return(slice_sample(.x, n = 1000))  # Sample 1000 rows if more than 1000
    }
  }) %>%
  ungroup()  # Remove grouping

# Restore row names
dt_sub <- dt_sub %>%
  column_to_rownames(var = "rowname")

# Check the result
head(dt_sub)

dt_sub = as.data.frame(dt_sub)
dt_sub = dt_sub[-which(dt_sub$sample==0|dt_sub$sample==1|
                         dt_sub$sample==11|dt_sub$sample==4),]





