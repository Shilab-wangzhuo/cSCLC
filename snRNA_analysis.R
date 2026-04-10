# Load necessary libraries
library(Seurat)
library(ggplot2)
library(dplyr)
sall <- readRDS("/mnt/storage/qinluo/zhuo/sall.rds")

# Extract unique batch information
batches <- unique(sall@meta.data$batch)

# Initialize an empty list to store results for each batch
all_results <- list()
# Run your analysis here (replace with your specific analysis steps)
finalref = readRDS("/mnt/storage/qinluo/zhuo/publicData1/finalref_lcnce.rds")
#finalref = readRDS("finalref.rds")
finalref$cell_type [which(finalref$cell_type=="epithelial cell")]="LCNCE"

finalref = subset(finalref,subset=cell_type!="Endothelial"&
                    cell_type!="DC")

finalref$cell_type[which(finalref$cell_type=="LUAD")] = "NSCLC"
finalref$cell_type[which(finalref$cell_type=="LUSC")] = "NSCLC"
finalref$cell_type[which(finalref$cell_type=="LCNCE")] = "NSCLC"

# Perform SCT normalization on the reference object
finalref <- SCTransform(finalref, verbose = FALSE)
clustercol=c("#eeac99", "#b8a9c9", "orange", "#77AADD", 
             "#117755",  
             "#ffef96", "#99CCBB", "#c1946a", "#82b74b",  
             "#034f84","#b8a9c9", "#F5DEB3", "#7e4a35","#778822") 

sall = NormalizeData(sall)
#sall = FindVariableFeatures(sall,nfeatures=3000)
#sall = ScaleData(sall)
# Loop over each batch, subset the data, and run the analysis
for (x in batches) {
  print(x)
  # Subset the Seurat object for the current batch
  batch_data <- subset(sall, subset = batch == x)
  
 # batch_data= SCTransform(batch_data)
  # Re-run the FindTransferAnchors function
  anchors <- FindTransferAnchors(
    reference = finalref,
    query = batch_data,
    normalization.method = "SCT"
  )
  batch_data = NormalizeData(batch_data)
  batch_data = FindVariableFeatures(batch_data,nfeatures=3000)
  batch_data = ScaleData(batch_data)
  batch_data = RunPCA(batch_data)
  
  batch_data <- TransferData(
    anchorset = anchors,
    query=batch_data,
    refdata = finalref$cell_type, # Replace "cell_type" with the correct metadata column name
    prediction.assay = TRUE,
    weight.reduction = batch_data[["pca"]],
    dims = 1:20 # Adjust the number of dimensions as needed
  )
  
 
 # batch_data = subset(batch_data,subset=predicted.id.score>.7)
  # Store the results in the list
  all_results[[x]] <- batch_data@meta.data
}

# Combine all the metadata into a single data frame
combined_metadata <- do.call(rbind, all_results)
saveRDS(combined_metadata,"combined_metadata.rds")
# Load the stringr package
library(stringr)

combined_metadata = combined_metadata[which(combined_metadata$predicted.id.score>.5),]
# Create a summary table of cell composition across all batches

combined_metadata$batch = str_extract(combined_metadata$batch, "\\d{3}(?=_matrix)")

cell_composition <- combined_metadata %>%
  group_by(batch, predicted.id) %>%
  summarise(count = n()) %>%
  ungroup()

colors <- c("Fibroblast"="#3e4444", 
            "Macrophage"="red",
            "Neutrophil"="#8DA0CB",
            "NSCLC"="#82b74b",
            "Plasma cell"="orange", 
            "SCLC-A"="#c1946a",
            "SCLC-N"="#1F78B4", 
            "SCLC-P"="#E5C494", 
            "T cell"="#7e4a35",
            "B cell"="royalblue")

# Plot the bar plot
ggplot(cell_composition, aes(x = batch, y = count, 
                             fill = predicted.id)) +
  geom_bar(stat = "identity", position = "fill") +
  labs(title = "Cell Composition Across All Batches",
       x = "Sample",
       y = "Proportion",
       fill = "Cell Type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_fill_manual(values= colors)

#########################################################

# Combine all the metadata into a single data frame
combined_metadata <- do.call(rbind, all_results)
combined_metadata = combined_metadata[which(combined_metadata$predicted.id.score>.5),]

for (batch in batches2) {
  print(batch)
  # Subset the Seurat object for the current batch
  batch_data <- subset(sall, subset = batch == batch)
  
  cells = rownames(combined_metadata)[which(combined_metadata$batch==batch)]
  cells =  sub(".*\\.", "", cells)
  batch_data <- subset(batch_data, cells = cells)
  
  batch_data$cell = combined_metadata$predicted.id[which(combined_metadata$batch==batch)]
  batch_data = FindVariableFeatures(batch_data,nfeatures = 3000)
  batch_data = subset(batch_data,subset=cell!="Plasma cell")
  batch_data = RunPCA(batch_data)
  batch_data = RunUMAP(batch_data,dims=1:30, n.neighbors =10 )
  
  p1=DimPlot(batch_data,group.by = "cell",cols=colors)+coord_fixed()+
    theme_void()
  
  dc = as.data.frame(table(batch_data$cell))
  dc = dc[order(dc$Freq,decreasing = T),]
  cell_rm = as.character(dc$Var1[1])
  batch_data2 = subset(batch_data,subset=cell!=cell_rm)
  if(dim(batch_data2)[2]<30) {
    pdf(paste0(batch,"_Dimplot.pdf"),width = 8,height = 4)
    print(p1)
    dev.off()
  }else{
    batch_data2 = RunPCA(batch_data2)
    batch_data2 = RunUMAP(batch_data2,dims=1:30, n.neighbors = 10)
    p2=DimPlot(batch_data2,group.by = "cell",cols=colors)+coord_fixed()+
    theme_void()
    pdf(paste0(batch,"_Dimplot.pdf"),width = 8,height = 4)
    print(p1+p2)
    dev.off()
  }
  
  
  batch_data3 <- subset(batch_data, subset = cell=="SCLC-N"| cell=="SCLC-A"|
                         cell=="SCLC-P")
  cells_rankings <- AUCell_buildRankings(batch_data@assays$RNA@counts, 
                                         nCores=1, plotStats=F)

  cells_AUC1 <- AUCell_calcAUC(g2,cells_rankings)
  auc = cells_AUC1@assays@data$AUC
  batch_data3$g2AUC = auc
  
  cells_AUC2 <- AUCell_calcAUC(g3,cells_rankings)
  auc2 = cells_AUC2@assays@data$AUC
  batch_data3$g3AUC = auc2
  
  cells_AUC3<- AUCell_calcAUC(g5,cells_rankings)
  auc3 = cells_AUC3@assays@data$AUC
  batch_data3$g5AUC = auc3
  
  batch_data3$TSCLC_LUAD= batch_data3$g2AUC/batch_data3$g3AUC 
  batch_data3$TSCLC_SCLC = batch_data3$g5AUC/batch_data3$g3AUC
  dt = data.frame(TSCLC_LUAD= batch_data3$g2AUC/batch_data3$g3AUC ,
                  TSCLC_SCLC=batch_data3$g5AUC/batch_data3$g3AUC ,
                  Cell=colnames(batch_data3))
  # Reshape the data for plotting
  library(reshape2)
  dt = dt[order(dt$TSCLC_LUAD), ]
  data_melted <- reshape2::melt(dt, id.vars = "Cell", 
                                variable.name = "Metric",
                                value.name = "Value")
  library(ggplot2)
  # Create a barplot
  p=ggplot(data_melted, aes(x = reorder(Cell, Value), 
                          y = Value, color = Metric,
                          group = Metric)) +
    geom_line() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + # Hide x-axis labels for clarity
    labs(title = batch, 
         x = "Cell (Sorted)", y = "AUC Value")+ylim(c(0,1))+
    scale_color_manual(values = c("orange","black"))
  pdf(paste0(batch,"_AUC1.pdf"),width = 8,height = 4)
  print(p)
  dev.off()
  pp=  FeaturePlot(batch_data3,features = c("TSCLC_LUAD","TSCLC_SCLC"),cols=c("white","#d64161"))+
    coord_fixed()
  pdf(paste0(batch,"_AUC2.pdf"),width = 8,height = 4)
  print(pp)
  dev.off()
    
}


g2 <- readRDS("/mnt/storage/qinluo/zhuo/TSCLC>LUAD.rds")
#g3 <- readRDS("/mnt/storage/qinluo/zhuo/SCLC_>_LUAD.rds")
g5  <- readRDS("/mnt/storage/qinluo/zhuo/TSCLC_>_SCLC.rds")


###
library(AUCell)
batch_data3 = subset(finalref,subset=cell_type=="SCLC-P"|
                       cell_type=="SCLC-N"| cell_type=="SCLC-P")
counts = as.matrix(batch_data3@assays$SCT@counts)
rownames(counts) = rownames(batch_data3)
colnames(counts) = colnames(batch_data3)

cells_rankings <- AUCell_buildRankings( counts, 
                                       nCores=1, plotStats=F)
cells_AUC1 <- AUCell_calcAUC(g2,cells_rankings)
auc = cells_AUC1@assays@data$AUC
batch_data3$g2AUC = auc

cells_AUC2 <- AUCell_calcAUC(g3,cells_rankings)
auc2 = cells_AUC2@assays@data$AUC
batch_data3$g3AUC = auc2

cells_AUC5 <- AUCell_calcAUC(g5,cells_rankings)
auc5= cells_AUC5@assays@data$AUC
batch_data3$g5AUC = auc3
dt = data.frame(TSCLC_LUAD= batch_data3$g2AUC/batch_data3$g3AUC ,
                TSCLC_SCLC=batch_data3$g5AUC/batch_data3$g3AUC,
                Cell=colnames(batch_data3))
# Reshape the data for plotting
library(reshape2)
dt = dt[order(dt$TSCLC_LUAD), ]
data_melted <- reshape2::melt(dt, id.vars = "Cell", 
                              variable.name = "Metric",
                              value.name = "Value")
library(ggplot2)
# Create a barplot
pfinalref = ggplot(data_melted, aes(x = reorder(Cell, Value), 
                        y = Value, color = Metric,
                        group = Metric)) +
  geom_line() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + # Hide x-axis labels for clarity
  labs(#title = "Line Plot of g2AUC and g3AUC Sorted by g2AUC", 
    x = "Cell (Sorted)", y = "AUC Value")+ylim(c(0,1))+
  scale_color_manual(values = c("orange","black"))

pdf("sn_s151_Signautre.pdf",width=5,height = 2)
ggplot(data_melted, aes(x = reorder(Cell, Value), y = Value,
                        fill = Metric, group = Metric)) +
  geom_area(alpha = 0.6) + # Use geom_area to fill under the lines
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1), # Add the frame
    panel.background = element_blank() # Make the background transparent if desired
  ) +
  labs(#title = "Area Plot of g2AUC and g3AUC Sorted by g2AUC",
    x = "Cell (Sorted)", y = "AUC Value")+
  scale_fill_manual(values = c("orange","black","grey"))+ylim(c(0,.2))
dev.off()
#


######
