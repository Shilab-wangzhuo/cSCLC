.libPaths( c( .libPaths(),c( "/usr/software/miniconda3/lib/R/library",
                             "/home/qinluo/R/x86_64-conda-linux-gnu-library/4.2",
                             "/mnt/data/qinluo/software/miniconda3/lib/R/library")) )
options(warn=-1)
library(Seurat)
setwd("/mnt/storage/qinluo/zhuo/publicData1")

library(ggplot2)
library(dplyr)
library(Seurat)
library(infercnv)
clustercol=c("#eeac99", "#77AADD", #"#b2c2bf", "#AA4455",
             "#DDCC66", 
             "#99CCBB", "#bd5734","#2288BB",  
             "#87bdd8",  "#c1946a", "#50394c","#f4e1d2", "#DD99AA", 
             "#b8a9c9", "#F5DEB3")

clustercol2=c( "orange",  "#117755","#777711",  
               "#c1946a", "#82b74b", "#034f84", 
               "#7e4a35", "#113370",  "#33BBDD", "#116666", 
               "#11AA77", "#2288BB", "#999922", "#AABB33", "#DDCC66",
               "#332222", "#AA4455", "#BB2277", "#44AA22", "#BB2244", "#DD99AA")

########plot merger
merged_seurat = readRDS("merged_seurat.rds")
merged_seurat = JoinLayers(merged_seurat,assay = "Spatial")
count_matrix = GetAssayData(merged_seurat,assay="Spatial",slot="count")

FinalAnnotationsForExport <- data.frame(
  Barcode=colnames(count_matrix),
  Histology=merged_seurat$seurat_clusters)

genes = read.table("gene_order_file.tsv")
genes = intersect(rownames(count_matrix),
                  genes[,1])

Breast_ENSBMLID_Counts = count_matrix[genes,]
Breast_ENSBMLID_Counts = as.data.frame(Breast_ENSBMLID_Counts)

write.table(Breast_ENSBMLID_Counts, "merger_spatial_counts_30k.tsv", sep = "\t")
write.table(FinalAnnotationsForExport, "merger_spatial_Annotations_30k.tsv", 
            sep = "\t",
            quote = FALSE, 
            col.names = FALSE, 
            row.names = FALSE)

options("Seurat.object.assay.version" = "v3")

Userguide_10xBreastCancer_infCNV <- infercnv::CreateInfercnvObject(
  raw_counts_matrix="merger_spatial_counts_30k.tsv", 
  gene_order_file="./gene_order_file.tsv",
  annotations_file="merger_spatial_Annotations_30k.tsv",
  delim="\t",
  ref_group_names="1")
saveRDS(Userguide_10xBreastCancer_infCNV,
        "merger_spatial_inferCNVObject")

.libPaths( c( .libPaths(),c( "/usr/software/miniconda3/lib/R/library",
                             "/home/qinluo/R/x86_64-conda-linux-gnu-library/4.2",
                             "/mnt/data/qinluo/software/miniconda3/lib/R/library")) )
options(warn=-1)
library(Seurat)
setwd("/mnt/storage/qinluo/zhuo/publicData1")

library(ggplot2)
library(dplyr)
library(Seurat)
library(infercnv)
options(expressions = 50000)
options("Seurat.object.assay.version" = "v3")



Userguide_10xBreastCancer_infCNV = infercnv::run(Userguide_10xBreastCancer_infCNV, 
                                                 cutoff=0.1, #(see infercnv::run documentation)
                                                 out_dir=paste0( "./InferCNVrun_","spatial_merger"), 
                                                 scale_data =T,  
                                                 cluster_by_groups=T, #unsupervised analysis
                                                 HMM = T, # turns gradient filtering on
                                                 HMM_type =  "i3",
                                                 window_length=150,
                                                 denoise=TRUE) #denoising applies noise reduction for the plot 







