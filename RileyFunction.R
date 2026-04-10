
.libPaths( c( .libPaths(),c( "/usr/software/miniconda3/lib/R/library",
                             "/home/qinluo/R/x86_64-conda-linux-gnu-library/4.2",
                             "/mnt/data/qinluo/software/miniconda3/lib/R/library")) )
options(warn=-1)
setwd("/mnt/storage/qinluo/zhuo/")
library(spatstat.explore)
library(spatstat.geom)
library(dplyr)
library(Seurat)

s076 = readRDS("s076.rds")
s173 = readRDS("s173.rds")
s396 = readRDS("s396.rds")
s354 = readRDS("./S354/s354.rds")
s151 = readRDS("/mnt/storage/qinluo/zhuo/s151.rds")
s074 = readRDS("./S6-074/s074.rds")
spaL = list(s151,s354,s396,
            s173,s076,s074)


setwd("/mnt/storage/qinluo/zhuo/publicData1")
#s354res = readRDS("resList_s151.rds")
s354res = readRDS("resList_s354.rds")


spa = spaL[[2]]
spa@meta.data = cbind(spa@meta.data,s354res$mat)
# Extract spatialObjecttial coordinates from the Seurat object
spatial_coords <- as.data.frame(spa@images$sampleid@
                                  coordinates[, c("row", "col")])
spa$SCLC <- spa$`SCLC-A` + spa$`SCLC-N`
spa$I = spa$`T cell` + spa$Neutrophil+ spa$`Plasma cell`+spa$Macrophage
# Combine spatial coordinates with cell type proportions
spatial_data <- cbind(spatial_coords, spa$LUSC,spa$SCLC,spa$I,
                      spa$Fibroblast )
colnames(spatial_data) = c("imagerow","imagecol",
                           "lusc_proportion","sclc_proportion",
                           "infiltrated_cell_proportion",
                           "Fibroblast")
# Subset the spots where LUAD or SCLC is present
# Define thresholds for "presence" (e.g., proportion > 0.2)

lusc_spots <- spatial_data %>% filter(lusc_proportion > 0.1)
sclc_spots <- spatial_data %>% filter(sclc_proportion > 0.7392)

i_spots <- spatial_data %>% filter(infiltrated_cell_proportion > 0.3)
fibro_spots <- spatial_data %>% filter(Fibroblast > 0.3)
#plasma_spots <- spatial_data %>% filter(Plasma > 0.3)

# Define the window (observation area) for Ripley's analysis
xrange <- range(spatial_data$imagerow)
yrange <- range(spatial_data$imagecol)
window <- owin(xrange, yrange)


# Create point patterns for LUAD and SCLC spots
luad_ppp <- ppp(lusc_spots$imagerow, lusc_spots$imagecol, window = window)
sclc_ppp <- ppp(sclc_spots$imagerow, sclc_spots$imagecol, window = window)
t_ppp <- ppp(i_spots$imagerow, i_spots$imagecol, window = window)
fibro_ppp <- ppp(fibro_spots$imagerow, fibro_spots$imagecol, window = window)

# Calculate Ripley's K-function for LUAD
K_luad <- Kest(luad_ppp, correction = "Ripley")
# Calculate Ripley's K-function for SCLC
K_sclc <- Kest(sclc_ppp, correction = "Ripley")
k_t <- Kest(t_ppp, correction = "Ripley")
k_fibro <- Kest(fibro_ppp, correction = "Ripley")

# Create a multi-type point pattern by combining LUAD and SCLC
multi_ppp <- superimpose(LUAD = luad_ppp, SCLC = sclc_ppp, W = window)
luad_t_ppp <- superimpose(LUAD = luad_ppp, Immune_cells = t_ppp, W = window)
sclc_t_ppp <- superimpose(SCLC = sclc_ppp, Immune_cells = t_ppp, W = window)

luad_fibro_ppp <- superimpose(LUAD = luad_ppp, Fibroblast = fibro_ppp, W = window)
sclc_fibro_ppp <- superimpose(SCLC = sclc_ppp, Fibroblast = fibro_ppp, W = window)


# Calculate the cross K-function to see how LUAD and SCLC interact spatially
K_cross <- Kcross(multi_ppp, i = "LUAD", j = "SCLC", 
                  correction = "Ripley")
K_cross_luad_t <- Kcross(luad_t_ppp, i = "LUAD", j = "Immune_cells", 
                         correction = "Ripley")
K_cross_sclc_t <- Kcross(sclc_t_ppp, i = "SCLC", j = "Immune_cells", 
                         correction = "Ripley")
K_cross_luad_fibro <- Kcross(luad_fibro_ppp, i = "LUAD", j = "Fibroblast", 
                             correction = "Ripley")
K_cross_sclc_fibro<- Kcross(sclc_fibro_ppp, i = "SCLC", j = "Fibroblast", 
                            correction = "Ripley")





# Extract r, theo and K(r) values for each dataset
r_sclc_t <- K_cross_sclc_t$r # only need one r value since all the r values are the same
theo_sclc_t <- K_cross_sclc_t$theo #also only needs one theo (background value)
#k values
K_r_sclc_t <- K_cross_sclc_t$iso  # The 'iso' correction, replace if you used another correction
K_r_luad_t <- K_cross_luad_t$iso  # The 'iso' correction, replace if you used another correction
K_r_luad_fibro <- K_cross_luad_fibro$iso  # The 'iso' correction
K_r_sclc_fibro <- K_cross_sclc_fibro$iso  # The 'iso' correction
K_r_luad_sclc <- K_cross$iso   # The 'iso' correction


pdf("P11_Cell_Pattern_all.pdf",width=4,height=4)
# Create an empty plot with appropriate labels and ranges
plot(r_sclc_t, type = "l", col = "blue", lty = 1, lwd = 0, 
     xlab = "Distance (r)", ylab = "K(r)", 
     main = "P11 Cell Pattern (Ripley's K-Cross Function)", 
     xlim = range(r_sclc_t),
     ylim = range(K_r_sclc_t, K_r_luad_fibro, K_r_sclc_fibro))

# Add the lines
lines(r_sclc_t, theo_sclc_t, col = "black", lty = 3, lwd = 1.5)

lines(r_sclc_t, K_r_luad_t, col = "#c83349", lty = 1, lwd = 1.5)
lines(r_sclc_t, K_r_sclc_t, col = "#f9d5e5", lty =1, lwd = 1.5)
lines(r_sclc_t, K_r_luad_fibro, col = "#034f84", lty = 1, lwd = 1.5)
lines(r_sclc_t, K_r_sclc_fibro, col = "#87bdd8", lty =1, lwd = 1.5)
lines(r_sclc_t, K_r_luad_sclc, col = "black", lty =1, lwd = 1.5)

legend("topleft", legend = c("Background","LUSC vs Immune_cells", 
                             "SCLC vs Immune_cells",
                             "LUSC vs Fibroblast", 
                             "SCLC vs Fibroblast"),
       col = c("black","#c83349","#f9d5e5", "#034f84", "#87bdd8"), cex = .5,
       lty = c(3, 1,1,1, 1), lwd = 1)

dev.off()


####
s354$cell = 0
s354$cell [which(spa$Fibroblast>.1)] =1
legend_plot <- ggplot() +
  geom_point(aes(x = 1, y = 2, color = "Fibroblast"), size = 4) +
  geom_point(aes(x = 1, y = 1, color = "None"), size = 4) +
  scale_color_manual(name = "Legend", values = c("Fibroblast" = "black", "None" = "white")) +
  theme_void() +  # Remove background and axes
  theme(legend.position = "right", legend.title = element_text(size = 12), legend.text = element_text(size = 10))

p1=SpatialFeaturePlot(s354,features="cell",image.alpha=0,
                   pt.size.factor=2,
                   stroke=.1)+scale_fill_gradientn(colors =
              colorRampPalette(c( "white","white","black"))(100))+
  theme(legend.position = "none") +  # Remove the original legend
  ggtitle("Fibroblast") +  # Add a title in the middle
  theme(plot.title = element_text(hjust = 0.5, vjust = -0.5, size = 16))  # Center and adjust title position

s354$cell2 = 0
s354$cell2 [which(spa$LUSC>.1)] = 1
p2=SpatialFeaturePlot(s354,features="cell2",image.alpha=0,
                   pt.size.factor=2,
                   stroke=.1)+scale_fill_gradientn(colors =
             colorRampPalette(c( "white","white","black"))(100))+
  theme(legend.position = "none") +  # Remove the original legend
  ggtitle("Fibroblast") +  # Add a title in the middle
  theme(plot.title = element_text(hjust = 0.5, vjust = -0.5, size = 16))  # Center and adjust title position


s354$cell3 = 0
s354$cell3 [which(spa$SCLC>0.6003)] = 1
p3=SpatialFeaturePlot(s354,features="cell3",image.alpha=0,
                      pt.size.factor=2,
                      stroke=.1)+scale_fill_gradientn(colors =
                      colorRampPalette(c( "white","white","black"))(100))+
  theme(legend.position = "none") +  # Remove the original legend
  ggtitle("Fibroblast") +  # Add a title in the middle
  theme(plot.title = element_text(hjust = 0.5, vjust = -0.5, size = 16))  # Center and adjust title position


s354$cell4 = 0
s354$cell4 [which(spa$Macrophage >0.2)] = 1
p5=SpatialFeaturePlot(s354,features="cell4",image.alpha=0,
                      pt.size.factor=2,
                      stroke=.1)+scale_fill_gradientn(colors =
                                                        colorRampPalette(c( "white","white","black"))(100))+
  theme(legend.position = "none") +  # Remove the original legend
  ggtitle("Fibroblast") +  # Add a title in the middle
  theme(plot.title = element_text(hjust = 0.5, vjust = -0.5, size = 16))  # Center and adjust title position


p1+p2+p3+p5+legend_plot



immune_spots = list(i_spots,macro_spots,t_spots,fibro_spots)
names = c("Infiltrated", "Macrophage","T cell","Fibroblast")
for (i in c(1:4)) { 
  print(i)
# Create point patterns for LUAD and SCLC spots
luad_ppp <- ppp(lusc_spots$imagerow, lusc_spots$imagecol, window = window)
sclc_ppp <- ppp(sclc_spots$imagerow, sclc_spots$imagecol, window = window)
t_ppp <- ppp(immune_spots[[i]]$imagerow, immune_spots[[i]]$imagecol, window = window)

# Create a multi-type point pattern by combining LUAD and SCLC
multi_ppp <- superimpose(LUAD = luad_ppp, SCLC = sclc_ppp, W = window)
luad_t_ppp <- superimpose(LUAD = luad_ppp, T_cell = t_ppp, W = window)
sclc_t_ppp <- superimpose(SCLC = sclc_ppp, T_cell = t_ppp, W = window)

# Calculate the cross K-function to see how LUAD and SCLC interact spatially
K_cross <- Kcross(multi_ppp, i = "LUAD", j = "SCLC", 
                  correction = "Ripley")
K_cross_luad_t <- Kcross(luad_t_ppp, i = "LUAD", j = "T_cell", 
                         correction = "Ripley")
K_cross_sclc_t <- Kcross(sclc_t_ppp, i = "SCLC", j = "T_cell", 
                         correction = "Ripley")
# Extract r, theo and K(r) values for each dataset
r_sclc_t <- K_cross_sclc_t$r # only need one r value since all the r values are the same
theo_sclc_t <- K_cross_sclc_t$theo #also only needs one theo (background value)
#k values
K_r_sclc_t <- K_cross_sclc_t$iso  # The 'iso' correction, replace if you used another correction
K_r_luad_t <- K_cross_luad_t$iso  # The 'iso' correction, replace if you used another correction
K_r_sclc_lusc <- K_cross$iso 

pdf(paste0("P11_Cell_Pattern_",names[i],".pdf"),width=4,height=4)
# Create an empty plot with appropriate labels and ranges
plot(r_sclc_t, type = "l", col = "blue", lty = 1, lwd = 0, 
     xlab = "Distance (r)", ylab = "K(r)", 
     main = paste0("P11_Cell_Pattern_",names[i]), 
     xlim = range(r_sclc_t),
     ylim = range(K_r_sclc_t, K_r_luad_t,K_r_sclc_lusc))

# Add the lines
lines(r_sclc_t, theo_sclc_t, col = "black", lty = 3, lwd = 1)

lines(r_sclc_t, K_r_luad_t, col = "#c83349", lty = 1, lwd = 1)
lines(r_sclc_t, K_r_sclc_t, col = "#f9d5e5", lty =1, lwd = 1)
lines(r_sclc_t, K_r_sclc_lusc, col = "#034f84", lty =1, lwd = 1)

legend("topleft", legend = c("Background",paste0("LUSC vs ",names[i]), 
                             paste0("SCLC vs ",names[i]),"SCLC vs LUSC"),
       col = c("black","#c83349","#f9d5e5","#034f84"), cex = .5,
       lty = c(3, 1,1,1, 1), lwd = 1)
dev.off()
}



