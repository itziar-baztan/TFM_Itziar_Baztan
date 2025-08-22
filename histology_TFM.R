# HISTOLOGY

library(Seurat)
library(ggplot2)

regions_df <- read.csv("Brain_regions/Regiones.csv", sep = ";", stringsAsFactors = FALSE)

rownames(regions_df) <- regions_df$Barcode
regions_df$Barcode <- NULL

# filtra el data frame solo a los barcodes que estén en el objeto Seurat
regions_df <- regions_df[colnames(DMG_pons), , drop = FALSE]

# agrega la columna "region" como metadato
DMG_pons <- AddMetaData(DMG_pons, metadata = regions_df["Regiones"])

pdf("Brain_regions/Brain_regions_spatial_plot.pdf", width = 10, height = 10, pointsize = 12,family = "sans", bg = "white")
SpatialDimPlot(DMG_pons,group.by = 
                 "Regiones") + theme(legend.position = "right")
dev.off()