# FUNCIONES scRNA-seq
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(patchwork)
library(SingleCellExperiment)
library(scDblFinder)
library(clustree)
library(celldex)
library(SingleR)

# Function to generate violin plots
plot_vln <- function(objs, title) {
  wrap_plots(lapply(objs, function(x)
    VlnPlot(x, features = c("nFeature_RNA","nCount_RNA","percent.mt"), ncol = 3)
  )) + plot_annotation(
    title = title,
    theme = theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
  )
}

# Function to filter doublets
filter_doublets <- function(seurat_obj) {
  # Seurat object to SingleCellExperiment
  sce <- as.SingleCellExperiment(seurat_obj)
  # Detect doublets
  sce <- scDblFinder(sce)
  # Doublet information to Seurat metadata
  seurat_obj$predicted_doublets <- colData(sce)$scDblFinder.class
  return(seurat_obj)
}

# Function to plot doublets vs singlets
plot_doublet_summary <- function(seurat_list, pdf_file = "dobletes_singletes.pdf") {
  
  # Extract metadata and combine into a single dataframe
  df_list <- lapply(names(seurat_list), function(nm) {
    df <- seurat_list[[nm]]@meta.data
    df$orig.ident[df$orig.ident == "singleCell_tme"] <- "TME"
    df$orig.ident[df$orig.ident == "singleCell_tumor"] <- "Tumor"
    df$sample <- nm 
    df
  })
  
  df_combined <- do.call(rbind, df_list)
  df_combined <- as.data.frame(df_combined)
  
  # Summarise counts and percentages by sample and doublet class
  df_summary <- df_combined %>%
    mutate(
      orig.ident = as.character(orig.ident),
      scDblFinder.class = as.character(scDblFinder.class)
    ) %>%
    group_by(orig.ident, scDblFinder.class) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(orig.ident) %>%
    mutate(percent = n / sum(n) * 100)
  
  # Plot
  pdf(pdf_file, width = 10, height = 10, pointsize = 12, family = "sans", bg = "white")
  
  p <- ggplot(df_summary, aes(x = orig.ident, y = percent, fill = scDblFinder.class)) +
    geom_bar(stat = "identity", position = "fill") +
    scale_fill_manual(values = c("singlet" = "#8da0cb", "doublet" = "#fc8d62")) +
    labs(x = "Sample", y = "Proportion", fill = "Cell Type",
         title = "Proportion of Singlets and Doublets by Sample") +
    theme_minimal(base_size = 16) +
    theme(axis.title = element_text(size = 18),
          axis.text = element_text(size = 14),
          plot.title = element_text(hjust = 0.5, size = 20, face = "bold"))
  
  print(p)
  dev.off()
  
  return(df_summary)
}

# Function to integrate objects
integrate_seurat_objects <- function(seurat_list, dims = 1:30) {
  set.seed(123)
  # Normalize and find variable features
  seurat_list <- lapply(seurat_list, function(obj) {
    obj <- NormalizeData(obj)
    obj <- FindVariableFeatures(obj)
    return(obj)
  })
  # Select integration features
  features <- SelectIntegrationFeatures(object.list = seurat_list)
  # Find integration anchors
  anchors <- FindIntegrationAnchors(object.list = seurat_list,
                                    anchor.features = features)
  # Integrate data
  seurat_integrated <- IntegrateData(anchorset = anchors, dims = dims)
  
  return(seurat_integrated)
}

# Function to save as pdf
save_pdf <- function(filename, printed_plot, w = 10, h = 10) {
  pdf(filename, width = w, height = h, pointsize = 12, family = "sans", bg = "white")
  print(printed_plot)
  dev.off()
}

# Function to plot by known markers (manual annotation)
plot_celltype <- function(celltype, features, reduction = "tsne") {
  pdf(paste0("single_cell_notacion_0.5/", celltype, " feature.pdf"),
      width = 10, height = 10, pointsize = 12, family = "sans", bg = "white")
  
  print(DotPlot(clustering_int, features = features, dot.scale = 15))
  print(FeaturePlot(clustering_int, features = features,
                    order = TRUE, pt.size = 0.1, reduction = reduction))
  dev.off()
}

# Function to plot automatic annotation
plot_labels <- function(object, label_col, file, text_size = 8) {
  object <- SetIdent(object, value = label_col)
  pdf(file, width = 10, height = 10, pointsize = 12, family = "sans", bg = "white")
  p <- DimPlot(object, reduction = "tsne", group.by = label_col)
  p <- LabelClusters(plot = p, id = label_col, repel = TRUE, size = text_size) +
    theme(legend.position = "none",
          legend.text = element_text(size = text_size * 2.5),
          legend.title = element_text(size = 16))
  print(p)
  dev.off()
}