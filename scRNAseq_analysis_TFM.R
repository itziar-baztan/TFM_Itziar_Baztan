## scRNA-seq ANALYSIS


# Libraries ---------------------------------------------------------------
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
# library(DropletUtils) #conseguir objetos STRIDE

setwd("C:/Users/itzia/OneDrive/Desktop/METODOS_COMPUTACIONALES/TFM/Single_cell_analysis")
# Load the functions
source("funciones/funciones_scrnaseq.R")

# Load data ---------------------------------------------------------------

tumor_data <- Read10X_h5(filename = "C:/Users/itzia/OneDrive/Desktop/METODOS_COMPUTACIONALES/TFM/Single_cell_analysis/single_cell_reyes/scRNAseq_26C-7_d30/TUMOR2_withGFP_polyA/filtered_feature_bc_matrix.h5")
tme_data <- Read10X_h5(filename = "C:/Users/itzia/OneDrive/Desktop/METODOS_COMPUTACIONALES/TFM/Single_cell_analysis/single_cell_reyes/scRNAseq_26C-7_d30/TME1_withGFP_polyA/filtered_feature_bc_matrix.h5")

data_names <- c("tumor", "tme")

data_list <- list(tumor_data, tme_data)
names(data_list) <- data_names

rm(tumor_data, tme_data)


# Preprocessing -----------------------------------------------------------
## QC filtering -----------------------------------------------------------

# Create Seurat objects and calculate %mt and %rb
seurat_list <- lapply(data_list, function(x) {
  obj <- CreateSeuratObject(x)
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^mt-")
  obj[["percent.rb"]] <- PercentageFeatureSet(obj, pattern = "^Rps|^Rpl")
  obj
})

# Filter cells
seurat_filt <- lapply(seurat_list, function(obj) {
  subset(obj, subset = nFeature_RNA > 500 & percent.mt < 20)
})

# Table with number of cells before and after filtering
cell_counts <- data.frame(
  sample   = names(seurat_list),
  before   = sapply(seurat_list, ncol),
  after    = sapply(seurat_filt, ncol)
)
print(cell_counts)

# Plot before and after filtering
pdf("single_cell_graficos/VlnPlot Tumor_vs_TME_before.pdf", width = 10, height = 10)
plot_vln(seurat_list, "Tumor vs TME before filtering")
dev.off()

pdf("single_cell_graficos/VlnPlot Tumor_vs_TME_after.pdf", width = 10, height = 10)
plot_vln(seurat_filt, "Tumor vs TME after filtering")
dev.off()

rm(seurat_list)

# Save objects
# saveRDS(seurat_filt[[1]], file = "objetos_R/tumor_filt.rds")
# saveRDS(seurat_filt[[2]], file = "objetos_R/tme_filt.rds")

## Remove doublets ---------------------------------------------------------
# Apply function
seurat_filt <- lapply(seurat_filt, filter_doublets)
names(seurat_filt) <- data_names

# Save objects
# saveRDS(seurat_filt[[1]], file = "objetos_R/tumor_filt_doublets.rds")
# saveRDS(seurat_filt[[2]], file = "objetos_R/tme_filt_doublets.rds")

# Apply function
doublet_summary <- plot_doublet_summary(seurat_filt, pdf_file = "single_cell_graficos/dobletes_singletes.pdf")


# Normalization and integration -------------------------------------------
seurat_integrated <- integrate_seurat_objects(seurat_filt)

# Save object
# saveRDS(seurat_integrated, file = "objetos_R/seurat_integrated.rds")


# Scale by cell cycle -----------------------------------------------------
# Cell cycle genes
s.genes <- c("Mcm5", "Pcna", "Tyms", "Fen1", "Mcm7", "Mcm4", "Rrm1", "Ung", "Gins2", "Mcm6", 
             "Cdca7", "Dtl", "Prim1", "Uhrf1", "Cenpu", "Hells", "Rfc2", "Polr1b", "Nasp", "Rad51ap1", 
             "Gmnn", "Wdr76", "Slbp", "Ccne2", "Ubr7", "Msh2", "Rad51", "Rrm2", "Cdc45", "Cdc6", 
             "Exo1", "Tipin", "Dscc1", "Blm", "Casp8ap2", "Usp1", "Clspn", "Pola1", "Chaf1b", "Mrpl36", 
             "E2f8")
g2m.genes <- c("Hmgb2", "Cdk1", "Nusap1", "Ube2c", "Birc5", "Tpx2", "Top2a", "Ndc80", "Cks2", "Nuf2", 
               "Cks1b", "Mki67", "Tmpo", "Cenpf", "Tacc3", "Pimreg", "Smc4", "Ccnb2", "Ckap2l", "Ckap2", 
               "Aurkb", "Bub1", "Kif11", "Anp32e", "Tubb4b", "Gtse1", "Kif20b", "Hjurp", "Cdca3", "Jpt1", 
               "Cdc20", "Ttk", "Cdc25c", "Kif2c", "Rangap1", "Ncapd2", "Dlgap5", "Cdca2", "Cdca8", 
               "Ect2", "Kif23", "Hmmr", "Aurka", "Psrc1", "Anln", "Lbr", "Ckap5", "Cenpe", "Ctcf", 
               "Nek2", "G2e3", "Gas2l3", "Cbx5", "Cenpa")

# Cell cycle scoring and regression
set.seed(123)
seurat_integrated <- seurat_integrated %>%
  CellCycleScoring(s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE) %>%
  ScaleData(vars.to.regress = c("S.Score", "G2M.Score"))

# Save object
# saveRDS(seurat_integrated, file = "objetos_R/seurat_integrated_scaled.rds")


# Dimensionality reduction ------------------------------------------------
# Set integrated assay as default
DefaultAssay(seurat_integrated) <- "integrated"

# PCA
set.seed(123)
seurat_integrated <- RunPCA(seurat_integrated, features = VariableFeatures(seurat_integrated))

# Elbow plot
save_pdf("single_cell_graficos/ElbowPlot_PCA.pdf", ElbowPlot(seurat_integrated, ndims = 30))

# t-SNE
set.seed(123)
seurat_integrated <- RunTSNE(seurat_integrated, dims = 1:30)

# Set identity to original sample
Idents(seurat_integrated) <- "orig.ident"

# t-SNE plots
save_pdf("single_cell_graficos/tSNE_by_sample.pdf", DimPlot(seurat_integrated, reduction = "tsne"))
save_pdf("single_cell_graficos/tSNE_by_phase.pdf", DimPlot(seurat_integrated, reduction = "tsne", group.by = "Phase"))


# Clustering --------------------------------------------------------------
# Clustering workflow
set.seed(123)
clustering_int <- seurat_integrated %>%
  FindNeighbors(dims = 1:30) %>%
  FindClusters(resolution = seq(0.1, 1, 0.1))

# Save object
# saveRDS(clustering_int, file = "objetos_R/clustering_int.rds")

# Clustree plot
save_pdf("single_cell_graficos/clustree.pdf",
         clustree(clustering_int@meta.data, prefix = "integrated_snn_res."))

# Keep only chosen resolution (0.5)
clustering_int@meta.data <- clustering_int@meta.data[
  , !grepl("integrated_snn_res\\.(?!0\\.5)", colnames(clustering_int@meta.data), perl = TRUE)
]
clustering_int <- SetIdent(clustering_int, value = "Final_Labels")

# t-SNE with chosen resolution
set.seed(123)
clustering_int <- RunTSNE(clustering_int, dims = 1:30)
save_pdf("single_cell_graficos/t_SNE_clusters_0.5.pdf",
         DimPlot(clustering_int, reduction = "tsne", label = TRUE, pt.size = 1,
                 label.size = 12, repel = TRUE) +
           ggplot2::ggtitle("t-SNE con resolución 0.5") +
           theme(legend.position = "none",
                 legend.text = element_text(size = 21),
                 legend.title = element_text(size = 16),
                 axis.title.x = element_text(size = 18),
                 axis.title.y = element_text(size = 18),
                 axis.text.x = element_text(size = 14),
                 axis.text.y = element_text(size = 14))
)

# Find markers and top genes per cluster
cluster_markers <- FindAllMarkers(clustering_int,
                                  logfc.threshold = 0.25,
                                  test.use = "wilcox",
                                  min.pct = 0.1,
                                  only.pos = TRUE)
top_markers <- cluster_markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)

# Dotplot of top markers
save_pdf("single_cell_graficos/top_cluster_markers_dotplot.pdf",
         DotPlot(clustering_int, features = top_markers$gene) + RotatedAxis())


# Manual annotation -------------------------------------------------------
## Based on known markers --------------------------------------------------
# Define marker lists
markers <- list(
  tumoral_cells = c("GFP"),
  immune_cells  = c("Ptprc"),
  B_cells       = c("Cd19","Ms4a1","Cd79a","Igkc","Igha","Ighm","Fcgr1","Fcmr"),
  T_cells_CD4-CD8 = c("Cd3e","Cd4","Cd8a"),
  eritrocytes           = c("Hbb-bs","Hba-a1","Hba-a2","Hbb-bt"),
  Tregs = c("Cd3e","Cd4","Foxp3","Izumo1r","Tnfrsf4","Icos"),
  NK_cells      = c("Klrk1","Ncr1","Klrd1","Prf1"),
  NK-NKT = c("Klrk1","Ncr1","Klrd1","Prf1","Cd3e"),
  myeloid_cells = c("Itgam","Csf1r","Ccr2"),
  monocytes     = c("Ly6c1","Ly6c2","Spn","Ccr7","Cacnb3","Nudt17","Bcl2l14"),
  cDC2          = c("Itgax","Sirpa","Irf4"),
  cDC1          = c("Flt3","Irf8","Xcr1","Clec9a","H2-Ab1","H2-Aa","H2-Eb1","H2-Eb2","Cd40","Ly75","Cd24a"),
  macrophages     = c("Adgre1", "Itga4","Grk3", "P2ry12", "Gpr34", "Slco2b1", "Ptgs1", "Cx3cr1","Itgal"),
  ILCs            = c("Il5", "Il9", "Il13", "Il17", "Il22", "Infg", "Tnf","Ifng","Eomes","Rora","Gata3","Rorc","Runx1"),
  astrocytes      = c("Gpc5", "Gfap", "Aqp4", "Apoe", "Wdr17", "Plpp3","Mertk"),
  endothelial     = c("Flt1", "Slco1a4", "Adgrl4", "Slc2a1", "Klf2", "Mecom"),
  microglia       = c("Inpp5d", "Hexb", "Tgfbr1", "C1qa", "Ctss", "C1qb","Tmem119", "Adgre1"),
  fibroblast      = c("Slc4a4","Kcnma1","Bnc2","Cemip","Bicc1","Flrt2"),
  neuron          = c("Celf2", "Ptprd", "Arpp21", "Sv2b", "Pcsk2", "Phactr1","Grip1", "Galntl6", "Gad1", "Gad2", "Dlx6os1", "Cntnap2"),
  oligodendrocytes= c("Plp1", "Mbp", "St18", "Prr5l", "Mobp", "Mal"),
  pericyte        = c("Atp13a5", "Vtn", "Cald1", "Ebf1", "Abcc9",  "Pdgfrb"),
  granulocytes    = c("Ly6g"),
  ependymal       = c("Foxj1", "Pifo", "Dynlrb2","Clec9a"),
  NSCs            = c("Gfap","Sox2","Blbp","Pax6","Dcx","Psa-ncam","Neun","Tubb3","Map2")
)

# Run loop
lapply(names(markers), function(ct) {
  plot_celltype(ct, markers[[ct]])
})


# Automatic annotation ----------------------------------------------------
# Load references
immgen_ref <- celldex::ImmGenData()
mouse_ref  <- celldex::MouseRNAseqData()

# Integrate layers
clustering_int <- IntegrateLayers(
  object = clustering_int, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "integrated.harmony",
  verbose = TRUE, assay = "RNA"
)

# Re-join and convert layers
clustering_int[["RNA"]] <- as(JoinLayers(clustering_int[["RNA"]]), "Assay")
clustering_single_cell  <- as.SingleCellExperiment(clustering_int)

# SingleR annotation
annotations <- list(
  Mouse_Main = SingleR(test = clustering_single_cell, ref = mouse_ref,  labels = mouse_ref$label.main),
  Mouse_Fine = SingleR(test = clustering_single_cell, ref = mouse_ref,  labels = mouse_ref$label.fine),
  ImmGen     = SingleR(test = clustering_single_cell, ref = immgen_ref, labels = immgen_ref$label.main)
)

# Add metadata
for (nm in names(annotations)) {
  clustering_int <- AddMetaData(clustering_int, 
                                metadata = annotations[[nm]]$pruned.labels, 
                                col.name = paste0(nm, "_Labels"))
}

plot_labels(clustering_int, "Mouse_Main_Labels", "Figuras/Mouse_Main_Labels.pdf", text_size = 8)
plot_labels(clustering_int, "Mouse_Fine_Labels", "Figuras/Mouse_Fine_Labels.pdf", text_size = 7)
plot_labels(clustering_int, "ImmGen_Labels",     "Figuras/ImmGen_Labels.pdf",     text_size = 8)


# Rename clusters ---------------------------------------------------------
# Set resolution and mapping
Idents(clustering_int) <- "integrated_snn_res.0.5"

cluster_map <- c(
  "19" = "Oligodendrocytes", "18" = "Astrocytes", "17" = "B Cells",
  "16" = "ILCs", "15" = "Endothelial cells", "14" = "TAMs",
  "13" = "NSCs", "12" = "Tumoral cells", "11" = "T cells",
  "10" = "TAMs", "9" = "NK cells", "8" = "Neurons",
  "7" = "Dendritic cells", "6" = "Monocytes", "5" = "Ependymal cells",
  "4" = "T cells", "3" = "TAMs", "2" = "TAMs", "1" = "TAMs", "0" = "TAMs"
)

# Keep original labels and rename
clustering_int$Original_Labels <- Idents(clustering_int)
clustering_int <- RenameIdents(clustering_int, cluster_map)
clustering_int$Renamed_Labels <- Idents(clustering_int)

# Plot
save_pdf("Figuras/tSNE_clusters_renamed_res0.5.pdf", 
         DimPlot(clustering_int, reduction = "tsne", label = TRUE, pt.size = 1, 
                 label.size = 9, repel = TRUE) +
           ggtitle("t-SNE renamed (resolution 0.5)") +
           theme(
             legend.text  = element_text(size = 21),
             legend.title = element_text(size = 16)
           )
)


# TAMs Reclustering -------------------------------------------------------
# Define TAMs
clustering_int <- SetIdent(clustering_int, value = "integrated_snn_res.0.5")
TAMs <- subset(clustering_int, idents = c(0,1,2,3,10,14))

## Preprocessing -----------------------------------------------------------
set.seed(123)
TAMs <- TAMs %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA(features = VariableFeatures(.)) %>%
  RunTSNE(dims = 1:30) %>%
  FindNeighbors(dims = 1:30) %>%
  FindClusters(resolution = 0.5)

# General t-SNE 
tsne_plot <- DimPlot(TAMs, reduction = "tsne", label = TRUE, pt.size = 1, label.size = 9, repel = TRUE) +
  ggtitle("t-SNE TAMs 0.5") +
  theme(legend.text = element_text(size = 21),
        legend.title = element_text(size = 16))

save_pdf("C:/Users/itzia/OneDrive/Desktop/METODOS_COMPUTACIONALES/TFM/versiones_TFM/Figuras/t-SNE TAMs.pdf", tsne_plot)


## Manual annotation -----------------------------------------------------------
DefaultAssay(TAMs) <- "RNA"
TAMs <- SetIdent(TAMs, value = "RNA_snn_res.0.5")

# All cell types and its genes
TAMs_markers <- list(
  tumoral_cells = c("GFP"),
  microglia = c("Inpp5d", "Hexb", "Tgfbr1", "C1qa", "Ctss", "C1qb","Tmem119", "Adgre1"),
  macrophages = c("Adgre1", "Itga4","Grk3", "P2ry12", "Gpr34", "Slco2b1", "Ptgs1", "Cx3cr1","Itgal"),
  macrophages = c("P2ry12","P2ry13", "Mrc1","Lyve1", "Fcgr1", "Cd68", "Adgre1", "Csf1r"),
  neuron = c("Celf2", "Ptprd", "Arpp21", "Sv2b", "Pcsk2", "Phactr1","Grip1", "Galntl6", "Gad1", "Gad2", "Dlx6os1", "Cntnap2"),
  astrocytes = c("Gpc5", "Gfap", "Aqp4", "Apoe", "Wdr17", "Plpp3","Mertk"),
  oligodendrocytes = c("Plp1", "Mbp", "St18", "Prr5l", "Mobp", "Mal"),
  NSCs = c("Gfap","Sox2","Blbp","Pax6","Dcx","Psa-ncam","Neun","Tubb3","Map2")
)

lapply(names(TAMs_markers), function(ct) {
  plot_celltype(ct, TAMs_markers[[ct]])
})

## Automatic annotation -----------------------------------------------------------
TAMs_single_cell <- as.SingleCellExperiment(as(TAMs[["RNA"]], "Assay"))

# SingleR annotations
annotations <- list(
  ImmGen     = SingleR(test = TAMs_single_cell, ref = immgen_ref, labels = immgen_ref$label.main),
  Mouse_Main = SingleR(test = TAMs_single_cell, ref = mouse_ref,  labels = mouse_ref$label.main),
  Mouse_Fine = SingleR(test = TAMs_single_cell, ref = mouse_ref,  labels = mouse_ref$label.fine)
)

# Add metadata
for(nm in names(annotations)){
  TAMs <- AddMetaData(TAMs, metadata = annotations[[nm]]$pruned.labels, col.name = paste0(nm, "_Labels"))
}

# Plots
plot_labels(TAMs, "ImmGen_Labels",     "single_cell_notacion_0.5/ImmGen_Labels_TAMs.pdf", text_size = 8)
plot_labels(TAMs, "Mouse_Main_Labels", "single_cell_notacion_0.5/Mouse_Main_Labels_TAMs.pdf", text_size = 8)
plot_labels(TAMs, "Mouse_Fine_Labels", "single_cell_notacion_0.5/Mouse_Fine_Labels_TAMs.pdf", text_size = 7)

# Identify markers per cluster
cluster_markers_TAMs <- FindAllMarkers(
  TAMs, logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.1, only.pos = TRUE
)

top_markers_TAMs <- cluster_markers_TAMs %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)


## Rename TAMs -------------------------------------------------------------
# Set cluster identities
Idents(TAMs) <- "RNA_snn_res.0.5"

# New names
new_TAMs_names <- c(
  "0" = "Microglia", "1" = "Macrophages", "2" = "Microglia",
  "3" = "Macrophages", "4" = "Macrophages", "5" = "Macrophages",
  "6" = "Microglia", "7" = "Macrophages", "8" = "Tumoral cells",
  "9" = "Macrophages"
)

# Rename clusters
TAMs <- RenameIdents(TAMs, new_TAMs_names)

# Save updated labels in metadata
TAMs$Updated_Labels <- Idents(TAMs)

# t-SNE plot with renamed clusters
tsne_plot <- DimPlot(TAMs, reduction = "tsne", label = TRUE, pt.size = 1, label.size = 9, repel = TRUE) +
  ggtitle("t-SNE RENAMED TAMs 0.5") +
  theme(
    legend.text = element_text(size = 21),
    legend.title = element_text(size = 16)
  )

save_pdf(
  "C:/Users/itzia/OneDrive/Desktop/METODOS_COMPUTACIONALES/TFM/versiones_TFM/Figuras/t_SNE_RENAMED_TAMs_0.5.pdf",
  tsne_plot
)


# Rename final clusters ---------------------------------------------------
# Updated_Labels and copy from TAMs
clustering_int$Updated_Labels <- as.character(TAMs$Updated_Labels)
clustering_int$Original_Labels <- as.character(clustering_int$Original_Labels)

# Final_Labels: prefer Updated_Labels, fallback to Original_Labels
clustering_int$Final_Labels <- coalesce(clustering_int$Updated_Labels, clustering_int$Original_Labels)

# Set identities to Final_Labels
Idents(clustering_int) <- clustering_int$Final_Labels

# Create t-SNE plot with combined annotations
final_tsne <- DimPlot(clustering_int, reduction = "tsne", label = TRUE, group.by = "Final_Labels",
                      label.size = 8, repel = TRUE) +
  ggtitle("Clusters with Combined Annotations") +
  theme(
    legend.text = element_text(size = 21),
    legend.title = element_text(size = 16)
  )

save_pdf("C:/Users/itzia/OneDrive/Desktop/METODOS_COMPUTACIONALES/TFM/versiones_TFM/Figuras/t_SNE_FINAL.pdf",
         final_tsne)

# Save annotated object
# saveRDS(clustering_int, file = "objetos_R/seurat.anotado.rds")


# Re-annotate tumoral cells -----------------------------------------------
# Add GFP module score
seurat.anotado <- AddModuleScore(object = seurat.anotado,
                                 features = "GFP",
                                 name = "GFP_Score")

# GFP expression per cell type
gfp_summary <- seurat.anotado@meta.data %>%
  group_by(Final_Labels) %>%
  summarise(
    Min = min(GFP_Score1, na.rm = TRUE),
    Q1 = quantile(GFP_Score1, 0.25, na.rm = TRUE),
    Median = median(GFP_Score1, na.rm = TRUE),
    Q3 = quantile(GFP_Score1, 0.75, na.rm = TRUE),
    Max = max(GFP_Score1, na.rm = TRUE),
    Mean = mean(GFP_Score1, na.rm = TRUE)
  )

# Jitter plot of GFP expression by cell type
gfp_plot <- ggplot(seurat.anotado@meta.data, aes(x = Final_Labels, y = GFP_Score1, color = Final_Labels)) +
  geom_jitter(width = 0.5, alpha = 0.6) +  # Spread points for better visualization
  theme_minimal() +
  labs(title = "Tumoral cell distribution according to GFP",
       x = "Cell type",
       y = "GFP") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 17),
    legend.position = "none",
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.text.y = element_text(size = 14)
  )

save_pdf("C:/Users/itzia/OneDrive/Desktop/METODOS_COMPUTACIONALES/TFM/versiones_TFM/Figuras/expresion_GFP.pdf",
         gfp_plot)

# Update: assign "Tumoral cells" for neurons/NSCs with high GFP
seurat.anotado@meta.data <- seurat.anotado@meta.data %>%
  mutate(Final_Labels = ifelse(Final_Labels %in% c("Neurons", "NSCs") & GFP_Score1 > 0.4,
                               "Tumoral cells",
                               Final_Labels))





