# Spatial Transcriptomics analysis


# Libraries ---------------------------------------------------------------
library(Seurat)
library(SpatialExperiment)
library(ggplot2)
library(scales)
library(patchwork)
library(dplyr)
library(CARD)
library(pbmcapply)
library(SingleCellExperiment)
library(MuSiC)
library(STdeconvolve)

# Load data ---------------------------------------------------------------
data_dir = "Visium_d30_porta5/outs/"
h5_mat_name = "filtered_feature_bc_matrix.h5"

DMG_pons <- Load10X_Spatial(
  data.dir = data_dir, filename = h5_mat_name, assay = "Spatial",
  slice = "Pons_30_porta5", filter.matrix = TRUE, to.upper = FALSE)


# Preprocessing -----------------------------------------------------------
## QC filtering -----------------------------------------------------------
# Create SpatialExperiment from Seurat object
DMG_pons <- SpatialExperiment(
  assays = list(counts = GetAssayData(DMG_pons, "Spatial", "counts")),
  colData = DMG_pons@meta.data,
  spatialCoords = as.matrix(GetTissueCoordinates(DMG_pons)[, c("x","y")])
)

plotSpots(DMG_pons, in_tissue = NULL)  # visualize all spots
dim(DMG_pons)  # check object dimensions


### Filtering by library size -----------------------------------------------
# Add library size to colData
colData(DMG_pons)$library_size <- colSums(assay(DMG_pons, "counts"))

# Histogram + density plot of library sizes
libsize_plot <- colData(DMG_pons) %>% 
  as.data.frame() %>% 
  ggplot(aes(x = library_size)) +
  geom_histogram(aes(y = after_stat(density)), colour = "black", fill = "grey") +
  geom_density(alpha = 0.5, adjust = 1.0, fill = "#A0CBE8", colour = "#4E79A7") +
  scale_x_continuous(breaks = pretty_breaks(10)) +
  scale_y_continuous(breaks = pretty_breaks(10)) +
  xlab("Library size") + ylab("Density") +
  theme_classic() +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 16))

save_pdf("histogram_library_sizes.pdf", libsize_plot)

# Apply QC threshold for library size
qc_threshold <- 700
colData(DMG_pons)$qc_lib_size <- colData(DMG_pons)$library_size < qc_threshold

# Plot spatial distribution of low-QC spots
lowQC_plot <- plotSpotQC(
  DMG_pons,
  plot_type = "spot",
  x_metric = "library_size",
  annotate = "qc_lib_size",
  point_size = 3
) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black"), name = "Low QC") +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 16))

save_pdf("spots_low_library_size.pdf", lowQC_plot)

### Filtering by number of genes -----------------------------------------------
# Plot expressed genes per spot
plot_density_hist(DMG_pons@colData, "nFeature_Spatial", 
                  "Genes expressed in each spot", 
                  "density_histogram_expressed_genes.pdf")

# Apply QC threshold for expressed genes
qc_detected <- colData(DMG_pons)$nFeature_Spatial < 500
colData(DMG_pons)$qc_detected <- qc_detected
message("Number of low-QC spots by detected genes: ", sum(qc_detected))

# Plot spatial pattern for low expressed genes
plot_qc_spots(DMG_pons, "library_size", "qc_detected", 
              "spots_low_expressed_genes.pdf")

# Combine library size and detected genes QC
discard <- colData(DMG_pons)$qc_lib_size | colData(DMG_pons)$qc_detected
colData(DMG_pons)$discard <- discard
message("Number of total discarded spots: ", sum(discard))

# Plot spatial pattern for all discarded spots
plot_qc_spots(DMG_pons, "library_size", "discard", "spots_discarded.pdf", point_size = 1.5)


# Normalization and dimensionality reduction ------------------------------
DMG_pons <- SCTransform(DMG_pons, assay = "Spatial", verbose = FALSE)
#Dimensionality reduction, clustering, and visualization
DMG_pons <- RunPCA(DMG_pons, assay = "SCT", verbose = FALSE)
DMG_pons <- FindNeighbors(DMG_pons, reduction = "pca", dims = 1:30)
DMG_pons <- FindClusters(DMG_pons, verbose = FALSE)
DMG_pons <- RunTSNE(DMG_pons, reduction = "pca", dims = 1:30)


# Seurat Deconvolution ----------------------------------------------------
seurat.anotado <- readRDS("C:/Users/itzia/OneDrive/Desktop/METODOS_COMPUTACIONALES/TFM/Single_cell_analysis/objetos_R/seurat.anotado_tumoral.rds")

# Preprocess reference Seurat object
seurat.anotado <- preprocess_seurat(seurat.anotado)

# Transfer labels from reference to query
DMG_pons <- transfer_labels(seurat.anotado, DMG_pons, common.features)

# Plot and save spatial predictions
plot_seurat_predictions(
  DMG_pons,
  features = c("Tumoral cells","Macrophages","NSCs","Neurons"),
  filename = "seurat_deconv_BUENO.pdf"
)


# CARD Deconvolution ------------------------------------------------------
# CARD object
CARD_obj <- prepare_CARD(seurat.anotado, DMG_pons)

# CARD deconvolution
CARD_obj <- CARD_deconvolution(CARD_obj)
print(CARD_obj@Proportion_CARD[1:2,])

# Filter low proportions
pro_filtered <- CARD_obj@Proportion_CARD %>% 
  as.data.frame() %>% 
  filter_proportions(threshold = 0.05)

# Plot and save pie chart
plot_CARD_pie(
  pro = pro_filtered,
  location = CARD_obj@spatial_location,
  filename = "CARD_deconv_comparacion.pdf",
  colors = c("#FFD92F","#FCCDE5","#D9D9D9","#377EB8","#BFEFFF",
             "#4B0082","#8B184B","#D2691E","#FF0000","#DEB887",
             "#E7298A","#66A61E"),
  radius = 9
)


# STdeconvolve Deconvolution ---------------------------------------------
# Run STdeconvolve
st_results <- run_STdeconvolve(
  spatial_obj = DMG_pons,
  sc_meta = seurat.anotado@meta.data,
  ct_var = "Final_Labels",
  K_range = 9:15
)

# Visualization colors
my_colors <- c(
  "#E6194B", "#3CB44B", "#0082C8", "#F58231", "#911EB4",
  "#46F0F0", "#F032E6", "#D2F53C", "#FABEBE", "#008080",
  "#AA6E28", "#FFFAC8", "#800000", "#A9A9A9", "#000000"
)

# Plot
plot_STdeconvolve(
  deconProp = st_results$deconProp,
  spatial_coords = GetTissueCoordinates(DMG_pons),
  cell_types = st_results$cell_types,
  filename = "STdeconvolve/STdeconv_con_annot.pdf",
  colors = my_colors
)

## Annotate STdeconvolve clusters ------------------------------------------

# Proxy matrices
proxy <- create_proxy(cell_types_filtered, DMG_cd_matrix)

# Correlation Beta (transcriptional)
cor_beta <- compute_correlation(
  m1 = deconGexp, 
  m2 = t(proxy$ProxyGexp), 
  type = "b"
)
save_corr_plot(
  mat = cor_beta,
  filename = "STdeconvolve/Transcriptional_correlation.pdf",
  colLabs = "Deconvolved cell-types",
  rowLabs = "Ground truth cell-types",
  title = "Transcriptional correlation"
)

# Correlation Theta (proportional)
cor_theta <- compute_correlation(
  m1 = deconProp, 
  m2 = proxy$ProxyTheta, 
  type = "t"
)
save_corr_plot(
  mat = cor_theta,
  filename = "STdeconvolve/Proportional_correlation.pdf",
  colLabs = "Deconvolved cell-types",
  rowLabs = "Ground truth cell-types",
  title = "Proportional correlation"
)

# Ordered pairs based on best matches
pairs <- lsatPairs(t(cor_theta))
m <- t(cor_theta)[pairs$rowix, pairs$colsix]

correlationPlot(
  mat = t(m),
  colLabs = "Deconvolved cell-types",
  rowLabs = "Ground truth cell-types",
  title = "Ordered proportional correlation",
  annotation = TRUE
) + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0))





