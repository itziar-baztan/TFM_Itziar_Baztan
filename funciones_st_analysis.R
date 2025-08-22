# FUNCIONES ST ANALYSIS
library(Seurat)
library(dplyr)
library(ggplot2)
library(scales)
library(ggspavis)
library(scater)
library(scran)
library(igraph)
library(pheatmap)
library(ggExtra)
library(CARD)
library(STdeconvolve)
library(CARD)
library(pbmcapply)
library(SingleCellExperiment)
library(MuSiC)
library(SpatialExperiment)
library(STdeconvolve)


# Function to save PDF plots
save_pdf <- function(filename, plot, w = 10, h = 10) {
  pdf(filename, width = w, height = h, pointsize = 12, family = "sans", bg = "white")
  print(plot)
  dev.off()
}

# Function to create density + histogram plots
plot_density_hist <- function(data, metric, xlabel, filename) {
  p <- data %>% 
    as.data.frame() %>% 
    ggplot(aes(x = !!sym(metric))) +
    geom_histogram(aes(y = after_stat(density)), colour = "black", fill = "grey") +
    geom_density(alpha = 0.5, adjust = 1.0, fill = "#A0CBE8", colour = "#4E79A7") +
    scale_x_continuous(breaks = scales::pretty_breaks(10)) +
    scale_y_continuous(breaks = scales::pretty_breaks(10)) +
    xlab(xlabel) + ylab("Density") +
    theme_classic() +
    theme(axis.title = element_text(size = 20),
          axis.text = element_text(size = 16))
  
  save_pdf(filename, p)
}

# Function to plot spatial QC spots
plot_qc_spots <- function(obj, metric, annotate_col, filename, point_size = 3) {
  p <- plotSpotQC(obj,
                  plot_type = "spot",
                  x_metric = metric,
                  annotate = annotate_col,
                  point_size = point_size) +
    scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black"), name = "Low QC") +
    theme(axis.title = element_text(size = 20),
          axis.text = element_text(size = 16))
  
  save_pdf(filename, p)
}


# FUNCIONES SEURAT --------------------------------------------------------

# Function to perform SCT, PCA, and t-SNE in a single pipeline
preprocess_seurat <- function(seurat_obj, assay = "RNA", ncells = 3000, dims = 1:30) {
  seurat_obj %>%
    SCTransform(assay = assay, ncells = ncells, verbose = TRUE) %>%
    RunPCA(verbose = TRUE) %>%
    RunTSNE(dims = dims)
}

# Function to transfer labels from reference to query
transfer_labels <- function(reference, query, features, dims = 1:30, ref_labels_col = "Final_Labels") {
  anchors <- FindTransferAnchors(
    reference = reference,
    query = query,
    normalization.method = "SCT",
    features = features,
    verbose = TRUE
  )
  
  predictions <- TransferData(
    anchorset = anchors,
    refdata = reference[[ref_labels_col]],
    prediction.assay = TRUE,
    weight.reduction = query[["pca"]],
    dims = dims
  )
  
  query[["predictions"]] <- predictions
  DefaultAssay(query) <- "predictions"
  return(query)
}

# Function to plot spatial feature predictions and save as PDF
plot_seurat_predictions <- function(seurat_obj, features, filename, pt.size.factor = 2.3, ncol = 2) {
  p <- SpatialFeaturePlot(
    seurat_obj,
    features = features,
    pt.size.factor = pt.size.factor,
    ncol = ncol,
    crop = TRUE,
    alpha = c(0.5, 1),
    image.alpha = 0.3
  )
  save_pdf(filename, p)
}


# FUNCIONES CARD ----------------------------------------------------------
# Function to prepare CARD object
prepare_CARD <- function(sc_obj, spatial_obj, ct_var = "Final_Labels", sample_var = "orig.ident",
                         minCountGene = 100, minCountSpot = 5) {
  sc_count <- as.matrix(GetAssayData(sc_obj, slot = "counts", assay = "RNA"))
  sc_meta <- sc_obj@meta.data
  spatial_count <- GetAssayData(spatial_obj, slot = "counts", assay = "Spatial")
  
  location <- as.data.frame(spatial_obj@images[[1]]@boundaries$centroids@coords)
  colnames(location) <- c("x", "y")
  rownames(location) <- spatial_obj@images[[1]]@boundaries$centroids@cells
  
  ct_select <- as.character(unique(sc_meta[[ct_var]]))
  
  createCARDObject(
    sc_count = sc_count,
    sc_meta = sc_meta,
    spatial_count = spatial_count,
    spatial_location = location,
    ct.varname = ct_var,
    ct.select = ct_select,
    sample.varname = sample_var,
    minCountGene = minCountGene,
    minCountSpot = minCountSpot
  )
}

# Function to filter low proportions and remove empty columns
filter_proportions <- function(pro, threshold = 0.05) {
  pro[pro <= threshold] <- NA
  pro <- pro[, colSums(!is.na(pro)) > 0, drop = FALSE]
  pro[is.na(pro)] <- 0
  return(pro)
}

# Function to plot CARD pie chart and save PDF
plot_CARD_pie <- function(pro, location, filename, colors = NULL, radius = 9) {
  pdf(filename, width = 10, height = 10, pointsize = 12, family = "sans", bg = "white")
  if (is.null(colors)) colors <- NULL # use package defaults if not provided
  p <- CARD.visualize.pie(proportion = pro, spatial_location = location, colors = colors, radius = radius)
  print(p)
  dev.off()
}


# FUNCIONES STdeconvolve --------------------------------------------------
# Prepare corpus
prepare_corpus <- function(count_matrix, removeAbove = 0.99, removeBelow = 0.001,
                           alpha = 0.05, nTopOD = 500) {
  corpus <- restrictCorpus(count_matrix,
                           removeAbove = removeAbove,
                           removeBelow = removeBelow,
                           alpha = alpha,
                           nTopOD = nTopOD)
  t_corpus <- t(as.matrix(corpus))
  t_corpus[rowSums(t_corpus) > 0, ]
}

# Fit LDA models and select optimal
fit_optimal_lda <- function(corpus, K_range = 9:15, opt_metric = "min") {
  ldas <- fitLDA(corpus, Ks = K_range)
  optimalModel(models = ldas, opt = opt_metric)
}

# Run STdeconvolve pipeline
run_STdeconvolve <- function(spatial_obj, sc_meta, ct_var = "Final_Labels",
                             K_range = 9:15, perc_filt = 0.05, betaScale = 1000) {
  # Prepare matrices
  cd_matrix <- spatial_obj[["Spatial"]]@layers$counts
  dimnames(cd_matrix) <- list(rownames(spatial_obj), colnames(spatial_obj))
  
  cell_types <- as.factor(sc_meta[[ct_var]])
  names(cell_types) <- colnames(spatial_obj)
  
  # Corpus + LDA
  clean_corpus <- prepare_corpus(cd_matrix)
  opt_lda <- fit_optimal_lda(clean_corpus, K_range = K_range)
  
  # Extract proportions and gene expression
  results <- getBetaTheta(opt_lda, perc.filt = perc_filt, betaScale = betaScale)
  
  list(
    deconProp = results$theta,
    deconGexp = results$beta,
    cell_types = cell_types
  )
}

# Visualization
plot_STdeconvolve <- function(deconProp, spatial_coords, cell_types, filename, colors, r = 10, lwd = 0.7) {
  pdf(filename, width = 10, height = 10, pointsize = 12,
      family = "sans", bg = "white")
  p <- vizAllTopics(
    deconProp, spatial_coords,
    groups = cell_types,
    group_cols = colors,
    r = r,
    lwd = lwd
  )
  print(p)
  dev.off()
}

## ANOTAR STdeconvolve -----------------------------------------------------
# Proxy theta and proxy gene expression
create_proxy <- function(cell_types, count_matrix) {
  ProxyTheta <- model.matrix(~ 0 + cell_types)
  rownames(ProxyTheta) <- names(cell_types)
  colnames(ProxyTheta) <- gsub("^cell_types", "", colnames(ProxyTheta))
  
  ProxyGexp <- count_matrix %*% ProxyTheta
  list(ProxyTheta = ProxyTheta, ProxyGexp = ProxyGexp)
}

# Correlation matrix 
compute_correlation <- function(m1, m2, type = c("b", "t")) {
  corMtx <- getCorrMtx(m1 = as.matrix(m1), m2 = as.matrix(m2), type = type)
  rownames(corMtx) <- paste0("decon_", seq(nrow(corMtx)))
  corMtx
}

# Save plot
save_corr_plot <- function(mat, filename, colLabs, rowLabs, title) {
  pdf(filename, width = 10, height = 10, pointsize = 12,
      family = "sans", bg = "white")
  p <- correlationPlot(
    mat = mat,
    colLabs = colLabs,
    rowLabs = rowLabs,
    title = title,
    annotation = TRUE
  ) + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0))
  print(p)
  dev.off()
}








