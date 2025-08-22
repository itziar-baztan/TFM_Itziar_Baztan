# STATISTIC

# Libraries ---------------------------------------------------------------
library(dplyr)
library(tidyr)
library(purrr)
library(car)
library(crayon)
library(lme4) 

setwd("C:/Users/itzia/OneDrive/Desktop/METODOS_COMPUTACIONALES/TFM/Spatial_analysis")
# Load the functions
source("funciones/funciones_statistic_analysis.R")

# Data preparation --------------------------------------------------------
df_combined <- bind_rows(seurat_prop, card_prop, c2l_prop, stride_prop) %>%
  mutate(spot_id = c(rownames(seurat_prop), rownames(card_prop),
                     rownames(c2l_prop), rownames(stride_prop)))

DMG_pons$Regiones[is.na(DMG_pons$Regiones) | DMG_pons$Regiones == ""] <- "Tumor"

spots_tumor  <- rownames(DMG_pons@meta.data[DMG_pons$Regiones %in% c("Tumor", "Peritumoral area"), ])
spots_normal <- rownames(DMG_pons@meta.data[DMG_pons$Regiones %in% "Normal brain", ])

df_tumor  <- prepare_df(df_combined, spots_tumor)
df_normal <- prepare_df(df_combined, spots_normal)

datasets <- list(Tumor = df_tumor, Normal = df_normal)

# Model fitting -----------------------------------------------------------
modelo_tumor  <- lmer(proportion ~ method * cell_type + (1 | spot_id), data = df_tumor)
modelo_normal <- lmer(proportion ~ method * cell_type + (1 | spot_id), data = df_normal)

# Run assumption tests ----------------------------------------------------
run_assumption_tests(df_tumor, modelo_tumor, "Tumor")
run_assumption_tests(df_normal, modelo_normal, "Normal")

# Friedman tests ----------------------------------------------------------
friedman_results <- map(datasets, run_friedman)
print(friedman_results)

# Post-hoc tests ----------------------------------------------------------
posthoc_results <- map(datasets, run_posthoc)
print(posthoc_results)

# Print p-values for a given cell type

extract_and_print <- function(posthoc_res, cell_type) {
  p_valores <- posthoc_res$p_values[[which(posthoc_res$cell_type == cell_type)]]
  
  pval_matrix <- apply(p_valores, c(1, 2), function(x) {
    if (is.na(x)) return(NA)
    formatC(as.numeric(x), format = "e", digits = 4)
  })
  
  rownames(pval_matrix) <- c("CELL2LOCATION", "SEURAT", "STRIDE")
  colnames(pval_matrix) <- c("CARD", "CELL2LOCATION", "SEURAT")
  
  print_with_color(pval_matrix)
}

extract_and_print(posthoc_results$Tumor, "tumoral_cells")
extract_and_print(posthoc_results$Normal, "tumoral_cells")



