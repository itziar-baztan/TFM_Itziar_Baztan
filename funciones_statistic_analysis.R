# FUNCIONES ESTADÍSTICA

# Prepare dataset
prepare_df <- function(df, spots) {
  df %>%
    filter(spot_id %in% spots) %>%
    pivot_longer(-c("method", "spot_id"),
                 names_to = "cell_type",
                 values_to = "proportion")
}

# Normality and homogeneity of variance tests
run_assumption_tests <- function(df, modelo, name) {
  cat("\n--- Assumption tests for", name, "---\n")
  set.seed(123)
  res <- residuals(modelo)
  print(shapiro.test(sample(res, 5000)))   # Shapiro-Wilk
  print(leveneTest(proportion ~ method * cell_type, data = df))  # Homogeneity of variance
}

# Friedman test for each cell type
run_friedman <- function(df) {
  df %>%
    group_by(cell_type) %>%
    group_split() %>%
    map_dfr(function(subdf) {
      ct <- unique(subdf$cell_type)
      if (n_distinct(subdf$method) == 4 && n_distinct(subdf$spot_id) > 1) {
        test <- tryCatch(
          friedman.test(proportion ~ method | spot_id, data = subdf),
          error = function(e) NULL
        )
        if (!is.null(test)) {
          return(tibble(cell_type = ct,
                        statistic = test$statistic,
                        p_value = test$p.value))
        }
      }
      tibble(cell_type = ct, statistic = NA, p_value = NA)
    })
}

# Wilcoxon post-hoc tests
run_posthoc <- function(df) {
  df %>%
    group_by(cell_type) %>%
    summarise(
      p_values = list(
        pairwise.wilcox.test(
          proportion, method,
          paired = TRUE,
          p.adjust.method = "bonferroni",
          exact = FALSE
        )$p.value
      ),
      .groups = "drop"
    )
}

# Print matrix with colors
print_with_color <- function(mat) {
  cat("              ", paste(sprintf("%-15s", colnames(mat)), collapse = ""), "\n")
  for (i in 1:nrow(mat)) {
    cat(sprintf("%-15s", rownames(mat)[i]))
    for (j in 1:ncol(mat)) {
      val <- mat[i, j]
      if (is.na(val)) {
        cat(sprintf("%-15s", "NA"))
      } else if (as.numeric(val) < 0.05) {
        cat(red(sprintf("%-15s", val)))
      } else {
        cat(sprintf("%-15s", val))
      }
    }
    cat("\n")
  }
}