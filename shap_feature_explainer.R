#!/usr/bin/env Rscript

# SHAP analysis for XGBoost models using SHAPforxgboost
# Installs the package if missing, then fits an example model and generates SHAP plots.

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

# Helper to ensure package is installed
ensure_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = 'https://cloud.r-project.org')
  }
}

# Ensure SHAPforxgboost is installed (from CRAN if available, fallback to GitHub)
if (!requireNamespace('SHAPforxgboost', quietly = TRUE)) {
  try({
    install.packages('SHAPforxgboost', repos = 'https://cloud.r-project.org')
  }, silent = TRUE)
}

if (!requireNamespace('SHAPforxgboost', quietly = TRUE)) {
  ensure_pkg('devtools')
  devtools::install_github('liuyanguu/SHAPforxgboost')
}

suppressPackageStartupMessages({
  library(SHAPforxgboost)
  library(xgboost)
  library(ggplot2)
})

# CLI
opt_list <- list(
  make_option('--csv', type='character', default=NA, help='Path to CSV with features/outcome'),
  make_option('--outcome', type='character', default=NA, help='Outcome column name for regression'),
  make_option('--outdir', type='character', default='output/shap', help='Output directory'),
  make_option('--rounds', type='integer', default=200, help='Training rounds (nrounds)'),
  make_option('--eta', type='double', default=0.02),
  make_option('--max_depth', type='integer', default=10),
  make_option('--gamma', type='double', default=0.01),
  make_option('--subsample', type='double', default=0.98),
  make_option('--colsample_bytree', type='double', default=0.86)
)
opt <- parse_args(OptionParser(option_list = opt_list))

outdir <- opt$outdir
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

if (!is.na(opt$csv) && !is.na(opt$outcome)) {
  # Use provided dataset
  dt <- data.table::fread(opt$csv)
  if (!(opt$outcome %in% names(dt))) {
    stop('Outcome column not found: ', opt$outcome)
  }
  y_var <- opt$outcome
  X <- as.matrix(dt[, setdiff(names(dt), y_var), with=FALSE])
  y <- as.matrix(dt[[y_var]])
} else {
  # Fall back to built-in example from SHAPforxgboost vignette-like data
  set.seed(1)
  n <- 1000
  dt <- data.table::data.table(
    dayint = sample(1:365, n, replace = TRUE),
    Column_WV = rnorm(n),
    AOT_Uncertainty = runif(n),
    diffcwv = rnorm(n) + 0.5 * rnorm(n)
  )
  y_var <- 'diffcwv'
  X <- as.matrix(dt[, setdiff(names(dt), y_var), with=FALSE])
  y <- as.matrix(dt[[y_var]])
}

params <- list(
  objective = 'reg:squarederror',
  eta = opt$eta,
  max_depth = opt$max_depth,
  gamma = opt$gamma,
  subsample = opt$subsample,
  colsample_bytree = opt$colsample_bytree
)

mod <- xgboost::xgboost(
  data = X,
  label = y,
  params = params,
  nrounds = opt$rounds,
  verbose = FALSE,
  early_stopping_rounds = 8
)

# Compute SHAP values
shap_values <- shap.values(xgb_model = mod, X_train = X)
# Ranked features by mean |SHAP|
ranked <- data.frame(feature = names(shap_values$mean_shap_score),
                     mean_abs_shap = as.numeric(shap_values$mean_shap_score))
ranked <- ranked[order(-ranked$mean_abs_shap), ]
write.csv(ranked, file = file.path(outdir, 'shap_ranked_features.csv'), row.names = FALSE)

# Long format
data_long <- shap.prep(xgb_model = mod, X_train = X)
# Alternative using given shap_contrib
# data_long <- shap.prep(shap_contrib = shap_values$shap_score, X_train = X)

# Summary plot
png(file.path(outdir, 'shap_summary.png'), width = 1000, height = 800)
print(shap.plot.summary(data_long))
dev.off()

# Lighter summary for preview
png(file.path(outdir, 'shap_summary_light.png'), width = 1000, height = 800)
print(shap.plot.summary(data_long, x_bound = 1.2, dilute = 10))
dev.off()

# Dependence plot for the top feature
top_feat <- ranked$feature[1]
png(file.path(outdir, paste0('shap_dependence_', top_feat, '.png')), width = 900, height = 700)
print(shap.plot.dependence(data_long = data_long, x = top_feat))
dev.off()

# If at least two features, color by second feature
if (nrow(ranked) >= 2) {
  color_feat <- ranked$feature[2]
  png(file.path(outdir, paste0('shap_dependence_', top_feat, '_color_', color_feat, '.png')), width = 900, height = 700)
  print(shap.plot.dependence(data_long = data_long, x = top_feat, color_feature = color_feat))
  dev.off()
}

# Interaction (can be slow)
try({
  shap_int <- shap.prep.interaction(xgb_mod = mod, X_train = X)
  if (nrow(ranked) >= 2) {
    png(file.path(outdir, paste0('shap_interaction_', ranked$feature[1], '_', ranked$feature[2], '.png')), width = 1000, height = 800)
    print(shap.plot.dependence(data_long = data_long, data_int = shap_int, x = ranked$feature[1], y = ranked$feature[2], color_feature = ranked$feature[2]))
    dev.off()
  }
}, silent = TRUE)

# Force plots (stacked)
plot_data <- shap.prep.stack.data(shap_contrib = shap_values$shap_score, top_n = min(4, ncol(X)), n_groups = 6)

png(file.path(outdir, 'shap_force_plot.png'), width = 1600, height = 900)
print(shap.plot.force_plot(plot_data, zoom_in_location = 500, y_parent_limit = c(-1,1)))
dev.off()

png(file.path(outdir, 'shap_force_plot_bygroup.png'), width = 1600, height = 900)
print(shap.plot.force_plot_bygroup(plot_data))
dev.off()

cat('✅ SHAP analysis complete. Outputs in: ', outdir, '\n')
