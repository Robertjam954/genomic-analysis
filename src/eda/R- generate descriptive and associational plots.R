#!/usr/bin/env Rscript
"""
Produce generate exploratory (outcome-related) plots and tables (Table 1 and Table 2) from cleaned data.

This script expects the cleaned data RDS produced by `scripts/prepare_descriptive_vars.R`
and will write plots and tables to `output/descriptive` by default.

Usage:
  Rscript scripts/analysis_tables_and_plots.R --cleaned output/descriptive/merged_genie_cleaned.rds --outdir output/descriptive
"""

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(ggplot2)
  library(naniar)
  library(GGally)
  library(gtsummary)
  library(broom)
  library(readr)
  library(forcats)
})

option_list <- list(
  make_option(c("--cleaned"), type = "character", default = "output/descriptive/merged_genie_cleaned.rds",
              help = "Path to cleaned RDS produced earlier"),
  make_option(c("--outdir"), type = "character", default = "output/descriptive",
              help = "Output directory")
)

opt <- parse_args(OptionParser(option_list = option_list))
cleaned_path <- opt$cleaned
outdir <- opt$outdir
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

if (!file.exists(cleaned_path)) {
  stop('Cleaned data not found at ', cleaned_path, '. Run scripts/prepare_descriptive_vars.R first.')
}

df <- readRDS(cleaned_path)
message('Loaded cleaned data: rows=', nrow(df), ' cols=', ncol(df))

# 1) Missingness visualization (naniar)
png(file.path(outdir, 'missingness_plot.png'), width = 1200, height = 800)
try({
  vis_miss(df) + ggtitle('Missingness heatmap')
}, silent = TRUE)
dev.off()

# 2) Numeric correlation matrix (GGally) for numeric vars
num_df <- df %>% select(where(is.numeric))
if (ncol(num_df) >= 2) {
  p_corr <- try(ggcorr(num_df, label = TRUE), silent = TRUE)
  if (!inherits(p_corr, 'try-error')) {
    ggsave(file.path(outdir, 'numeric_correlation.png'), plot = p_corr, width = 10, height = 8)
  }
}

# Variant allele count distribution and boxplots (cohort-wide and by DS)
va_col <- c('variant_allele_count', 'variant_allele_freq', 'mutant_allele_count')
va_col <- va_col[va_col %in% names(df)][1]
if (!is.null(va_col)) {
  df <- df %>% mutate(!!sym(va_col) := as.numeric(.data[[va_col]]))
  # cohort-wide histogram + boxplot layout
  png(file.path(outdir, 'variant_allele_count_layout.png'), width = 1200, height = 800)
  layout(matrix(c(1,1,2,3), nrow = 2, byrow = TRUE))
  hist(df[[va_col]], breaks = 30, border = NA, col = rgb(0.1,0.8,0.3,0.5), xlab = va_col, main = 'Distribution')
  boxplot(df[[va_col]], xlab = va_col, col = rgb(0.8,0.8,0.3,0.5))
  # by DS
  if ('dfs_status' %in% names(df)) {
    boxplot(split(df[[va_col]], df$dfs_status), xlab = 'DFS status', col = c('lightblue','pink'))
  }
  dev.off()

  # Quartile boxplot + density by dfs_status
  if ('dfs_status' %in% names(df)) {
    p <- ggplot(df, aes(x = dfs_status, y = .data[[va_col]], fill = dfs_status)) +
      geom_boxplot(alpha = 0.6) + theme_minimal() + xlab('DFS status') + ylab(va_col)
    ggsave(file.path(outdir, 'variant_by_dfs_boxplot.png'), p, width = 8, height = 6)
  }
}

# Probability of DS=1 by demographic variables (age_cat, race, ethnicity)
demo_vars <- c('age_cat','race','ethnicity')
demo_vars <- demo_vars[demo_vars %in% names(df)]
if (length(demo_vars) > 0 && 'DS_status_bin' %in% names(df)) {
  for (v in demo_vars) {
    p <- df %>% filter(!is.na(.data[[v]])) %>%
      group_by(.data[[v]]) %>%
      summarise(n = n(), pct_ds1 = 100 * mean(DS_status_bin, na.rm = TRUE)) %>%
      ggplot(aes(x = reorder(.data[[v]], pct_ds1), y = pct_ds1, fill = pct_ds1)) +
      geom_col() + coord_flip() + theme_minimal() + xlab(v) + ylab('Probability of DS=1 (%)')
    ggsave(file.path(outdir, paste0('ds1_by_', v, '.png')), p, width = 8, height = 6)
  }
}

# Stacked bar for subtype by DFS (fill = subtype, position = fill)
if ('subtype' %in% names(df) && 'dfs_status' %in% names(df)) {
  p <- ggplot(df) + geom_bar(aes(x = dfs_status, fill = subtype), position = 'fill') + theme_minimal() + ylab('Proportion')
  ggsave(file.path(outdir, 'subtype_by_dfs_stacked.png'), p, width = 8, height = 6)
}

# Gene prevalence plot: prevalence in DS=1 (left) and DS=0 (right)
if (file.exists(file.path(outdir, 'gene_prevalence.csv')) && 'dfs_status' %in% names(df)) {
  gp <- read.csv(file.path(outdir, 'gene_prevalence.csv'))
  # compute prevalence by status
  hugo_col <- intersect(c('hugo_symbol','hugo_symbol_list','hugo'), names(df))[1]
  if (!is.na(hugo_col)) {
    longg <- df %>% select(id = row_number(), dfs_status, !!sym(hugo_col)) %>%
      filter(!is.na(.data[[hugo_col]])) %>% mutate(tmp = str_replace_all(.data[[hugo_col]], ',', ';')) %>%
      separate_rows(tmp, sep = ';') %>% mutate(gene = str_trim(tmp)) %>% distinct(id, dfs_status, gene)
    prev_by_status <- longg %>% group_by(dfs_status, gene) %>% summarise(n = n()) %>% group_by(dfs_status) %>% mutate(pct = 100 * n / sum(n))
    # wide table for plotting
    wide <- prev_by_status %>% pivot_wider(names_from = dfs_status, values_from = pct, values_fill = 0)
    # plot top genes side-by-side (DS=1 left, DS=0 right)
    top_genes <- readLines(file.path(outdir, 'top5_genes.txt'))
    wide_top <- wide %>% filter(gene %in% top_genes) %>% replace(is.na(.), 0)
    p <- wide_top %>% pivot_longer(-gene, names_to = 'dfs_status', values_to = 'pct') %>%
      ggplot(aes(x = reorder(gene, pct), y = pct, fill = dfs_status)) + geom_col(position = 'dodge') + coord_flip() + theme_minimal()
    ggsave(file.path(outdir, 'top_genes_by_dfs.png'), p, width = 8, height = 6)
  }
}

# Histogram of predicted probability (example) by mutation covariates if a model exists
pred_file <- file.path(outdir, 'multivar_glm_top5_vs_ds.csv')
if (file.exists(pred_file)) {
  # if glm exists, create predicted probabilities and facet by top genes
  tbl <- read.csv(pred_file)
  # we can't recreate predictions from coefficients here without refitting; skip
}

# Use fct_count-like summaries for Table 1 variables
if (length(table1_vars_present) > 0) {
  fct_counts <- map_dfr(table1_vars_present, function(v) {
    if (v %in% names(df)) {
      tibble(variable = v, counts = list(df %>% count(!!sym(v))))
    } else tibble(variable = v, counts = list(NULL))
  })
  saveRDS(fct_counts, file.path(outdir, 'table1_factor_counts.rds'))
}

# 3) Simple regression diagnostics example
# Ensure there is a binary DFS status coded as factor with levels (Recurrence..., Disease Free) or 0/1
if ('dfs_status' %in% names(df)) {
  # create a numeric binary if needed
  df <- df %>% mutate(dfs_status_bin = case_when(
    dfs_status %in% c('Recurrence or Progression of Disease') ~ 1,
    dfs_status %in% c('Disease Free') ~ 0,
    TRUE ~ as.numeric(NA)
  ))
}

# Example linear model: substitute with appropriate variables
if ('top_5_gene_status' %in% names(df) && 'dfs_status_bin' %in% names(df)) {
  fit <- try(lm(dfs_status_bin ~ top_5_gene_status, data = df), silent = TRUE)
  if (!inherits(fit, 'try-error')) {
    sink(file.path(outdir, 'lm_top5_vs_dfs_summary.txt'))
    print(summary(fit))
    sink()
    # diagnostic plots
    png(file.path(outdir, 'lm_diagnostics.png'), width = 1200, height = 1000)
    par(mfrow = c(2,2))
    plot(fit)
    dev.off()
  }
}

# 4) Tables
# Table 1: descriptive characteristics. Choose a standard set of variables (modify as needed)
table1_vars <- c('age', 'age_cat', 'sex', 'race', 'ethnicity', 'subtype',
                 'ajcc_pathologic_tumor_stage', 'tbl_score_quartile', 'buffa_hypoxia_quartile',
                 'aneuploidy_score_quartile', 'mantis_bin')
table1_vars_present <- table1_vars[table1_vars %in% names(df)]

if (length(table1_vars_present) > 0) {
  # If dfs_status exists, stratify by it
  by <- if ('dfs_status' %in% names(df)) 'dfs_status' else NULL
  tbl1 <- df %>% select(all_of(unique(c(by, table1_vars_present)))) %>%
    gtsummary::tbl_summary(by = by, missing = 'no') %>%
    gtsummary::add_overall() %>%
    gtsummary::modify_header(label = '**Variable**')

  # save as html and csv
  try(gtsave(tbl1, filename = file.path(outdir, 'table1.html')))
  write.csv(tbl1$table_body %>% as.data.frame(), file.path(outdir, 'table1_table_body.csv'), row.names = FALSE)
  message('Wrote Table 1 outputs to ', outdir)
}

# Table 2: association between top_5_gene_status and DFS_status (adjusted)
if ('top_5_gene_status' %in% names(df) && 'dfs_status_bin' %in% names(df)) {
  # select confounders present
  confs <- c('age', 'age_cat', 'sex', 'race', 'subtype', 'ajcc_pathologic_tumor_stage')
  confs_present <- confs[confs %in% names(df)]
  form <- as.formula(paste('dfs_status_bin ~ top_5_gene_status', paste(confs_present, collapse = ' + '), sep = ifelse(length(confs_present)>0, ' + ', '')))
  fit_glm <- try(glm(form, data = df, family = binomial(link = 'logit')), silent = TRUE)
  if (!inherits(fit_glm, 'try-error')) {
    tbl2 <- broom::tidy(fit_glm, exponentiate = TRUE, conf.int = TRUE)
    write.csv(tbl2, file.path(outdir, 'table2_glm_top5_vs_dfs.csv'), row.names = FALSE)
    message('Wrote Table 2 (glm) to ', outdir, '/table2_glm_top5_vs_dfs.csv')
  }
}

message('analysis_tables_and_plots.R finished; outputs in ', outdir)

