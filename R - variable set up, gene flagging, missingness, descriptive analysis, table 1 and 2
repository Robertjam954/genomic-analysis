#!/usr/bin/env Rscript

# Comprehensive Genomic Data Analysis Script
# Combines descriptive variables preparation, mutation analysis, and visualization
# Creates cleaned datasets with properly coded factors and comprehensive analysis outputs

suppressPackageStartupMessages({
  library(optparse)
  library(readxl)
  library(readr)
  library(dplyr)
  library(stringr)
  library(janitor)
  library(skimr)
  library(purrr)
  library(forcats)
  library(tidyr)
  library(tidyverse)
  library(gtsummary)
  library(knitr)
  library(kableExtra)
  library(tfrmt)
  library(naniar)      # For missing data visualization
  library(GGally)      # For correlation matrix
  library(ggplot2)     # For plots
  library(broom)       # For tidy model outputs
  library(viridisLite) # For color palettes
  library(scales)      # For formatting
})

# ===== Configuration and Setup =====
message("📂 Working directory: ", getwd())

# Create output directory
outdir <- "output_descriptive_vars"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

# ===== Data Loading =====
# Load main dataset
df_raw <- read.csv("~/Research/Projects/genomics_brain_mets_genie_bpc/tcga_gene_alteration_dfs_os_xgb/merged_tcga_df.csv", na = c('NA', ''))
colnames_original <- colnames(df_raw)
message("📊 Loaded dataset with ", nrow(df_raw), " rows and ", ncol(df_raw), " columns")

# Load breast cancer dataset if available
if (file.exists("complete_genie_bpc_dataset_deid.xlsx")) {
  bc_data <- read_xlsx("complete_genie_bpc_dataset_deid.xlsx")
  message("📊 Loaded breast cancer dataset with ", nrow(bc_data), " rows and ", ncol(bc_data), " columns")
} else {
  message("⚠️ Breast cancer dataset not found, skipping BC-specific analysis")
  bc_data <- NULL
}

# ===== Data Cleaning and Preparation =====
# Clean column names
df <- clean_names(df_raw) %>% glimpse()

# ===== Factor Creation and Coding =====

# AJCC stage: Convert to ordered factor
ajcc_levels <- c("Stage I", "Stage II", "Stage III", "Stage IV")
df <- df %>%
  mutate(
    ajcc_pathologic_tumor_stage = as.character(ajcc_pathologic_tumor_stage),
    ajcc_pathologic_tumor_stage = case_when(
      ajcc_pathologic_tumor_stage == "1" ~ "Stage I",
      ajcc_pathologic_tumor_stage == "2" ~ "Stage II",
      ajcc_pathologic_tumor_stage == "3" ~ "Stage III",
      ajcc_pathologic_tumor_stage == "4" ~ "Stage IV",
      TRUE ~ NA_character_
    ),
    ajcc_pathologic_tumor_stage = factor(ajcc_pathologic_tumor_stage, levels = ajcc_levels, ordered = TRUE)
  )

# Subtype: Create ordered factor
subtype_levels <- c("Luminal A", "Luminal B", "HER-2 enriched", "Basal-type (TNBC)")
df <- df %>%
  mutate(
    subtype = as.character(subtype),
    subtype = case_when(
      subtype == "1" ~ "Luminal A",
      subtype == "2" ~ "Luminal B",
      subtype == "3" ~ "HER-2 enriched",
      subtype == "4" ~ "Basal-type (TNBC)",
      TRUE ~ NA_character_
    ),
    subtype = factor(subtype, levels = subtype_levels, ordered = TRUE)
  )

# Age category: Create ordered factor
age_levels <- c('<50', '50-75', 'Over 75')
if ('age_cat' %in% names(df)) {
  df <- df %>%
    mutate(age_cat = factor(age_cat, levels = age_levels, ordered = TRUE))
}

# Race/Ethnicity factors
if ('race' %in% names(df)) {
  race_levels <- c('White or Caucasian', 'Black or African American', 'Asian', 'Indian or Alaska Native')
  df <- df %>% mutate(race = factor(race, levels = race_levels))
}

if ('ethnicity' %in% names(df)) {
  eth_levels <- c('Non-Hispanic', 'Hispanic or Latino')
  df <- df %>% mutate(ethnicity = factor(ethnicity, levels = eth_levels))
}

# Sex: normalize and create binary indicator
sex_candidates <- intersect(tolower(names(df)), c('sex', 'gender'))
sex_col <- if (length(sex_candidates) > 0) names(df)[which(tolower(names(df)) %in% sex_candidates)[1]] else NULL

if (!is.null(sex_col)) {
  message('Using sex column: ', sex_col)
  df <- df %>%
    mutate(
      sex = case_when(
        tolower(as.character(.data[[sex_col]])) %in% c('female','f') ~ 'Female',
        tolower(as.character(.data[[sex_col]])) %in% c('male','m') ~ 'Male',
        TRUE ~ NA_character_
      ),
      sex = factor(sex, levels = c('Female','Male')),
      sex_bin = case_when(
        sex == 'Male' ~ 1L,
        sex == 'Female' ~ 0L,
        TRUE ~ NA_integer_
      )
    )
}

# DFS status
dfs_status_col <- NULL
for (candidate in c('dfs_status', 'disease_free_status', 'recurrence_status')) {
  if (candidate %in% names(df)) {
    dfs_status_col <- candidate
    break
  }
}

if (!is.null(dfs_status_col)) {
  dfs_levels <- c('Recurrence or Progression of Disease', 'Disease Free')
  df <- df %>% 
    mutate(
      dfs_status = case_when(
        .data[[dfs_status_col]] == 1 ~ dfs_levels[1],
        .data[[dfs_status_col]] == 0 ~ dfs_levels[2],
        TRUE ~ as.character(.data[[dfs_status_col]])
      ),
      dfs_status = factor(dfs_status, levels = dfs_levels),
      # Create binary and factor versions for plotting
      dfs_raw = .data[[dfs_status_col]],
      dfs_bin = case_when(
        is.numeric(dfs_raw) ~ as.integer(dfs_raw),
        grepl("1", as.character(dfs_raw)) ~ 1L,
        grepl("yes|recurrence|progression", tolower(as.character(dfs_raw))) ~ 1L,
        TRUE ~ 0L
      ),
      dfs_f = factor(dfs_bin, levels = c(0,1), labels = c("DFS = 0", "DFS = 1"))
    )
}

# MANTIS/MSI factor
if ('mantis_bin' %in% names(df)) {
  mantis_levels <- c("MSI-H", "MSI-L")
  df <- df %>%
    mutate(
      mantis_bin = case_when(
        mantis_bin == "MSI-H" ~ ">0.4",
        mantis_bin == "MSI-L" ~ "<0.4",
        TRUE ~ NA_character_
      ),
      mantis_bin = factor(mantis_bin, levels = mantis_levels)
    )
}

# TBL score quartiles
if ('tbl_score' %in% names(df)) {
  df <- df %>% mutate(tbl_score = as.numeric(tbl_score))
  tbl_quantile <- quantile(df$tbl_score, probs = c(0, .25, .5, .75, 1), na.rm = TRUE)
  df <- df %>%
    mutate(
      tbl_score_quartile = cut(tbl_score, breaks = tbl_quantile, include.lowest = TRUE, labels = c('Q1','Q2','Q3','Q4'))
    )
}

# BUFFA hypoxia quartiles
if ('buffa_hypoxia_score' %in% names(df)) {
  df <- df %>% mutate(buffa_hypoxia_score = as.numeric(buffa_hypoxia_score))
  qh <- quantile(df$buffa_hypoxia_score, probs = c(0, .25, .5, .75, 1), na.rm = TRUE)
  df <- df %>% 
    mutate(buffa_hypoxia_quartile = cut(buffa_hypoxia_score, breaks = qh, include.lowest = TRUE, labels = c('Q1','Q2','Q3','Q4')))
}

# Aneuploidy quartiles
aneuploidy_col <- NULL
for (candidate in c('aneuploidy_score', 'aneuploidy_score')) {
  if (candidate %in% names(df)) {
    aneuploidy_col <- candidate
    break
  }
}

if (!is.null(aneuploidy_col)) {
  df <- df %>% mutate(!!sym(aneuploidy_col) := as.numeric(.data[[aneuploidy_col]]))
  qa <- quantile(df[[aneuploidy_col]], probs = c(0, .25, .5, .75, 1), na.rm = TRUE)
  df <- df %>% 
    mutate(aneuploidy_score_quartile = cut(.data[[aneuploidy_col]], breaks = qa, include.lowest = TRUE, labels = c('Q1','Q2','Q3','Q4')))
}

# Variant type
if ('variant_type' %in% names(df)) {
  variant_type_levels <- c('SNP','INS','ONP','DEL')
  df <- df %>% mutate(variant_type = factor(variant_type, levels = variant_type_levels))
}

# Variant allele count quartiles
var_acol <- NULL
for (candidate in c('variant_allele_count', 'variant_allele_freq', 'mutant_allele_count', 'mutant_allele_freq')) {
  if (candidate %in% names(df)) { 
    var_acol <- candidate
    break 
  }
}

if (!is.null(var_acol)) {
  df <- df %>% mutate(!!sym(var_acol) := as.numeric(.data[[var_acol]]))
  vq <- quantile(df[[var_acol]], probs = c(0, .25, .5, .75, 1), na.rm = TRUE)
  df <- df %>% 
    mutate(variant_allele_count_quartile = cut(.data[[var_acol]], breaks = vq, include.lowest = TRUE, labels = c('Q1','Q2','Q3','Q4')))
  message('Created variant_allele_count_quartile from ', var_acol)
}

# ===== Mutation Counts and Top-Gene Flags =====

# Find Hugo symbol column
hugo_col <- NULL
for (candidate in c('hugo_symbol', 'Hugo_Symbol', 'gene', 'gene_symbol')) {
  if (candidate %in% names(df)) {
    hugo_col <- candidate
    break
  }
}

if (!is.null(hugo_col)) {
  # Normalize Hugo symbol column to character
  df <- df %>%
    mutate(hugo_symbol = as.character(.data[[hugo_col]]))
  
  # Per-sample mutation count
  df <- df %>%
    mutate(
      mutation_list = case_when(
        is.na(hugo_symbol) ~ NA_character_,
        TRUE ~ str_replace_all(hugo_symbol, "\\s+", "")
      ),
      mutation_count = case_when(
        is.na(mutation_list) | mutation_list == "" ~ 0L,
        TRUE ~ as.integer(str_count(mutation_list, ";") + str_count(mutation_list, ",") + 1)
      )
    )
  
  # More robust splitting: split on ; or ,
  long_genes <- df %>%
    group_by(merge_key = row_number(), mutation_list) %>%
    filter(!is.na(mutation_list) & mutation_list != '') %>%
    mutate(mutation_list = str_replace_all(mutation_list, ',', ';')) %>%
    separate_rows(mutation_list, sep = ';') %>%
    mutate(gene = na_if(mutation_list, '')) %>%
    filter(!is.na(gene)) %>%
    mutate(gene = str_trim(gene))
  
  # Compute gene prevalence (per unique sample)
  gene_prev <- long_genes %>%
    distinct(merge_key, gene) %>%
    count(gene, name = 'n_samples') %>%
    mutate(n_total = nrow(df), prevalence = n_samples / n_total) %>%
    arrange(desc(prevalence))
  
  # Top 10 and top 5 genes
  top10_genes <- gene_prev %>% slice_head(n = 10) %>% pull(gene)
  top5_genes <- gene_prev %>% slice_head(n = 5) %>% pull(gene)
  
  # Create presence flags for top10 genes
  if (length(top10_genes) > 0) {
    for (g in top10_genes) {
      colname <- paste0('G__', make.names(g))
      df[[colname]] <- as.integer(str_detect(df[["hugo_symbol"]] %||% '', fixed(g)))
    }
  }
  
  # Combined top-5 flag
  if (length(top5_genes) > 0) {
    df <- df %>% 
      mutate(top_5_gene_status = as.integer(
        reduce(map(top5_genes, ~ str_detect(.data[["hugo_symbol"]] %||% '', fixed(.x))), `|`, .init = FALSE)
      ))
  } else {
    df <- df %>% mutate(top_5_gene_status = 0L)
  }
  
  # Save gene prevalence and top lists
  write.csv(gene_prev, file = file.path(outdir, 'gene_prevalence.csv'), row.names = FALSE)
  writeLines(top10_genes, file = file.path(outdir, 'top10_genes.txt'))
  writeLines(top5_genes, file = file.path(outdir, 'top5_genes.txt'))
  
  # Group by DFS status and summarize top_5 presence
  if (!is.null(dfs_status_col)) {
    top5_by_status <- df %>%
      group_by(.data[[dfs_status_col]]) %>%
      summarise(
        n = n(),
        top5_count = sum(top_5_gene_status, na.rm = TRUE),
        top5_pct = 100 * top5_count / n
      ) %>%
      arrange(desc(top5_pct))
    
    write.csv(top5_by_status, file = file.path(outdir, 'top5_by_status.csv'), row.names = FALSE)
    message('Wrote top5_by_status.csv')
  }
}

# ===== Missingness Analysis =====
missing_summary <- df %>%
  summarise(across(everything(), ~ sum(is.na(.)))) %>%
  pivot_longer(everything(), names_to = 'variable', values_to = 'n_missing') %>%
  mutate(n_total = nrow(df), prop_missing = n_missing / n_total)

write.csv(missing_summary, file = file.path(outdir, 'missing_summary.csv'), row.names = FALSE)
message('Wrote missing_summary.csv')

# ===== Skim Summary =====
sk <- skim(df)
write.csv(as.data.frame(sk), file = file.path(outdir, 'skim_summary.csv'), row.names = FALSE)
message('Wrote skim summary (csv)')

# ===== Visualization Functions =====

# Theme for plots
theme_msk <- function(base = 10) {
  theme_minimal(base_size = base) +
    theme(
      panel.grid.minor = element_blank(),
      panel.spacing = unit(0.6, "lines"),
      strip.text = element_text(face = "bold"),
      legend.position = "right"
    )
}

# Safe factor with explicit unknown level "."
as_factor_with_dot <- function(x, lvl = NULL, ordered = FALSE) {
  x <- as.character(x)
  x[is.na(x) | trimws(x) == ""] <- "."
  if (is.null(lvl)) return(factor(x, ordered = ordered))
  factor(x, levels = c(lvl, if (!"." %in% lvl) "."), ordered = ordered)
}

# ===== Exploratory Analysis =====

# Clean numeric subset for correlation analysis
numeric_vars <- df %>% select(where(is.numeric))

# Correlation matrix plot
if (ncol(numeric_vars) > 1) {
  png(file.path(outdir, 'correlation_matrix.png'), width = 800, height = 800)
  print(ggcorr(numeric_vars, label = TRUE))
  dev.off()
}

# Visualize missingness
png(file.path(outdir, 'missingness_plot.png'), width = 800, height = 600)
print(gg_miss_var(df))
dev.off()

# Missingness pattern
png(file.path(outdir, 'missingness_pattern.png'), width = 800, height = 600)
print(vis_miss(df))
dev.off()

# ===== Statistical Analysis =====

if (!is.null(dfs_status_col) && 'top_5_gene_status' %in% names(df)) {
  
  # Univariate linear regression
  fit_uni <- lm(dfs_bin ~ top_5_gene_status, data = df)
  summary(fit_uni)
  tidy(fit_uni) %>% write.csv(file.path(outdir, 'univariate_lm.csv'), row.names = FALSE)
  
  # Multivariate regression with possible confounders
  confounders <- intersect(names(df), c('sex_bin', 'age_cat', 'ajcc_pathologic_tumor_stage', 
                                        'buffa_hypoxia_quartile', 'aneuploidy_score_quartile'))
  
  if (length(confounders) > 0) {
    formula_str <- paste("dfs_bin ~ top_5_gene_status +", paste(confounders, collapse = " + "))
    fit_multi <- lm(as.formula(formula_str), data = df)
    summary(fit_multi)
    tidy(fit_multi) %>% write.csv(file.path(outdir, 'multivariable_lm.csv'), row.names = FALSE)
    
    # Residual diagnostics
    png(file.path(outdir, 'residuals_diagnostics.png'), width = 800, height = 800)
    par(mfrow = c(2, 2))
    plot(fit_multi)
    dev.off()
  }
}

# ===== Advanced Visualizations =====

# Variant Allele Count Analysis
if (!is.null(var_acol) && var_acol %in% names(df)) {
  # Histogram of variant allele count
  p_hist_vac <- ggplot(df, aes(x = .data[[var_acol]])) +
    geom_histogram(bins = 30, alpha = 0.8, fill = "skyblue") +
    labs(x = "Variant Allele Count", y = "Frequency", title = "Distribution of Variant Allele Count") +
    theme_msk()
  
  ggsave(file.path(outdir, 'plot_vac_hist.png'), p_hist_vac, width = 8, height = 6, dpi = 300)
  
  # Boxplot by DFS status if available
  if ('dfs_f' %in% names(df)) {
    p_box_vac <- ggplot(df, aes(x = dfs_f, y = .data[[var_acol]], fill = dfs_f)) +
      geom_boxplot(outlier.alpha = 0.4, show.legend = FALSE) +
      labs(x = "DFS Status", y = "Variant Allele Count", title = "Variant Allele Count by DFS Status") +
      theme_msk()
    
    ggsave(file.path(outdir, 'plot_vac_by_dfs.png'), p_box_vac, width = 8, height = 6, dpi = 300)
  }
}

# Probability plots by demographics
plot_grouped_prop <- function(data, group_var, out_file) {
  if (!'dfs_bin' %in% names(data) || !group_var %in% names(data)) return()
  
  p <- data %>%
    group_by(.data[[group_var]]) %>%
    summarise(
      n = n(),
      prop_dfs1 = mean(dfs_bin == 1, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    filter(!is.na(.data[[group_var]])) %>%
    mutate(group_factor = fct_reorder(as.factor(.data[[group_var]]), prop_dfs1, .desc = TRUE)) %>%
    ggplot(aes(x = group_factor, y = prop_dfs1)) +
    geom_col(fill = "steelblue", alpha = 0.7) +
    geom_text(aes(label = percent(prop_dfs1, accuracy = 0.1)), vjust = -0.3, size = 3) +
    labs(x = NULL, y = "P(DFS = 1)", title = paste("Probability of DFS = 1 by", group_var)) +
    scale_y_continuous(labels = percent_format()) +
    theme_msk() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(file.path(outdir, out_file), p, width = 8, height = 6, dpi = 300)
}

# Create demographic plots if variables exist
if ('dfs_bin' %in% names(df)) {
  demo_vars <- intersect(names(df), c('race', 'ethnicity', 'subtype', 'sex'))
  for (var in demo_vars) {
    plot_grouped_prop(df, var, paste0('plot_prob_dfs1_by_', var, '.png'))
  }
}

# Gene prevalence plots
if (exists('gene_prev') && nrow(gene_prev) > 0 && 'dfs_f' %in% names(df)) {
  # Create long format for gene analysis
  if (!is.null(hugo_col) && exists('long_genes')) {
    gene_dfs_data <- long_genes %>%
      left_join(df %>% select(merge_key = row_number(), dfs_f), by = "merge_key") %>%
      group_by(gene, dfs_f) %>%
      summarise(n_with_gene = n(), .groups = 'drop') %>%
      left_join(
        df %>% group_by(dfs_f) %>% summarise(total_n = n(), .groups = 'drop'),
        by = "dfs_f"
      ) %>%
      mutate(prevalence = n_with_gene / total_n) %>%
      filter(gene %in% top10_genes)
    
    # Plot by DFS status
    for (dfs_level in unique(gene_dfs_data$dfs_f)) {
      if (is.na(dfs_level)) next
      
      p <- gene_dfs_data %>%
        filter(dfs_f == dfs_level) %>%
        arrange(prevalence) %>%
        mutate(gene = factor(gene, levels = gene)) %>%
        ggplot(aes(x = gene, y = prevalence)) +
        geom_segment(aes(xend = gene, yend = 0), color = "gray50") +
        geom_point(size = 3, color = "orange") +
        coord_flip() +
        labs(x = NULL, y = "Prevalence", title = paste("Gene Mutation Prevalence —", dfs_level)) +
        theme_msk()
      
      ggsave(file.path(outdir, paste0('plot_gene_prev_', gsub("[^A-Za-z0-9]", "_", dfs_level), '.png')), 
             p, width = 8, height = 6, dpi = 300)
    }
  }
}

# ===== Breast Cancer Specific Analysis =====
if (!is.null(bc_data)) {
  message("🔬 Processing breast cancer specific analysis...")
  
  # Create derived variables for BC data
  bc_data <- bc_data %>%
    mutate(
      # Age category
      age_category = case_when(
        age_at_diagnosis < 50 ~ "<50",
        age_at_diagnosis >= 50 & age_at_diagnosis <= 70 ~ "50-70",
        age_at_diagnosis > 70 ~ ">70",
        TRUE ~ NA_character_
      ),
      age_category = factor(age_category, levels = c("<50", "50-70", ">70")),
      
      # Race factor
      race = fct_explicit_na(factor(race, levels = c("White", "Black or African American", "Asian", "Other")), 
                            na_level = "Unknown"),
      
      # Histology factor
      histology = fct_explicit_na(factor(histology, levels = c("Invasive Ductal Carcinoma",
                                                               "Invasive Lobular Carcinoma",
                                                               "Mixed Histology",
                                                               "Other")), 
                                  na_level = "Unknown"),
      
      # Subtype factors
      subtype_initial = fct_explicit_na(factor(subtype_initial, levels = c("HR+/HER2-",
                                                                           "HR+/HER2+",
                                                                           "HR-/HER2+",
                                                                           "HR-/HER2-")), 
                                       na_level = "Unknown"),
      
      # Grade factor
      grade = factor(grade, levels = c(1, 2, 3),
                    labels = c("Well differentiated", "Moderately differentiated", "Poorly differentiated")),
      grade = fct_explicit_na(grade, na_level = "Unknown")
    )
  
  # Table 1: All Breast Cancer Patients
  total_n <- nrow(bc_data)
  
  table1_data <- list()
  
  # Age category distribution
  table1_data$age <- bc_data %>%
    count(age_category) %>%
    mutate(
      percent = 100 * n / total_n,
      variable = "Age at Diagnosis (years)",
      category = as.character(age_category)
    )
  
  # Race distribution
  table1_data$race <- bc_data %>%
    count(race) %>%
    mutate(
      percent = 100 * n / total_n,
      variable = "Race/Ethnicity",
      category = as.character(race)
    )
  
  # Histology distribution
  table1_data$histology <- bc_data %>%
    count(histology) %>%
    mutate(
      percent = 100 * n / total_n,
      variable = "Histology",
      category = as.character(histology)
    )
  
  # Subtype at initial diagnosis distribution
  table1_data$subtype_init <- bc_data %>%
    count(subtype_initial) %>%
    mutate(
      percent = 100 * n / total_n,
      variable = "Subtype at Initial Diagnosis",
      category = as.character(subtype_initial)
    )
  
  # Grade distribution
  table1_data$grade <- bc_data %>%
    count(grade) %>%
    mutate(
      percent = 100 * n / total_n,
      variable = "Grade",
      category = as.character(grade)
    )
  
  # Combine all summary rows for Table 1
  table1_df <- bind_rows(table1_data$age, table1_data$race, table1_data$histology,
                         table1_data$subtype_init, table1_data$grade) %>%
    arrange(factor(variable, levels = c("Age at Diagnosis (years)", "Race/Ethnicity",
                                        "Histology", "Subtype at Initial Diagnosis", "Grade"))) %>%
    mutate(
      percent = round(percent, 1),
      stat = paste0(n, " (", percent, "%)")
    )
  
  write.csv(table1_df, file = file.path(outdir, 'bc_table1_demographics.csv'), row.names = FALSE)
  
  # Table 2: Metastatic Cohort (if time_to_met exists)
  if ('time_to_met' %in% names(bc_data)) {
    mets_data <- bc_data %>% filter(!is.na(time_to_met))
    mets_n <- nrow(mets_data)
    
    if (mets_n > 0) {
      # Similar structure for metastatic cohort
      table2_data <- list()
      
      table2_data$age <- mets_data %>%
        count(age_category) %>%
        mutate(percent = 100 * n / mets_n, variable = "Age at Diagnosis (years)", category = as.character(age_category))
      
      table2_data$race <- mets_data %>%
        count(race) %>%
        mutate(percent = 100 * n / mets_n, variable = "Race/Ethnicity", category = as.character(race))
      
      # Add time to metastasis summary
      mean_time <- mean(mets_data$time_to_met, na.rm = TRUE)
      sd_time <- sd(mets_data$time_to_met, na.rm = TRUE)
      
      table2_time_row <- tibble(
        variable = "Time from Initial Dx to Metastasis (months)",
        category = "Mean (SD)",
        n = NA,
        percent = NA,
        stat = paste0(round(mean_time, 1), " (", round(sd_time, 1), ")")
      )
      
      table2_df <- bind_rows(table2_data$age, table2_data$race, table2_time_row) %>%
        mutate(
          percent = round(percent, 1),
          stat = if_else(is.na(stat), paste0(n, " (", percent, "%)"), stat)
        )
      
      write.csv(table2_df, file = file.path(outdir, 'bc_table2_metastatic.csv'), row.names = FALSE)
    }
  }
}

# ===== Final Data Save =====
# Save cleaned data
cleaned_path <- file.path(outdir, 'merged_tcga_cleaned.rds')
saveRDS(df, cleaned_path)
message('💾 Saved cleaned data to: ', cleaned_path)

# Summary message
message("✅ Comprehensive genomic analysis complete!")
message("📂 All outputs saved to: ", outdir)
message("🔍 Key outputs:")
message("   - Cleaned dataset: merged_tcga_cleaned.rds")
message("   - Gene prevalence: gene_prevalence.csv")
message("   - Missing data analysis: missing_summary.csv")
message("   - Statistical models: univariate_lm.csv, multivariable_lm.csv")
message("   - Visualizations: Various PNG files")
if (!is.null(bc_data)) {
  message("   - Breast cancer tables: bc_table1_demographics.csv, bc_table2_metastatic.csv")
}

message("🎉 Analysis pipeline completed successfully!")
