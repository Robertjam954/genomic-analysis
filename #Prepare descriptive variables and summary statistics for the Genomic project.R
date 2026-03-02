#!/usr/bin/env Rscript

#Prepare descriptive variables and summary statistics for the Genomic project.

#Creates a cleaned dataset with properly coded factors
#(AJCC stage as ordered factor) and writes summary CSVs into an output folder.

#Usage:
#Rscript scripts/prepare_descriptive_vars.
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
})

df_raw <- read.csv("~/Research/Projects/genomics_brain_mets_genie_bpc/tcga_gene_alteration_dfs_os_xgb/merged_tcga_df.csv", na = c('NA', ''))
colnames <- colnames(df_raw)
head(df_raw)

# Print working directory to help debug path issues
message("📂 Working directory: ", getwd())

# snake case col names
df <- clean_names(df_raw) |> glimpse()

# AJCC stage: Convert to ordered factor with labels as type character to match label.
  # Assume values 1-4 map to these levels (user-specified mapping)
  ajcc_levels <- c("Stage I", "Stage II", "Stage III", "Stage IV")
  df <- df %>%
    mutate(
      ajcc_pathologic_tumor_stage = as.character(ajcc_pathologic_tumor_stage),
      ajcc_pathologic_tumor_stage = case_when(
        ajcc_pathologic_tumor_stage == "1" ~ "Stage I",
        ajcc_pathologic_tumor_stage == "2" ~ "Stage II",
        ajcc_pathologic_tumor_stage == "3" ~ "Stage III",
        ajcc_pathologic_tumor_stage == "4" ~ "Stage IV",
        TRUE ~ NA_character_),
      ajcc_pathologic_tumor_stage = factor(ajcc_pathologic_tumor_stage, levels = ajcc_levels, ordered = TRUE)
      )

#Subtype: expect a column named subtype or subtype_label
  subtype_levels <- c("Luminal A", "Luminal B", "HER-2 enriched", "Basal-type (TNBC)")

  df <- df |>
    mutate(
      subtype = as.character(subtype),
      subtype = case_when(
           subtype == "1" ~  "Luminal A",
           subtype == "2" ~  "Luminal B",
           subtype == "3" ~  "HER-2 enriched",
           subtype == "4" ~  "Basal-type (TNBC)",
           TRUE ~ NA_character_),
      subtype = factor(subtype, levels = subtype_levels, ordered = TRUE)
    )

# Age category: if present as numeric/character, make ordered factor
age_levels <- c('<50', '50-75', 'Over 75')
age_cat <- factor(df_raw$age_cat, levels = age_levels, ordered = TRUE)

# Race / Ethnicity
if ('race' %in% names(df)) {
  race_levels <- c('White or Caucasian', 'Black or African American', 'Asian', 'Indian or Alaska Native')
  df <- df %>% mutate(race = factor(race, levels = race_levels))
}
if ('ethnicity' %in% names(df)) {
  eth_levels <- c('Non-Hispanic', 'Hispanic or Latino')
  df <- df %>% mutate(ethnicity = factor(ethnicity, levels = eth_levels))
}

# Sex: normalize and create binary indicator sex_bin (0=Female, 1=Male)
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
if ('dfs_status' %in% names(df)) {
  dfs_levels <- c('Recurrence or Progression of Disease', 'Disease Free')
  # assume 1 = recurrence/progression, 0 = disease free
  df <- df %>% mutate(dfs_status = case_when(
    dfs_status == 1 ~ dfs_levels[1],
    dfs_status == 0 ~ dfs_levels[2],
    TRUE ~ as.character(dfs_status)
  ), dfs_status = factor(dfs_status, levels = dfs_levels))
}

# MANTIS / MSI factor
mantis_levels <- c("MSI-H", "MSI-L")

df <- df |>
  mutate(
    mantis_bin = case_when(
      mantis_bin == "MSI-H" ~ ">0.4",
      mantis_bin == "MSI-L" ~ "<0.4",
      TRUE ~ NA_character_
    ),
    mantis_bin = factor(mantis_bin, levels = mantis_levels)
  )

# TBL score quartiles (use tbl_quantile variable)
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
  df <- df %>% mutate(buffa_hypoxia_quartile = cut(buffa_hypoxia_score, breaks = qh, include.lowest = TRUE, labels = c('Q1','Q2','Q3','Q4')))
}

# Aneuploidy quartiles
if ('aneuploidy_score' %in% names(df) || 'aneuploidy_score' %in% names(df)) {
  an_col <- ifelse('aneuploidy_score' %in% names(df), 'aneuploidy_score', 'aneuploidy_score')
  df <- df %>% mutate(!!sym(an_col) := as.numeric(.data[[an_col]]))
  qa <- quantile(df[[an_col]], probs = c(0, .25, .5, .75, 1), na.rm = TRUE)
  df <- df %>% mutate(aneuploidy_score_quartile = cut(.data[[an_col]], breaks = qa, include.lowest = TRUE, labels = c('Q1','Q2','Q3','Q4')))
}

# Variant type
if ('variant_type' %in% names(df)) {
  variant_type_levels <- c('SNP','INS','ONP','DEL')
  df <- df %>% mutate(variant_type = factor(variant_type, levels = variant_type_levels))
}

# Variant allele count quartiles (create if variant_allele_count exists)
var_acol <- NULL
for (candidate in c('variant_allele_count', 'variant_allele_freq', 'mutant_allele_count', 'mutant_allele_freq')) {
  if (candidate %in% names(df)) { var_acol <- candidate; break }
}
if (!is.null(var_acol)) {
  df <- df %>% mutate(!!sym(var_acol) := as.numeric(.data[[var_acol]]))
  vq <- quantile(df[[var_acol]], probs = c(0, .25, .5, .75, 1), na.rm = TRUE)
  df <- df %>% mutate(variant_allele_count_quartile = cut(.data[[var_acol]], breaks = vq, include.lowest = TRUE, labels = c('Q1','Q2','Q3','Q4')))
  message('Created variant_allele_count_quartile from ', var_acol)
}

# ----- Mutation counts and top-gene flags -----
# Normalize Hugo symbol column to character
df <- df %>%
  mutate(hugo_symbol = as.character(.data[["hugo_symbol"]]))

# Per-sample mutation count
df <- df %>%
  mutate(
    mutation_list = case_when(
      is.na(.data[["hugo_symbol"]]) ~ NA_character_,
      TRUE ~ str_replace_all(.data[["hugo_symbol"]], "\\s+", "")
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

  # Create presence flags for top10 genes and a combined top_5_gene_status
  if (length(top10_genes) > 0) {
    for (g in top10_genes) {
      colname <- paste0('G__', make.names(g))
      df[[colname]] <- as.integer(str_detect(df[["hugo_symbol"]] %||% '', fixed(g)))
    }
  }

  # Combined top-5 flag
  if (length(top5_genes) > 0) {
    df <- df %>% mutate(top_5_gene_status = as.integer(
      reduce(map(top5_genes, ~ str_detect(.data[["hugo_symbol"]] %||% '', fixed(.x))), `|`, .init = FALSE)
    ))
  } else {
    df <- df %>% mutate(top_5_gene_status = 0L)
  }

  # Save gene prevalence and top lists
  write.csv(gene_prev, file = 'gene_prevalence.csv', row.names = FALSE)
  writeLines(top10_genes, file = 'top10_genes.txt')
  writeLines(top5_genes, file = 'top5_genes.txt')

  # Group by DFS/DS status and summarize top_5 presence
    top5_by_status <- df %>%
      group_by(.data[[dfs_status]]) %>%
      summarise(n = n(),
                top5_count = sum(top_5_gene_status, na.rm = TRUE),
                top5_pct = 100 * top5_count / n) %>%
      arrange(desc(top5_pct))
    write.csv(top5_by_status, file = 'top5_by_status.csv', row.names = FALSE)
    message('Wrote top5_by_status.csv')

# Create an output of missingness per variable
missing_summary <- df %>%
  summarise(across(everything(), ~ sum(is.na(.)))) %>%
  pivot_longer(everything(), names_to = 'variable', values_to = 'n_missing') %>%
  mutate(n_total = nrow(df), prop_missing = n_missing / n_total)

write.csv(missing_summary, file = 'missing_summary.csv', row.names = FALSE)
message('Wrote missing_summary.csv')

# Skim and save a compact version
sk <- skim(df)
write.csv(as.data.frame(sk), file = file.path(outdir, 'skim_summary.csv'), row.names = FALSE)
message('Wrote skim summary (csv)')

# Save cleaned data (compressed)
cleaned_path <- file.path(outdir, 'merged_tcga_cleaned.rds')
saveRDS(df, cleaned_path)

packages
Library(gtsummary)
library(knitr)
library(kableExtra)

## visualize missingness, add proportion NA to visual
vis_dat(df1_clean)
vis_misss(df1_clean)

#table 1\
\resetSession()
library(tidyverse)
library(tfrmt)

# Get data
data("cadsl", package = "random.cdisc.data")

# Number of unique subjects per ARM
big_n <- cadsl |>
  dplyr::group_by(ARM) |>
  dplyr::summarize(
    N = dplyr::n_distinct(USUBJID)
  )

# Join big_n with adsl
adsl_with_n <- cadsl |>
  dplyr::left_join(big_n, by = "ARM")

# Explore column: AGE
age_stats <-
  adsl_with_n |>
  group_by(ARM) |>
  reframe(
    n = n_distinct(USUBJID),
    Mean = mean(AGE),
    SD = sd(AGE),
    Median = median(AGE),
    Min = min(AGE),
    Max = max(AGE)
  ) |>
  pivot_longer(
    c("n", "Mean", "SD", "Median", "Min", "Max")
  ) |>
  mutate(
    group = "Age (years)",
    label = case_when(name == "Mean" ~ "Mean (SD)",
                      name == "SD" ~ "Mean (SD)",
                      name == "Min" ~ "Min - Max",
                      name == "Max" ~ "Min - Max",
                      TRUE ~ name)
  )

sex_n <-
  adsl_with_n |>
  group_by(ARM, SEX) |>
  reframe(
    n = n(),
    pct = (n/N)*100
  ) |>
  distinct() |>
  pivot_longer(
    c("n", "pct")
  ) |>
  rename(
    label = SEX
  ) |>
  mutate(
    group = "Sex"
  )
