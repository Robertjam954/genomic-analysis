# ===== Setup =====
library(tidyverse)
library(forcats)
library(scales)
library(viridisLite)

# ---- Path to your CSV (Windows-friendly forward slashes) ----
data_path <- "C:/Users/jamesr4/OneDrive - Memorial Sloan Kettering Cancer Center/Documents/Research/Projects/genomics_brain_mets_genie_bpc/tcga_gene_alteration_dfs_os_xgb/merged_tcga_df.csv"

df <- readr::read_csv(data_path, show_col_types = FALSE)

# ===== Helpers =====
theme_msk <- function(base = 10) {
  theme_minimal(base_size = base) +
    theme(
      panel.grid.minor = element_blank(),
      panel.spacing    = unit(0.6, "lines"),
      strip.text       = element_text(face = "bold"),
      legend.position  = "right"
    )
}

# Safe factor with explicit unknown level "."
as_factor_with_dot <- function(x, lvl = NULL, ordered = FALSE) {
  x <- as.character(x)
  x[is.na(x) | trimws(x) == ""] <- "."
  if (is.null(lvl)) return(factor(x, ordered = ordered))
  factor(x, levels = c(lvl, if (!"." %in% lvl) "."), ordered = ordered)
}

# ===== Column mappings (EDIT HERE if your names differ) =====
# Variant allele count (numeric)
col_vac   <- "variant_allele_count"    # e.g., numeric; change if needed
# DFS status: 0/1 (or character). Will create a labeled factor for plotting.
col_dfs   <- "dfs_status"              # values like 0/1 or "DFS0/DFS1"
# Demographics for grouped bars:
col_race  <- "race"                    # factor-like
col_eth   <- "ethnicity"               # 0/1 or "Hispanic/Not Hispanic"
col_subty <- "subtype"                 # e.g., ER+/HER2-, etc.
# For gene prevalence segment plot:
col_gene  <- "gene"                    # gene symbol or name
col_prev  <- "prevalence"              # numeric prevalence (0-100 or 0-1)

# For probability faceted histograms by mutation covariate:
# A long-frame with columns: 'text' (covariate label) & 'value' (probability %)
# If you don’t have it long already, we’ll build a demo from wide columns (EDIT HERE).
prob_columns <- c()  # e.g., c("prob_tp53", "prob_pik3ca", ...)

# ===== Light coercions / standardizations =====
# DFS factor
df <- df %>%
  mutate(
    dfs_raw = !!sym(col_dfs),
    dfs_bin = case_when(
      is.numeric(dfs_raw) ~ as.integer(dfs_raw),
      grepl("1", as.character(dfs_raw)) ~ 1L,
      grepl("yes|recurrence|progression", tolower(as.character(dfs_raw))) ~ 1L,
      TRUE ~ 0L
    ),
    dfs_f   = factor(dfs_bin, levels = c(0,1),
                     labels = c("DFS = 0", "DFS = 1"))
  )

# Ethnicity recode to explicit labels (per your dictionary)
if (col_eth %in% names(df)) {
  df <- df %>%
    mutate(
      !!sym(col_eth) := case_when(
        is.numeric(!!sym(col_eth)) ~ !!sym(col_eth),
        grepl("hisp", tolower(as.character(!!sym(col_eth)))) ~ 1,
        grepl("non",  tolower(as.character(!!sym(col_eth)))) ~ 0,
        TRUE ~ NA
      ),
      ethnicity_f = fct_explicit_na(
        factor(ifelse(!!sym(col_eth) == 1, "Hispanic",
                      ifelse(!!sym(col_eth) == 0, "Not Hispanic", ".")),
               levels = c("Not Hispanic","Hispanic","."))
      )
    )
}

# Race & Subtype factors w/ "." for unknown
if (col_race %in% names(df))   df <- df %>% mutate(race_f   = as_factor_with_dot(!!sym(col_race)))
if (col_subty %in% names(df))  df <- df %>% mutate(subtype_f= as_factor_with_dot(!!sym(col_subty)))
if (col_vac %in% names(df))    df <- df %>% mutate(vac = suppressWarnings(as.numeric(!!sym(col_vac))))

# ====== 1) Variant Allele Count: Histogram + Boxplots layout ======
# Layout analogous to your base R example, but in ggplot (cleaner + saveable).
# If you must use base R layout, see the base version further down.

p_hist_vac <- ggplot(df, aes(x = vac)) +
  geom_histogram(bins = 30, alpha = 0.8) +
  labs(x = "Variant Allele Count", y = "Frequency", title = "Distribution of Variant Allele Count") +
  theme_msk()

p_box_vac_all <- ggplot(df, aes(y = vac)) +
  geom_boxplot(outlier.alpha = 0.4) +
  coord_flip() +
  labs(x = "", y = "Variant Allele Count", title = "Overall") +
  theme_msk()

p_box_vac_by_dfs <- ggplot(df, aes(x = dfs_f, y = vac, fill = dfs_f)) +
  geom_boxplot(outlier.alpha = 0.4, show.legend = FALSE) +
  coord_flip() +
  labs(x = "", y = "Variant Allele Count", title = "By DFS Status") +
  theme_msk()

# Arrange (patchwork is convenient; no extra theming imposed)
suppressWarnings({
  if (requireNamespace("patchwork", quietly = TRUE)) {
    library(patchwork)
    p_vac_layout <- p_hist_vac / (p_box_vac_all | p_box_vac_by_dfs)
    ggsave("plot_vac_hist_box.png", p_vac_layout, width = 10, height = 8, dpi = 300)
  } else {
    # Fallback: save each panel separately
    ggsave("plot_vac_hist.png", p_hist_vac, width = 7, height = 4, dpi = 300)
    ggsave("plot_vac_box_overall.png", p_box_vac_all, width = 5, height = 3.5, dpi = 300)
    ggsave("plot_vac_box_by_dfs.png", p_box_vac_by_dfs, width = 5, height = 3.5, dpi = 300)
  }
})

# ----- OPTIONAL: strict base R layout (close to your snippet) -----
# if (interactive()) {
#   a <- df$vac[!is.na(df$vac)]
#   b <- df$vac[!is.na(df$vac) & df$dfs_bin == 1]
#   nf <- layout(matrix(c(1,1,2,3), nrow=2, byrow=TRUE))
#   hist(a, breaks=30, border=FALSE, col=rgb(0.1,0.8,0.3,0.5),
#        xlab="Variant Allele Count", main="")
#   boxplot(a, xlab="Overall", col=rgb(0.8,0.8,0.3,0.5), las=2)
#   boxplot(b, xlab="DFS = 1", col=rgb(0.4,0.2,0.3,0.5), las=2)
# }

# ====== 2) Probability of DFS = 1 across demographics (grouped bars) ======
# Here “probability” is proportion of DFS=1 per group.
plot_grouped_prop <- function(data, group_var, out_file) {
  data %>%
    group_by({{ group_var }}) %>%
    summarise(
      n = n(),
      prop_dfs1 = mean(dfs_bin == 1, na.rm = TRUE)
    ) %>%
    mutate({{ group_var }} := fct_reorder(as.factor({{ group_var }}), prop_dfs1, .desc = TRUE)) %>%
    ggplot(aes(x = {{ group_var }}, y = prop_dfs1)) +
    geom_col() +
    geom_text(aes(label = percent(prop_dfs1, accuracy = 0.1)), vjust = -0.3, size = 3) +
    labs(x = NULL, y = "P(DFS = 1)", title = "Probability of DFS = 1 by Group") +
    scale_y_continuous(labels = percent_format()) +
    theme_msk() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) -> p

  ggsave(out_file, p, width = 7.5, height = 4.5, dpi = 300)
}

if ("race_f" %in% names(df))      plot_grouped_prop(df, race_f, "plot_prob_dfs1_by_race.png")
if ("ethnicity_f" %in% names(df)) plot_grouped_prop(df, ethnicity_f, "plot_prob_dfs1_by_ethnicity.png")
if ("subtype_f" %in% names(df))   plot_grouped_prop(df, subtype_f, "plot_prob_dfs1_by_subtype.png")

# ====== 3) Gene mutation prevalence — DFS=1 on left vs DFS=0 on right ======
# Expect a long frame with columns: gene (col_gene), prevalence (col_prev), dfs_f
# If you have raw mutation counts, compute prevalence first (EDIT if needed).

if (all(c(col_gene, col_prev) %in% names(df))) {
  gene_prev_df <- df %>%
    select(all_of(c(col_gene, col_prev, "dfs_f"))) %>%
    filter(!is.na(!!sym(col_gene)), !is.na(!!sym(col_prev)))

  # Order genes by prevalence within each DFS group, then plot with segments + points
  for (g in levels(gene_prev_df$dfs_f)) {
    tmp <- gene_prev_df %>% filter(dfs_f == g) %>%
      arrange(!!sym(col_prev)) %>%
      mutate(gene_ord = factor(!!sym(col_gene), levels = !!sym(col_gene)))
    p <- ggplot(tmp, aes(x = gene_ord, y = !!sym(col_prev))) +
      geom_segment(aes(xend = gene_ord, yend = 0)) +
      geom_point(size = 2.8, color = "orange") +
      coord_flip() +
      theme_bw() +
      labs(x = NULL, y = "Prevalence", title = paste("Gene Mutation Prevalence —", g)) +
      theme_msk()
    ggsave(paste0("plot_gene_prev_", g, ".png"), p, width = 6.5, height = 6.5, dpi = 300)
  }
}

# ====== 4) Stacked bar: DFS 1 vs 0 across a categorical x (e.g., subtype) ======
if ("subtype_f" %in% names(df)) {
  p_stack <- ggplot(df, aes(x = subtype_f, fill = dfs_f)) +
    geom_bar(position = "fill") +
    scale_y_continuous(labels = percent_format()) +
    labs(x = "Subtype", y = "Share", fill = "DFS", title = "DFS = 1 vs 0 by Subtype") +
    theme_msk() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave("plot_dfs_stack_by_subtype.png", p_stack, width = 7.5, height = 4.8, dpi = 300)
}

# ====== 5) Faceted histograms of P(DFS=1) across mutation covariates ======
# If you already have a long frame with columns: text (facet label) & value (probability %),
# set prob_columns <- character(0) and instead supply df_long below.
if (length(prob_columns) > 0 && all(prob_columns %in% names(df))) {
  df_long <- df %>%
    select(all_of(prob_columns)) %>%
    pivot_longer(everything(), names_to = "text", values_to = "value") %>%
    mutate(text  = fct_reorder(text, value, .fun = median, .desc = FALSE),
           value = as.numeric(value))

  p_fac <- df_long %>%
    ggplot(aes(x = value, color = text, fill = text)) +
    geom_histogram(alpha = 0.6, binwidth = 5) +
    scale_fill_viridis_d() +
    scale_color_viridis_d() +
    theme_msk() +
    theme(
      legend.position = "none",
      strip.text.x    = element_text(size = 8)
    ) +
    xlab("") +
    ylab("Assigned Probability (%)") +
    facet_wrap(~ text, scales = "free_y")
  ggsave("plot_prob_dfs1_faceted.png", p_fac, width = 10, height = 8, dpi = 300)
}

message("All plots saved in the working directory:")
message(getwd())
