# Exploratory Correlation and Association Analysis Script

# Load required libraries
library(naniar)      # For missing data visualization
library(GGally)      # For correlation matrix
library(ggplot2)     # For plots
library(broom)       # For tidy model outputs
library(dplyr)

# Output directory
if (!dir.exists("../output_descriptive_vars")) dir.create("../output_descriptive_vars", recursive = TRUE)

# Clean numeric subset for correlation analysis
numeric_vars <- df %>% select(where(is.numeric))

# Correlation matrix plot
png("../output_descriptive_vars/correlation_matrix.png", width = 800, height = 800)
ggcorr(numeric_vars, label = TRUE)
dev.off()

# Visualize missingness
png("../output_descriptive_vars/missingness_plot.png", width = 800, height = 600)
gg_miss_var(df)
dev.off()

# Univariate linear regression (e.g., DS_Status ~ top_5_gene_status)
fit_uni <- lm(dfs_status ~ top_5_gene_status, data = df)
summary(fit_uni)
tidy(fit_uni) %>% write.csv("../output_descriptive_vars/univariate_lm.csv", row.names = FALSE)

# Multivariate regression with possible confounders
fit_multi <- lm(dfs_status ~ top_5_gene_status + sex_bin + age_cat + ajcc_pathologic_tumor_stage + buffa_hypoxia_quartile + aneuploidy_score_quartile, data = df)
summary(fit_multi)
tidy(fit_multi) %>% write.csv("../output_descriptive_vars/multivariable_lm.csv", row.names = FALSE)

# Scatterplot: DS_Status vs Top 5 gene status, faceted by confounders
ggplot(df, aes(x = top_5_gene_status, y = as.numeric(dfs_status))) +
  geom_jitter(width = 0.2, height = 0.1) +
  facet_wrap(~ sex + age_cat) +
  labs(title = "DFS Status vs Top 5 Gene Status by Sex and Age", y = "DFS Status (0/1)", x = "Top 5 Gene Mutated") +
  ggsave("../output_descriptive_vars/dfs_vs_top5_facet.png", width = 10, height = 6)

# Residual diagnostics
png("../output_descriptive_vars/residuals_diagnostics.png", width = 800, height = 800)
par(mfrow = c(2, 2))
plot(fit_multi)
dev.off()

message("✅ Exploratory correlation and regression analysis complete. Output written to output_descriptive_vars folder.")
