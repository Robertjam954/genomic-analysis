# Time-to-Event Survival Analysis

# Install required packages (only once needed)
# install.packages(c("survival", "flexsurv", "survminer", "ggplot2", "readxl"))

# Load libraries
library(survival)
library(flexsurv)
library(survminer)
library(tidyverse)
library(readxl)

# Read in GENIE dataset (adjust path if needed)
df <- read_excel("genie_breast_data.xlsx")

# Ensure necessary variables are properly formatted
df <- df %>% 
  mutate(
    DFS_MONTHS = as.numeric(DFS_MONTHS),
    DFS_STATUS = as.numeric(DFS_STATUS),
    PFS_MONTHS = as.numeric(PFS_MONTHS),
    PFS_STATUS = as.numeric(PFS_STATUS),
    AGE_CAT = factor(AGE_CAT, levels = c("<50", "50-75", "Over 75"), ordered = TRUE),
    SUBTYPE = factor(SUBTYPE, ordered = TRUE),
    AJCC_STAGE = factor(AJCC_STAGE, ordered = TRUE)
  )

# Create a DFS survival object
surv_obj_dfs <- Surv(time = df$DFS_MONTHS, event = df$DFS_STATUS)

# Kaplan-Meier curve for DFS
km_dfs <- survfit(surv_obj_dfs ~ 1)
ggsurvplot(km_dfs, conf.int = TRUE, ggtheme = theme_minimal(), title = "DFS Kaplan-Meier Estimate")

# Compare parametric fits for DFS
dfs_fits <- list(
  weibull = flexsurvreg(surv_obj_dfs ~ 1, data = df, dist = "weibull"),
  lognormal = flexsurvreg(surv_obj_dfs ~ 1, data = df, dist = "lognormal"),
  loglogistic = flexsurvreg(surv_obj_dfs ~ 1, data = df, dist = "loglogistic")
)

# Compare AICs
sapply(dfs_fits, AIC)

# Create a PFS survival object
surv_obj_pfs <- Surv(df$PFS_MONTHS, df$PFS_STATUS)

# Parametric AFT model for PFS with covariates
aft_model <- flexsurvreg(
  surv_obj_pfs ~ AGE_CAT + SUBTYPE + AJCC_STAGE + TMB_NONSYNONYMOUS + ADI,
  data = df,
  dist = "weibull"
)

summary(aft_model)
