# Executive Summary & Code Map

## Project Overview

**genomic-analysis-scripts** is a clinical-genomic research pipeline for brain metastasis survival analysis. It ingests AACR GENIE BPC, TCGA, and MSK-IMPACT data, performs de-identification, exploratory analysis, survival modeling (Kaplan-Meier, Cox PH, Fine-Gray, XGBoost AFT), and generates manuscript-ready figures and tables.

**Domain:** Breast cancer brain metastasis — gene-level alteration indicators, microsatellite instability, hypoxia scores, aneuploidy, and time-to-event endpoints (DFS, OS).

**Languages:** Python 3.9+ and R 4.x, with Jupyter notebooks for prototyping.

---

## Directory Structure

```
genomic-analysis-scripts/
├── notebooks/                  Analysis notebooks (Jupyter)
├── data/
│   ├── raw/                    cBioPortal API outputs (clinical, mutation, profile CSVs)
│   ├── processed/              Non-PHI analysis outputs
│   ├── features/               PCA components, cluster assignments
│   └── splits/                 Train/test index definitions
├── datasets_analysis_dictionary/  Primary merged dataset (merged_genie.xlsx/.csv)
├── docs/                       Documentation
├── conferences/                Conference abstract drafts
├── reports/
│   ├── figures/                Exported figures (PNG/PDF)
│   └── tables/                 Exported tables (CSV/Excel)
├── references/                 Papers, cheat sheets, external resource index
├── src/                        Source code (see Code Map below)
│   ├── collection/             Data acquisition
│   ├── cleaning/               De-identification, variable extraction
│   ├── eda/                    Exploratory data analysis
│   ├── modeling/               Survival models, ML
│   └── reports/                Figure and table generation scripts
│       ├── figures/
│       └── tables/
├── eval/                       Evaluation schemas, missingness classification
├── models/                     Model configs (XGBoost, Cox PH params)
├── experiments/                Run tracking (run_id, commit, results)
└── tools/                      Utility scripts
```

---

## Code Map — Every File, What It Does, and What It Produces

### Root Files

| File | Purpose | Outputs |
|---|---|---|
| `README.md` | Project overview, directory tree, component listings, quick-start guide | — |
| `LICENCE` | MIT license | — |
| `.env.example` | Template for `PROJECT_ROOT`, `DATA_PRIVATE_DIR`, API keys | Copy to `.env` |
| `.gitignore` | Excludes PHI data, large files, env, caches; whitelists `data/processed/**/*.csv` | — |
| `.gitattributes` | Git line-ending and diff config | — |
| `pyproject.toml` | Project metadata, Python dependencies, build config, `[dev]` extras | — |
| `requirements.txt` | Pip-compatible dependency list: pandas, numpy, matplotlib, seaborn, lifelines, xgboost, openpyxl, statsmodels, scikit-learn | — |

### `src/collection/` — Data Acquisition

| File | Language | What It Does | Inputs | Outputs |
|---|---|---|---|---|
| `fetch_molecular_data.py` | Python | cBioPortal REST API client. Fetches molecular profiles, mutation data, clinical attributes. Includes retry logic and CLI. | cBioPortal API (no key needed) | `data/raw/<study>_clinical.csv`, `data/raw/<study>_profiles.csv`, `data/raw/<profile>_mutations.csv`, `data/raw/molecular_data.csv` |
| `build_merged_dataset.py` | Python | Pivots long-format clinical data to wide, builds `G__*` gene indicators from mutation data, merges into single analysis-ready dataset. | `data/raw/brca_tcga_clinical.csv`, `data/raw/brca_tcga_mutations_mutations.csv` | `datasets_analysis_dictionary/merged_genie.xlsx`, `datasets_analysis_dictionary/merged_genie.csv` |

### `src/cleaning/` — De-identification & Variable Preparation

| File | Language | What It Does | Inputs | Outputs |
|---|---|---|---|---|
| `deidentify.py` | Python | Full de-identification engine: HMAC date shifting, MRN hashing, ZIP→ZIP3, age >89→90+, PHI regex scrubbing in free text, provenance JSON. | Raw CSV with PHI, salt files | De-identified CSV, `provenance.json`, optional `linkage_map.csv` |
| `deidentify_redcap.py` | Python | Same engine as above, specifically tuned for REDCap-style exports. | REDCap export CSV, `--salt-file`, `--hash-salt-file` | De-identified CSV, provenance JSON, optional linkage map |
| `extract variables from dataset.py` | Python | Extracts and validates variables from GENIE BPC metastatic dataset against a data dictionary. 50+ clinical/genomic variables including treatment, metastatic site, and therapy fields. | `genie_bpc_metastatic_dataset_latest.csv/.xlsx`, `genie_data_dic.xlsx` | `output/extracted_subset.csv`, `output/qa_report.csv`, `output/missing_variables_suggestions.csv` |
| `extract variables from dataset.R` | R | Reusable extractor/validator: coerces types (date, binary, numeric, factor), validates against allowed value sets, suggests data sources for missing variables. | Any CSV/Excel + data dictionary Excel | `extracted_subset.csv`, `qa_report.csv`, `missing_variables_suggestions.csv` |
| `extract_impact_variables.py` | R (mislabeled .py) | Calls the R extractor for MSK-IMPACT data. Includes helper functions: `snake()`, `coerce_to_class()`, `parse_allowed()`, `validate_column()`, `suggest_source()`. | `filtered_goel_impact_deid.csv`, `goel_data_dictionary_from_snakecase.xlsx` | Same as above |
| `prepare_descriptive_vars.r` | R | Core cleaning: snake_case normalization, AJCC stage → ordered factor, subtype mapping (Luminal A/B, HER-2, TNBC), sex normalization, DFS/MANTIS/TBL/hypoxia/aneuploidy variable coding, gene prevalence flags (`G__*`), top-5 gene status. | `merged_tcga_df.csv` | `merged_tcga_cleaned.rds`, `gene_prevalence.csv`, `top10_genes.txt`, `top5_genes.txt`, `top5_by_status.csv`, `missing_summary.csv`, `skim_summary.csv` |
| `date_handling_template.R` | R | Reference template for date/time handling with lubridate. Covers string→Date parsing, component extraction, arithmetic, intervals, durations, time zones. | — (template) | — (educational) |

### `src/eda/` — Exploratory Data Analysis

| File | Language | What It Does | Inputs | Outputs |
|---|---|---|---|---|
| `variable_setup_gene_flagging_plots.py` | Python | Missingness analysis, patient characteristics table, gene-level descriptive statistics, prevalence bar plots, KM curves by gene, correlation heatmaps. | `merged_genie.xlsx` (or candidates) | Missingness CSV, patient characteristics table, gene prevalence plots, KM PNGs, correlation heatmap |
| `comprehensive_descriptive_analysis.R` | R | Comprehensive analysis: factor coding (AJCC, subtype, race, sex, DFS, MANTIS), gene flags, Table 1/Table 2 generation via `gtsummary`, missingness visualization via `naniar`. | `merged_tcga_df.csv` | `output_descriptive_vars/` — cleaned RDS, skim summary, missing summary, gene prevalence CSVs |
| `R- generate descriptive and associational plots.R` | R | Produces outcome-related plots and tables from cleaned data. Table 1, Table 2, mosaic plots, boxplots by DFS status, subtype distributions. | `merged_genie_cleaned.rds` | `output/descriptive/` — Table 1/2 CSVs, association plots (PNGs) |
| `R - generate plots for missing_associational_correlation_matrix_regression_scatterfacetedbyconfounders.R` | R | Correlation matrix, missingness visualization, univariate/multivariate regressions, scatter plots faceted by confounders, interaction terms. | `df` object (expects upstream script to have run) | `output_descriptive_vars/` — `correlation_matrix.png`, `missingness_plot.png`, `univariate_lm.csv`, multivariate regression CSVs, interaction plots |

### `src/modeling/` — Survival Analysis & Machine Learning

| File | Language | What It Does | Inputs | Outputs |
|---|---|---|---|---|
| `survival_analysis_km_cox_finegray_xgb_aft.py` | Python | **Main survival pipeline.** Kaplan-Meier (overall + stratified by SUBTYPE, MANTIS_BIN), univariate + multivariate Cox PH, Fine-Gray competing risks, XGBoost AFT. Auto-detects time/event columns. | `merged_genie.xlsx`, optional `per_gene_cs_cox_sksurv_bootstrap.csv` | `output/survival/` — KM PNGs, Cox PH summary CSVs, Fine-Gray CSVs, `xgb/xgb_aft_model.json`, C-index report |
| `Py- survival and time to event models.py` | Python | Risk modeling for OS and PFS. KM curves by subtype, multivariable Cox PH, gene-level Cox tests with FDR correction. | Merged dataset | KM PNGs, Cox PH CSVs, per-gene HR/p-value/q-value CSV |
| `Py - xgb_aft_prediction_1020.py` | Python | XGBoost AFT model training and prediction. Drops TBL/quartile features, one-hot encodes categoricals, fits AFT with `survival:aft` objective, reports concordance index. | Pre-cleaned CSV with `survival_time` and `event` columns | Trained XGBoost model, C-index metric, predictions CSV |
| `xgb_aft_shap_summary.py` | Python | SHAP-like feature contribution analysis for trained XGBoost AFT model. Loads saved model JSON, computes `pred_contribs`, ranks features by mean |contrib|. | `merged_genie.xlsx`, `xgb_aft_model.json`, optional per-gene CSV | `output/survival/xgb/` — `shap_contribs.csv`, `shap_feature_ranking.csv`, `shap_topN_bar.png` |
| `shap analysis and plot setup.R` | R | SHAP analysis using `SHAPforxgboost` R package. Installs package if missing, fits XGBoost model, generates SHAP summary/dependence/interaction plots. | Dataset + XGBoost model | SHAP summary plot PNG, dependence plots, interaction plots |
| `time_event_analysis_aft.R` | R | Parametric AFT survival analysis. Compares Weibull, log-normal, and log-logistic distributions. KM curves for DFS, model comparison by AIC/BIC. | `genie_breast_data.xlsx` | KM plot, AFT model summaries, AIC/BIC comparison |

### `src/reports/figures/` — Figure Generation

| File | Language | What It Does | Inputs | Outputs |
|---|---|---|---|---|
| `export and assemble_multipanel_figure_for_genetic_variables.py` | Python | Assembles a multi-panel results figure from existing PNGs (covariate grid, KM OS/PFS, gene prevalence, molecular markers, Fine-Gray volcano). | PNGs in `output/` subdirectories | `output/figures/results_figure.png`, `results_figure.pdf` |
| `export and generate genetic mut prevalence plots.py` | Python | Computes top mutated genes, plots per-subtype mutation prevalence bar charts. | Merged dataset with `Hugo_Symbol_list` | `output/gene_mutation_plots/` — prevalence bar PNGs |
| `export gene and outcome plots.R` | R | Generates gene-outcome visualizations: faceted bar plots by gene and DFS status, prevalence by subtype, MSI/hypoxia distributions. Uses MSK-themed `ggplot2`. | `merged_tcga_df.csv` | Gene-outcome PNGs in `output/` |
| `R - export and generate shap analysis and plot generation.R` | R | SHAP visualization using `SHAPforxgboost`. Fits XGBoost, generates SHAP summary/dependence/force/interaction plots. | Dataset | SHAP PNGs in output directory |

### `src/reports/tables/` — Table Generation

| File | Language | What It Does | Inputs | Outputs |
|---|---|---|---|---|
| `table_1_and_2.R` | R | Generates manuscript Table 1 (all patients: age, race, histology, subtype, grade) and Table 2 (metastatic patients: subtype, met sites, therapies) using `tfrmt` formatting. | `complete_genie_bpc_dataset_deid.xlsx` | Formatted Table 1 and Table 2 |
| `BBB regulation table.py` | Python | Reads and reformats a blood-brain barrier regulation table from Excel. Identifies header rows, parses layer structure. **Note:** references `caas_jupyter_tools` (notebook environment). | `table_1_BBB_regulation.xlsx` | Formatted BBB regulation table |

### `notebooks/` — Jupyter Notebooks

| File | Description | Inputs | Outputs |
|---|---|---|---|
| `xgb_aft_1020.ipynb` | Original XGBoost AFT exploration notebook | Merged dataset | In-notebook analysis |
| `02_dimensionality_reduction_pca.ipynb` | PCA on `G__*` gene indicators — scree plot, biplot, variance analysis, PC export. Adapted from SalvatoreRa/tutorial. | `merged_genie.xlsx` | `reports/figures/pca_scree_plot.png`, `pca_scatter_pc1_pc2.png`, `data/features/pca_gene_features.csv` |
| `03_clustering_gene_features.ipynb` | K-Means, hierarchical, DBSCAN, GMM clustering. Elbow/silhouette analysis, dendrogram, method comparison. Adapted from SalvatoreRa/tutorial. | `data/features/pca_gene_features.csv` or raw gene columns | `reports/figures/clustering_*.png`, `data/features/cluster_assignments.csv` |
| `04_regression_gene_expression.ipynb` | Linear + logistic regression baselines. Coefficient plots, ROC curves, confusion matrices. Adapted from SalvatoreRa/tutorial. | `merged_genie.xlsx` | `reports/figures/linreg_top_coefficients.png`, `logreg_roc_confusion.png`, `data/processed/baseline_regression_metrics.csv`, `data/splits/train_indices.csv`, `test_indices.csv` |

### `eval/` — Evaluation

| File | Purpose | Inputs | Outputs |
|---|---|---|---|
| `missing_data_classification.py` | Classifies column missingness as MCAR/MAR/MNAR. Implements Little's MCAR χ² test, point-biserial correlations (MAR), KS distribution tests (MNAR). CLI and library. | Any CSV/Excel dataset | `data/processed/missingness_classification.csv`, `mar_correlations.csv`, `mnar_distribution_tests.csv` |
| `metric_definitions.md` | Documents all metrics used across the pipeline: C-index, log-rank, HR, AIC/BIC, SHAP, R², AUC, silhouette, Little's test. Includes experiment tracking YAML schema. | — (documentation) | — |

### `docs/` — Documentation

| File | Purpose |
|---|---|
| `executive_summary.md` | This file — full technical summary and code map |
| `dataset_metadata.md` | YAML front-matter dataset specification: schema (25+ columns), preprocessing steps, provenance, access controls, cBioPortal API usage |
| `README_survival_AFT_pipeline.md` | Detailed survival analysis pipeline documentation: inputs, column expectations, model parameters, outputs |
| `RESTRUCTURING_NOTES.md` | Records duplicate files and known issues found during repository restructuring |

### `references/`

| File | Purpose |
|---|---|
| `REFERENCES_INDEX.md` | Curated index of local files + external resources from SalvatoreRa/tutorial: genomic series, tabular learning, missing data, AI in medicine, graph ML |
| `CheatSheet-Python-5_-Functions-and-Tricks.pdf` | Python functions quick reference |
| `CheatSheet-Python-5_-Functions-and-Tricks-2.pdf` | Python functions quick reference (v2) |
| `Finxter_WorldsMostDensePythonCheatSheet.pdf` | Dense Python cheat sheet |
| `data wrangling python.pdf` | Python data wrangling reference |
| `Tidyverse_Cheat_Sheet-2.pdf` | R Tidyverse cheat sheet |
| `Mathematical modeling of the metastatic process.pdf` | Metastatic process modeling paper |
| `Table_AFT_Model_Assumptions.pdf` | AFT model assumptions reference |
| `Excel Quick Reference Cheat Sheet.pages` | Excel quick reference |

### `tools/`

| File | Purpose |
|---|---|
| `watch_repo.sh` | Bash script for macOS notifications when remote repo has new commits. Polls via `git fetch` at configurable interval. Requires `terminal-notifier`. |

---

## Pipeline Execution Order

```
1. DATA ACQUISITION
   src/collection/fetch_molecular_data.py   → data/raw/*.csv
   src/collection/build_merged_dataset.py   → datasets_analysis_dictionary/merged_genie.xlsx

2. DE-IDENTIFICATION (if working with raw PHI)
   src/cleaning/deidentify.py               → de-identified CSV + provenance JSON
   src/cleaning/deidentify_redcap.py        → de-identified REDCap CSV

3. VARIABLE EXTRACTION & CLEANING
   src/cleaning/extract variables from dataset.py/.R → extracted subset + QA report
   src/cleaning/extract_impact_variables.py          → IMPACT variables + QA
   src/cleaning/prepare_descriptive_vars.r           → cleaned RDS + gene prevalence + missingness

4. EXPLORATORY DATA ANALYSIS
   src/eda/variable_setup_gene_flagging_plots.py     → descriptive tables, gene plots, KM curves
   src/eda/comprehensive_descriptive_analysis.R      → Table 1/2 drafts, factor summaries
   src/eda/R- generate descriptive and associational plots.R → association plots
   src/eda/R - generate plots for missing_...R       → correlation matrix, regression, scatter plots
   eval/missing_data_classification.py               → MCAR/MAR/MNAR classification CSVs

5. DIMENSIONALITY REDUCTION & CLUSTERING (Notebooks)
   notebooks/02_dimensionality_reduction_pca.ipynb   → PCA features, scree/scatter plots
   notebooks/03_clustering_gene_features.ipynb       → cluster assignments, dendrograms
   notebooks/04_regression_gene_expression.ipynb     → baseline metrics, ROC, coefficients

6. SURVIVAL MODELING
   src/modeling/survival_analysis_km_cox_finegray_xgb_aft.py → KM, Cox PH, Fine-Gray, XGBoost AFT
   src/modeling/Py- survival and time to event models.py     → OS/PFS risk models, gene-level Cox
   src/modeling/Py - xgb_aft_prediction_1020.py              → XGBoost AFT predictions
   src/modeling/time_event_analysis_aft.R                    → Parametric AFT (Weibull, log-normal)

7. EXPLAINABILITY
   src/modeling/xgb_aft_shap_summary.py              → SHAP feature ranking + bar plot
   src/modeling/shap analysis and plot setup.R        → SHAP summary/dependence/interaction plots

8. REPORT GENERATION
   src/reports/figures/export and assemble_multipanel_figure_for_genetic_variables.py → multi-panel figure
   src/reports/figures/export and generate genetic mut prevalence plots.py            → gene prevalence plots
   src/reports/figures/export gene and outcome plots.R                                → gene-outcome plots
   src/reports/figures/R - export and generate shap analysis and plot generation.R    → SHAP plots
   src/reports/tables/table_1_and_2.R                                                → manuscript tables
   src/reports/tables/BBB regulation table.py                                        → BBB table
```

---

## Key Outputs Summary

| Category | Files Produced | Location |
|---|---|---|
| **Raw API data** | `*_clinical.csv`, `*_mutations.csv`, `*_profiles.csv` | `data/raw/` |
| **Merged dataset** | `merged_genie.xlsx`, `merged_genie.csv` | `datasets_analysis_dictionary/` |
| **Cleaned data** | `merged_tcga_cleaned.rds`, extracted subsets | `output/`, `data/processed/` |
| **QA reports** | `qa_report.csv`, `missing_variables_suggestions.csv` | `output/` |
| **Missingness** | `missingness_classification.csv`, `mar_correlations.csv`, `mnar_distribution_tests.csv` | `data/processed/` |
| **Gene prevalence** | `gene_prevalence.csv`, `top10_genes.txt`, `top5_genes.txt` | Working directory |
| **PCA features** | `pca_gene_features.csv` | `data/features/` |
| **Cluster assignments** | `cluster_assignments.csv` | `data/features/` |
| **Regression metrics** | `baseline_regression_metrics.csv` | `data/processed/` |
| **Train/test splits** | `train_indices.csv`, `test_indices.csv` | `data/splits/` |
| **KM curves** | `kaplan_os.png`, `kaplan_pfs.png`, stratified KM PNGs | `output/survival/`, `reports/figures/` |
| **Cox PH summaries** | Univariate/multivariate CSV summaries, per-gene HR/p/q | `output/survival/` |
| **XGBoost AFT model** | `xgb_aft_model.json`, C-index report | `output/survival/xgb/` |
| **SHAP analysis** | `shap_contribs.csv`, `shap_feature_ranking.csv`, `shap_topN_bar.png` | `output/survival/xgb/` |
| **Multi-panel figure** | `results_figure.png`, `results_figure.pdf` | `output/figures/` |
| **Manuscript tables** | Table 1, Table 2 | `output/` |
| **Descriptive plots** | Correlation matrix, missingness, scatter, bar charts | `output_descriptive_vars/`, `reports/figures/` |

---

## Known Issues

1. **`src/cleaning/extract_impact_variables.py`** — File extension is `.py` but content is R code.

2. **`src/reports/tables/BBB regulation table.py`** — References `caas_jupyter_tools` (notebook-only); `py` on line 7 is a syntax error.

3. **Hardcoded Windows paths** in several scripts (e.g., `C:\Users\jamesr4\...`). The pipeline uses candidate-path fallback logic to handle this on other machines.

## Cleanup Performed

The following duplicate/broken files were removed:
- `src/modeling/run survival models and export summary csv.py` — duplicate of `survival_analysis_km_cox_finegray_xgb_aft.py`
- `src/eda/R- descriptive variables preparation, mutation analysis, and visualization.R` — duplicate of `comprehensive_descriptive_analysis.R`
- `src/eda/Py_missing data_descriptive_analysis_genetic_descriptive_tables_plots.py` — duplicate of `variable_setup_gene_flagging_plots.py`
- `src/collection/import.py` — broken (R syntax in Python); replaced by `fetch_molecular_data.py` + `build_merged_dataset.py`
- `#Prepare descriptive variables and summary statistics for the Genomic project.R` — stale root copy of `src/cleaning/prepare_descriptive_vars.r`
- `.Rhistory` — session artifact

---

## Dependencies

### Python
```
pandas>=2.0, numpy>=1.23, matplotlib>=3.6, seaborn>=0.12,
lifelines>=0.27, xgboost>=1.7, openpyxl>=3.1, statsmodels>=0.14,
scikit-learn>=1.3, requests, scipy, Pillow
```

### R
```
tidyverse, readxl, readr, dplyr, stringr, janitor, skimr, purrr, forcats,
tidyr, gtsummary, knitr, kableExtra, tfrmt, naniar, GGally, ggplot2,
broom, viridisLite, scales, optparse, survival, flexsurv, survminer,
xgboost, SHAPforxgboost, data.table, lubridate
```

---

## Data Privacy

- **PHI data** resides exclusively in `DATA_PRIVATE_DIR` (configured in `.env`, never committed).
- **De-identification** uses HMAC-based date shifting and irreversible ID hashing.
- **Salt files** (`*_salt.txt`, `linkage_map.csv`) are gitignored.
- Only non-PHI processed data in `data/processed/`, `data/features/`, `data/splits/` is committed.
