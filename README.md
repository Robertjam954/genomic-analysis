# Genomic Analysis Scripts

A collection of scripts for genomic data analysis, including data processing, exploratory analysis, survival modeling, and visualization.

## Project Organization

This repository follows a structured data science workflow:

```
├── README.md                      <- The top-level README for developers using this project
├── requirements.txt               <- Python package dependencies
│
├── data/                          <- Data directories (not tracked by default)
│   ├── external                   <- Data from third party sources
│   ├── interim                    <- Intermediate data that has been transformed
│   ├── processed                  <- The final, canonical data sets for modeling
│   └── raw                        <- The original, immutable data dump
│
├── notebooks/                     <- Jupyter notebooks for exploration and prototyping
│
├── references/                    <- Data dictionaries, manuals, cheat sheets, and other explanatory materials
│
└── src/                           <- Source code for this project
    │
    ├── data collection and processing/   <- Data collection, cleaning, and processing code
    │   ├── collection/                    <- Scripts to gather data from databases, APIs, etc.
    │   └── cleaning/                      <- Cleaning & processing routines (encoding, missing values, etc.)
    │
    ├── data reports/                <- Generated analysis outputs and exportable tables
    │   ├── reports/                  <- Rendered analysis outputs (figures, visualizations)
    │   └── tables/                   <- CSV/TSV/Excel tables and summary tables
    │
    ├── exploratory data analysis/   <- EDA scripts for understanding data patterns and relationships
    │
    └── modeling/                    <- Statistical and machine learning models
        ├── Survival analysis scripts (Cox PH, Kaplan-Meier, Fine-Gray)
        ├── XGBoost AFT models
        └── SHAP analysis and explainability
```

## Key Components

### Data Collection and Processing

**Collection:**
- `import.py` - Import data from various sources

**Cleaning:**
- `deidentify.py` / `deidentify_redcap.py` - De-identification scripts for protected health information
- `extract variables from dataset.py` / `.R` - Extract and validate variables from datasets
- `extract_impact_variables.py` - Extract genomic variables from IMPACT data
- `prepare_descriptive_vars.r` - Prepare and clean variables for analysis
- `date_handling_template.R` - Template for standardized date handling

### Exploratory Data Analysis

- `Py_missing data_descriptive_analysis_genetic_descriptive_tables_plots.py` - Missingness analysis and descriptive statistics
- `variable_setup_gene_flagging_plots.py` - Gene indicator setup and visualization
- `comprehensive_descriptive_analysis.R` - Comprehensive descriptive statistics and Table 1/2 generation
- `R- descriptive variables preparation, mutation analysis, and visualization.R` - Mutation analysis
- `R- generate descriptive and associational plots.R` - Descriptive visualizations
- `R - generate plots for missing_associational_correlation_matrix_regression_scatterfacetedbyconfounders.R` - Complex associational plots

### Modeling

**Survival Analysis:**
- `survival_analysis_km_cox_finegray_xgb_aft.py` - Comprehensive survival pipeline (Kaplan-Meier, Cox PH, Fine-Gray competing risks, XGBoost AFT)
- `Py- survival and time to event models.py` - Time-to-event modeling
- `run survival models and export summary csv.py` - Batch survival model execution
- `time_event_analysis_aft.R` - AFT models in R

**Machine Learning:**
- `Py - xgb_aft_prediction_1020.py` - XGBoost AFT predictions
- `xgb_aft_shap_summary.py` - SHAP analysis for XGBoost models
- `shap analysis and plot setup.R` - SHAP visualization in R

### Data Reports

**Reports (Figures/Visualizations):**
- `export and assemble_multipanel_figure_for_genetic_variables.py` - Multi-panel figures
- `export and generate genetic mut prevalence plots.py` - Mutation prevalence plots
- `export gene and outcome plots.R` - Gene-outcome visualizations
- `R - export and generate shap analysis and plot generation.R` - SHAP plots

**Tables:**
- `table_1_and_2.R` - Manuscript tables
- `BBB regulation table.py` - Blood-brain barrier regulation tables

## Survival and ML Analysis Pipeline

See `README_survival_AFT_pipeline.md` for detailed documentation on the survival analysis and AFT modeling pipeline.

### Quick Start

```bash
# Install Python dependencies
pip install -r requirements.txt

# Example: Run survival analysis
python src/modeling/survival_analysis_km_cox_finegray_xgb_aft.py \
  --xlsx datasets_analysis_dictionary/merged_genie.xlsx

# Example: Extract variables from dataset
python "src/data collection and processing/cleaning/extract variables from dataset.py" \
  --data /path/to/data.csv \
  --dict /path/to/dictionary.xlsx \
  --outdir output
```

## References

The `references/` directory contains:
- Python cheat sheets and quick references
- R Tidyverse cheat sheets
- Mathematical modeling references
- AFT model assumptions documentation
- Excel quick reference

## Notes

- Data directories (`data/`) should typically be excluded from version control
- Keep sensitive data secure and use de-identification scripts before sharing
- Column names are automatically normalized in many scripts (snake_case)
- Most scripts support both CSV and Excel input formats
