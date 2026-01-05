# Genomic Analysis Scripts

A collection of genomic analysis scripts for survival analysis, statistical modeling, and exploratory data analysis, organized following data science best practices.

## Project Organization

This project follows the [Cookiecutter Data Science](https://cookiecutter-data-science.drivendata.org) template structure:

```
├── README.md                      <- The top-level README for developers using this project
├── requirements.txt               <- Python package dependencies
│
├── notebooks/                     <- Jupyter notebooks for exploration and prototyping
│
├── references/                    <- Data dictionaries, manuals, and other explanatory materials
│   ├── README_survival_AFT_pipeline.md  <- Documentation for survival analysis pipeline
│   └── *.pdf                      <- Cheat sheets and reference materials
│
└── src/                           <- Source code for this project
    │
    ├── data collection and processing/   <- Data collection, cleaning, and processing code
    │   ├── collection/            <- Gather data from databases, APIs, and other sources
    │   │   ├── extract variables from dataset.py
    │   │   ├── extract variables from dataset.R
    │   │   ├── extract_impact_variables.py
    │   │   └── import.py
    │   │
    │   └── cleaning/              <- Cleaning & processing routines
    │       ├── deidentify.py      <- De-identification scripts for data privacy
    │       ├── deidentify_redcap.py
    │       ├── date_handling_template.R
    │       └── prepare_descriptive_vars.r
    │
    ├── data reports/              <- Generated analysis outputs and exportable tables
    │   ├── reports/               <- Rendered analysis outputs (HTML, PDF, LaTeX, etc.)
    │   └── tables/                <- CSV/TSV/Excel tables and summary tables
    │       ├── table_1_and_2.R    <- Generate descriptive tables
    │       └── BBB regulation table.py
    │
    ├── exploratory data analysis/ <- EDA: use statistics and visualizations to understand data
    │   ├── Py_missing data_descriptive_analysis_genetic_descriptive_tables_plots.py
    │   ├── R- descriptive variables preparation, mutation analysis, and visualization.R
    │   ├── R- generate descriptive and associational plots.R
    │   └── export and generate genetic mut prevalence plots.py
    │
    └── modeling/                  <- Statistical and machine learning models
        ├── run survival models and export summary csv.py  <- Main survival analysis pipeline
        ├── Py- survival and time to event models.py
        ├── Py - xgb_aft_prediction_1020.py              <- XGBoost AFT models
        ├── xgb_aft_shap_summary.py                      <- SHAP explainability
        ├── time_event_analysis_aft.R
        └── shap analysis and plot setup.R
```

## Getting Started

### Installation

Install Python dependencies:

```bash
pip install -r requirements.txt
```

### Usage

Scripts are organized by function:

- **Data Collection**: Scripts in `src/data collection and processing/collection/` extract and import data
- **Data Cleaning**: Scripts in `src/data collection and processing/cleaning/` handle de-identification and preprocessing
- **Exploratory Analysis**: Scripts in `src/exploratory data analysis/` perform statistical summaries and visualizations
- **Modeling**: Scripts in `src/modeling/` run survival analyses, Cox models, XGBoost AFT, and SHAP explanations
- **Reports**: Scripts in `src/data reports/tables/` generate publication-ready tables

For detailed information about the survival analysis pipeline, see [references/README_survival_AFT_pipeline.md](references/README_survival_AFT_pipeline.md).

## Key Features

- **Survival Analysis**: Kaplan-Meier curves, Cox proportional hazards models, Fine-Gray competing risks
- **Machine Learning**: XGBoost Accelerated Failure Time (AFT) models with SHAP explainability
- **Data Privacy**: De-identification scripts for PHI/PII removal
- **Visualization**: Comprehensive plotting for genetic mutations, survival curves, and SHAP summaries

## Requirements

See `requirements.txt` for Python package dependencies. R scripts require:
- tidyverse
- survival
- gtsummary
- SHAPforxgboost

## License

See LICENSE file for details.
