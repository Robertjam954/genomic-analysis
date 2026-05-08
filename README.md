# Genomic Analysis Scripts

> Clinical-genomic research pipeline for brain metastasis survival analysis — from raw data ingestion through manuscript-ready figures and tables.

[![Python](https://img.shields.io/badge/Python-3.10%2B-3776AB?logo=python&logoColor=white)](https://www.python.org/)
[![R](https://img.shields.io/badge/R-4.x-276DC3?logo=r&logoColor=white)](https://www.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENCE)

---

## About

**genomic-analysis-scripts** is a clinical-genomic analysis repository focused on brain metastasis outcomes in breast cancer. It combines data engineering, exploratory analysis, survival modeling, and publication-oriented reporting in one reproducible workflow.

---

## Overview

**genomic-analysis-scripts** ingests AACR GENIE BPC, TCGA, and MSK-IMPACT data and runs a reproducible end-to-end pipeline:

| Stage | What happens |
|---|---|
| **Collection** | Pull molecular profiles, mutations, and clinical data from the cBioPortal REST API |
| **Cleaning** | De-identify PHI, extract and standardize clinical variables |
| **EDA** | Gene indicator setup, missingness analysis, descriptive statistics |
| **Modeling** | Kaplan-Meier, Cox PH, Fine-Gray competing risks, XGBoost AFT |
| **Reports** | Manuscript-ready figures (PNG/PDF) and tables (CSV/Excel) |

**Domain:** Breast cancer brain metastasis — gene-level alteration indicators, microsatellite instability, hypoxia scores, aneuploidy, and time-to-event endpoints (DFS, OS).

---

## Table of Contents

1. [About](#about)
2. [Project Structure](#project-structure)
3. [Quick Start](#quick-start)
4. [Source Code](#source-code)
5. [Notebooks](#notebooks)
6. [Evaluation](#evaluation)
7. [Documentation](#documentation)
8. [References](#references)
9. [Notes](#notes)

---

## Project Structure

```
genomic-analysis-scripts/
│
├── src/                        Source code (Python & R)
│   ├── collection/             Data acquisition from cBioPortal API
│   ├── cleaning/               De-identification and variable extraction
│   ├── eda/                    Exploratory data analysis
│   ├── modeling/               Survival models and machine learning
│   └── reports/
│       ├── figures/            Figure-generation scripts
│       └── tables/             Table-generation scripts
│
├── notebooks/                  Jupyter notebooks for exploration and prototyping
│
├── data/                       Non-PHI data (generated locally, not committed)
│   ├── processed/              Analysis outputs (metrics, classification results)
│   ├── features/               Derived feature sets (gene indicators, embeddings)
│   └── splits/                 Train/test definitions
│
├── reports/                    Generated outputs (not committed)
│   ├── figures/                Exported figures (PNG/PDF)
│   └── tables/                 Exported tables (CSV/Excel)
│
├── docs/                       Documentation
│   ├── executive_summary.md    Technical overview and code map
│   ├── dataset_metadata.md     Dataset specification
│   ├── README_survival_AFT_pipeline.md  Pipeline deep-dive
│   ├── manuscript/             Project outline and methods
│   └── manuscript_components/  Abstract, appendix, supplementary materials
│
├── eval/                       Evaluation schemas and metric definitions
├── models/                     Model configurations (XGBoost, Cox PH params)
├── experiments/                Run tracking (run_id, commit hash, results)
├── conferences/                Conference abstract drafts and submissions
├── references/                 Papers, cheat sheets, and external resource index
├── tools/
│   ├── colab/                  Cloud notebook helpers
│   └── watch_repo.sh           macOS notifications for remote repo changes
│
├── .env.example                Path and API key template (copy → .env, never commit)
├── pyproject.toml              Project metadata and dependencies
├── requirements.txt            Pip-compatible dependency list
└── LICENCE                     MIT

LOCAL ONLY — never committed:
<DATA_PRIVATE_DIR>/             Set in .env
├── raw/                        Source data with PHI
├── deidentified/               Redacted outputs + case_document_mapping.csv
├── extracted_text/             Per-case de-identified text files
└── extracted_text_comparison/  Method comparison text outputs
```

---

## Quick Start

```bash
# 1. Install Python dependencies
pip install -r requirements.txt

# Or install with dev tools (Jupyter, Black, Ruff)
pip install -e ".[dev]"

# 2. Configure environment
cp .env.example .env
# Edit .env to set DATA_PRIVATE_DIR and any API keys

# 3. Run survival analysis
python src/modeling/survival_analysis_km_cox_finegray_xgb_aft.py \
  --xlsx /path/to/merged_genie.xlsx

# 4. Extract variables from dataset
python "src/cleaning/extract variables from dataset.py" \
  --data /path/to/data.csv \
  --dict /path/to/dictionary.xlsx \
  --outdir output
```

---

## Source Code

### Data Collection — `src/collection/`

| Script | Description |
|---|---|
| `fetch_molecular_data.py` | cBioPortal REST API client — molecular profiles, mutations, clinical data |
| `build_merged_dataset.py` | Assembles clinical + mutation API data into merged analysis dataset |

### Cleaning — `src/cleaning/`

| Script | Description |
|---|---|
| `deidentify.py` | De-identification for PHI |
| `deidentify_redcap.py` | REDCap-specific de-identification |
| `extract variables from dataset.py / .R` | Extract and validate variables from datasets |
| `extract_impact_variables.py` | Extract genomic variables from IMPACT data |
| `prepare_descriptive_vars.r` | Prepare and clean variables for analysis |
| `date_handling_template.R` | Standardized date handling template |

### Exploratory Data Analysis — `src/eda/`

| Script | Description |
|---|---|
| `variable_setup_gene_flagging_plots.py` | Gene indicator setup, missingness analysis, visualization |
| `comprehensive_descriptive_analysis.R` | Descriptive statistics and Table 1/2 generation |
| `R- generate descriptive and associational plots.R` | Descriptive visualizations |
| `R - generate plots for missing_associational...R` | Complex associational and correlation plots |

### Modeling — `src/modeling/`

| Script | Description |
|---|---|
| `survival_analysis_km_cox_finegray_xgb_aft.py` | Full survival pipeline (KM, Cox PH, Fine-Gray, XGBoost AFT) |
| `Py- survival and time to event models.py` | Time-to-event modeling |
| `time_event_analysis_aft.R` | AFT models in R |
| `Py - xgb_aft_prediction_1020.py` | XGBoost AFT predictions |
| `xgb_aft_shap_summary.py` | SHAP analysis for XGBoost models |
| `shap analysis and plot setup.R` | SHAP visualization in R |

### Report Generation — `src/reports/`

| Script | Description |
|---|---|
| `figures/export and assemble_multipanel_figure_for_genetic_variables.py` | Multi-panel figures |
| `figures/export and generate genetic mut prevalence plots.py` | Mutation prevalence plots |
| `figures/export gene and outcome plots.R` | Gene-outcome visualizations |
| `figures/R - export and generate shap analysis and plot generation.R` | SHAP plots |
| `tables/table_1_and_2.R` | Manuscript tables |
| `tables/BBB regulation table.py` | Blood-brain barrier regulation tables |

---

## Notebooks

Ordered analysis notebooks in `notebooks/`:

| # | Notebook | Description |
|---|---|---|
| 01 | `xgb_aft_1020.ipynb` | XGBoost AFT exploration (original) |
| 02 | `02_dimensionality_reduction_pca.ipynb` | PCA on gene indicators — scree plot, biplot, variance analysis |
| 03 | `03_clustering_gene_features.ipynb` | K-Means, hierarchical, DBSCAN, GMM clustering on gene features |
| 04 | `04_regression_gene_expression.ipynb` | Linear and logistic regression baselines before survival modeling |

> Notebooks 02–04 adapted from [SalvatoreRa/tutorial](https://github.com/SalvatoreRa/tutorial) genomic series (Apache-2.0).

---

## Evaluation

See [`eval/metric_definitions.md`](eval/metric_definitions.md) for all metrics used across the pipeline (C-index, SHAP, AUC, silhouette, etc.) and the experiment tracking schema.

```bash
# Classify column missingness as MCAR / MAR / MNAR
python eval/missing_data_classification.py \
  --input /path/to/merged_genie.xlsx \
  --outdir data/processed
```

---

## Documentation

| Document | Description |
|---|---|
| [`docs/executive_summary.md`](docs/executive_summary.md) | Full technical overview and code map |
| [`docs/dataset_metadata.md`](docs/dataset_metadata.md) | YAML dataset specification |
| [`docs/README_survival_AFT_pipeline.md`](docs/README_survival_AFT_pipeline.md) | Survival analysis and AFT pipeline deep-dive |

---

## References

See [`references/REFERENCES_INDEX.md`](references/REFERENCES_INDEX.md) for the full curated index, including:

**Local files:**
- Python and R (Tidyverse) cheat sheets
- Mathematical modeling of the metastatic process
- AFT model assumptions documentation

**External (linked) — from [SalvatoreRa/tutorial](https://github.com/SalvatoreRa/tutorial):**
- Genomic series: PCA, clustering, regression on gene expression data
- Tabular learning: tree-vs-DL benchmarks, missing data, synthetic data, KANs
- AI in medicine: scGPT, Med-PaLM, ClinicalGPT, LLMs for gene editing
- Graph ML: NetworkX, iGraph, graph visualization

---

## Notes

- `reports/` and `data/` outputs are **not committed** — generate them locally by running the pipeline
- Private/PHI data lives in `DATA_PRIVATE_DIR` (set in `.env`, never committed)
- Column names are automatically normalized (snake_case) across most scripts
- Most scripts support both CSV and Excel input formats
