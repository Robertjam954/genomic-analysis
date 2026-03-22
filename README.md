# Genomic Analysis Scripts

A collection of scripts for genomic data analysis, including data processing, exploratory analysis, survival modeling, and visualization.

## Project Organization

```
genomic-analysis-scripts/                          ← PROJECT_ROOT (GitHub)
├── notebooks/              Jupyter notebooks for exploration and prototyping
├── data/
│   ├── processed/          Non-PHI CSVs: metrics, analysis outputs                    ← committed
│   ├── features/           Derived feature sets (gene indicators, embeddings)          ← committed
│   └── splits/             Train/test definitions                                     ← committed
├── docs/
│   ├── executive_summary.md    Full technical executive summary + code map
│   ├── dataset_metadata.md     YAML front-matter dataset specification
│   ├── manuscript/             Project outline, methods documentation
│   └── manuscript_components/  Abstract, appendix, supplementary methods, cover letter
├── conferences/            Conference abstract drafts and submissions
├── reports/
│   ├── figures/            Exported figures for manuscript/presentations               ← committed
│   └── tables/             Exported tables for manuscript/presentations                ← committed
├── references/             Academic papers, cheat sheets, model references
├── src/                    Source code
│   ├── collection/         Scripts to gather data from databases, APIs, etc.
│   ├── cleaning/           Cleaning & processing (de-identification, variable extraction)
│   ├── eda/                Exploratory data analysis scripts
│   ├── modeling/           Survival analysis, XGBoost AFT, SHAP explainability
│   └── reports/            Scripts that generate figures and tables
│       ├── figures/        Figure-generation scripts
│       └── tables/         Table-generation scripts
├── eval/                   Evaluation schemas and metric definitions
├── models/                 Model configurations (XGBoost, Cox PH params)
├── experiments/            Run tracking (run_id, commit hash, results)
├── tools/
│   ├── colab/              Cloud notebook helpers
│   └── watch_repo.sh       Local macOS notifications for remote repo changes
├── .env.example            Path variables + API key template (copy to .env, never commit .env)
├── pyproject.toml          Project metadata + dependencies
├── requirements.txt        Pip-compatible dependency list
└── LICENCE                 MIT

LOCAL ONLY (never committed):
<DATA_PRIVATE_DIR>/                   ← set in .env
├── raw/                    Source data with PHI
├── deidentified/           Redacted outputs + case_document_mapping.csv
├── extracted_text/         Per-case deidentified text files
└── extracted_text_comparison/  Method comparison text outputs
```

## Key Components

### Data Collection and Processing (`src/collection/`, `src/cleaning/`)

**Collection:**
- `fetch_molecular_data.py` — cBioPortal REST API client for molecular profiles, mutations, and clinical data
- `build_merged_dataset.py` — Assemble clinical + mutation API data into merged analysis dataset

**Cleaning:**
- `deidentify.py` / `deidentify_redcap.py` — De-identification scripts for PHI
- `extract variables from dataset.py` / `.R` — Extract and validate variables from datasets
- `extract_impact_variables.py` — Extract genomic variables from IMPACT data
- `prepare_descriptive_vars.r` — Prepare and clean variables for analysis
- `date_handling_template.R` — Template for standardized date handling

### Exploratory Data Analysis (`src/eda/`)

- `variable_setup_gene_flagging_plots.py` — Gene indicator setup, missingness analysis, and visualization
- `comprehensive_descriptive_analysis.R` — Comprehensive descriptive statistics and Table 1/2 generation
- `R- generate descriptive and associational plots.R` — Descriptive visualizations
- `R - generate plots for missing_associational_correlation_matrix_regression_scatterfacetedbyconfounders.R` — Complex associational plots

### Modeling (`src/modeling/`)

**Survival Analysis:**
- `survival_analysis_km_cox_finegray_xgb_aft.py` — Comprehensive survival pipeline (Kaplan-Meier, Cox PH, Fine-Gray, XGBoost AFT)
- `Py- survival and time to event models.py` — Time-to-event modeling
- `time_event_analysis_aft.R` — AFT models in R

**Machine Learning:**
- `Py - xgb_aft_prediction_1020.py` — XGBoost AFT predictions
- `xgb_aft_shap_summary.py` — SHAP analysis for XGBoost models
- `shap analysis and plot setup.R` — SHAP visualization in R

### Report Generation (`src/reports/`)

**Figures** (`src/reports/figures/`):
- `export and assemble_multipanel_figure_for_genetic_variables.py` — Multi-panel figures
- `export and generate genetic mut prevalence plots.py` — Mutation prevalence plots
- `export gene and outcome plots.R` — Gene-outcome visualizations
- `R - export and generate shap analysis and plot generation.R` — SHAP plots

**Tables** (`src/reports/tables/`):
- `table_1_and_2.R` — Manuscript tables
- `BBB regulation table.py` — Blood-brain barrier regulation tables

## Notebooks

Ordered analysis notebooks in `notebooks/`:

| # | Notebook | Description |
|---|---|---|
| 01 | `xgb_aft_1020.ipynb` | XGBoost AFT exploration (original) |
| 02 | `02_dimensionality_reduction_pca.ipynb` | PCA on `G__*` gene indicators — scree plot, biplot, variance analysis |
| 03 | `03_clustering_gene_features.ipynb` | K-Means, hierarchical, DBSCAN, GMM clustering on gene features |
| 04 | `04_regression_gene_expression.ipynb` | Linear & logistic regression baselines before survival modeling |

Notebooks 02–04 adapted from [SalvatoreRa/tutorial](https://github.com/SalvatoreRa/tutorial) genomic series (Apache-2.0).

## Evaluation (`eval/`)

- [`metric_definitions.md`](eval/metric_definitions.md) — All metrics used across the pipeline (C-index, SHAP, AUC, silhouette, etc.) with experiment tracking schema
- [`missing_data_classification.py`](eval/missing_data_classification.py) — Classify column missingness as MCAR / MAR / MNAR via Little's test, point-biserial correlations, and KS distribution tests

```bash
# Run missingness classification
python eval/missing_data_classification.py \
  --input datasets_analysis_dictionary/merged_genie.xlsx \
  --outdir data/processed
```

## Survival and ML Analysis Pipeline

See [`docs/README_survival_AFT_pipeline.md`](docs/README_survival_AFT_pipeline.md) for detailed documentation on the survival analysis and AFT modeling pipeline.

## Quick Start

```bash
# Install Python dependencies
pip install -r requirements.txt

# Or install with dev tools
pip install -e ".[dev]"

# Copy environment template
cp .env.example .env
# Edit .env to set DATA_PRIVATE_DIR and any API keys

# Example: Run survival analysis
python src/modeling/survival_analysis_km_cox_finegray_xgb_aft.py \
  --xlsx datasets_analysis_dictionary/merged_genie.xlsx

# Example: Extract variables from dataset
python "src/cleaning/extract variables from dataset.py" \
  --data /path/to/data.csv \
  --dict /path/to/dictionary.xlsx \
  --outdir output
```

## References

See [`references/REFERENCES_INDEX.md`](references/REFERENCES_INDEX.md) for a full curated index, including:

**Local files:**
- Python cheat sheets and quick references
- R Tidyverse cheat sheets
- Mathematical modeling of the metastatic process
- AFT model assumptions documentation
- Excel quick reference

**External (linked) — from [SalvatoreRa/tutorial](https://github.com/SalvatoreRa/tutorial):**
- Genomic series: PCA, clustering, regression on gene expression data
- Tabular learning: tree-vs-DL benchmarks, missing data (MAR/MCAR/MNAR), synthetic data, KANs
- AI in medicine: scGPT, Med-PaLM, ClinicalGPT, LLMs for gene editing
- Graph ML: NetworkX, iGraph, graph visualization (potential gene-interaction extension)

## Notes

- Data in `data/processed/`, `data/features/`, and `data/splits/` is committed (non-PHI only)
- Private/PHI data lives in `DATA_PRIVATE_DIR` (set in `.env`, never committed)
- Column names are automatically normalized in many scripts (snake_case)
- Most scripts support both CSV and Excel input formats
