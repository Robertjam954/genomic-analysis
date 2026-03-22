# Dataset metadata (YAML front-matter)
dataset_name: "merged_genie_tcga"
dataset_id: "merged_genie_tcga_v1"
dataset_version: "1.0"
description: |
  Merged GENIE BPC / TCGA clinical-genomic dataset for brain metastasis survival
  analysis. Contains de-identified patient demographics, tumor staging, molecular
  markers (MSI, TMB, hypoxia scores, aneuploidy), gene-level alteration indicators
  (G__* columns derived from HUGO symbols), and time-to-event endpoints (DFS, OS).
  Used across the full analysis pipeline: EDA, Cox PH, Fine-Gray competing risks,
  and XGBoost AFT survival modeling.
source: "AACR GENIE BPC / TCGA via cBioPortal; MSK-IMPACT DMP"
source_url: "https://www.cbioportal.org/"
date_acquired: "2024-12-15"
license: "internal-research-use"
access_restrictions: "controlled"
sensitive: true
data_root: "data/raw"
files:
  - path: "datasets_analysis_dictionary/merged_genie.xlsx"
    filename: "merged_genie.xlsx"
    format: "xlsx"
    size_bytes: null
    checksum: null
    rows: null
  - path: "datasets_analysis_dictionary/combined_patient_sample_hypoxia_data.xlsx"
    filename: "combined_patient_sample_hypoxia_data.xlsx"
    format: "xlsx"
    size_bytes: null
    checksum: null
    rows: null
    sheets: ["data_clinical_patient", "data_clinical_supp_hypoxia", "data_clinical_sample"]
  - path: "data/raw/merged_tcga_df.csv"
    filename: "merged_tcga_df.csv"
    format: "csv"
    size_bytes: null
    checksum: null
    rows: null
  - path: "data/raw/filtered_goel_impact_deid.csv"
    filename: "filtered_goel_impact_deid.csv"
    format: "csv"
    size_bytes: null
    checksum: null
    rows: null
  - path: "data/raw/goel_data_dictionary_from_snakecase.xlsx"
    filename: "goel_data_dictionary_from_snakecase.xlsx"
    format: "xlsx"
    size_bytes: null
    checksum: null
    rows: null
schema: |
  # ── Identifiers ──────────────────────────────────────────────────────
  - name: patient_id
    type: string
    description: "De-identified patient identifier (hashed from original MRN)"
    nullable: false
  - name: sample_id
    type: string
    description: "Sequencing sample identifier"
    nullable: true

  # ── Demographics ─────────────────────────────────────────────────────
  - name: age
    type: integer
    description: "Age at diagnosis (years); >89 binned as 90+"
    nullable: true
  - name: age_cat
    type: categorical
    description: "Age category: <50, 50-75, Over 75 (ordered factor)"
    nullable: true
  - name: sex
    type: categorical
    description: "Patient sex: Female, Male"
    nullable: true
  - name: sex_bin
    type: integer
    description: "Binary sex indicator (0=Female, 1=Male)"
    nullable: true
  - name: race
    type: categorical
    description: "Race: White or Caucasian, Black or African American, Asian, Indian or Alaska Native"
    nullable: true
  - name: ethnicity
    type: categorical
    description: "Ethnicity: Non-Hispanic, Hispanic or Latino"
    nullable: true

  # ── Staging & Pathology ──────────────────────────────────────────────
  - name: ajcc_pathologic_tumor_stage
    type: categorical
    description: "AJCC stage (ordered factor): Stage I–IV"
    nullable: true
  - name: subtype
    type: categorical
    description: "Molecular subtype (ordered): Luminal A, Luminal B, HER-2 enriched, Basal-type (TNBC)"
    nullable: true
  - name: variant_type
    type: categorical
    description: "Variant classification: SNP, INS, ONP, DEL"
    nullable: true

  # ── Genomic / Molecular Markers ─────────────────────────────────────
  - name: hugo_symbol
    type: string
    description: "HUGO gene symbol(s) for detected mutations (semicolon- or comma-delimited)"
    nullable: true
  - name: G__*
    type: integer
    description: "Binary gene indicator flags (1=mutation present, 0=absent); derived from hugo_symbol for top-prevalence genes"
    nullable: false
  - name: mutation_count
    type: integer
    description: "Number of distinct gene mutations per sample"
    nullable: true
  - name: top_5_gene_status
    type: integer
    description: "Binary flag: 1 if any of the top-5 most prevalent genes are mutated"
    nullable: true
  - name: mantis_bin
    type: categorical
    description: "Microsatellite instability status: MSI-H (>0.4), MSI-L (<0.4)"
    nullable: true
  - name: tbl_score
    type: float
    description: "Tumor burden / TMB-like continuous score"
    nullable: true
  - name: tbl_score_quartile
    type: categorical
    description: "TBL score quartile (Q1–Q4)"
    nullable: true
  - name: buffa_hypoxia_score
    type: float
    description: "Buffa hypoxia gene signature score (continuous)"
    nullable: true
  - name: buffa_hypoxia_quartile
    type: categorical
    description: "Buffa hypoxia score quartile (Q1–Q4)"
    nullable: true
  - name: aneuploidy_score
    type: float
    description: "Aneuploidy score (continuous)"
    nullable: true
  - name: aneuploidy_score_quartile
    type: categorical
    description: "Aneuploidy score quartile (Q1–Q4)"
    nullable: true
  - name: variant_allele_count
    type: float
    description: "Variant allele count or frequency"
    nullable: true
  - name: variant_allele_count_quartile
    type: categorical
    description: "Variant allele count quartile (Q1–Q4)"
    nullable: true

  # ── Survival / Time-to-Event Endpoints ──────────────────────────────
  - name: dfs_months
    type: float
    description: "Disease-free survival time in months"
    nullable: true
  - name: dfs_status
    type: categorical
    description: "DFS event indicator: Recurrence or Progression of Disease, Disease Free"
    nullable: true
  - name: os_months
    type: float
    description: "Overall survival time in months"
    nullable: true
  - name: os_status
    type: integer
    description: "Overall survival event indicator (1=deceased, 0=alive/censored)"
    nullable: true
key_identifiers: ["patient_id", "sample_id"]
primary_key: "patient_id"
sample_size:
  observations: null        # fill after data load
  units: "patients"
date_range:
  start: "2010-01-01"
  end: "2023-12-31"
update_frequency: "as-needed"
preprocessing_notes: |
  - Column names normalized to snake_case via janitor::clean_names()
  - AJCC stage mapped from numeric (1–4) to ordered factor (Stage I–IV)
  - Subtype mapped from numeric (1–4) to ordered factor (Luminal A/B, HER-2, TNBC)
  - Sex normalized from free text to Female/Male; binary indicator sex_bin created
  - DFS status mapped from 0/1 to labeled factor
  - MANTIS score binarized to MSI-H / MSI-L categories
  - TBL, hypoxia, aneuploidy, variant allele count → quartile bins (Q1–Q4)
  - Gene indicators (G__*) derived from hugo_symbol using top-10 prevalence genes
  - Per-sample mutation_count computed by counting delimited hugo symbols
  - De-identification: date shifting (HMAC-based), MRN hashing, ZIP→ZIP3, age >89→90+
  - PHI scrubbing in free-text fields (regex removal of emails, phones, SSNs, URLs)
missing_value_strategy: |
  - Numeric: NaN preserved for statistical analysis; impute_median for XGBoost AFT
  - Categorical: NA preserved as factor level; mode imputation when required
  - Gene indicators (G__*): NaN → 0 assumption (absence of mutation)
  - See eval/missing_data_classification.py for MCAR/MAR/MNAR diagnostics
transformations: |
  - Continuous scores (TBL, hypoxia, aneuploidy) → quartile categorical variables
  - Hugo symbol string → binary gene indicator matrix (G__*)
  - Age → ordered categorical (age_cat)
  - Date columns → HMAC-shifted dates (de-identification)
  - Log-transform of DFS/OS months explored in AFT modeling
provenance: |
  - created_by: "src/cleaning/prepare_descriptive_vars.r"
  - data_source_import: "src/collection/import.py"
  - deidentification: "src/cleaning/deidentify.py, src/cleaning/deidentify_redcap.py"
  - variable_extraction: "src/cleaning/extract_impact_variables.py"
  - gene_flagging: "src/eda/variable_setup_gene_flagging_plots.py"
  - cbioportal_api: "src/collection/fetch_molecular_data.py"
  - created_on: "2024-12-15"
  - commit: null
contact:
  name: "Robert James"
  email: ""
  role: "analyst"
manifest_file: null
notes: |
  - Sensitive clinical variables present; dataset access controlled.
  - Raw files with PHI reside in DATA_PRIVATE_DIR (see .env); never committed.
  - merged_genie.xlsx is the primary analysis file; other files are upstream sources.
  - Gene indicator columns are prefixed G__ followed by make.names(gene) for R compatibility.
  - Competing event for Fine-Gray model: OS_STATUS used as competing risk to DFS.
  - cBioPortal API functions (src/collection/fetch_molecular_data.py) can refresh molecular data.
  - See references/REFERENCES_INDEX.md for external resource links.
tags: ["clinical", "genomics", "brain-metastasis", "survival", "genie-bpc", "tcga", "msk-impact"]
---

# Dataset README / Summary

## Quick summary
- Name: **merged_genie_tcga**
- Version: **1.0**
- Source: **AACR GENIE BPC / TCGA via cBioPortal; MSK-IMPACT DMP**
- Acquired: **2024-12-15**
- Size: **TBD** rows (fill after data load)

## Files included
- `datasets_analysis_dictionary/merged_genie.xlsx` — Primary merged clinical-genomic dataset (main analysis input)
- `datasets_analysis_dictionary/combined_patient_sample_hypoxia_data.xlsx` — Multi-sheet workbook (patient, hypoxia, sample data)
- `data/raw/merged_tcga_df.csv` — Merged TCGA dataframe with gene alteration, DFS, OS, and XGBoost features
- `data/raw/filtered_goel_impact_deid.csv` — De-identified Goel Lab MSK-IMPACT cohort
- `data/raw/goel_data_dictionary_from_snakecase.xlsx` — Data dictionary for variable extraction/validation

## Schema (high level)

### Identifiers
- `patient_id` — string — De-identified patient ID (hashed)
- `sample_id` — string — Sequencing sample ID

### Demographics
- `age` — integer — Age at diagnosis (>89 binned as 90+)
- `age_cat` — categorical — Age category: <50, 50-75, Over 75
- `sex` / `sex_bin` — categorical / integer — Sex (Female/Male; 0/1)
- `race` — categorical — Race (4 levels)
- `ethnicity` — categorical — Ethnicity (2 levels)

### Staging & Pathology
- `ajcc_pathologic_tumor_stage` — ordered factor — Stage I–IV
- `subtype` — ordered factor — Luminal A, Luminal B, HER-2 enriched, Basal-type (TNBC)
- `variant_type` — categorical — SNP, INS, ONP, DEL

### Genomic / Molecular Markers
- `hugo_symbol` — string — Gene symbol(s) for detected mutations
- `G__*` — integer — Binary gene indicator flags (one per top-prevalence gene)
- `mutation_count` — integer — Number of distinct gene mutations per sample
- `top_5_gene_status` — integer — Any top-5 gene mutated (0/1)
- `mantis_bin` — categorical — MSI status: MSI-H / MSI-L
- `tbl_score` / `tbl_score_quartile` — float / categorical — Tumor burden score and quartile
- `buffa_hypoxia_score` / `buffa_hypoxia_quartile` — float / categorical — Hypoxia signature
- `aneuploidy_score` / `aneuploidy_score_quartile` — float / categorical — Aneuploidy score
- `variant_allele_count` / `variant_allele_count_quartile` — float / categorical — Variant allele count

### Survival Endpoints
- `dfs_months` — float — Disease-free survival (months)
- `dfs_status` — categorical — DFS event (Recurrence/Disease Free)
- `os_months` — float — Overall survival (months)
- `os_status` — integer — OS event (1=deceased, 0=censored)

## Key fields and identifiers
- Primary key: `patient_id`
- ID columns: `patient_id`, `sample_id`

## Preprocessing and transformations

| Task | Description | Tool |
|---|---|---|
| Snake-case normalization | Normalize all column names to snake_case | `janitor::clean_names()` in R |
| De-identification | Date shifting (HMAC), MRN hashing, ZIP→ZIP3, age cap, PHI scrub | `src/cleaning/deidentify.py` |
| Variable extraction & validation | Extract columns per data dictionary, coerce types, validate ranges | `src/cleaning/extract_impact_variables.py` |
| Factor encoding | Encode nominal (race, sex) and ordered (AJCC stage, subtype) factors | `src/cleaning/prepare_descriptive_vars.r` |
| Gene indicator derivation | Parse hugo_symbol → binary `G__*` columns for top-prevalence genes | `src/cleaning/prepare_descriptive_vars.r` |
| Quartile binning | Continuous scores (TBL, hypoxia, aneuploidy) → Q1–Q4 categories | `src/cleaning/prepare_descriptive_vars.r` |
| Missing value handling | Preserve NaN/NA for statistics; impute median/mode for ML models | `eval/missing_data_classification.py` |
| Feature scaling | StandardScaler for regression baselines | `notebooks/04_regression_gene_expression.ipynb` |

Output: Clean, well-structured dataset ready for survival analysis and ML modeling.

## Missingness and quality notes
- Run `python eval/missing_data_classification.py --input <dataset>` for per-column MCAR/MAR/MNAR classification
- Gene indicators (`G__*`): 0% missing by construction (NaN → 0)
- Key clinical variables (age, sex, stage): typically <5% missing
- Survival endpoints (dfs_months, os_months): censoring rate varies by endpoint
- Known issues: `import.py` contains syntax errors (R syntax in Python); see `docs/RESTRUCTURING_NOTES.md`

## Provenance and reproducibility
- Import script: `src/collection/import.py`
- cBioPortal API fetch: `src/collection/fetch_molecular_data.py`
- De-identification: `src/cleaning/deidentify.py`
- Variable preparation: `src/cleaning/prepare_descriptive_vars.r`
- Gene flagging: `src/eda/variable_setup_gene_flagging_plots.py`
- Repo commit: *(fill with `git rev-parse --short HEAD`)*

## Access and licensing
- License: **Internal research use**
- Access restrictions: **Controlled** — raw data contains PHI; stored in `DATA_PRIVATE_DIR` only

## Contact
- Robert James

---

# How to use this template
1. Copy this file into the dataset folder or repository root and fill the YAML front-matter.
2. Keep the `manifest_file` and `checksum` up to date when files change.
3. Link the `commit` field to the code that produces the dataset so others can reproduce it.

# Refreshing molecular data from cBioPortal
Use `src/collection/fetch_molecular_data.py` to pull updated molecular profiles:

```python
from src.collection.fetch_molecular_data import fetch_molecular_data

# Fetch BRCA1/BRCA2 expression across BRCA TCGA
df = fetch_molecular_data(
    entrez_gene_ids=["672", "675"],
    molecular_profile_ids=["brca_tcga_mrna"],
)
```

See the script docstrings for full parameter documentation.
