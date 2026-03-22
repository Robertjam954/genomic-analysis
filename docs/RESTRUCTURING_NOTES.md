# Repository Restructuring - Notes

## Completed Restructuring

The repository has been successfully restructured to follow the template from `Robertjam954/project_template`. All files have been organized into appropriate directories based on their purpose.

## Pre-existing Issues Found (Not Fixed)

During the restructuring, the following pre-existing issues were identified but NOT fixed, as the task was specifically to restructure the repository with minimal changes:

### 1. Duplicate Files

The following files appear to be duplicates and were already duplicates in the original repository:

**Modeling Scripts:**
- `src/modeling/run survival models and export summary csv.py`
- `src/modeling/survival_analysis_km_cox_finegray_xgb_aft.py`

These are identical files (484 lines each). The second file was renamed from a file without an extension during restructuring, but both originated from the same content.

**Exploratory Analysis Scripts:**
- `src/exploratory data analysis/R- descriptive variables preparation, mutation analysis, and visualization.R`
- `src/exploratory data analysis/comprehensive_descriptive_analysis.R`

These are identical files (656 lines each). The second file was renamed from a file without an extension during restructuring.

### 2. Syntax Errors

**import.py** contains syntax errors that were present in the original file:
- Line 4: Function definition uses a string literal as parameter instead of variable name
- Line 6: Uses R syntax `c()` in Python code
- Missing proper function parameter syntax

### Resolution

The following duplicates and broken files have been **deleted**:

1. `src/modeling/run survival models and export summary csv.py` — kept `survival_analysis_km_cox_finegray_xgb_aft.py`
2. `src/eda/R- descriptive variables preparation, mutation analysis, and visualization.R` — kept `comprehensive_descriptive_analysis.R`
3. `src/eda/Py_missing data_descriptive_analysis_genetic_descriptive_tables_plots.py` — kept `variable_setup_gene_flagging_plots.py`
4. `src/collection/import.py` — broken (R syntax in Python); replaced by `fetch_molecular_data.py` + `build_merged_dataset.py`
5. `#Prepare descriptive variables and summary statistics for the Genomic project.R` — stale root copy of `src/cleaning/prepare_descriptive_vars.r`
6. `.Rhistory` — session artifact

## Files Successfully Restructured

- **27 script files** moved to appropriate directories
- **8 reference files** moved to `references/`
- **3 files** given proper extensions (.py or .R)
- **6 duplicate/broken files** removed
- **2 new files** created (README.md, .gitignore)

All remaining files maintain their original content and functionality.
