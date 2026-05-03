# Evaluation Schemas & Metric Definitions

Metrics and evaluation criteria used across the genomic analysis pipeline.

---

## Survival Analysis Metrics

| Metric | Description | Script(s) |
|---|---|---|
| **Concordance Index (C-index)** | Proportion of concordant pairs in time-to-event prediction. Range [0, 1]; 0.5 = random. | `src/modeling/survival_analysis_km_cox_finegray_xgb_aft.py` |
| **Log-rank p-value** | Tests equality of survival curves across groups (Kaplan-Meier). | `src/modeling/survival_analysis_km_cox_finegray_xgb_aft.py` |
| **Hazard Ratio (HR)** | Exponentiated Cox PH coefficient; HR > 1 = increased risk. | `src/modeling/survival_analysis_km_cox_finegray_xgb_aft.py` |
| **Schoenfeld residuals** | Test proportional hazards assumption; significant p → time-varying effect. | `src/modeling/survival_analysis_km_cox_finegray_xgb_aft.py` |
| **AIC / BIC** | Model selection for AFT distributions (Weibull, log-normal, log-logistic). | `src/modeling/time_event_analysis_aft.R` |

## XGBoost AFT Metrics

| Metric | Description | Script(s) |
|---|---|---|
| **AFT-nloglik** | Negative log-likelihood under the AFT survival model. | `src/modeling/Py - xgb_aft_prediction_1020.py` |
| **RMSE (predicted vs observed time)** | Root mean squared error on predicted survival time (uncensored subset). | `src/modeling/Py - xgb_aft_prediction_1020.py` |
| **SHAP feature importance** | Mean |SHAP value| per feature; identifies drivers of predicted survival time. | `src/modeling/xgb_aft_shap_summary.py` |

## Baseline Regression Metrics (Notebook 04)

| Metric | Description | Notebook |
|---|---|---|
| **R²** | Proportion of variance explained by linear regression. | `notebooks/04_regression_gene_expression.ipynb` |
| **MSE** | Mean squared error for continuous outcome prediction. | `notebooks/04_regression_gene_expression.ipynb` |
| **AUC-ROC** | Area under receiver operating characteristic; discrimination for binary outcome. | `notebooks/04_regression_gene_expression.ipynb` |
| **Accuracy** | Proportion of correct binary predictions. | `notebooks/04_regression_gene_expression.ipynb` |

## Clustering Quality Metrics (Notebook 03)

| Metric | Description | Notebook |
|---|---|---|
| **Silhouette Score** | Mean ratio of within-cluster vs between-cluster distance. Range [-1, 1]; higher = better separation. | `notebooks/03_clustering_gene_features.ipynb` |
| **Inertia (WCSS)** | Within-cluster sum of squares (K-Means). Used in elbow plots. | `notebooks/03_clustering_gene_features.ipynb` |

## Missing Data Evaluation

| Metric | Description | Script |
|---|---|---|
| **Little's MCAR χ² test** | Global test for MCAR; p > 0.05 → fail to reject MCAR. | `eval/missing_data_classification.py` |
| **Point-biserial r** | Correlation between missingness indicator and observed columns (MAR detection). | `eval/missing_data_classification.py` |
| **KS statistic** | Kolmogorov-Smirnov test comparing tail vs middle distributions (MNAR heuristic). | `eval/missing_data_classification.py` |

---

## Experiment Tracking Schema

Each run logged in `experiments/` should capture:

```yaml
run_id: <uuid>
timestamp: <ISO 8601>
commit_hash: <git short hash>
script: <path to script>
dataset: <path or name>
parameters:
  model: <model name>
  hyperparameters: {}
  features: <list or count>
  train_n: <int>
  test_n: <int>
results:
  primary_metric: <value>
  secondary_metrics: {}
notes: <free text>
```
