#!/usr/bin/env python3
"""
Missing Data Classification — MAR / MCAR / MNAR Detection
==========================================================
Adapted from SalvatoreRa/tutorial (Apache-2.0):
  - machine learning/scripts/MAR.py
  - machine learning/scripts/MCAR.py
  - machine learning/scripts/MNAR.py

Provides utilities to:
  1. Run Little's MCAR test to distinguish MCAR from MAR/MNAR.
  2. Compute missingness indicator correlations (MAR detection).
  3. Compare observed vs missing value distributions (MNAR detection).
  4. Generate a per-column missingness classification report.

Usage:
    python eval/missing_data_classification.py \
        --input datasets_analysis_dictionary/merged_genie.xlsx \
        --outdir data/processed

References:
    Little, R. J. A. (1988). A Test of Missing Completely at Random for
    Multivariate Data with Missing Values. JASA, 83(404), 1198-1202.
"""

import argparse
import os
import sys
from itertools import combinations

import numpy as np
import pandas as pd
from scipy import stats


# ── Little's MCAR Test ───────────────────────────────────────────────────────

def littles_mcar_test(df: pd.DataFrame) -> dict:
    """
    Simplified Little's MCAR test.

    Compares group means of observations sharing the same missingness pattern
    against the overall (EM-estimated) means. Under H0 the data are MCAR.

    Returns dict with test_statistic, df, p_value.
    """
    numeric = df.select_dtypes(include=[np.number])
    if numeric.empty:
        return {"test_statistic": np.nan, "df": 0, "p_value": np.nan,
                "conclusion": "No numeric columns"}

    # Missingness patterns (binary mask → tuple key)
    mask = numeric.isna()
    patterns = mask.apply(tuple, axis=1)
    unique_patterns = patterns.unique()

    if len(unique_patterns) <= 1:
        return {"test_statistic": 0.0, "df": 0, "p_value": 1.0,
                "conclusion": "Single missingness pattern → trivially MCAR"}

    grand_mean = numeric.mean()
    grand_cov = numeric.cov()

    d2_sum = 0.0
    dof = 0

    for pat in unique_patterns:
        group = numeric.loc[patterns == pat]
        n_j = len(group)
        if n_j < 2:
            continue
        observed_cols = [c for c, m in zip(numeric.columns, pat) if not m]
        if not observed_cols:
            continue

        group_mean = group[observed_cols].mean()
        diff = group_mean - grand_mean[observed_cols]

        sub_cov = grand_cov.loc[observed_cols, observed_cols]
        try:
            inv_cov = np.linalg.pinv(sub_cov.values)
        except np.linalg.LinAlgError:
            continue

        d2 = float(n_j * diff.values @ inv_cov @ diff.values)
        d2_sum += d2
        dof += len(observed_cols)

    dof -= len(numeric.columns)
    dof = max(dof, 1)

    p_value = 1.0 - stats.chi2.cdf(d2_sum, dof)

    conclusion = "MCAR (fail to reject H0)" if p_value > 0.05 else "Not MCAR (reject H0 at α=0.05)"

    return {
        "test_statistic": round(d2_sum, 4),
        "df": dof,
        "p_value": round(p_value, 6),
        "conclusion": conclusion,
    }


# ── MAR Detection ────────────────────────────────────────────────────────────

def mar_indicator_correlations(df: pd.DataFrame, threshold: float = 0.1) -> pd.DataFrame:
    """
    For each column with missing values, correlate its missingness indicator
    (1=missing, 0=observed) against every other observed column.

    High correlations suggest the missingness depends on other observed
    variables → MAR.

    Returns a DataFrame of (column, correlated_with, pearson_r, p_value).
    """
    numeric = df.select_dtypes(include=[np.number])
    cols_with_na = [c for c in numeric.columns if numeric[c].isna().any()]

    records = []
    for col in cols_with_na:
        indicator = numeric[col].isna().astype(int)
        for other in numeric.columns:
            if other == col:
                continue
            valid = numeric[other].notna() & True  # keep all
            if valid.sum() < 10:
                continue
            x = indicator[valid]
            y = numeric[other][valid].fillna(0)
            if x.std() == 0 or y.std() == 0:
                continue
            r, p = stats.pointbiserialr(x, y)
            if abs(r) >= threshold:
                records.append({
                    "column": col,
                    "correlated_with": other,
                    "pearson_r": round(r, 4),
                    "p_value": round(p, 6),
                })

    return pd.DataFrame(records).sort_values("pearson_r", key=abs, ascending=False)


# ── MNAR Detection ───────────────────────────────────────────────────────────

def mnar_distribution_test(df: pd.DataFrame) -> pd.DataFrame:
    """
    Heuristic MNAR check: for each column with missingness, compare the
    distribution of *other* columns between the missing-group and
    observed-group using a KS test. If no other column explains the
    missingness well (low MAR signal), but the column's own distribution
    is skewed near the missing boundary, MNAR is suspected.

    Returns a DataFrame of (column, pct_missing, ks_stat, ks_pvalue, mnar_flag).
    """
    numeric = df.select_dtypes(include=[np.number])
    cols_with_na = [c for c in numeric.columns if numeric[c].isna().any()]

    records = []
    for col in cols_with_na:
        pct = numeric[col].isna().mean() * 100
        observed = numeric[col].dropna()

        # Compare observed values near the tails: if missingness is
        # concentrated at extremes of the observed range, MNAR is plausible.
        if len(observed) < 20:
            records.append({
                "column": col,
                "pct_missing": round(pct, 2),
                "ks_stat": np.nan,
                "ks_pvalue": np.nan,
                "mnar_flag": "insufficient data",
            })
            continue

        q10, q90 = observed.quantile([0.10, 0.90])
        lower_tail = observed[observed <= q10]
        upper_tail = observed[observed >= q90]
        middle = observed[(observed > q10) & (observed < q90)]

        if len(middle) < 5 or (len(lower_tail) + len(upper_tail)) < 5:
            ks_stat, ks_p = np.nan, np.nan
        else:
            tails = pd.concat([lower_tail, upper_tail])
            ks_stat, ks_p = stats.ks_2samp(tails, middle)

        mnar_flag = "possible MNAR" if (ks_p is not np.nan and ks_p < 0.05) else "unlikely"

        records.append({
            "column": col,
            "pct_missing": round(pct, 2),
            "ks_stat": round(ks_stat, 4) if not np.isnan(ks_stat) else np.nan,
            "ks_pvalue": round(ks_p, 6) if not np.isnan(ks_p) else np.nan,
            "mnar_flag": mnar_flag,
        })

    return pd.DataFrame(records)


# ── Unified Report ───────────────────────────────────────────────────────────

def classify_missingness(df: pd.DataFrame, mar_threshold: float = 0.1) -> pd.DataFrame:
    """
    Produce a per-column missingness classification summary.
    """
    mcar_result = littles_mcar_test(df)
    mar_df = mar_indicator_correlations(df, threshold=mar_threshold)
    mnar_df = mnar_distribution_test(df)

    numeric = df.select_dtypes(include=[np.number])
    cols_with_na = [c for c in numeric.columns if numeric[c].isna().any()]

    # Columns with strong MAR signal
    mar_cols = set(mar_df["column"].unique()) if not mar_df.empty else set()
    mnar_cols = set(mnar_df.loc[mnar_df["mnar_flag"] == "possible MNAR", "column"]) if not mnar_df.empty else set()

    records = []
    for col in cols_with_na:
        pct = numeric[col].isna().mean() * 100
        if col in mar_cols and col in mnar_cols:
            classification = "MAR + possible MNAR"
        elif col in mar_cols:
            classification = "MAR"
        elif col in mnar_cols:
            classification = "possible MNAR"
        elif mcar_result["p_value"] > 0.05:
            classification = "MCAR"
        else:
            classification = "indeterminate"

        records.append({
            "column": col,
            "pct_missing": round(pct, 2),
            "classification": classification,
        })

    report = pd.DataFrame(records).sort_values("pct_missing", ascending=False)
    return report, mcar_result, mar_df, mnar_df


# ── CLI ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="Classify missingness patterns (MCAR/MAR/MNAR)")
    parser.add_argument("--input", required=True, help="Path to CSV or Excel dataset")
    parser.add_argument("--outdir", default="data/processed", help="Output directory")
    parser.add_argument("--mar-threshold", type=float, default=0.1,
                        help="Minimum |r| to flag a MAR correlation (default 0.1)")
    parser.add_argument("--sheet", default=None, help="Excel sheet name (optional)")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # Load
    ext = os.path.splitext(args.input)[1].lower()
    if ext in (".xlsx", ".xls"):
        df = pd.read_excel(args.input, sheet_name=args.sheet)
    else:
        df = pd.read_csv(args.input)

    print(f"Loaded {args.input}: {df.shape[0]} rows × {df.shape[1]} cols")

    report, mcar, mar_df, mnar_df = classify_missingness(df, mar_threshold=args.mar_threshold)

    # Print summary
    print("\n── Little's MCAR Test ──────────────────────────────")
    for k, v in mcar.items():
        print(f"  {k}: {v}")

    print(f"\n── Missingness Report ({len(report)} columns with missing data) ──")
    if not report.empty:
        print(report.to_string(index=False))

    # Save
    report.to_csv(os.path.join(args.outdir, "missingness_classification.csv"), index=False)
    if not mar_df.empty:
        mar_df.to_csv(os.path.join(args.outdir, "mar_correlations.csv"), index=False)
    if not mnar_df.empty:
        mnar_df.to_csv(os.path.join(args.outdir, "mnar_distribution_tests.csv"), index=False)

    print(f"\nOutputs saved to {args.outdir}/")


if __name__ == "__main__":
    main()
