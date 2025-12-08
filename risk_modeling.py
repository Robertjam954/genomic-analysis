#!/usr/bin/env python3
"""
Risk modeling for OS and PFS.

Produces:
- Kaplan-Meier curves by subtype (KM PNGs) and log-rank test results
- Multivariable Cox PH models for OS and PFS (CSV summaries)
- Gene-level feature selection: per-gene Cox PH tests (HR, p) and FDR q-values

Notes:
- Requires `lifelines` and `statsmodels` for modeling and multiple-testing correction.
"""
import os
import argparse
import warnings
from collections import Counter

import pandas as pd
import numpy as np

try:
    from lifelines import KaplanMeierFitter, CoxPHFitter
    from lifelines.statistics import logrank_test
except Exception as e:
    raise ImportError("This script requires the 'lifelines' package. Install with `pip install lifelines` or use your conda env.")

try:
    from statsmodels.stats.multitest import multipletests
except Exception:
    # provide a lightweight Benjamini-Hochberg implementation as a fallback
    def multipletests(pvals, alpha=0.05, method='fdr_bh'):
        """Return (reject, pvals_corrected, alphacSidak, alphacBonf) similar to statsmodels API
        Only supports 'fdr_bh' method here."""
        import numpy as _np
        p = _np.asarray(pvals)
        n = p.size
        sort_idx = _np.argsort(p)
        sort_p = p[sort_idx]
        bh = _np.empty(n, dtype=float)
        bh[sort_idx] = _np.minimum.accumulate((sort_p * n) / (_np.arange(1, n + 1)))
        qvals = _np.minimum(bh, 1.0)
        reject = qvals <= alpha
        return reject, qvals, None, None

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def ensure_dir(path):
    os.makedirs(path, exist_ok=True)


def parse_event_column(s):
    # Accept many encodings: numeric 1/0, '1:DECEASED', '0:LIVING', True/False, DEAD/ALIVE, etc.
    # Return a boolean Series (True = event occurred)
    s = s.copy()
    # If numeric dtype, treat >0 as event
    try:
        if pd.api.types.is_numeric_dtype(s):
            return (s.fillna(0).astype(float) > 0).astype(bool)
    except Exception:
        pass

    s = s.astype(str).fillna("").str.strip().str.upper()
    # handle '0:LIVING' or '1:DECEASED' by taking part before ':' if present
    has_colon = s.str.contains(":")
    if has_colon.any():
        left = s.where(~has_colon, s.str.split(":", n=1).str[0])
    else:
        left = s

    # Map common true values
    true_vals = set(["1", "TRUE", "T", "Y", "YES", "DECEASED", "DEAD", "D"])    
    false_vals = set(["0", "FALSE", "F", "N", "NO", "ALIVE", "LIVING", "ALIVE;LIVING"]) 
    mapped = left.isin(list(true_vals))
    # If explicit false tokens present, ensure False
    mapped = mapped & (~left.isin(list(false_vals))) | (left.isin(list(true_vals)))
    # As a fallback, look for keywords inside the original string
    mapped = mapped | s.str.contains("DECEASED|DEAD|DIED")
    mapped = mapped & ~s.str.contains("ALIVE|LIVING")
    mapped = mapped.fillna(False).astype(bool)
    return mapped


def fit_km_by_subtype(df, time_col, event_col, group_col, out_png, out_csv):
    kmf = KaplanMeierFitter()
    # coerce types for lifelines
    df = df.copy()
    df[time_col] = pd.to_numeric(df[time_col], errors="coerce")
    df[event_col] = parse_event_column(df[event_col]).astype(int)
    groups = df[group_col].fillna("<missing>").unique()
    plt.figure(figsize=(8, 6))
    results = []
    MIN_GROUP = 5
    for g in sorted(groups):
        sub = df[df[group_col].fillna("<missing>") == g].copy()
        # drop rows with missing time or event
        sub = sub.dropna(subset=[time_col, event_col])
        if len(sub) < MIN_GROUP:
            # skip very small groups
            continue
        # lifelines expects numeric durations and boolean/int event
        try:
            kmf.fit(sub[time_col].astype(float), event_observed=sub[event_col].astype(int), label=str(g))
            kmf.plot_survival_function(ci_show=True)
            results.append((g, len(sub)))
        except Exception as e:
            print(f"KM fit failed for group {g}: {e}")
            continue
    plt.title(f"Kaplan-Meier by {group_col}")
    plt.xlabel("Time (months)")
    plt.ylabel("Survival probability")
    plt.tight_layout()
    plt.savefig(out_png)
    plt.close()

    # pairwise log-rank: for brevity, do overall test between groups (logrank with groups combined)
    # lifelines' multigroup logrank can be approximated by pairwise tests; we'll compute pairwise and save table
    rows = []
    groups = sorted(df[group_col].fillna("<missing>").unique())
    for i in range(len(groups)):
        for j in range(i + 1, len(groups)):
            g1 = groups[i]
            g2 = groups[j]
            s1 = df[df[group_col].fillna("<missing>") == g1]
            s2 = df[df[group_col].fillna("<missing>") == g2]
            if len(s1) < 5 or len(s2) < 5:
                continue
            s1 = s1.dropna(subset=[time_col, event_col])
            s2 = s2.dropna(subset=[time_col, event_col])
            if len(s1) < 5 or len(s2) < 5:
                continue
            try:
                lr = logrank_test(s1[time_col], s2[time_col], event_observed_A=s1[event_col].astype(int), event_observed_B=s2[event_col].astype(int))
                rows.append({"group1": g1, "group2": g2, "p_value": lr.p_value})
            except Exception as e:
                print(f"Log-rank test failed for {g1} vs {g2}: {e}")
                continue
            rows.append({"group1": g1, "group2": g2, "p_value": lr.p_value})
    if len(rows) > 0:
        pd.DataFrame(rows).to_csv(out_csv, index=False)


def prepare_df_for_cox(df, time_col, event_col, covariates):
    sub = df[[time_col, event_col] + covariates].copy()
    # convert event to boolean
    sub[event_col] = parse_event_column(sub[event_col]).astype(int)
    # one-hot encode categorical covariates
    sub = pd.get_dummies(sub, columns=[c for c in covariates if sub[c].dtype == object or sub[c].dtype == "category"], drop_first=True)
    # convert numeric-like strings to numeric where possible
    for c in sub.columns:
        if c in [time_col, event_col]:
            continue
        if sub[c].dtype == object:
            try:
                sub[c] = pd.to_numeric(sub[c], errors="raise")
            except Exception:
                pass
    sub = sub.dropna()
    return sub


def fit_cox(df, time_col, event_col, covariates, out_csv):
    sub = prepare_df_for_cox(df, time_col, event_col, covariates)
    if len(sub) < 20:
        warnings.warn(f"Too few rows ({len(sub)}) to fit CoxPH reliably. Skipping.")
        return None
    cph = CoxPHFitter()
    try:
        cph.fit(sub, duration_col=time_col, event_col=event_col)
    except Exception as e:
        print("CoxPH fit failed (retrying with L2 penalizer):", e)
        try:
            cph = CoxPHFitter(penalizer=0.1)
            cph.fit(sub, duration_col=time_col, event_col=event_col)
        except Exception as e2:
            print("CoxPH penalized fit also failed:", e2)
            return None
    summary = cph.summary.reset_index().rename(columns={"index": "covariate"})
    summary.to_csv(out_csv, index=False)
    return summary


def build_gene_matrix(df, gene_list_col, top_k=50, sep=";"):
    # explode gene list and build binary matrix for top_k genes
    s = df[gene_list_col].fillna("").astype(str)
    all_items = []
    per_row = []
    for v in s:
        items = [it.strip() for it in v.split(sep) if it.strip()]
        per_row.append(items)
        all_items.extend(items)
    counts = Counter(all_items)
    top = [g for g, _ in counts.most_common(top_k)]
    mat = pd.DataFrame(0, index=df.index, columns=top)
    for i, items in enumerate(per_row):
        for it in items:
            if it in mat.columns:
                mat.at[i, it] = 1
    return mat, counts


def gene_feature_tests(df, time_col, event_col, covariates, gene_mat, out_csv_prefix):
    rows = []
    for gene in gene_mat.columns:
        tmp = df[[time_col, event_col] + covariates].copy()
        tmp[gene] = gene_mat[gene].values
        tmp[event_col] = parse_event_column(tmp[event_col]).astype(int)
        tmp = pd.get_dummies(tmp, drop_first=True)
        tmp = tmp.dropna()
        if tmp[gene].sum() < 5:
            continue
        cph = CoxPHFitter()
        try:
            cph.fit(tmp, duration_col=time_col, event_col=event_col)
            s = cph.summary.loc[[gene]].reset_index()
            hr = float(np.exp(s.loc[0, "coef"]))
            p = float(s.loc[0, "p"])
            rows.append({"gene": gene, "hr": hr, "p": p, "n_mut": int(gene_mat[gene].sum())})
        except Exception as e:
            # skip genes that cause model failures
            continue
    if len(rows) == 0:
        return None
    gdf = pd.DataFrame(rows)
    # FDR correction
    rej, qvals, _, _ = multipletests(gdf["p"].values, alpha=0.05, method="fdr_bh")
    gdf["q"] = qvals
    gdf.to_csv(out_csv_prefix + ".csv", index=False)
    return gdf


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--input", default="output/merged_genie_with_mutations.csv")
    p.add_argument("--outdir", default="output/risk_models")
    p.add_argument("--top_genes", type=int, default=50)
    args = p.parse_args()

    ensure_dir(args.outdir)
    df = pd.read_csv(args.input, dtype=str)
    # normalize numeric columns
    for c in ["OS_MONTHS", "PFS_MONTHS"]:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")

    # default covariates (subset that exist in data)
    candidate_covs = [
        "AGE_CAT",
        "SEX",
        "RACE",
        "ETHNICITY",
        "SUBTYPE",
        "BUFFA_HYPOXIA_SCORE",
        "HISTORY_NEOADJUVANT_TRTYN",
        "AJCC_PATHOLOGIC_TUMOR_STAGE",
        "TMB_NONSYNONYMOUS",
    ]
    covs = [c for c in candidate_covs if c in df.columns]
    print("Using covariates:", covs)

    # --- OS analyses ---
    if "OS_MONTHS" in df.columns and "OS_STATUS" in df.columns:
        print("Running OS KM and CoxPH")
        # KM by SUBTYPE
        if "SUBTYPE" in df.columns:
            fit_km_by_subtype(df, "OS_MONTHS", "OS_STATUS", "SUBTYPE",
                              os.path.join(args.outdir, "OS_KM_by_SUBTYPE.png"),
                              os.path.join(args.outdir, "OS_logrank_pairs.csv"))
        # CoxPH
        os_cox = fit_cox(df, "OS_MONTHS", "OS_STATUS", covs, os.path.join(args.outdir, "OS_cox_summary.csv"))
    else:
        print("OS columns missing; skipping OS analyses")

    # --- PFS analyses ---
    if "PFS_MONTHS" in df.columns and "PFS_STATUS" in df.columns:
        print("Running PFS KM and CoxPH")
        if "SUBTYPE" in df.columns:
            fit_km_by_subtype(df, "PFS_MONTHS", "PFS_STATUS", "SUBTYPE",
                              os.path.join(args.outdir, "PFS_KM_by_SUBTYPE.png"),
                              os.path.join(args.outdir, "PFS_logrank_pairs.csv"))
        pfs_cox = fit_cox(df, "PFS_MONTHS", "PFS_STATUS", covs, os.path.join(args.outdir, "PFS_cox_summary.csv"))
    else:
        print("PFS columns missing; skipping PFS analyses")

    # Gene-level feature selection using Hugo_Symbol_list (if present)
    if "Hugo_Symbol_list" in df.columns:
        print("Building gene matrix for top genes")
        gene_mat, counts = build_gene_matrix(df, "Hugo_Symbol_list", top_k=args.top_genes)
        # Cox PH gene tests for OS and PFS
        if "OS_MONTHS" in df.columns and "OS_STATUS" in df.columns:
            g_os = gene_feature_tests(df, "OS_MONTHS", "OS_STATUS", covs, gene_mat, os.path.join(args.outdir, "genes_OS"))
        else:
            g_os = None
        if "PFS_MONTHS" in df.columns and "PFS_STATUS" in df.columns:
            g_pfs = gene_feature_tests(df, "PFS_MONTHS", "PFS_STATUS", covs, gene_mat, os.path.join(args.outdir, "genes_PFS"))
        else:
            g_pfs = None
        # combine selections
        rows = []
        genes = set()
        if g_os is not None:
            genes.update(g_os["gene"].tolist())
        if g_pfs is not None:
            genes.update(g_pfs["gene"].tolist())
        for g in sorted(genes):
            r = {"gene": g}
            if g_os is not None and g in g_os["gene"].values:
                tmp = g_os[g_os["gene"] == g].iloc[0]
                r.update({"os_hr": tmp["hr"], "os_p": tmp["p"], "os_q": tmp.get("q", np.nan)})
            if g_pfs is not None and g in g_pfs["gene"].values:
                tmp = g_pfs[g_pfs["gene"] == g].iloc[0]
                r.update({"pfs_hr": tmp["hr"], "pfs_p": tmp["p"], "pfs_q": tmp.get("q", np.nan)})
            rows.append(r)
        if rows:
            pd.DataFrame(rows).to_csv(os.path.join(args.outdir, "gene_selection_summary.csv"), index=False)

    print("Risk modeling complete. Results saved to", args.outdir)


if __name__ == "__main__":
    main()
