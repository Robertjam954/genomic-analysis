#!/usr/bin/env python3
import argparse
import os
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


FALLBACK_ENTREZ_MAP = {
    5290: {
        "hugo_symbol": "PIK3CA",
        "gene_name": "phosphatidylinositol-4,5-bisphosphate 3-kinase catalytic subunit alpha",
        "gene_function": "PI3K pathway kinase / oncogenic signaling",
        "locus_type": "gene with protein product",
    },
    7157: {
        "hugo_symbol": "TP53",
        "gene_name": "tumor protein p53",
        "gene_function": "tumor suppressor / DNA damage response",
        "locus_type": "gene with protein product",
    },
    999: {
        "hugo_symbol": "CDH1",
        "gene_name": "cadherin 1",
        "gene_function": "cell adhesion / epithelial integrity",
        "locus_type": "gene with protein product",
    },
    2064: {
        "hugo_symbol": "ERBB2",
        "gene_name": "erb-b2 receptor tyrosine kinase 2",
        "gene_function": "receptor tyrosine kinase signaling",
        "locus_type": "gene with protein product",
    },
    675: {
        "hugo_symbol": "BRCA2",
        "gene_name": "BRCA2 DNA repair associated",
        "gene_function": "homologous recombination DNA repair",
        "locus_type": "gene with protein product",
    },
    672: {
        "hugo_symbol": "BRCA1",
        "gene_name": "BRCA1 DNA repair associated",
        "gene_function": "homologous recombination DNA repair",
        "locus_type": "gene with protein product",
    },
    3845: {
        "hugo_symbol": "KRAS",
        "gene_name": "KRAS proto-oncogene, GTPase",
        "gene_function": "RAS/MAPK signaling",
        "locus_type": "gene with protein product",
    },
    1956: {
        "hugo_symbol": "EGFR",
        "gene_name": "epidermal growth factor receptor",
        "gene_function": "receptor tyrosine kinase signaling",
        "locus_type": "gene with protein product",
    },
    4893: {
        "hugo_symbol": "NRAS",
        "gene_name": "NRAS proto-oncogene, GTPase",
        "gene_function": "RAS/MAPK signaling",
        "locus_type": "gene with protein product",
    },
}


def ensure_dir(path: str) -> None:
    os.makedirs(path, exist_ok=True)


def to_num(s: pd.Series) -> pd.Series:
    return pd.to_numeric(s, errors="coerce")


def detect_patient_id_col(df: pd.DataFrame) -> str:
    candidates = [
        "patientId",
        "PATIENTID",
        "PATIENT_ID",
        "bcr_patient_barcode",
        "BCR_PATIENT_BARCODE",
    ]
    col_lut = {c.lower(): c for c in df.columns}
    for cand in candidates:
        if cand.lower() in col_lut:
            return col_lut[cand.lower()]
    return "_row_id"


def parse_entrez_from_gene_col(col: str):
    m = re.match(r"^G__([0-9]+)$", str(col))
    if not m:
        return None
    return int(m.group(1))


def load_hugo_map(hugo_xlsx: str) -> pd.DataFrame:
    hugo = pd.read_excel(hugo_xlsx, engine="openpyxl")
    req = ["symbol", "name", "gene_group", "locus_type", "entrez_id"]
    for c in req:
        if c not in hugo.columns:
            raise ValueError(f"Missing required HUGO column: {c}")

    hugo = hugo.copy()
    hugo["entrez_id"] = to_num(hugo["entrez_id"])
    hugo = hugo.dropna(subset=["entrez_id"]) 
    hugo["entrez_id"] = hugo["entrez_id"].astype(int)

    hugo = hugo.drop_duplicates(subset=["entrez_id"], keep="first")
    return hugo[["entrez_id", "symbol", "name", "gene_group", "locus_type"]]


def build_gene_mapping(df: pd.DataFrame, hugo_map: pd.DataFrame) -> pd.DataFrame:
    gene_cols = [c for c in df.columns if str(c).startswith("G__")]
    if not gene_cols:
        raise ValueError("No G__ gene indicator columns found in merged dataset.")

    mut_counts = df[gene_cols].apply(to_num).fillna(0).sum(axis=0).sort_values(ascending=False)
    n_patients = len(df)

    rows = []
    for col, n_mut in mut_counts.items():
        entrez = parse_entrez_from_gene_col(col)
        if entrez is None:
            continue

        match = hugo_map[hugo_map["entrez_id"] == entrez]
        if len(match) > 0:
            m = match.iloc[0]
            symbol = m["symbol"]
            gene_name = m["name"]
            gene_group = m["gene_group"]
            locus_type = m["locus_type"]
            mapping_source = "hugo_xlsx"
            function_param = gene_group if pd.notna(gene_group) and str(gene_group).strip() else locus_type
        elif entrez in FALLBACK_ENTREZ_MAP:
            fb = FALLBACK_ENTREZ_MAP[entrez]
            symbol = fb["hugo_symbol"]
            gene_name = fb["gene_name"]
            function_param = fb["gene_function"]
            locus_type = fb["locus_type"]
            mapping_source = "fallback_top_genes"
        else:
            symbol = np.nan
            gene_name = np.nan
            function_param = np.nan
            locus_type = np.nan
            mapping_source = "unmapped"

        rows.append({
            "gene_col": col,
            "entrez_id": entrez,
            "hugo_symbol": symbol,
            "gene_name": gene_name,
            "gene_function": function_param,
            "locus_type": locus_type,
            "mapping_source": mapping_source,
            "n_mutated": int(n_mut),
            "mutation_prevalence_pct": float((n_mut / n_patients) * 100.0),
        })

    out = pd.DataFrame(rows).sort_values("n_mutated", ascending=False).reset_index(drop=True)
    return out


def add_dfs_stats(df: pd.DataFrame, mapping_df: pd.DataFrame) -> pd.DataFrame:
    out = mapping_df.copy()

    if "DFS_STATUS" not in df.columns:
        out["dfs_available"] = False
        out["dfs_event_rate_mutated_pct"] = np.nan
        out["dfs_event_rate_non_mutated_pct"] = np.nan
        return out

    dfs = to_num(df["DFS_STATUS"]) 

    out_rows = []
    for _, r in out.iterrows():
        gcol = r["gene_col"]
        g = to_num(df[gcol]).fillna(0)

        mut_mask = g > 0
        non_mask = g <= 0

        mut_dfs = dfs[mut_mask]
        non_dfs = dfs[non_mask]

        mut_rate = (mut_dfs > 0).mean() * 100.0 if mut_dfs.notna().sum() > 0 else np.nan
        non_rate = (non_dfs > 0).mean() * 100.0 if non_dfs.notna().sum() > 0 else np.nan

        row = r.to_dict()
        row["dfs_available"] = True
        row["dfs_event_rate_mutated_pct"] = float(mut_rate) if pd.notna(mut_rate) else np.nan
        row["dfs_event_rate_non_mutated_pct"] = float(non_rate) if pd.notna(non_rate) else np.nan
        out_rows.append(row)

    return pd.DataFrame(out_rows)


def make_top5_plot(top5_df: pd.DataFrame, out_png: str):
    plt.figure(figsize=(10, 6))
    labels = [
        f"{sym} ({int(eid)})" if pd.notna(sym) else f"Entrez {int(eid)}"
        for sym, eid in zip(top5_df["hugo_symbol"], top5_df["entrez_id"])
    ]

    plt.bar(labels, top5_df["mutation_prevalence_pct"].values)
    plt.ylabel("Mutation prevalence (%)")
    plt.title("Top 5 Mutated Genes (mapped to HUGO)")
    plt.xticks(rotation=20, ha="right")
    plt.tight_layout()
    plt.savefig(out_png, dpi=250)
    plt.close()


def main():
    p = argparse.ArgumentParser(description="Map top mutated genes to HUGO + function and run EDA")
    p.add_argument("--merged-csv", default="datasets_analysis_dictionary/merged_genie.csv")
    p.add_argument("--hugo-xlsx", required=True)
    p.add_argument("--outdir", default="output/eda")
    args = p.parse_args()

    ensure_dir(args.outdir)

    df = pd.read_csv(args.merged_csv)
    pid_col = detect_patient_id_col(df)
    if pid_col == "_row_id":
        df[pid_col] = np.arange(len(df))

    hugo_map = load_hugo_map(args.hugo_xlsx)
    full_map = build_gene_mapping(df, hugo_map)
    full_map = add_dfs_stats(df, full_map)

    full_map_csv = os.path.join(args.outdir, "gene_mapping_with_function.csv")
    full_map.to_csv(full_map_csv, index=False)

    top5 = full_map.head(5).copy()
    top5_csv = os.path.join(args.outdir, "top5_mutated_genes_with_function.csv")
    top5.to_csv(top5_csv, index=False)

    long_rows = []
    for _, r in top5.iterrows():
        gcol = r["gene_col"]
        for i in range(len(df)):
            long_rows.append({
                "patient_id": df.iloc[i][pid_col],
                "gene_col": gcol,
                "entrez_id": r["entrez_id"],
                "hugo_symbol": r["hugo_symbol"],
                "gene_function": r["gene_function"],
                "mutated": int(to_num(pd.Series([df.iloc[i][gcol]])).fillna(0).iloc[0] > 0),
                "dfs_status": to_num(pd.Series([df.iloc[i]["DFS_STATUS"]])).iloc[0] if "DFS_STATUS" in df.columns else np.nan,
            })

    long_df = pd.DataFrame(long_rows)
    long_csv = os.path.join(args.outdir, "top5_gene_function_patient_long.csv")
    long_df.to_csv(long_csv, index=False)

    plot_png = os.path.join(args.outdir, "top5_mutated_genes_prevalence.png")
    make_top5_plot(top5, plot_png)

    print("Wrote", full_map_csv)
    print("Wrote", top5_csv)
    print("Wrote", long_csv)
    print("Wrote", plot_png)


if __name__ == "__main__":
    main()
