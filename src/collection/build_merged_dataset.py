#!/usr/bin/env python3
"""
build_merged_dataset.py
-----------------------
Assemble clinical + mutation data fetched from cBioPortal into a single
wide-format dataset that matches the schema expected by the analysis pipeline.

Usage:
    python src/collection/build_merged_dataset.py \
        --clinical data/raw/brca_tcga_clinical.csv \
        --mutations data/raw/brca_tcga_mutations_mutations.csv \
        --outdir datasets_analysis_dictionary
"""

import argparse
import os

import numpy as np
import pandas as pd


def pivot_clinical(clin_path: str) -> pd.DataFrame:
    """Pivot long-format clinical data to one row per patient."""
    clin = pd.read_csv(clin_path)
    # Pivot: patientId × clinicalAttributeId → value
    wide = clin.pivot_table(
        index="patientId", columns="clinicalAttributeId",
        values="value", aggfunc="first"
    ).reset_index()
    wide.columns.name = None
    return wide


def build_gene_indicators(mut_path: str) -> pd.DataFrame:
    """Build per-patient binary gene indicator matrix from mutation data."""
    mut = pd.read_csv(mut_path)

    # Get patient-gene pairs
    if "patientId" not in mut.columns and "sampleId" in mut.columns:
        # Extract patient ID from sample ID (TCGA convention: first 12 chars)
        mut["patientId"] = mut["sampleId"].str[:12]

    gene_col = "hugoGeneSymbol" if "hugoGeneSymbol" in mut.columns else "hugo_symbol"
    if gene_col not in mut.columns:
        # Try to find any gene column
        for c in mut.columns:
            if "gene" in c.lower() or "hugo" in c.lower():
                gene_col = c
                break

    patient_genes = mut[["patientId", gene_col]].drop_duplicates()
    patient_genes["value"] = 1

    # Pivot to wide: one column per gene
    gene_matrix = patient_genes.pivot_table(
        index="patientId", columns=gene_col,
        values="value", fill_value=0
    ).reset_index()
    gene_matrix.columns.name = None

    # Rename gene columns to G__<gene>
    rename = {c: f"G__{c}" for c in gene_matrix.columns if c != "patientId"}
    gene_matrix = gene_matrix.rename(columns=rename)

    return gene_matrix


def build_mutation_metadata(mut_path: str) -> pd.DataFrame:
    """Extract per-patient mutation count and hugo_symbol list."""
    mut = pd.read_csv(mut_path)

    if "patientId" not in mut.columns and "sampleId" in mut.columns:
        mut["patientId"] = mut["sampleId"].str[:12]

    gene_col = "hugoGeneSymbol" if "hugoGeneSymbol" in mut.columns else "hugo_symbol"
    for c in mut.columns:
        if "gene" in c.lower() or "hugo" in c.lower():
            gene_col = c
            break

    meta = mut.groupby("patientId").agg(
        hugo_symbol=(gene_col, lambda x: ";".join(sorted(x.astype(str).unique()))),
        mutation_count=(gene_col, "nunique"),
    ).reset_index()

    # Variant type if available
    vtype_col = None
    for c in ["variantType", "variant_type", "mutationType"]:
        if c in mut.columns:
            vtype_col = c
            break
    if vtype_col:
        vtype = mut.groupby("patientId")[vtype_col].first().reset_index()
        vtype.columns = ["patientId", "variant_type"]
        meta = meta.merge(vtype, on="patientId", how="left")

    return meta


def main():
    parser = argparse.ArgumentParser(description="Build merged analysis dataset")
    parser.add_argument("--clinical", required=True, help="Path to clinical CSV")
    parser.add_argument("--mutations", required=True, help="Path to mutations CSV")
    parser.add_argument("--outdir", default="datasets_analysis_dictionary")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # Pivot clinical
    print("Pivoting clinical data...")
    df = pivot_clinical(args.clinical)
    print(f"  {len(df)} patients, {len(df.columns)} attributes")

    # Numeric conversions
    for col in ["AGE", "DFS_MONTHS", "OS_MONTHS", "PFS_MONTHS"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    # DFS_STATUS → binary (1 = recurrence, 0 = disease free)
    if "DFS_STATUS" in df.columns:
        df["DFS_STATUS_LABEL"] = df["DFS_STATUS"]
        df["DFS_STATUS"] = df["DFS_STATUS"].map(
            lambda x: 1 if "recur" in str(x).lower() or "progress" in str(x).lower()
            else (0 if "disease free" in str(x).lower() or "diseasefree" in str(x).lower()
                  else np.nan)
        )

    # OS_STATUS → binary
    if "OS_STATUS" in df.columns:
        df["OS_STATUS_LABEL"] = df["OS_STATUS"]
        df["OS_STATUS"] = df["OS_STATUS"].map(
            lambda x: 1 if "deceased" in str(x).lower() or "dead" in str(x).lower()
            else (0 if "living" in str(x).lower() or "alive" in str(x).lower()
                  else np.nan)
        )

    # Build gene indicators
    print("Building gene indicators...")
    gene_df = build_gene_indicators(args.mutations)
    print(f"  {len(gene_df)} patients with mutations, {len(gene_df.columns) - 1} genes")

    # Mutation metadata
    mut_meta = build_mutation_metadata(args.mutations)

    # Merge
    print("Merging...")
    merged = df.merge(gene_df, on="patientId", how="left")
    merged = merged.merge(mut_meta, on="patientId", how="left")

    # Fill gene indicator NaN with 0 (no mutation detected)
    gene_cols = [c for c in merged.columns if c.startswith("G__")]
    merged[gene_cols] = merged[gene_cols].fillna(0).astype(int)
    merged["mutation_count"] = merged["mutation_count"].fillna(0).astype(int)

    # Top-5 gene flag
    if len(gene_cols) >= 5:
        top5 = gene_cols[:5]
        merged["top_5_gene_status"] = merged[top5].max(axis=1)
    elif gene_cols:
        merged["top_5_gene_status"] = merged[gene_cols].max(axis=1)
    else:
        merged["top_5_gene_status"] = 0

    # SUBTYPE from clinical if available
    subtype_col = None
    for c in ["SUBTYPE", "CANCER_TYPE_DETAILED", "HISTOLOGICAL_DIAGNOSIS"]:
        if c in merged.columns:
            subtype_col = c
            break
    if subtype_col and subtype_col != "SUBTYPE":
        merged["SUBTYPE"] = merged[subtype_col]

    # Save
    out_xlsx = os.path.join(args.outdir, "merged_genie.xlsx")
    out_csv = os.path.join(args.outdir, "merged_genie.csv")
    merged.to_excel(out_xlsx, index=False)
    merged.to_csv(out_csv, index=False)

    print(f"\nMerged dataset: {merged.shape[0]} rows × {merged.shape[1]} columns")
    print(f"Gene indicator columns: {len(gene_cols)}")
    print(f"Saved → {out_xlsx}")
    print(f"Saved → {out_csv}")

    # Summary
    print("\n── Column summary ──")
    for col in ["AGE", "DFS_MONTHS", "DFS_STATUS", "OS_MONTHS", "OS_STATUS",
                "SUBTYPE", "AJCC_PATHOLOGIC_TUMOR_STAGE",
                "hugo_symbol", "mutation_count", "top_5_gene_status"] + gene_cols[:5]:
        if col in merged.columns:
            nn = merged[col].notna().sum()
            print(f"  {col}: {nn}/{len(merged)} non-null")


if __name__ == "__main__":
    main()
