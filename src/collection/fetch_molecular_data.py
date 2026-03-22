#!/usr/bin/env python3
"""
fetch_molecular_data.py
-----------------------
cBioPortal REST API client for retrieving molecular/genomic data.

Wraps the cBioPortal public API (https://www.cbioportal.org/api) to fetch
molecular profiles, mutation data, and clinical attributes. Returns results
as pandas DataFrames for downstream analysis.

Usage (CLI):
    python src/collection/fetch_molecular_data.py \
        --genes 672,675 \
        --profiles brca_tcga_mrna \
        --outdir data/raw

Usage (library):
    from src.collection.fetch_molecular_data import fetch_molecular_data
    df = fetch_molecular_data(
        entrez_gene_ids=["672", "675"],
        molecular_profile_ids=["brca_tcga_mrna"],
    )

Notes:
    - The public cBioPortal API does not require an API key.
    - For institutional instances, set CBIOPORTAL_URL in your .env file.
    - Rate-limit-friendly: adds small delays between paginated requests.
"""

import argparse
import os
import sys
import time
from typing import List, Optional

import pandas as pd
import requests

# ── Configuration ────────────────────────────────────────────────────────────

BASE_URL = os.environ.get("CBIOPORTAL_URL", "https://www.cbioportal.org/api")
REQUEST_TIMEOUT = 30  # seconds
RETRY_DELAY = 1.0     # seconds between retries


# ── Helpers ──────────────────────────────────────────────────────────────────

def _get(endpoint: str, params: Optional[dict] = None) -> requests.Response:
    """GET request with retry."""
    url = f"{BASE_URL}{endpoint}"
    for attempt in range(3):
        try:
            resp = requests.get(url, params=params, timeout=REQUEST_TIMEOUT)
            resp.raise_for_status()
            return resp
        except requests.exceptions.RequestException as e:
            if attempt == 2:
                raise
            time.sleep(RETRY_DELAY * (attempt + 1))
    raise RuntimeError("Unreachable")


def _post(endpoint: str, params: Optional[dict] = None,
          json_body: Optional[dict] = None) -> requests.Response:
    """POST request with retry."""
    url = f"{BASE_URL}{endpoint}"
    for attempt in range(3):
        try:
            resp = requests.post(url, params=params, json=json_body,
                                 timeout=REQUEST_TIMEOUT)
            resp.raise_for_status()
            return resp
        except requests.exceptions.RequestException as e:
            if attempt == 2:
                raise
            time.sleep(RETRY_DELAY * (attempt + 1))
    raise RuntimeError("Unreachable")


def process_response(response: requests.Response, error_msg: str) -> pd.DataFrame:
    """Convert JSON response to DataFrame or raise with context."""
    if response.status_code != 200:
        raise RuntimeError(f"{error_msg} (HTTP {response.status_code}): {response.text[:500]}")
    data = response.json()
    if not data:
        return pd.DataFrame()
    return pd.DataFrame(data)


# ── Public API Functions ─────────────────────────────────────────────────────

def fetch_molecular_data(
    entrez_gene_ids: Optional[List[str]] = None,
    molecular_profile_ids: Optional[List[str]] = None,
    sample_molecular_identifiers: Optional[List[dict]] = None,
    projection: str = "SUMMARY",
) -> pd.DataFrame:
    """
    Fetch molecular data from cBioPortal.

    :param entrez_gene_ids: List of Entrez Gene IDs (e.g., ["672", "675"]).
    :param molecular_profile_ids: List of MolecularProfile IDs
        (e.g., ["brca_tcga_mrna", "acc_tcga_rna_seq_v2_mrna"]).
    :param sample_molecular_identifiers: List of dicts, each with
        "molecular_profile_id" and "sample_ids" keys. Example::

            [{"molecular_profile_id": "brca_tcga_mrna",
              "sample_ids": ["TCGA-AR-A1AR-01", "TCGA-BH-A1EO-01"]}]

    :param projection: Detail level — "DETAILED", "ID", "META", or "SUMMARY".
    :returns: DataFrame of molecular data.
    """
    endpoint = "/molecular-data/fetch"
    params = {"projection": projection}
    body = {}

    if entrez_gene_ids:
        body["entrezGeneIds"] = [int(g) for g in entrez_gene_ids]

    if molecular_profile_ids:
        body["molecularProfileIds"] = molecular_profile_ids

    if sample_molecular_identifiers:
        body["sampleMolecularIdentifiers"] = []
        for item in sample_molecular_identifiers:
            profile_id = item["molecular_profile_id"]
            for sid in item["sample_ids"]:
                body["sampleMolecularIdentifiers"].append({
                    "molecularProfileId": profile_id,
                    "sampleId": sid,
                })

    resp = _post(endpoint, params=params, json_body=body)
    return process_response(resp, "Failed to fetch molecular data.")


def get_all_molecular_data_in_molecular_profile(
    molecular_profile_id: str,
    sample_list_id: str,
    entrez_gene_id: str,
    projection: str = "SUMMARY",
) -> pd.DataFrame:
    """
    Get all molecular data in a molecular profile for a specific gene.

    :param molecular_profile_id: e.g. "acc_tcga_rna_seq_v2_mrna"
    :param sample_list_id: e.g. "acc_tcga_all"
    :param entrez_gene_id: e.g. "1"
    :param projection: Detail level.
    :returns: DataFrame of molecular data for the specified gene.
    """
    endpoint = f"/molecular-profiles/{molecular_profile_id}/molecular-data"
    params = {
        "entrezGeneId": entrez_gene_id,
        "projection": projection,
        "sampleListId": sample_list_id,
    }
    resp = _get(endpoint, params=params)
    return process_response(resp, "Failed to get molecular data in molecular profile.")


def fetch_all_molecular_data_in_molecular_profile(
    molecular_profile_id: str,
    entrez_gene_ids: Optional[List[str]] = None,
    sample_ids: Optional[List[str]] = None,
    sample_list_id: Optional[str] = None,
    projection: str = "SUMMARY",
) -> pd.DataFrame:
    """
    Fetch molecular data in a molecular profile for a list of genes.

    :param molecular_profile_id: e.g. "acc_tcga_rna_seq_v2_mrna"
    :param entrez_gene_ids: List of Entrez Gene IDs (e.g., ["672", "675"]).
    :param sample_ids: List of Sample IDs (e.g., ["TCGA-AR-A1AR-01"]).
    :param sample_list_id: Sample List ID (e.g., "brca_tcga_all").
    :param projection: Detail level.
    :returns: DataFrame of molecular data for the specified genes.
    """
    endpoint = f"/molecular-profiles/{molecular_profile_id}/molecular-data/fetch"
    params = {"projection": projection}
    body = {}

    if entrez_gene_ids:
        body["entrezGeneIds"] = [int(g) for g in entrez_gene_ids]

    if sample_ids:
        body["sampleIds"] = sample_ids

    if sample_list_id:
        body["sampleListId"] = sample_list_id

    resp = _post(endpoint, params=params, json_body=body)
    return process_response(
        resp, "Failed to fetch molecular data in molecular profile."
    )


def get_mutations(
    molecular_profile_id: str,
    sample_list_id: Optional[str] = None,
    sample_ids: Optional[List[str]] = None,
    entrez_gene_ids: Optional[List[str]] = None,
    projection: str = "DETAILED",
) -> pd.DataFrame:
    """
    Fetch mutation data for a study's mutation profile.

    :param molecular_profile_id: e.g. "brca_tcga_mutations"
    :param sample_list_id: e.g. "brca_tcga_all"
    :param sample_ids: Explicit sample IDs (alternative to sample_list_id).
    :param entrez_gene_ids: Filter to specific genes.
    :param projection: Detail level.
    :returns: DataFrame of mutations.
    """
    endpoint = f"/molecular-profiles/{molecular_profile_id}/mutations/fetch"
    params = {"projection": projection}
    body = {}

    if sample_list_id:
        body["sampleListId"] = sample_list_id
    if sample_ids:
        body["sampleIds"] = sample_ids
    if entrez_gene_ids:
        body["entrezGeneIds"] = [int(g) for g in entrez_gene_ids]

    resp = _post(endpoint, params=params, json_body=body)
    return process_response(resp, "Failed to fetch mutations.")


def get_clinical_data(
    study_id: str,
    attribute_ids: Optional[List[str]] = None,
    clinical_data_type: str = "PATIENT",
    projection: str = "SUMMARY",
) -> pd.DataFrame:
    """
    Fetch clinical data for a study.

    :param study_id: e.g. "brca_tcga"
    :param attribute_ids: Specific clinical attribute IDs to fetch.
    :param clinical_data_type: "PATIENT" or "SAMPLE".
    :param projection: Detail level.
    :returns: DataFrame of clinical data.
    """
    endpoint = f"/studies/{study_id}/clinical-data"
    params = {
        "clinicalDataType": clinical_data_type,
        "projection": projection,
    }
    if attribute_ids:
        params["attributeId"] = attribute_ids

    resp = _get(endpoint, params=params)
    return process_response(resp, "Failed to fetch clinical data.")


def list_molecular_profiles(study_id: str) -> pd.DataFrame:
    """
    List all molecular profiles available for a study.

    :param study_id: e.g. "brca_tcga"
    :returns: DataFrame with molecularProfileId, name, description, etc.
    """
    endpoint = f"/studies/{study_id}/molecular-profiles"
    resp = _get(endpoint)
    return process_response(resp, "Failed to list molecular profiles.")


def list_sample_lists(study_id: str) -> pd.DataFrame:
    """
    List all sample lists available for a study.

    :param study_id: e.g. "brca_tcga"
    :returns: DataFrame with sampleListId, name, description, sampleCount.
    """
    endpoint = f"/studies/{study_id}/sample-lists"
    resp = _get(endpoint)
    return process_response(resp, "Failed to list sample lists.")


# ── CLI ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Fetch molecular data from cBioPortal API"
    )
    parser.add_argument("--genes", type=str, default=None,
                        help="Comma-separated Entrez Gene IDs (e.g. 672,675)")
    parser.add_argument("--profiles", type=str, default=None,
                        help="Comma-separated Molecular Profile IDs")
    parser.add_argument("--study", type=str, default=None,
                        help="Study ID (e.g. brca_tcga) — lists profiles if no genes given")
    parser.add_argument("--sample-list", type=str, default=None,
                        help="Sample List ID (e.g. brca_tcga_all)")
    parser.add_argument("--mutations", action="store_true",
                        help="Fetch mutations instead of expression data")
    parser.add_argument("--clinical", action="store_true",
                        help="Fetch clinical data for --study")
    parser.add_argument("--outdir", type=str, default="data/raw",
                        help="Output directory for saved CSVs")
    parser.add_argument("--projection", type=str, default="SUMMARY",
                        choices=["SUMMARY", "DETAILED", "ID", "META"])
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # List profiles for a study
    if args.study and not args.genes and not args.clinical and not args.mutations:
        print(f"Listing molecular profiles for study: {args.study}")
        df = list_molecular_profiles(args.study)
        out = os.path.join(args.outdir, f"{args.study}_profiles.csv")
        df.to_csv(out, index=False)
        print(df[["molecularProfileId", "name", "molecularAlterationType"]].to_string(index=False))
        print(f"\nSaved → {out}")
        return

    # Clinical data
    if args.clinical and args.study:
        print(f"Fetching clinical data for study: {args.study}")
        df = get_clinical_data(args.study, projection=args.projection)
        out = os.path.join(args.outdir, f"{args.study}_clinical.csv")
        df.to_csv(out, index=False)
        print(f"Fetched {len(df)} records → {out}")
        return

    # Molecular / mutation data
    gene_ids = args.genes.split(",") if args.genes else None
    profile_ids = args.profiles.split(",") if args.profiles else None

    if args.mutations and profile_ids:
        print(f"Fetching mutations from {profile_ids[0]}")
        df = get_mutations(
            molecular_profile_id=profile_ids[0],
            sample_list_id=args.sample_list,
            entrez_gene_ids=gene_ids,
            projection=args.projection,
        )
        out = os.path.join(args.outdir, f"{profile_ids[0]}_mutations.csv")
    elif profile_ids and gene_ids:
        if len(profile_ids) == 1:
            print(f"Fetching molecular data from {profile_ids[0]} for genes: {gene_ids}")
            df = fetch_all_molecular_data_in_molecular_profile(
                molecular_profile_id=profile_ids[0],
                entrez_gene_ids=gene_ids,
                sample_list_id=args.sample_list,
                projection=args.projection,
            )
        else:
            print(f"Fetching molecular data across profiles: {profile_ids}")
            df = fetch_molecular_data(
                entrez_gene_ids=gene_ids,
                molecular_profile_ids=profile_ids,
                projection=args.projection,
            )
        out = os.path.join(args.outdir, "molecular_data.csv")
    else:
        parser.print_help()
        sys.exit(1)

    df.to_csv(out, index=False)
    print(f"Fetched {len(df)} records → {out}")


if __name__ == "__main__":
    main()
