# ================================================================
# extract_items_from_dataset.py — Reusable extractor/validator
# ================================================================

import pandas as pd
import numpy as np
import re
import os
from datetime import datetime
import argparse
from pathlib import Path


def main():
    # Use absolute paths based on the current file location
    base_dir = os.path.dirname(os.path.abspath(__file__))
    
    # Try CSV first (faster), fallback to Excel if needed
    csv_path = os.path.join(base_dir, "datasets_analysis_dictionary", "genie_bpc_metastatic_dataset_latest.csv")
    excel_path = r"C:\Users\jamesr4\OneDrive - Memorial Sloan Kettering Cancer Center\Documents\Research\Projects\genomics_brain_mets_genie_bpc\datasets_analysis_dictionary\genie_bpc_metastatic_dataset_latest.xlsx"
    
    if os.path.exists(csv_path):
        print(f"Using CSV file for faster loading: {csv_path}")
        data_path = csv_path
    else:
        print(f"CSV not found, using Excel file: {excel_path}")
        data_path = excel_path
        
    dict_path = os.path.join(base_dir, "datasets_analysis_dictionary", "genie_data_dic.xlsx")
    out_dir = os.path.join(base_dir, "datasets_analysis_dictionary", "output")
    vars_list = [
        "Age at Primary Diagnosis", "Center", "Primary Race", "Cancer Type", "Oncotree Code", "Sex", 
        "Cancer Type Detailed", "Age at Metastatic Diagnosis", "AJCC Stage", "Overall", "Grade", 
        "HER2 Status", "Histology\nEthnicity Category\nHormone Receptor Status", 
        "Discordance receptor status with primary", "Number of Samples Per Patient", "Sample Type", 
        "Sites of Distant Metastasis at the Time of Cancer Diagnosis (Stage IV Patients)", 
        "Sequencing Method", "Age at sequencing", "Mutation Count", "Fraction Genome Altered", 
        "Was a biopsy performed of metastatic site?", "Met Site: Bone \nOnly", "Met Site: Bone", 
        "Met Site: Brain", "Met Site: Lung", "Met Site: Lymph Node", "Met Site: Liver", 
        "Met Site: Multiple Sites", "Met Site: Soft Tissue", "Met Site: Visceral", "Met Site: Other Site", 
        "AKT Inhibitor Overall", "AKT Mutation Status", "Aromatase Inhibitor Overall", 
        "CDK4/6 Inhibitor Overall", "CDK4/6 Censoring Status (any line)", "Chemotherapy Overall", 
        "Chemotherapy received in LRR Treatment?", "Chemotherapy received in Primary Disease Treatment", 
        "C\nChemotherapy in first-line treatment?", "Chemotherapy in second-line treatment?", 
        "Chemotherapy lines received in Metastatic Treatment", "Overall Endocrine Therapy Sensitivity", 
        "Endocrine Therapy in first-line treatment?", "Endocrine Therapy in second-line treatment?", 
        "First-line treatment Censoring Status\nOverall Fulvestrant Censoring status (any line)", 
        "Fulvestrant \nHas this Patient experienced a Loco-Regional Recurrence (LRR)?", 
        "Total number of therapies received in metastatic disease treatment", 
        "Interval between sequencing and metastatis daignosis", "mTOR Inhibitor Overall", 
        "Overall mTOR Censoring Status (any line) (Months)", "Tamoxifen Overall", 
        "Therapeutic clinical trial overall", "Ovarian function suppression Overall", 
        "Overall Survival (Months)", "Overall Survival Status"
    ]
    
    process_dataset(
        data_path=data_path,
        dict_path=dict_path,
        out_dir=out_dir,
        vars_list=vars_list
    )


def snake_case(text):
    """Convert text to snake_case"""
    text = str(text).strip().lower()
    text = re.sub(r'[^a-z0-9]+', '_', text)
    text = re.sub(r'^_|_$', '', text)
    return text


def coerce_to_class(x, target_class):
    """Coerce data to target class"""
    if target_class is None or pd.isna(target_class) or target_class == "":
        return {"vec": x, "note": "no target class"}
    
    tc = str(target_class).lower()
    
    if "date" in tc:
        if pd.api.types.is_datetime64_any_dtype(x):
            return {"vec": x, "note": "already Date"}
        if pd.api.types.is_numeric_dtype(x):
            # Convert Excel serial dates
            return {"vec": pd.to_datetime(x, origin='1899-12-30', unit='D'), "note": "serial -> Date"}
        try:
            v = pd.to_datetime(x, errors='coerce')
            return {"vec": v, "note": "parsed -> Date"}
        except:
            return {"vec": x, "note": "date parse failed"}
    
    if "binary" in tc:
        if pd.api.types.is_numeric_dtype(x):
            v = x.copy()
        else:
            v = x.astype(str).str.lower()
            v = v.map({
                '1': 1, 'yes': 1, 'y': 1, 'true': 1, 't': 1,
                '0': 0, 'no': 0, 'n': 0, 'false': 0, 'f': 0
            })
            v = pd.to_numeric(v, errors='coerce')
        return {"vec": v.astype('Int64'), "note": "coerced to 0/1"}
    
    if "logical" in tc:
        v = x.astype(str).str.lower()
        v = v.map({
            'true': True, 't': True, '1': True, 'yes': True, 'y': True,
            'false': False, 'f': False, '0': False, 'no': False, 'n': False
        })
        return {"vec": v, "note": "to logical"}
    
    if "numeric" in tc:
        return {"vec": pd.to_numeric(x, errors='coerce'), "note": "to numeric"}
    
    if "factor" in tc:
        return {"vec": x.astype(str), "note": "kept char (factor validated by levels)"}
    
    if "character" in tc:
        return {"vec": x.astype(str), "note": "to character"}
    
    return {"vec": x, "note": "no rule"}


def parse_allowed(txt):
    """Parse allowed values specification"""
    if pd.isna(txt) or txt == "":
        return {"kind": "any", "values": None, "range": None, "note": "none"}
    
    txt = str(txt).strip()
    
    if "mm/dd/yyyy" in txt.lower():
        return {"kind": "date", "values": None, "range": None, "note": "date"}
    
    if txt.lower().startswith("range:"):
        # Extract numeric range
        rng = re.sub(r'[^0-9eE+\.\-]', ' ', txt)
        nums = [float(x) for x in rng.split() if x and not pd.isna(pd.to_numeric(x, errors='coerce'))]
        if len(nums) >= 2:
            return {"kind": "range", "values": None, "range": [min(nums[:2]), max(nums[:2])], "note": "range"}
    
    if "," in txt:
        vals = [v.strip() for v in txt.split(",")]
        return {"kind": "set", "values": list(set(vals)), "range": None, "note": "set"}
    
    return {"kind": "any", "values": None, "range": None, "note": "none"}


def validate_column(vec, target_class, allowed_spec):
    """Validate column values against specifications"""
    n = len(vec)
    invalid_idx = np.zeros(n, dtype=bool)
    
    if target_class is not None and not pd.isna(target_class) and "binary" in str(target_class).lower():
        invalid_idx = ~(pd.isna(vec) | vec.isin([0, 1]))
    elif allowed_spec["kind"] == "set":
        invalid_idx = ~(pd.isna(vec) | vec.astype(str).str.strip().isin(allowed_spec["values"]))
    elif allowed_spec["kind"] == "range":
        vnum = pd.to_numeric(vec, errors='coerce')
        range_min, range_max = allowed_spec["range"]
        invalid_idx = ~(pd.isna(vnum) | ((vnum >= range_min) & (vnum <= range_max)))
    elif allowed_spec["kind"] == "date":
        try:
            vdate = pd.to_datetime(vec, errors='coerce')
            invalid_idx = pd.isna(vdate) & ~pd.isna(vec)
        except:
            invalid_idx = np.ones(n, dtype=bool)
    
    invalid_rows = np.where(invalid_idx)[0][:20].tolist()  # First 20 invalid rows
    
    return {
        "n": n,
        "n_na": pd.isna(vec).sum(),
        "n_invalid": invalid_idx.sum(),
        "invalid_rows": invalid_rows
    }


def suggest_source(var_name):
    """Suggest data source based on variable name"""
    v = str(var_name).lower()
    
    if any(term in v for term in ["impact", "mutation", "gene", "variant", "oncoprint", "germline", "somatic", "msi", "tmb", "maf"]):
        return "Genomic → DMP/MSK-IMPACT (cBioPortal/OncoKB export; DMP LIMS)."
    if any(term in v for term in ["stage", "tnm", "pt", "pn", "pm", "clinical_stag", "pathologic"]):
        return "Staging → DMT or synoptic pathology (tumor board notes)."
    if any(term in v for term in ["adi", "deprivation", "zipcode", "census"]):
        return "Socioeconomic → ADI linkage (Neighborhood Atlas)."
    if any(term in v for term in ["death", "vital", "dod"]):
        return "Vital status → EHR demographics or registry linkage."
    if any(term in v for term in ["recurrence", "progression", "mets", "pfs", "dfs"]):
        return "Disease status → DMT; radiology/onc notes; registry abstractions."
    if "oncotype" in v:
        return "Oncotype DX → pathology/molecular reports; vendor portal."
    if any(term in v for term in ["insurance", "medicaid", "medicare", "private", "uninsured"]):
        return "Payer → registration/billing (EHR)."
    if any(term in v for term in ["lvi", "grade", "subtype", "er", "pr", "her2"]):
        return "Pathology → synoptic reports / LIS abstraction."
    
    return "Check EHR/registry; variable not found in current dataset."


def process_dataset(data_path, dict_path, out_dir="out", vars_list=None):
    """Main processing function"""
    
    # Create output directory
    os.makedirs(out_dir, exist_ok=True)
    
    # Load data
    print(f"Loading dataset from: {data_path}")
    if data_path.lower().endswith('.csv'):
        dat = pd.read_csv(data_path)
    else:
        print("Reading Excel file... this may take a moment for large files")
        # Try different approaches for large Excel files
        try:
            # First try with xlsxreader engine which is faster for large files
            dat = pd.read_excel(data_path, engine='openpyxl')
        except Exception as e:
            print(f"Error with openpyxl engine: {e}")
            print("Trying alternative approach...")
            # If that fails, try reading in chunks or converting to CSV first
            try:
                dat = pd.read_excel(data_path, engine='calamine')
            except:
                print("Could not read with calamine, trying default...")
                dat = pd.read_excel(data_path)
    
    print(f"Dataset loaded: {len(dat)} rows, {len(dat.columns)} columns")
    print(f"First 5 column names: {dat.columns[:5].tolist()}")
    
    # Convert column names to snake_case
    print("Converting column names to snake_case...")
    dat.columns = [snake_case(col) for col in dat.columns]
    
    # Load data dictionary
    print(f"Loading data dictionary from: {dict_path}")
    try:
        dd = pd.read_excel(dict_path, sheet_name=0)
        print(f"Data dictionary loaded: {len(dd)} entries")
        if len(dd) > 0 and 'variable' in dd.columns:
            dd['variable'] = dd['variable'].apply(snake_case)
        else:
            print("Warning: Data dictionary is empty or missing 'variable' column. Proceeding without dictionary validation.")
            dd = pd.DataFrame(columns=['variable', 'class', 'allowed_values'])
    except Exception as e:
        print(f"Error loading data dictionary: {e}")
        print("Proceeding without data dictionary validation.")
        dd = pd.DataFrame(columns=['variable', 'class', 'allowed_values'])
    
    # Process variable list
    if vars_list is None or len(vars_list) == 0:
        vars_list = dd['variable'].tolist()
    else:
        vars_list = [snake_case(v) for v in vars_list]
    
    # Find present and missing variables
    present = [v for v in vars_list if v in dat.columns]
    missing = [v for v in vars_list if v not in dat.columns]
    
    # Create dictionary mapping
    dd_map = dd[['variable', 'class', 'allowed_values']].drop_duplicates()
    dd_map['allowed_parsed'] = dd_map['allowed_values'].apply(parse_allowed)
    
    # Process each present variable
    logs = []
    for v in present:
        spec = dd_map[dd_map['variable'] == v]
        if not spec.empty:
            spec = spec.iloc[0]
            target_class = spec.get('class', None)
            allowed_spec = spec['allowed_parsed']
        else:
            target_class = None
            allowed_spec = {"kind": "any", "values": None, "range": None, "note": "none"}
        
        # Coerce data type
        coerced = coerce_to_class(dat[v], target_class)
        dat[v] = coerced["vec"]
        
        # Validate
        val = validate_column(dat[v], target_class, allowed_spec)
        val.update({
            "variable": v,
            "class": target_class,
            "coercion_note": coerced["note"],
            "allowed_note": allowed_spec["note"]
        })
        logs.append(val)
    
    # Create QA report
    qa_df = pd.DataFrame(logs)
    qa_df = qa_df[['variable', 'class', 'coercion_note', 'allowed_note', 'n', 'n_na', 'n_invalid', 'invalid_rows']]
    qa_df = qa_df.sort_values(['n_invalid', 'variable'], ascending=[False, True])
    
    # Create missing variables report
    missing_df = pd.DataFrame({
        'variable': missing,
        'suggestion': [suggest_source(v) for v in missing]
    })
    
    # Write outputs
    subset_path = os.path.join(out_dir, "extracted_subset.csv")
    qa_path = os.path.join(out_dir, "qa_report.csv")
    miss_path = os.path.join(out_dir, "missing_variables_suggestions.csv")
    
    dat[present].to_csv(subset_path, index=False, na_rep="")
    qa_df.to_csv(qa_path, index=False, na_rep="")
    missing_df.to_csv(miss_path, index=False, na_rep="")
    
    print(f"Wrote: {subset_path}")
    print(f"Wrote: {qa_path}")
    print(f"Wrote: {miss_path}")
    
    return {
        "subset": subset_path,
        "qa": qa_path,
        "missing": miss_path
    }


# CLI functionality
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract and validate variables from dataset")
    parser.add_argument("-d", "--data", type=str, help="Path to dataset (CSV/Excel)")
    parser.add_argument("-k", "--dict", type=str, help="Path to data dictionary (Excel)")
    parser.add_argument("-o", "--outdir", type=str, default="out", help="Output directory")
    parser.add_argument("-v", "--vars", type=str, default=None, 
                       help="Comma-separated variables to extract (snake_case). If omitted, uses all from dictionary.")
    
    args = parser.parse_args()
    
    # If no arguments provided, use the hardcoded defaults from main()
    if args.data is None and args.dict is None:
        print("No arguments provided, using default paths...")
        main()
    else:
        # Require both data and dict if either is provided
        if args.data is None or args.dict is None:
            parser.error("Both --data and --dict are required when using command line arguments")
        
        vars_vec = None
        if args.vars:
            vars_vec = [v.strip() for v in args.vars.split(",")]
        
        process_dataset(args.data, args.dict, out_dir=args.outdir, vars_list=vars_vec)

