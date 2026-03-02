import pandas as pd
import numpy as np

def import_data("C:\Users\jamesr4\OneDrive - Memorial Sloan Kettering Cancer Center\Documents\Research\Projects\genomics_brain_mets_genie_bpc\datasets_analysis_dictionary\combined_patient_sample_hypoxia_data.xlsx"):
    # Load the Excel file
    df = pd.read_excel("C:\Users\jamesr4\OneDrive - Memorial Sloan Kettering Cancer Center\Documents\Research\Projects\genomics_brain_mets_genie_bpc\datasets_analysis_dictionary\combined_patient_sample_hypoxia_data.xlsx", sheet_name=c("data_clinical_patient", "data_clinical_supp_hypoxia", "data_clinical_sample")

    # Display the first few rows of the dataframe
    print(df.head())

    return df