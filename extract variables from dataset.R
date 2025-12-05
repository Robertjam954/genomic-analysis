# ================================================================
# run_genie_extract.R — Reusable extractor/validator + GENIE run
# ================================================================

library(readr)
library(readxl)
library(dplyr)
library(stringr)
library(lubridate)
library(tidyr)
library(purrr)

# ---- helper operators & utilities ----

`%||%` <- function(a, b) if (is.null(a) || (length(a) == 1 && is.na(a))) b else a

snake <- function(x) {
  x |>
    str_trim() |>
    str_to_lower() |>
    str_replace_all("[^a-z0-9]+", "_") |>
    str_replace_all("^_|_$", "")
}

coerce_to_class <- function(x, target_class) {
  if (is.null(target_class) || is.na(target_class) || target_class == "") {
    return(list(vec = x, note = "no target class"))
  }
  tc <- tolower(target_class)
  
  if (grepl("date", tc)) {
    if (inherits(x, "Date")) return(list(vec = x, note = "already Date"))
    if (is.numeric(x)) return(list(vec = as.Date(x, origin = "1899-12-30"),
                                   note = "serial -> Date"))
    v <- suppressWarnings(mdy(x))
    if (all(is.na(v))) v <- suppressWarnings(ymd(x))
    return(list(vec = as.Date(v), note = "parsed mdy/ymd -> Date"))
  }
  
  if (grepl("binary", tc)) {
    v <- x
    if (!is.numeric(v)) {
      v <- dplyr::case_when(
        str_to_lower(as.character(x)) %in% c("1","yes","y","true","t") ~ 1,
        str_to_lower(as.character(x)) %in% c("0","no","n","false","f") ~ 0,
        TRUE ~ suppressWarnings(as.numeric(as.character(x)))
      )
    }
    return(list(vec = as.integer(v), note = "coerced to 0/1"))
  }
  
  if (grepl("logical", tc)) {
    v <- str_to_lower(as.character(x))
    v <- dplyr::case_when(
      v %in% c("true","t","1","yes","y") ~ TRUE,
      v %in% c("false","f","0","no","n") ~ FALSE,
      TRUE ~ NA
    )
    return(list(vec = as.logical(v), note = "to logical"))
  }
  
  if (grepl("numeric", tc)) {
    return(list(vec = suppressWarnings(as.numeric(as.character(x))),
                note = "to numeric"))
  }
  
  if (grepl("factor", tc)) {
    return(list(vec = as.character(x),
                note = "kept char (factor validated by levels)"))
  }
  
  if (grepl("character", tc)) {
    return(list(vec = as.character(x), note = "to character"))
  }
  
  list(vec = x, note = "no rule")
}

parse_allowed <- function(txt) {
  if (is.na(txt) || txt == "") {
    return(list(kind = "any", values = NULL, range = NULL, note = "none"))
  }
  
  txt <- str_trim(as.character(txt))
  
  if (grepl("^mm/dd/yyyy", str_to_lower(txt))) {
    return(list(kind = "date", values = NULL, range = NULL, note = "date"))
  }
  
  if (grepl("^range:", str_to_lower(txt))) {
    rng  <- str_replace_all(txt, "[^0-9eE+\\.-]", " ")
    nums <- suppressWarnings(as.numeric(unlist(str_split(rng, "\\s+"))))
    nums <- nums[!is.na(nums)]
    if (length(nums) >= 2) {
      return(list(kind = "range",
                  values = NULL,
                  range  = range(nums[1:2]),
                  note   = "range"))
    }
  }
  
  if (grepl(",", txt)) {
    vals <- str_split(txt, ",")[[1]] |> str_trim()
    return(list(kind = "set", values = unique(vals), range = NULL, note = "set"))
  }
  
  list(kind = "any", values = NULL, range = NULL, note = "none")
}

validate_column <- function(vec, target_class, allowed_spec) {
  n <- length(vec)
  invalid_idx <- rep(FALSE, n)
  
  if (!is.null(target_class) &&
      !is.na(target_class) &&
      grepl("binary", tolower(target_class))) {
    
    invalid_idx <- !(is.na(vec) | vec %in% c(0L, 1L))
    
  } else if (allowed_spec$kind == "set") {
    
    invalid_idx <- !(is.na(vec) |
                       str_trim(as.character(vec)) %in% allowed_spec$values)
    
  } else if (allowed_spec$kind == "range") {
    
    vnum <- suppressWarnings(as.numeric(vec))
    invalid_idx <- !(is.na(vnum) |
                       (vnum >= allowed_spec$range[1] &
                          vnum <= allowed_spec$range[2]))
    
  } else if (allowed_spec$kind == "date") {
    
    vdate  <- suppressWarnings(mdy(vec))
    vdate2 <- suppressWarnings(ymd(vec))
    ok <- !is.na(vdate) | !is.na(vdate2) | inherits(vec, "Date")
    invalid_idx <- !ok
  }
  
  tibble::tibble(
    n           = n,
    n_na        = sum(is.na(vec)),
    n_invalid   = sum(invalid_idx),
    invalid_rows = list(which(invalid_idx)[1:min(20, sum(invalid_idx))])
  )
}

suggest_source <- function(var_name) {
  v <- tolower(var_name)
  if (str_detect(v, "impact|mutation|gene|variant|oncoprint|germline|somatic|msi|tmb|maf"))
    return("Genomic → DMP/MSK-IMPACT (cBioPortal/OncoKB export; DMP LIMS).")
  if (str_detect(v, "stage|tnm|pt|pn|pm|clinical_stag|pathologic"))
    return("Staging → DMT or synoptic pathology (tumor board notes).")
  if (str_detect(v, "adi|deprivation|zipcode|census"))
    return("Socioeconomic → ADI linkage (Neighborhood Atlas).")
  if (str_detect(v, "death|vital|dod"))
    return("Vital status → EHR demographics or registry linkage.")
  if (str_detect(v, "recurrence|progression|mets|pfs|dfs"))
    return("Disease status → DMT; radiology/onc notes; registry abstractions.")
  if (str_detect(v, "oncotype"))
    return("Oncotype DX → pathology/molecular reports; vendor portal.")
  if (str_detect(v, "insurance|medicaid|medicare|private|uninsured"))
    return("Payer → registration/billing (EHR).")
  if (str_detect(v, "lvi|grade|subtype|er|pr|her2"))
    return("Pathology → synoptic reports / LIS abstraction.")
  "Check EHR/registry; variable not found in current dataset."
}

main <- function(data_path, dict_path, out_dir = "out", vars = NULL) {
  
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  if (grepl("\\.csv$", data_path, ignore.case = TRUE)) {
    dat <- suppressMessages(readr::read_csv(data_path, show_col_types = FALSE))
  } else {
    dat <- readxl::read_excel(data_path)
  }
  
  names(dat) <- snake(names(dat))
  
  dd <- readxl::read_excel(dict_path, sheet = 1) |>
    mutate(variable = snake(variable))
  
  if (is.null(vars) || length(vars) == 0) {
    vars_snake <- dd$variable
  } else {
    vars_snake <- snake(vars)
  }
  
  present <- intersect(vars_snake, names(dat))
  missing <- setdiff(vars_snake, names(dat))
  
  dd_map <- dd |>
    select(variable, class, allowed_values) |>
    distinct() |>
    mutate(allowed_parsed = purrr::map(allowed_values, parse_allowed))
  
  logs <- list()
  
  for (v in present) {
    spec <- dd_map |>
      filter(variable == v) |>
      slice(1)
    
    target_class <- spec$class %||% NA_character_
    allowed_spec <- spec$allowed_parsed[[1]] %||%
      list(kind = "any", values = NULL, range = NULL)
    
    coerced <- coerce_to_class(dat[[v]], target_class)
    dat[[v]] <- coerced$vec
    
    val <- validate_column(dat[[v]], target_class, allowed_spec)
    val$variable       <- v
    val$class          <- target_class
    val$coercion_note  <- coerced$note
    val$allowed_note   <- allowed_spec$note
    logs[[v]] <- val
  }
  
  qa <- bind_rows(logs) |>
    select(variable, class, coercion_note, allowed_note,
           n, n_na, n_invalid, invalid_rows) |>
    arrange(desc(n_invalid), variable)
  
  missing_df <- tibble::tibble(
    variable   = missing,
    suggestion = purrr::map_chr(missing, suggest_source)
  )
  
  subset_path <- file.path(out_dir, "extracted_subset_all_variables_of_interest_12_4.csv")
  qa_path     <- file.path(out_dir, "qa_report_all_variables_of_interest_12_4.csv")
  miss_path   <- file.path(out_dir, "missing_variables_suggestions_all_variables_of_interest_12_4.csv")
  
  readr::write_csv(dplyr::select(dat, dplyr::any_of(present)), subset_path, na = "")
  readr::write_csv(qa, qa_path, na = "")
  readr::write_csv(missing_df, miss_path, na = "")
  
  message("Wrote: ", subset_path)
  message("Wrote: ", qa_path)
  message("Wrote: ", miss_path)
  
  invisible(list(subset = subset_path,
                 qa      = qa_path,
                 missing = miss_path))
}

# ================================================================
# Run for your GENIE dataset
# ================================================================

genie_data_path <- "C:\\\\Users\\\\jamesr4\\\\OneDrive - Memorial Sloan Kettering Cancer Center\\\\Documents\\\\Research\\\\Projects\\\\genomics_brain_mets_genie_bpc\\\\datasets_analysis_dictionary\\\\clean_genie_bpc_metastatic_dataset_latest.csv"
genie_dict_path <- "C:\\\\Users\\\\jamesr4\\\\OneDrive - Memorial Sloan Kettering Cancer Center\\\\Documents\\\\Research\\\\Projects\\\\genomics_brain_mets_genie_bpc\\\\datasets_analysis_dictionary\\\\genie_data_dic.xlsx"
genie_out_dir   <- "C:\\\\Users\\\\jamesr4\\\\OneDrive - Memorial Sloan Kettering Cancer Center\\\\Documents\\\\Research\\\\Projects\\\\genomics_brain_mets_genie_bpc\\\\datasets_analysis_dictionary\\\\output"

genie_vars <- c(
  # Core identifiers
  "study_id",
  "patient_id", 
  "sample_id",
  
  # Demographics
  "center",
  "sex",
  "year_of_birth",
  "primary_race",
  "secondary_race",
  "tertiary_race",
  "ethnicity_category",
  "region_of_patient_origin",
  
  # Oncology history
  "age_at_diagnosis",
  "age_at_initial_diagnosis_years",
  "was_a_biopsy_performed_of_metastatic_site",
  "oncotree_code",
  "cancer_type",
  "cancer_type_detailed",
  "tumor_registry_icd_o_3_behavior_code",
  "age_at_primary_diagnosis",
  "additional_breast_cancer_diagnoses",
  "other_oncotree_diagnosis",
  "other_oncotree_diagnosis_1",
  "other_oncotree_diagnosis_2",
  "other_oncotree_diagnosis_3",
  "number_of_other_invasive_cancer_diagnosis",
  "other_invasive_cancer_diagnosis",
  "number_of_cancers_any_type",
  "number_of_bpc_project_cancers_index_cancers",
  
  # Tumor characteristics
  "ki67_percent",
  "oncotype_dx_score",
  "ajcc_stage",
  "stage_at_diagnosis",
  "grade",
  "her2_status",
  "hormone_receptor_status",
  "er_pr_receptor_change",
  "discordance_receptor_status_with_primary",
  "her2_receptor_status_at_metastatic_diagnosis",
  "her2_receptor_status_at_initial_diagnosis",
  "hr_status_at_metastatic_diagnosis",
  "hr_status_at_initial_diagnosis",
  "site_of_sample_tested",
  "stage_at_initial_breast_cancer_diagnosis",
  "subtype_of_metastatic_diagnosis",
  "subtype_of_initial_diagnosis",
  "histology",
  "histology_category",
  
  # Adjuvant and neoadjuvant treatment
  "number_of_cancer_directed_drug_regimens_curated",
  "received_adjuvant_ai",
  "received_adjuvant_tam_therapy",
  "primary_diagnosis_radiation_therapy",
  "cdk4_6_inhibitor_overall",
  "akt_inhibitor_overall",
  "akt_mutation_status",
  "aromatase_inhibitor_overall",
  "chemotherapy_overall",
  "chemotherapy_received_in_lrr_treatment",
  "chemotherapy_received_in_primary_disease_treatment",
  "chemotherapy_in_first_line_treatment",
  "chemotherapy_in_second_line_treatment",
  "chemotherapy_lines_received_in_metastatic_treatment",
  "overall_endocrine_therapy_sensitivity",
  "endocrine_therapy_in_first_line_treatment",
  "endocrine_therapy_in_second_line_treatment",
  "number_of_endocrine_therapies_received_in_metastatic_treatment",
  "number_of_endocrine_therapies_received_in_lrr_treatment",
  "tamoxifen_overall",
  "therapeutic_clinical_trial_overall",
  "neoadjuvant_chemotherapy_or_radiation_therapy_before_pathologic_stage_diagnosis",
  
  # Tumor and sequence characteristics
  "sequence_assay_id",
  "age_at_which_sequencing_was_reported",
  "age_at_sample_collection_4",
  "age_at_sequencing",
  "sequenced_sample",
  "how_long_sample_was_collected_after_starting_first_line_systemic_therapy_30_days_30_days_treatment_naive",
  "sample_site",
  "number_of_samples_per_patient",
  "sample_type",
  "sequencing_method",
  
  # Genetic alterations
  "fraction_genome_altered",
  "mutation_count",
  "mismatch_repair_mmr_testing_at_time_of_sample_acquisition",
  "msi_h_test_result_at_time_of_sample_acquisition",
  "microsatellite_instability_msi_testing_at_time_of_sample_acquisition",
  "tumor_mutational_burden",
  "tmb_nonsynonymous",
  "x1p_19q_codeletion",
  "idh1_2_mutation",
  "integrated_histomolecular_group",
  "mgmt_methylation_status",
  
  # Survival outcomes
  "patients_age_of_death_in_days",
  "age_at_first_distant_metastasis_in_days",
  "time_from_initial_diagnosis_to_metastatic_diagnosis_months",
  "cause_of_death",
  "patients_vital_status",
  "vital_status",
  "overall_survival_months",
  "overall_survival_status",
  "patients_age_of_last_follow_up_in_days",
  "interval_in_days_from_dob_to_date_of_last_contact",
  "interval_in_days_from_dob_to_dod",
  "sample_class",
  "year_of_last_contact",
  "year_of_death",
  "pfs_i_from_diagnosis_status",
  "pfs_m_from_diagnosis_status",
  
  # LRR/metastatic time
  "age_at_metastatic_diagnosis_years",
  "lrr_date_interval_days",
  "lrr_radiation_therapy",
  "site_of_lrr",
  "distant_metastatic_disease_interval",
  
  # Metastatic sites
  "met_site_soft_tissue",
  "met_site_visceral",
  "sites_of_distant_metastasis_at_the_time_of_cancer_diagnosis_stage_iv_patients",
  "year_of_next_generation_sequencing",
  "distant_mets_adrenal",
  "distant_mets_bone",
  "distant_mets_brain",
  "distant_mets_liver",
  "distant_mets_lung",
  "distant_mets_lymph_nodes",
  "distant_mets_other",
  "distant_mets_pleura",
  "distant_mets_subcutaneous_tissue",
  "age_at_metastatic_diagnosis",
  "met_site_liver",
  "has_this_patient_experienced_a_loco_regional_recurrence_lrr",
  "met_site_lung",
  "met_site_lymph_node",
  "total_number_of_therapies_received_in_metastatic_disease_treatment",
  "interval_between_sequencing_and_metastatis_daignosis",
  "met_site_multiple_sites",
  "met_site_other_site",
  "met_site_bone_only",
  "met_site_bone",
  "met_site_brain",
  
  # Additional identifier
  "patient_id_in_supp_table_1"
)

main(
  data_path = genie_data_path,
  dict_path = genie_dict_path,
  out_dir   = genie_out_dir,
  vars      = genie_vars
)


