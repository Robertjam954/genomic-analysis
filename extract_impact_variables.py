# ================================================================
# extract_items_from_dataset.R — Reusable extractor/validator
# ================================================================


  library(readr); library(readxl); library(dplyr); library(stringr)
  library(lubridate); library(tidyr); library(purrr); library(optparse)


source("extract_items_from_dataset.R")
main(
  data_path = "filtered_goel_impact_deid.csv",
  dict_path = "goel_data_dictionary_from_snakecase.xlsx",
  out_dir   = "
  vars      = c("mrn","age_years","bmi","goel_subtype_cat_4","stage_anatomic_initial_cat_4","recurrence_bin")
)


`%||%` <- function(a, b) if (is.null(a) || (length(a) == 1 && is.na(a))) b else a

snake <- function(x) x |> str_trim() |> str_to_lower() |>
  str_replace_all("[^a-z0-9]+", "_") |> str_replace_all("^_|_$", "")

coerce_to_class <- function(x, target_class)
  if (is.null(target_class) || is.na(target_class) || target_class == "") return(list(vec=x, note="no target class"))
  tc <- tolower(target_class)

  if (grepl("date", tc)) 
    if (inherits(x, "Date")) return(list(vec=x, note="already Date"))
    if (is.numeric(x)) return(list(vec=as.Date(x, origin="1899-12-30"), note="serial -> Date"))
    v <- suppressWarnings(mdy(x)); if (all(is.na(v))) v <- suppressWarnings(ymd(x))
    return(list(vec=as.Date(v), note="parsed mdy/ymd -> Date"))
  
  if (grepl("binary", tc)) 
    v <- x
    if (!is.numeric(v)) 
      v <- case_when(
        str_to_lower(as.character(x)) %in% c("1","yes","y","true","t") ~ 1,
        str_to_lower(as.character(x)) %in% c("0","no","n","false","f") ~ 0,
        TRUE ~ suppressWarnings(as.numeric(as.character(x)))
      )
    
    return(list(vec=as.integer(v), note="coerced to 0/1"))
  
  if (grepl("logical", tc)) 
    v <- str_to_lower(as.character(x))
    v <- case_when(
      v %in% c("true","t","1","yes","y") ~ TRUE,
      v %in% c("false","f","0","no","n") ~ FALSE,
      TRUE ~ NA
    )
    return(list(vec=as.logical(v), note="to logical"))
  
  if (grepl("numeric", tc)) return(list(vec=suppressWarnings(as.numeric(as.character(x))), note="to numeric"))
  if (grepl("factor",  tc)) return(list(vec=as.character(x), note="kept char (factor validated by levels)"))
  if (grepl("character", tc)) return(list(vec=as.character(x), note="to character"))
  list(vec=x, note="no rule")


parse_allowed <- function(txt)
  if (is.na(txt) || txt == "") return(list(kind="any", values=NULL, range=NULL, note="none"))
  txt <- str_trim(as.character(txt))
  if (grepl("^mm/dd/yyyy", str_to_lower(txt))) return(list(kind="date", values=NULL, range=NULL, note="date"))
  if (grepl("^range:", str_to_lower(txt))) 
    rng <- str_replace_all(txt, "[^0-9eE+\\.-]", " ")
    nums <- suppressWarnings(as.numeric(unlist(str_split(rng, "\\s+"))))
    nums <- nums[!is.na(nums)]
    if (length(nums) >= 2) return(list(kind="range", values=NULL, range=range(nums[1:2]), note="range"))
  
  if (grepl(",", txt)) 
    vals <- str_split(txt, ",")[[1]] |> str_trim()
    return(list(kind="set", values=unique(vals), range=NULL, note="set"))
  
  list(kind="any", values=NULL, range=NULL, note="none")


validate_column <- function(vec, target_class, allowed_spec)
  n <- length(vec); invalid_idx <- rep(FALSE, n)
  if (!is.null(target_class) && !is.na(target_class) && grepl("binary", tolower(target_class))) 
    invalid_idx <- !(is.na(vec) | vec %in% c(0L,1L))
   else if (allowed_spec$kind == "set") 
    invalid_idx <- !(is.na(vec) | str_trim(as.character(vec)) %in% allowed_spec$values)
   else if (allowed_spec$kind == "range") 
    vnum <- suppressWarnings(as.numeric(vec))
    invalid_idx <- !(is.na(vnum) | (vnum >= allowed_spec$range[1] & vnum <= allowed_spec$range[2]))
   else if (allowed_spec$kind == "date") 
    vdate <- suppressWarnings(mdy(vec)); vdate2 <- suppressWarnings(ymd(vec))
    ok <- !is.na(vdate) | !is.na(vdate2) | inherits(vec, "Date")
    invalid_idx <- !ok
  
  tibble::tibble(n=n, n_na=sum(is.na(vec)), n_invalid=sum(invalid_idx),
                 invalid_rows=list(which(invalid_idx)[1:min(20,sum(invalid_idx))]))


suggest_source <- function(var_name)
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


main <- function(data_path, dict_path, out_dir="out", vars=NULL)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  dat <- if (grepl("\\.csv$", data_path, ignore.case=TRUE)) 
    suppressMessages(readr::read_csv(data_path, show_col_types = FALSE))
   else readxl::read_excel(data_path)
  names(dat) <- snake(names(dat))

  dd <- readxl::read_excel(dict_path, sheet=1) |> mutate(variable = snake(variable))

  vars <- if (is.null(vars) || length(vars)==0) dd$variable else snake(vars)

  present <- intersect(vars, names(dat))
  missing <- setdiff(vars, names(dat))

  dd_map <- dd |> select(variable, class, allowed_values) |> distinct() |>
    mutate(allowed_parsed = map(allowed_values, parse_allowed))

  logs <- list()
  for (v in present) 
    spec <- dd_map |> filter(variable==v) |> slice(1)
    target_class <- spec$class %||% NA_character_
    allowed_spec <- spec$allowed_parsed[[1]] %||% list(kind="any", values=NULL, range=NULL)

    coerced <- coerce_to_class(dat[[v]], target_class)
    dat[[v]] <- coerced$vec

    val <- validate_column(dat[[v]], target_class, allowed_spec)
    val$variable <- v; val$class <- target_class
    val$coercion_note <- coerced$note; val$allowed_note <- allowed_spec$note
    logs[[v]] <- val
  

  qa <- bind_rows(logs) |>
    select(variable, class, coercion_note, allowed_note, n, n_na, n_invalid, invalid_rows) |>
    arrange(desc(n_invalid), variable)

  missing_df <- tibble(variable = missing,
                       suggestion = map_chr(missing, suggest_source))

  subset_path <- file.path(out_dir, "extracted_subset.csv")
  qa_path     <- file.path(out_dir, "qa_report.csv")
  miss_path   <- file.path(out_dir, "missing_variables_suggestions.csv")

  readr::write_csv(select(dat, any_of(present)), subset_path, na = "")
  readr::write_csv(qa, qa_path, na = "")
  readr::write_csv(missing_df, miss_path, na = "")

  message("Wrote: ", subset_path)
  message("Wrote: ", qa_path)
  message("Wrote: ", miss_path)
  invisible(list(subset=subset_path, qa=qa_path, missing=miss_path))


# CLI
if (sys.nframe() == 0) 
  option_list <- list(
    optparse::make_option(c("-d","--data"), type="character", help="Path to dataset (CSV/Excel)"),
    optparse::make_option(c("-k","--dict"), type="character", help="Path to data dictionary (Excel)"),
    optparse::make_option(c("-o","--outdir"), type="character", default="out", help="Output directory"),
    optparse::make_option(c("-v","--vars"), type="character", default=NULL, help="Comma-separated variables to extract (snake_case). If omitted, uses all from dictionary.")
  )
  opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))
  if (is.null(opt$data) || is.null(opt$dict)) stop("Provide --data and --dict")
  vars_vec <- if (!is.null(opt$vars)) strsplit(opt$vars, ",")[[1]] |> trimws() else NULL
  main(opt$data, opt$dict, out_dir = opt$outdir, vars = vars_vec)

