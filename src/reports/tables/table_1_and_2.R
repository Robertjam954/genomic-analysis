# Load necessary libraries
library(tidyverse)   # for data manipulation (dplyr, ggplot2, etc.)
library(readxl)      # for reading Excel files
library(forcats)     # for factor level handling (comes with tidyverse, but explicitly noted)
library(gtsummary)   # for creating summary tables (if needed for checking)
library(tfrmt)       # for formatting tables

# Read the custom dataset (Excel file) into R
# (Ensure the working directory is set to where the file is located, or provide full path)
bc_data <- read_xlsx("complete_genie_bpc_dataset_deid.xlsx")

# Inspect column names (optional, to identify relevant variables)
names(bc_data)
# We expect columns for age, race, histology, initial subtype, metastatic subtype, grade, etc.
# Let's assume the dataset has the following relevant columns (adjust if actual names differ):
# age_at_diagnosis, race, histology, subtype_initial (initial diagnosis subtype), subtype_mets (metastatic subtype), grade, and time_to_met (months from initial diagnosis to metastasis).

# 1. Create derived and factor variables for Table 1 (All Breast Cancer Patients)
#    - Age category (<50, 50-70, >70)
#    - Race (ensure proper factor labeling/order)
#    - Histology (column L in dataset)
#    - Subtype of initial diagnosis (column AC)
#    - Grade (tumor grade: 1=well, 2=moderate, 3=poorly differentiated)

bc_data <- bc_data %>%
  # Create age category as a factor with specified breaks
  mutate(
    age_category = case_when(
      # Use the appropriate age column from your dataset. Assuming it's named "age_at_diagnosis".
      age_at_diagnosis < 50 ~ "<50",
      age_at_diagnosis >= 50 & age_at_diagnosis <= 70 ~ "50-70",
      age_at_diagnosis > 70 ~ ">70",
      TRUE ~ NA_character_
    ),
    age_category = factor(age_category, levels = c("<50", "50-70", ">70")),

    # Convert existing columns to factors with explicit level order and labels.
    # For Race, we define common categories and an "Unknown" if missing or not provided.
    race = fct_explicit_na(factor(race, levels = c("White", "Black or African American", "Asian", "Other")), na_level = "Unknown"),

    # Histology: define factor levels explicitly (example levels provided; adjust to actual data).
    histology = fct_explicit_na(factor(histology, levels = c("Invasive Ductal Carcinoma",
                                                             "Invasive Lobular Carcinoma",
                                                             "Mixed Histology",
                                                             "Other")), na_level = "Unknown"),

    # Subtype at initial diagnosis: e.g., hormone receptor (HR) / HER2 status combinations.
    subtype_initial = fct_explicit_na(factor(subtype_initial, levels = c("HR+/HER2-",
                                                                         "HR+/HER2+",
                                                                         "HR-/HER2+",
                                                                         "HR-/HER2-")), na_level = "Unknown"),

    # Grade: numeric 1,2,3 mapped to labels; treat missing as "Unknown"
    grade = factor(grade, levels = c(1, 2, 3),
                   labels = c("Well differentiated", "Moderately differentiated", "Poorly differentiated")),
    grade = fct_explicit_na(grade, na_level = "Unknown")
  )

# At this point, age_category, race, histology, subtype_initial, and grade are all factors
# with specified levels and labels (including "Unknown" for missing where applicable).

# 2. Construct Table 1: Demographic Characteristics of All Breast Cancer Patients.
# We will summarize the selected variables for the entire cohort.
# Using gtsummary for convenience (it will handle counts and percentages and formatting):
table1_summary <- bc_data %>%
  select(age_category, race, histology, subtype_initial, grade) %>%
  tbl_summary(
    label = list(
      age_category   ~ "Age at Diagnosis (years)",
      race           ~ "Race/Ethnicity",
      histology      ~ "Histology",
      subtype_initial~ "Subtype at Initial Diagnosis",
      grade          ~ "Grade"
    ),
    missing = "ifany",            # include missing ("Unknown") in counts if present
    statistic = list(all_categorical() ~ "{n} ({p}%)")  # show count and column percent for categories
  )
# (The gtsummary table1_summary is a formatted summary. We will further refine presentation with tfrmt next.)

# Alternatively, we can manually calculate counts and percentages and format with tfrmt for fine control:
# Calculate total N for reference
total_n <- nrow(bc_data)

# Summarize each categorical variable with counts and percentages
table1_data <- list()

# Age category distribution
table1_data$age <- bc_data %>%
  count(age_category) %>%
  mutate(
    percent = 100 * n / total_n,
    variable = "Age at Diagnosis (years)",
    category = as.character(age_category)
  )

# Race distribution
table1_data$race <- bc_data %>%
  count(race) %>%
  mutate(
    percent = 100 * n / total_n,
    variable = "Race/Ethnicity",
    category = as.character(race)
  )

# Histology distribution
table1_data$histology <- bc_data %>%
  count(histology) %>%
  mutate(
    percent = 100 * n / total_n,
    variable = "Histology",
    category = as.character(histology)
  )

# Subtype at initial diagnosis distribution
table1_data$subtype_init <- bc_data %>%
  count(subtype_initial) %>%
  mutate(
    percent = 100 * n / total_n,
    variable = "Subtype at Initial Diagnosis",
    category = as.character(subtype_initial)
  )

# Grade distribution
table1_data$grade <- bc_data %>%
  count(grade) %>%
  mutate(
    percent = 100 * n / total_n,
    variable = "Grade",
    category = as.character(grade)
  )

# Combine all summary rows into one data frame
table1_df <- bind_rows(table1_data$age, table1_data$race, table1_data$histology,
                       table1_data$subtype_init, table1_data$grade) %>%
  arrange(factor(variable, levels = c("Age at Diagnosis (years)", "Race/Ethnicity",
                                      "Histology", "Subtype at Initial Diagnosis", "Grade"))) %>%
  mutate(
    percent = round(percent, 1),  # round percentages to 1 decimal place
    stat = paste0(n, " (", percent, "%)")
  )
# table1_df now contains columns: variable (group), category (levels), and stat (formatted "n (percent%)").

# Use tfrmt to format Table 1 with grouped rows and indentation for categories
table1_tfrmt <- tfrmt(
  title = "Table 1. Demographic Characteristics of All Breast Cancer Patients",
  group = variable,      # grouping variable (e.g., "Age at Diagnosis", "Race/Ethnicity", etc.)
  label = category,      # label for sub-categories within each group
  value = stat,          # the statistic value to display (combined count and percent string)
  column = NULL,         # no separate column variable (all stats in one column for overall cohort)
  row_grp_plan = row_grp_plan(
    label_loc = element_row_grp_loc(location = "indented")  # combine group & label, indent sub-categories under group
  ),
  body_plan = body_plan(
    frmt_structure(group_val = ".default", label_val = ".default", frmt(""))  # default format (we will apply our stat string directly)
  )
)

# Generate a gt table from tfrmt specification and the summary data
table1_gt <- print_to_gt(table1_tfrmt, .data = table1_df)

# Adjust styling for proper indentation display (preserve spaces, left-align stub column)
table1_gt <- table1_gt %>%
  gt::tab_style(
    style = gt::cell_text(whitespace = "pre", align = "left"),
    locations = gt::cells_stub()  # apply to the stub (row label) column
  )

# (At this point, table1_gt is a gt table object ready for rendering, with indentation for categories.)

# 3. Construct Table 2: Metastasis Cohort Demographics.
# Filter data to only patients who had metastasis (i.e., those with a recorded time to metastasis).
mets_data <- bc_data %>% filter(!is.na(time_to_met))  # assuming `time_to_met` is NA for no metastasis

# Prepare factors in the metastasis subset similarly (in case any factors lost levels after filtering, refactor if needed)
mets_data <- mets_data %>%
  mutate(
    age_category = factor(age_category, levels = c("<50", "50-70", ">70")),  # reuse the same factor levels
    race = fct_explicit_na(factor(race, levels = c("White", "Black or African American", "Asian", "Other")), na_level = "Unknown"),
    histology = fct_explicit_na(factor(histology, levels = c("Invasive Ductal Carcinoma",
                                                             "Invasive Lobular Carcinoma",
                                                             "Mixed Histology",
                                                             "Other")), na_level = "Unknown"),
    subtype_mets = fct_explicit_na(factor(subtype_mets, levels = c("HR+/HER2-",
                                                                   "HR+/HER2+",
                                                                   "HR-/HER2+",
                                                                   "HR-/HER2-")), na_level = "Unknown"),
    grade = fct_explicit_na(factor(grade, levels = c("Well differentiated",
                                                     "Moderately differentiated",
                                                     "Poorly differentiated")), na_level = "Unknown")
  )
# (We assume subtype_mets is already present in the dataset (column AB).)

# Summarize Table 2 variables:
# Include age_category, race, histology, subtype_mets, grade (categorical like in Table1)
# and add "time_to_met (months)" as a continuous variable summary.
mets_n <- nrow(mets_data)  # number of patients in metastasis cohort

# Summaries for categorical variables in metastasis cohort
table2_data <- list()

table2_data$age <- mets_data %>%
  count(age_category) %>%
  mutate(percent = 100 * n / mets_n,
         variable = "Age at Diagnosis (years)",
         category = as.character(age_category))

table2_data$race <- mets_data %>%
  count(race) %>%
  mutate(percent = 100 * n / mets_n,
         variable = "Race/Ethnicity",
         category = as.character(race))

table2_data$histology <- mets_data %>%
  count(histology) %>%
  mutate(percent = 100 * n / mets_n,
         variable = "Histology",
         category = as.character(histology))

table2_data$subtype_mets <- mets_data %>%
  count(subtype_mets) %>%
  mutate(percent = 100 * n / mets_n,
         variable = "Subtype at Metastasis",
         category = as.character(subtype_mets))

table2_data$grade <- mets_data %>%
  count(grade) %>%
  mutate(percent = 100 * n / mets_n,
         variable = "Grade",
         category = as.character(grade))

# Summary for time from initial diagnosis to metastasis (continuous)
mean_time <- mean(mets_data$time_to_met, na.rm = TRUE)
sd_time   <- sd(mets_data$time_to_met, na.rm = TRUE)
# (Optionally, could calculate median and IQR as well if desired:
median_time <- median(mets_data$time_to_met, na.rm = TRUE)
iqr_time <- IQR(mets_data$time_to_met, na.rm = TRUE)  # IQR = Q3 - Q1
# We'll include mean (SD) as the summary statistic for time to metastasis.
)

# Create a row (or rows) for time to metastasis
table2_time_row <- tibble(
  variable = "Time from Initial Dx to Metastasis (months)",
  category = "Mean (SD)",  # label for the statistic presented
  n = NA,                  # not applicable for continuous summary (we'll use stat instead)
  percent = NA,
  stat = paste0(round(mean_time, 1), " (", round(sd_time, 1), ")")  # format mean (SD) with 1 decimal
)
# (If we wanted to include median (IQR) as another row, we could add:
# tibble(variable = "Time from Initial Dx to Metastasis (months)",
#        category = "Median (IQR)",
#        stat = paste0(round(median_time,1), " (", round(quantile(mets_data$time_to_met, 0.25, na.rm=TRUE),1), ", ",
#                      round(quantile(mets_data$time_to_met, 0.75, na.rm=TRUE),1), ")") )

# Combine all summary rows for Table 2
table2_df <- bind_rows(
  table2_data$age, table2_data$race, table2_data$histology,
  table2_data$subtype_mets, table2_data$grade,
  table2_time_row
) %>%
  # Order the variables in a logical sequence for presentation
  arrange(factor(variable, levels = c("Age at Diagnosis (years)", "Race/Ethnicity",
                                      "Histology", "Subtype at Metastasis", "Grade",
                                      "Time from Initial Dx to Metastasis (months)"))) %>%
  # Round percentages and create stat string for categorical rows
  mutate(
    percent = round(percent, 1),
    stat = if_else(is.na(stat), paste0(n, " (", percent, "%)"), stat)
  )
# (For the continuous variable row, stat was already defined as mean (SD); for categorical rows, we paste n and percent.)

# Use tfrmt to format Table 2 with similar structure (group = variable, label = category)
table2_tfrmt <- tfrmt(
  title = "Table 2. Demographic Characteristics of Metastatic Breast Cancer Cohort",
  group = variable,
  label = category,
  value = stat,
  column = NULL,
  row_grp_plan = row_grp_plan(
    label_loc = element_row_grp_loc(location = "indented")
  ),
  body_plan = body_plan(
    frmt_structure(group_val = ".default", label_val = ".default", frmt(""))
  )
)

# Generate gt table for Table 2
table2_gt <- print_to_gt(table2_tfrmt, .data = table2_df) %>%
  gt::tab_style(
    style = gt::cell_text(whitespace = "pre", align = "left"),
    locations = gt::cells_stub()
  )

# The resulting table1_gt and table2_gt are gt objects that can be printed (e.g., in RMarkdown or viewer)
# to display Table 1 and Table 2 with proper formatting.
# Each table has grouped row labels (e.g., "Age at Diagnosis", "Race/Ethnicity", etc.) and indented sub-categories with counts and percentages.
# Table 2 additionally includes the continuous variable "Time from Initial Dx to Metastasis (months)" with a mean (SD) summary.
