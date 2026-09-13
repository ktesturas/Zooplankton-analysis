# ============================================================
# 09_Add_Alter_Columns_To_TRUE_Raw_Database.R
# Purpose:
#   Read the TRUE raw zooplankton database directly from file
#   before any filtering from 00_Dependencies.R.
#
#   Then add columns from Alter et al. 2024:
#     Metric
#     Metric_Type
#     Unit
#     Old_Response_category
#
# Important:
#   This should keep the original raw row count, e.g. 4412 rows.
#   It does NOT use dat1 because dat1 is already filtered.
#   It does NOT calculate yi/vi.
#   It does NOT use ES_ID for matching.
#   It does NOT use Year for matching.
#   It uses sheet 1 of the Alter file.
#
# Outputs:
#   Outputs/Tables/09_Add_Alter_Metric_Type_To_Raw_First_Database/
#   Outputs/RDS/09_Add_Alter_Metric_Type_To_Raw_First_Database/
# ============================================================

suppressPackageStartupMessages(invisible(
  lapply(
    c(
      "dplyr", "readxl", "readr", "stringr", "tibble", "tidyr"
    ),
    library,
    character.only = TRUE
  )
))

# ============================================================
# File paths
# ============================================================

raw_db_path <- file.path(
  "Data",
  "zooplankton_database_with_duplicate_flag.xlsx"
)

alter_path <- "C:/Users/Kris Jypson Esturas/OneDrive - Macquarie University/Documents/2025/00 Plankton Work/2025/Plankton/00 Final/Plankton/Experimental/Zooplankton_experiment_database/Reference Meta-analysis/Alter_etal_2024_DataInverts_EDITED.csv.xlsx"

# ============================================================
# Output folders
# Use same folder as before, so outputs are overwritten
# ============================================================

script_name <- "09_Add_Alter_Metric_Type_To_Raw_First_Database"

table_dir <- file.path("Outputs", "Tables", script_name)
rds_dir   <- file.path("Outputs", "RDS", script_name)

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# Helper functions
# ============================================================

clean_text <- function(x) {
  x %>%
    as.character() %>%
    stringr::str_to_lower() %>%
    stringr::str_replace_all("[^a-z0-9]+", " ") %>%
    stringr::str_squish()
}

clean_author_key <- function(x) {
  x %>%
    clean_text() %>%
    stringr::str_replace_all("\\bet al\\b", "") %>%
    stringr::str_squish()
}

num_key <- function(x, digits = 6) {
  suppressWarnings(as.numeric(x)) %>%
    round(digits = digits)
}

make_row_key <- function(df, key_cols, key_name = "match_key") {
  
  missing_cols <- setdiff(key_cols, names(df))
  
  if (length(missing_cols) > 0) {
    stop(
      "Missing key columns in make_row_key(): ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  df %>%
    mutate(
      "{key_name}" := do.call(
        paste,
        c(
          dplyr::across(dplyr::all_of(key_cols)),
          sep = " | "
        )
      )
    )
}

run_matching_stage <- function(first_pool, alter_pool, key_cols, key_name, stage_name) {
  
  if (nrow(first_pool) == 0) {
    return(list(
      matched = first_pool %>% slice(0),
      unmatched_first = first_pool,
      remaining_alter = alter_pool
    ))
  }
  
  if (nrow(alter_pool) == 0) {
    return(list(
      matched = first_pool %>% slice(0),
      unmatched_first = first_pool,
      remaining_alter = alter_pool
    ))
  }
  
  first_keyed <- first_pool %>%
    make_row_key(key_cols, key_name = key_name) %>%
    group_by(.data[[key_name]]) %>%
    mutate(dup_id_match = row_number()) %>%
    ungroup()
  
  alter_keyed <- alter_pool %>%
    make_row_key(key_cols, key_name = key_name) %>%
    group_by(.data[[key_name]]) %>%
    mutate(dup_id_match = row_number()) %>%
    ungroup() %>%
    select(
      all_of(key_name),
      dup_id_match,
      row_id_alter,
      Alter_Metric = Metric,
      Alter_Metric_Type = `Metric Type`,
      Alter_Unit = Unit,
      Alter_Old_Response_category = Old_Response_category,
      Alter_LifeStage = LifeStage,
      Alter_Response_Biota = Response_Biota,
      Alter_ES_ID = ES_ID,
      Alter_Author = Author,
      Alter_Year = Year,
      Alter_Species = Species,
      Alter_Treatment = Treatment,
      Alter_Biota_C = Biota_C,
      Alter_Biota_T = Biota_T,
      Alter_Biota_C_sd = Biota_C_sd,
      Alter_Biota_T_sd = Biota_T_sd,
      Alter_N_C = N_C,
      Alter_N_T = N_T
    )
  
  joined <- first_keyed %>%
    left_join(
      alter_keyed,
      by = c(key_name, "dup_id_match")
    )
  
  matched <- joined %>%
    filter(!is.na(row_id_alter)) %>%
    mutate(
      Match_status = "matched",
      match_stage = stage_name
    )
  
  unmatched_first <- joined %>%
    filter(is.na(row_id_alter)) %>%
    select(all_of(names(first_pool)))
  
  remaining_alter <- alter_pool %>%
    anti_join(
      matched %>%
        distinct(row_id_alter),
      by = "row_id_alter"
    )
  
  list(
    matched = matched,
    unmatched_first = unmatched_first,
    remaining_alter = remaining_alter
  )
}

# ============================================================
# Read TRUE raw first database directly
# ============================================================

if (!file.exists(raw_db_path)) {
  stop("Raw database file not found. Check path:\n", raw_db_path)
}

raw_first <- readxl::read_excel(raw_db_path) %>%
  mutate(row_id_first = row_number())

cat("\n==== TRUE raw first database ====\n")
cat("Rows in raw_first:", nrow(raw_first), "\n")
cat("Columns in raw_first:", ncol(raw_first), "\n\n")

# This should be 4412 if that is the original row count
if (nrow(raw_first) != 4412) {
  warning(
    "raw_first does not have 4412 rows. It has ",
    nrow(raw_first),
    " rows. Check whether raw_db_path points to the correct raw file."
  )
}

# ============================================================
# Read Alter database, sheet 1
# ============================================================

cat("\n==== Loading Alter et al. 2024 database ====\n")

alter_sheets <- readxl::excel_sheets(alter_path)

cat("Available sheets:\n")
print(alter_sheets)

alter_dat <- readxl::read_excel(
  alter_path,
  sheet = 1
)

cat("\nRows in alter_dat:", nrow(alter_dat), "\n")
cat("Columns in alter_dat:", ncol(alter_dat), "\n\n")

# ============================================================
# Check required columns
# ============================================================

needed_first_cols <- c(
  "Meta_analysis_reference",
  "Author",
  "Species",
  "Treatment",
  "Response_Biota",
  "Biota_C",
  "Biota_T",
  "Biota_C_sd",
  "Biota_T_sd",
  "N_C",
  "N_T"
)

missing_first_cols <- setdiff(needed_first_cols, names(raw_first))

if (length(missing_first_cols) > 0) {
  stop(
    "The following columns are missing from raw_first: ",
    paste(missing_first_cols, collapse = ", ")
  )
}

needed_alter_cols <- c(
  "ES_ID",
  "Author",
  "Year",
  "Species",
  "LifeStage",
  "Treatment",
  "Response_Biota",
  "Metric",
  "Metric Type",
  "Unit",
  "Old_Response_category",
  "Biota_C",
  "Biota_T",
  "Biota_C_sd",
  "Biota_T_sd",
  "N_C",
  "N_T"
)

missing_alter_cols <- setdiff(needed_alter_cols, names(alter_dat))

if (length(missing_alter_cols) > 0) {
  stop(
    "The following columns are missing from alter_dat: ",
    paste(missing_alter_cols, collapse = ", ")
  )
}

# ============================================================
# Filter raw first database to Alter rows only for matching
# Final product will still be the full 4412-row raw database.
# ============================================================

first_target <- raw_first %>%
  mutate(
    Meta_analysis_reference_chr = as.character(Meta_analysis_reference)
  ) %>%
  filter(
    stringr::str_detect(
      Meta_analysis_reference_chr,
      stringr::regex("Alter|Alters", ignore_case = TRUE)
    )
  )

cat("\n==== Raw first database Alter target rows ====\n")
cat("Rows in first_target:", nrow(first_target), "\n\n")

if (nrow(first_target) == 0) {
  stop(
    "No Alter/Alters rows found in Meta_analysis_reference. ",
    "Check exact spelling in Meta_analysis_reference."
  )
}

# ============================================================
# Filter Alter database to larval/embryo/postlarval entries
# ============================================================

alter_target <- alter_dat %>%
  mutate(
    row_id_alter = row_number(),
    LifeStage_clean = stringr::str_to_upper(stringr::str_squish(as.character(LifeStage)))
  ) %>%
  filter(
    LifeStage_clean %in% c("LARVA", "EMBRYO", "POSTLARVA")
  )

cat("\n==== Alter target rows ====\n")
cat("Rows in alter_target:", nrow(alter_target), "\n\n")

# ============================================================
# Prepare matching fields in raw first database
# ============================================================

first_match <- first_target %>%
  mutate(
    author_key = clean_author_key(Author),
    species_key = clean_text(Species),
    response_key = clean_text(Response_Biota),
    treatment_key = clean_text(Treatment),
    
    Biota_C_key = num_key(Biota_C),
    Biota_T_key = num_key(Biota_T),
    Biota_C_sd_key = num_key(Biota_C_sd),
    Biota_T_sd_key = num_key(Biota_T_sd),
    N_C_key = num_key(N_C),
    N_T_key = num_key(N_T)
  )

# ============================================================
# Prepare matching fields in Alter database
# ============================================================

alter_match <- alter_target %>%
  mutate(
    author_key = clean_author_key(Author),
    species_key = clean_text(Species),
    response_key = clean_text(Response_Biota),
    treatment_key = clean_text(Treatment),
    
    Biota_C_key = num_key(Biota_C),
    Biota_T_key = num_key(Biota_T),
    Biota_C_sd_key = num_key(Biota_C_sd),
    Biota_T_sd_key = num_key(Biota_T_sd),
    N_C_key = num_key(N_C),
    N_T_key = num_key(N_T)
  )

# ============================================================
# Matching stages
# No ES_ID.
# No Year.
# ============================================================

# Stage 1: strict match including Response_Biota, means, SDs, and N
strict_key_cols <- c(
  "author_key",
  "species_key",
  "treatment_key",
  "response_key",
  "Biota_C_key",
  "Biota_T_key",
  "Biota_C_sd_key",
  "Biota_T_sd_key",
  "N_C_key",
  "N_T_key"
)

stage1 <- run_matching_stage(
  first_pool = first_match,
  alter_pool = alter_match,
  key_cols = strict_key_cols,
  key_name = "strict_match_key",
  stage_name = "matched_by_strict_key"
)

cat("\n==== Match Stage 1: strict key ====\n")
cat("Matched rows:", nrow(stage1$matched), "\n")
cat("Unmatched first rows:", nrow(stage1$unmatched_first), "\n")
cat("Remaining Alter rows:", nrow(stage1$remaining_alter), "\n\n")

# Stage 2: moderate match excluding Response_Biota, keeping means, SDs, and N
moderate_key_cols <- c(
  "author_key",
  "species_key",
  "treatment_key",
  "Biota_C_key",
  "Biota_T_key",
  "Biota_C_sd_key",
  "Biota_T_sd_key",
  "N_C_key",
  "N_T_key"
)

stage2 <- run_matching_stage(
  first_pool = stage1$unmatched_first,
  alter_pool = stage1$remaining_alter,
  key_cols = moderate_key_cols,
  key_name = "moderate_match_key",
  stage_name = "matched_by_moderate_key"
)

cat("\n==== Match Stage 2: moderate key ====\n")
cat("Matched rows:", nrow(stage2$matched), "\n")
cat("Unmatched first rows:", nrow(stage2$unmatched_first), "\n")
cat("Remaining Alter rows:", nrow(stage2$remaining_alter), "\n\n")

# Stage 3: loose match excluding Response_Biota and SDs
loose_key_cols <- c(
  "author_key",
  "species_key",
  "treatment_key",
  "Biota_C_key",
  "Biota_T_key",
  "N_C_key",
  "N_T_key"
)

stage3 <- run_matching_stage(
  first_pool = stage2$unmatched_first,
  alter_pool = stage2$remaining_alter,
  key_cols = loose_key_cols,
  key_name = "loose_match_key",
  stage_name = "matched_by_loose_key"
)

cat("\n==== Match Stage 3: loose key ====\n")
cat("Matched rows:", nrow(stage3$matched), "\n")
cat("Unmatched first rows:", nrow(stage3$unmatched_first), "\n")
cat("Remaining Alter rows:", nrow(stage3$remaining_alter), "\n\n")

# ============================================================
# Combine matched and unmatched target rows
# ============================================================

unmatched_first_final <- stage3$unmatched_first %>%
  mutate(
    row_id_alter = NA_integer_,
    Alter_Metric = NA_character_,
    Alter_Metric_Type = NA_character_,
    Alter_Unit = NA_character_,
    Alter_Old_Response_category = NA_character_,
    Match_status = "not_matched_in_alter_sheet",
    match_stage = NA_character_
  )

matched_target_all <- bind_rows(
  stage1$matched,
  stage2$matched,
  stage3$matched,
  unmatched_first_final
)

cat("\n==== Matching summary ====\n")
print(table(matched_target_all$Match_status, useNA = "ifany"))
print(table(matched_target_all$match_stage, useNA = "ifany"))

# ============================================================
# Add new columns to the full TRUE raw database
# ============================================================

metric_fields <- matched_target_all %>%
  select(
    row_id_first,
    Metric = Alter_Metric,
    Metric_Type = Alter_Metric_Type,
    Unit = Alter_Unit,
    Old_Response_category = Alter_Old_Response_category,
    Match_status,
    match_stage,
    Alter_row_id = row_id_alter
  )

raw_first_with_metric <- raw_first %>%
  left_join(
    metric_fields,
    by = "row_id_first"
  ) %>%
  mutate(
    Match_status = dplyr::case_when(
      !stringr::str_detect(
        as.character(Meta_analysis_reference),
        stringr::regex("Alter|Alters", ignore_case = TRUE)
      ) ~ "not_alter_reference",
      is.na(Match_status) ~ "alter_reference_but_no_metric_match",
      TRUE ~ Match_status
    )
  ) %>%
  relocate(
    Metric,
    Metric_Type,
    Unit,
    Old_Response_category,
    .after = Response_Biota
  ) %>%
  relocate(
    Match_status,
    match_stage,
    Alter_row_id,
    .after = Old_Response_category
  )

cat("\n==== Final raw database with added columns ====\n")
cat("Rows in raw_first_with_metric:", nrow(raw_first_with_metric), "\n")
cat("Columns in raw_first_with_metric:", ncol(raw_first_with_metric), "\n\n")

# ============================================================
# Alter rows not matched into raw first database
# ============================================================

alter_not_in_first <- alter_match %>%
  anti_join(
    matched_target_all %>%
      filter(!is.na(row_id_alter)) %>%
      distinct(row_id_alter),
    by = "row_id_alter"
  ) %>%
  select(
    row_id_alter,
    ES_ID,
    Author,
    Year,
    Species,
    LifeStage,
    Treatment,
    Response_Biota,
    Metric,
    `Metric Type`,
    Unit,
    Old_Response_category,
    Biota_C,
    Biota_T,
    Biota_C_sd,
    Biota_T_sd,
    N_C,
    N_T,
    everything()
  )

first_not_matched <- matched_target_all %>%
  filter(is.na(row_id_alter))

# ============================================================
# Save outputs
# ============================================================

readr::write_csv(
  raw_first_with_metric,
  file.path(table_dir, "raw_4412_database_with_alter_columns.csv")
)

readr::write_csv(
  matched_target_all,
  file.path(table_dir, "alter_metric_matching_details.csv")
)

readr::write_csv(
  alter_not_in_first,
  file.path(table_dir, "alter_rows_not_matched_to_raw_database.csv")
)

readr::write_csv(
  first_not_matched,
  file.path(table_dir, "raw_database_alter_rows_without_metric_match.csv")
)

saveRDS(
  raw_first_with_metric,
  file.path(rds_dir, "raw_4412_database_with_alter_columns.rds")
)

# ============================================================
# Console summary
# ============================================================

cat("\n==== Final summary ====\n")
cat("Rows in TRUE raw first database:", nrow(raw_first), "\n")
cat("Rows in final raw database with added columns:", nrow(raw_first_with_metric), "\n")
cat("Rows in first Alter target:", nrow(first_target), "\n")
cat("Rows in Alter target sheet:", nrow(alter_target), "\n")
cat("Matched first target rows:", sum(matched_target_all$Match_status == "matched", na.rm = TRUE), "\n")
cat("Unmatched first target rows:", nrow(first_not_matched), "\n")
cat("Alter rows not matched to raw database:", nrow(alter_not_in_first), "\n\n")

cat("Match status in full raw database:\n")
print(table(raw_first_with_metric$Match_status, useNA = "ifany"))

cat("\n==== Saved outputs ====\n")
cat("Tables folder:", table_dir, "\n")
cat("RDS folder   :", rds_dir, "\n\n")

cat(
  "Main output:",
  file.path(table_dir, "raw_4412_database_with_alter_columns.csv"),
  "\n"
)

# ============================================================
# View in RStudio
# ============================================================

View(raw_first_with_metric)
View(alter_not_in_first)
View(first_not_matched)