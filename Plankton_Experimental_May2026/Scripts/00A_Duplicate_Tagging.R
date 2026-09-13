# ============================================================
# 00A_Duplicate_Tagging.R
# Purpose:
#   Replicate the previous Leung et al. 2022 duplicate-tagging
#   logic, but keep everything in memory.
#
# What this script does:
#   1. Reads the current database.
#   2. Handles the extra grouping row in the new Excel file.
#   3. Compares Leung et al. 2022 rows against non-Leung rows.
#   4. Flags Leung rows as duplicates if:
#        - strict match: surname + year + species + treatment
#        - moderate match: surname + species + treatment
#          AND Biota_C is close / matching
#   5. Creates zoo_db_with_duplicate_flag for 00_Dependencies.R.
#
# Main object created:
#   zoo_db_with_duplicate_flag
#
# Other objects created:
#   raw_db_no_extra_row
#   candidate_matches_all
#   moderate_matches_biotac_flagged
#   duplicate_tagging_summary
#
# Important:
#   This script does NOT save a new database file.
# ============================================================

# -----------------------------
# Load packages
# -----------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(readxl)
  library(stringr)
  library(tibble)
})


# -----------------------------
# Database path
# -----------------------------

# For CSV:
# db_path <- file.path("Data", "zooplankton_database_with_duplicate_flag.csv")

# For Excel:
# db_path <- file.path("Data", "zooplankton_database_with_duplicate_flag.xlsx") # old file
# db_path <- file.path("Data", "Updated_zooplankton_database_25May2026.xlsx") # new file
db_path <- file.path("Data", "Updated_zooplankton_database_10Sept2026.xlsx") # newer file


# -----------------------------
# Check database path
# -----------------------------

if (!file.exists(db_path)) {
  stop("Database file not found. Check db_path:\n", db_path)
}


# ============================================================
# 01 Read database
# ============================================================

file_ext <- tools::file_ext(db_path) %>%
  stringr::str_to_lower()

# New database has one grouping row before the actual column names.
# Old database does not.

if (stringr::str_detect(basename(db_path), "Updated_zooplankton_database_10Sept2026")) {
  skip_rows <- 1
} else {
  skip_rows <- 0
}

if (file_ext == "csv") {
  
  zoo_db <- readr::read_csv(
    db_path,
    show_col_types = FALSE
  )
  
} else if (file_ext %in% c("xlsx", "xls")) {
  
  zoo_db <- readxl::read_excel(
    path = db_path,
    sheet = 1,
    skip = skip_rows
  )
  
} else {
  
  stop(
    "Unsupported database file type: .", file_ext,
    "\nUse .csv, .xlsx, or .xls."
  )
}

raw_db_no_extra_row <- zoo_db

cat("\nDatabase loaded for duplicate tagging.\n")
cat("Database path:", db_path, "\n")
cat("Skipped rows:", skip_rows, "\n")
cat("Rows in raw database:", nrow(zoo_db), "\n")
cat("Columns in raw database:", ncol(zoo_db), "\n\n")


# ============================================================
# 02 Helper functions
# ============================================================

clean_text <- function(x) {
  x %>%
    as.character() %>%
    str_to_lower() %>%
    str_replace_all("[[:punct:]]", " ") %>%
    str_squish()
}

extract_year_from_author <- function(x) {
  str_extract(as.character(x), "(19|20)\\d{2}")
}

extract_final_year <- function(author_col, year_col) {
  year_col_chr   <- str_squish(as.character(year_col))
  year_from_year <- str_extract(year_col_chr, "(19|20)\\d{2}")
  year_from_auth <- extract_year_from_author(author_col)
  coalesce(year_from_year, year_from_auth)
}

extract_first_author_surname <- function(x) {
  x <- clean_text(x)
  x <- str_replace(x, "\\bet al\\b.*$", "")
  x <- str_replace(x, ",.*$", "")
  x <- str_replace(x, "\\sand\\s.*$", "")
  x <- str_squish(x)
  word(x, 1)
}


# ============================================================
# 03 Ensure required columns exist
# ============================================================

required_cols <- c(
  "Include",
  "Meta_analysis_reference",
  "Author",
  "Year",
  "Species",
  "Treatment",
  "Biota_C",
  "Duplicate"
)

for (col in required_cols) {
  if (!col %in% names(zoo_db)) {
    zoo_db[[col]] <- NA_character_
  }
}


# ============================================================
# 04 Standardise key columns
# ============================================================

dat <- zoo_db %>%
  mutate(
    row_id = row_number(),
    
    Existing_Duplicate = str_squish(as.character(Duplicate)),
    
    Include = str_squish(as.character(Include)),
    Meta_analysis_reference = str_squish(as.character(Meta_analysis_reference)),
    Author = str_squish(as.character(Author)),
    Year_original = str_squish(as.character(Year)),
    Species = str_squish(as.character(Species)),
    
    Treatment = case_when(
      Treatment %in% c("T1", "T2", "T3") ~ "T",
      Treatment %in% c("TpH1", "TpH2", "T3pH1") ~ "TpH",
      TRUE ~ as.character(Treatment)
    ),
    Treatment = str_squish(as.character(Treatment))
  ) %>%
  filter(Include == "Yes") %>%
  mutate(
    is_leung = str_detect(
      Meta_analysis_reference,
      regex(
        "Leung\\s*et\\s*al\\.?\\s*,?\\s*2022|Leung\\s*et\\s*al\\.?\\s*2022",
        ignore_case = TRUE
      )
    ),
    Year_extracted_from_author = extract_year_from_author(Author),
    Year_final = extract_final_year(Author, Year_original),
    Species_key = clean_text(Species),
    Treatment_key = clean_text(Treatment),
    First_author_surname = extract_first_author_surname(Author)
  )


# ============================================================
# 05 Split Leung vs non-Leung rows
# ============================================================

leung_df <- dat %>%
  filter(is_leung) %>%
  transmute(
    leung_row_id = row_id,
    leung_reference = Meta_analysis_reference,
    leung_author = Author,
    leung_year_original = Year_original,
    leung_year_extracted_from_author = Year_extracted_from_author,
    leung_year_final = Year_final,
    leung_species = Species,
    leung_treatment = Treatment,
    leung_species_key = Species_key,
    leung_treatment_key = Treatment_key,
    leung_surname_key = First_author_surname
  )

nonleung_df <- dat %>%
  filter(!is_leung) %>%
  transmute(
    nonleung_row_id = row_id,
    nonleung_reference = Meta_analysis_reference,
    nonleung_author = Author,
    nonleung_year_original = Year_original,
    nonleung_year_extracted_from_author = Year_extracted_from_author,
    nonleung_year_final = Year_final,
    nonleung_species = Species,
    nonleung_treatment = Treatment,
    nonleung_species_key = Species_key,
    nonleung_treatment_key = Treatment_key,
    nonleung_surname_key = First_author_surname
  )

cat("Rows with Include == 'Yes':", nrow(dat), "\n")
cat("Leung rows:", nrow(leung_df), "\n")
cat("Non-Leung rows:", nrow(nonleung_df), "\n\n")


# ============================================================
# 06 Strict match
# surname + year + species + treatment
# ============================================================

strict_matches <- leung_df %>%
  inner_join(
    nonleung_df,
    by = c(
      "leung_surname_key" = "nonleung_surname_key",
      "leung_year_final" = "nonleung_year_final",
      "leung_species_key" = "nonleung_species_key",
      "leung_treatment_key" = "nonleung_treatment_key"
    )
  ) %>%
  mutate(
    match_type = "strict",
    match_author_surname = TRUE,
    match_year = TRUE,
    match_species = TRUE,
    match_treatment = TRUE,
    match_score = 8
  )


# ============================================================
# 07 Moderate match
# surname + species + treatment
# ============================================================

moderate_matches <- leung_df %>%
  inner_join(
    nonleung_df,
    by = c(
      "leung_surname_key" = "nonleung_surname_key",
      "leung_species_key" = "nonleung_species_key",
      "leung_treatment_key" = "nonleung_treatment_key"
    )
  ) %>%
  mutate(
    match_type = "moderate",
    match_author_surname = TRUE,
    match_year = leung_year_final == nonleung_year_final,
    match_species = TRUE,
    match_treatment = TRUE,
    match_score = 6 + ifelse(match_year, 2, 0)
  ) %>%
  filter(
    !(paste(leung_row_id, nonleung_row_id) %in%
        paste(strict_matches$leung_row_id, strict_matches$nonleung_row_id))
  )


# ============================================================
# 08 Loose match
# year + species + treatment
# ============================================================

loose_matches <- leung_df %>%
  inner_join(
    nonleung_df,
    by = c(
      "leung_year_final" = "nonleung_year_final",
      "leung_species_key" = "nonleung_species_key",
      "leung_treatment_key" = "nonleung_treatment_key"
    )
  ) %>%
  mutate(
    match_type = "loose",
    match_author_surname = leung_surname_key == nonleung_surname_key,
    match_year = TRUE,
    match_species = TRUE,
    match_treatment = TRUE,
    match_score = 6 + ifelse(match_author_surname, 2, 0)
  ) %>%
  filter(
    !(paste(leung_row_id, nonleung_row_id) %in%
        paste(strict_matches$leung_row_id, strict_matches$nonleung_row_id))
  ) %>%
  filter(
    !(paste(leung_row_id, nonleung_row_id) %in%
        paste(moderate_matches$leung_row_id, moderate_matches$nonleung_row_id))
  )


# ============================================================
# 09 Combine candidate matches
# ============================================================

candidate_matches_all <- bind_rows(
  strict_matches,
  moderate_matches,
  loose_matches
) %>%
  transmute(
    match_type,
    match_score,
    
    leung_row_id,
    nonleung_row_id,
    
    leung_reference,
    nonleung_reference,
    
    leung_author,
    nonleung_author,
    
    leung_year_original,
    nonleung_year_original,
    
    leung_year_extracted_from_author,
    nonleung_year_extracted_from_author,
    
    leung_year_final,
    nonleung_year_final,
    
    leung_species,
    nonleung_species,
    
    leung_treatment,
    nonleung_treatment,
    
    match_author_surname,
    match_year,
    match_species,
    match_treatment
  ) %>%
  arrange(
    factor(match_type, levels = c("strict", "moderate", "loose")),
    desc(match_score),
    leung_row_id,
    nonleung_row_id
  )


# ============================================================
# 10 Investigate moderate matches using Biota_C
# ============================================================

moderate_matches_only <- candidate_matches_all %>%
  filter(match_type == "moderate")

zoo_db_with_rowid <- zoo_db %>%
  mutate(row_id = row_number())

leung_biotac <- zoo_db_with_rowid %>%
  transmute(
    leung_row_id = row_id,
    leung_Biota_C = suppressWarnings(as.numeric(Biota_C))
  )

nonleung_biotac <- zoo_db_with_rowid %>%
  transmute(
    nonleung_row_id = row_id,
    nonleung_Biota_C = suppressWarnings(as.numeric(Biota_C))
  )

moderate_matches_biotac <- moderate_matches_only %>%
  left_join(leung_biotac, by = "leung_row_id") %>%
  left_join(nonleung_biotac, by = "nonleung_row_id") %>%
  mutate(
    biotac_both_nonmissing = !is.na(leung_Biota_C) & !is.na(nonleung_Biota_C),
    
    biotac_exact_match = biotac_both_nonmissing &
      leung_Biota_C == nonleung_Biota_C,
    
    biotac_match_2dp = biotac_both_nonmissing &
      round(leung_Biota_C, 2) == round(nonleung_Biota_C, 2),
    
    biotac_match_3dp = biotac_both_nonmissing &
      round(leung_Biota_C, 3) == round(nonleung_Biota_C, 3),
    
    biotac_abs_diff = if_else(
      biotac_both_nonmissing,
      abs(leung_Biota_C - nonleung_Biota_C),
      NA_real_
    ),
    
    biotac_close_tol_0_001 = biotac_both_nonmissing &
      biotac_abs_diff <= 0.001,
    
    biotac_close_tol_0_01 = biotac_both_nonmissing &
      biotac_abs_diff <= 0.01,
    
    biotac_close_tol_0_1 = biotac_both_nonmissing &
      biotac_abs_diff <= 0.1
  )

moderate_matches_biotac_flagged <- moderate_matches_biotac %>%
  mutate(
    biotac_match_flag = case_when(
      biotac_exact_match ~ "Exact match",
      biotac_match_3dp ~ "Match to 3 d.p.",
      biotac_match_2dp ~ "Match to 2 d.p.",
      biotac_close_tol_0_001 ~ "Within 0.001",
      biotac_close_tol_0_01 ~ "Within 0.01",
      biotac_close_tol_0_1 ~ "Within 0.1",
      TRUE ~ "No close match"
    )
  ) %>%
  arrange(
    factor(
      biotac_match_flag,
      levels = c(
        "Exact match",
        "Match to 3 d.p.",
        "Match to 2 d.p.",
        "Within 0.001",
        "Within 0.01",
        "Within 0.1",
        "No close match"
      )
    )
  )


# ============================================================
# 11 Identify Leung rows to flag
# ============================================================

strict_leung_ids <- candidate_matches_all %>%
  filter(match_type == "strict") %>%
  distinct(leung_row_id) %>%
  pull(leung_row_id)

moderate_close_leung_ids <- moderate_matches_biotac_flagged %>%
  filter(biotac_match_flag != "No close match") %>%
  distinct(leung_row_id) %>%
  pull(leung_row_id)

duplicate_leung_ids <- union(strict_leung_ids, moderate_close_leung_ids)


# ============================================================
# 12 Update duplicate flag
# Only mark Leung rows, replicating the previous workflow.
# ============================================================

zoo_db_with_duplicate_flag <- zoo_db %>%
  mutate(
    row_id = row_number(),
    Meta_analysis_reference = str_squish(as.character(Meta_analysis_reference)),
    is_leung = str_detect(
      Meta_analysis_reference,
      regex(
        "Leung\\s*et\\s*al\\.?\\s*,?\\s*2022|Leung\\s*et\\s*al\\.?\\s*2022",
        ignore_case = TRUE
      )
    ),
    Duplicate = if_else(
      is_leung & row_id %in% duplicate_leung_ids,
      "duplicate",
      ""
    )
  ) %>%
  select(-row_id, -is_leung)

# For downstream scripts
zoo_db <- zoo_db_with_duplicate_flag
raw_db_no_extra_row <- zoo_db_with_duplicate_flag


# ============================================================
# 13 Summaries
# ============================================================

summary_table <- tibble(
  category = c(
    "Included rows",
    "Leung rows",
    "Non-Leung rows",
    "All candidate matches",
    "Strict matches",
    "Moderate matches",
    "Loose matches"
  ),
  n = c(
    nrow(dat),
    nrow(leung_df),
    nrow(nonleung_df),
    nrow(candidate_matches_all),
    sum(candidate_matches_all$match_type == "strict"),
    sum(candidate_matches_all$match_type == "moderate"),
    sum(candidate_matches_all$match_type == "loose")
  )
)

strict_leung_summary <- candidate_matches_all %>%
  filter(match_type == "strict") %>%
  summarise(
    n_strict_match_pairs = n(),
    n_unique_leung_rows_strictly_matched = n_distinct(leung_row_id)
  )

moderate_followup_summary <- moderate_matches_biotac_flagged %>%
  count(biotac_match_flag, sort = TRUE) %>%
  mutate(
    pct_of_moderate_pairs = round((n / sum(n)) * 100, 2)
  )

duplicate_flag_summary <- tibble(
  category = c(
    "Unique Leung rows strictly matched",
    "Unique Leung rows moderately matched with Biota_C support",
    "Unique Leung rows flagged as duplicate total"
  ),
  n = c(
    length(unique(strict_leung_ids)),
    length(unique(moderate_close_leung_ids)),
    length(unique(duplicate_leung_ids))
  )
)

duplicate_tagging_summary <- tibble(
  category = c(
    "Rows in full database",
    "Rows with Include == 'Yes'",
    "Leung rows checked",
    "Non-Leung rows checked",
    "Duplicate rows tagged by replicated Leung logic"
  ),
  n = c(
    nrow(zoo_db_with_duplicate_flag),
    nrow(dat),
    nrow(leung_df),
    nrow(nonleung_df),
    sum(zoo_db_with_duplicate_flag$Duplicate == "duplicate", na.rm = TRUE)
  )
)


# ============================================================
# 14 Print summaries
# ============================================================

cat("\n==============================\n")
cat("LEUNG DUPLICATE MATCH SUMMARY\n")
cat("==============================\n\n")
print(summary_table)

cat("\n==============================\n")
cat("STRICT LEUNG MATCH SUMMARY\n")
cat("==============================\n\n")
print(strict_leung_summary)

cat("\n==============================\n")
cat("MODERATE BIOTA_C FOLLOW-UP SUMMARY\n")
cat("==============================\n\n")
print(moderate_followup_summary)

cat("\n==============================\n")
cat("DUPLICATE FLAG SUMMARY\n")
cat("==============================\n\n")
print(duplicate_flag_summary)

cat("\n==============================\n")
cat("FINAL DUPLICATE TAGGING SUMMARY\n")
cat("==============================\n\n")
print(duplicate_tagging_summary)

cat("\n==============================\n")
cat("DUPLICATE COLUMN CHECK\n")
cat("==============================\n\n")
print(table(zoo_db_with_duplicate_flag$Duplicate, useNA = "ifany"))

cat("\nDone. No new database file was saved.\n")
cat("Object created: zoo_db_with_duplicate_flag\n")
cat("Object created: raw_db_no_extra_row\n")