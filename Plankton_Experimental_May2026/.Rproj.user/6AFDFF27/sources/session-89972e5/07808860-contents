# ============================================================
# 00B_Dependencies.R
# shared database cleaning
# ============================================================

# db_path <- file.path("Data", "zooplankton_database_with_duplicate_flag.csv")
# db_path <- file.path("Data", "zooplankton_database_with_duplicate_flag.xlsx") # old file
db_path <- file.path("Data", "Updated_zooplankton_database_25May2026.xlsx") # new file

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(readxl)
  library(stringr)
  library(tibble)
  library(tidyr)
})

if (!file.exists(db_path)) {
  stop("Database file not found. Check db_path:\n", db_path)
}

out_dir <- file.path("Outputs", "Tables", "00B_Dependencies")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# re-tag duplicates
source(file.path("Scripts", "00A_Duplicate_Tagging.R"))

if (!exists("zoo_db_with_duplicate_flag")) {
  stop("zoo_db_with_duplicate_flag not found.")
}

if (!exists("raw_db_no_extra_row")) {
  stop("raw_db_no_extra_row not found.")
}

zoo_db <- zoo_db_with_duplicate_flag

cat("\nDatabase loaded from 00A_Duplicate_Tagging.R\n")
cat("Database path:", db_path, "\n")
cat("Rows:", nrow(zoo_db), "\n")
cat("Columns:", ncol(zoo_db), "\n\n")

cat("Duplicate summary:\n")
print(table(zoo_db$Duplicate, useNA = "ifany"))

cat("\nColumn names:\n")
print(names(zoo_db))

# clean treatment values and deltas
zoo_db <- zoo_db %>%
  mutate(
    across(
      c(C_value_CO2, C_value_Temp, T_value_CO2, T_value_Temp),
      ~ readr::parse_number(
        as.character(.),
        na = c("", "NA", "na", "Na", "n/a", "N/A")
      )
    ),
    delta_CO2 = T_value_CO2 - C_value_CO2,
    delta_T   = T_value_Temp - C_value_Temp
  ) %>%
  mutate(
    Treatment = stringr::str_squish(as.character(Treatment)),
    Treatment = case_when(
      Treatment %in% c("T1", "T2", "T3") ~ "T",
      Treatment %in% c("TpH", "TpH1", "TpH2", "T3pH1") ~ "TpH",
      TRUE ~ Treatment
    ),
    Duplicate = stringr::str_squish(as.character(Duplicate))
  )

cat("\nTreatment values after cleaning:\n")
print(table(zoo_db$Treatment, useNA = "ifany"))

# remove duplicates and odd treatment values
zoo_db <- zoo_db %>%
  filter(is.na(Duplicate) | tolower(Duplicate) != "duplicate") %>%
  filter(
    (Treatment == "pH"  & !is.na(delta_CO2) & delta_CO2 > 0 & delta_CO2 < 5000) |
      (Treatment == "T"   & !is.na(delta_T)   & delta_T > 0 & delta_T <= 10) |
      (Treatment == "TpH" & !is.na(delta_CO2) & !is.na(delta_T) &
         delta_CO2 > 0 & delta_CO2 < 5000 & delta_T > 0 & delta_T <= 10)
  )

cat("\nRows after filtering:", nrow(zoo_db), "\n")
cat("Treatment values after filtering:\n")
print(table(zoo_db$Treatment, useNA = "ifany"))

# cleaned text fields for IDs
zoo_db <- zoo_db %>%
  mutate(
    Author_clean = str_replace_all(
      str_squish(coalesce(as.character(Author), "NA")),
      "[^A-Za-z0-9]+",
      "_"
    ),
    Title_clean = str_replace_all(
      str_squish(coalesce(as.character(Title), "NA")),
      "[^A-Za-z0-9]+",
      "_"
    ),
    Species_clean = str_replace_all(
      str_squish(coalesce(as.character(Species), "NA")),
      "[^A-Za-z0-9]+",
      "_"
    ),
    Group_clean = str_replace_all(
      str_squish(coalesce(as.character(Group), "NA")),
      "[^A-Za-z0-9]+",
      "_"
    ),
    Year_clean = coalesce(as.character(Year), "NA"),
    
    C_value_CO2_num  = suppressWarnings(as.numeric(C_value_CO2)),
    C_value_Temp_num = suppressWarnings(as.numeric(C_value_Temp)),
    N_C_num          = suppressWarnings(as.numeric(N_C))
  )

# Study_ID
zoo_db <- zoo_db %>%
  mutate(
    study_key = paste(
      Author_clean,
      Year_clean,
      Title_clean,
      sep = "__"
    )
  ) %>%
  group_by(study_key) %>%
  mutate(
    Study_ID = paste0("ZOO_STUDY_", sprintf("%04d", cur_group_id()))
  ) %>%
  ungroup()

cat("\nUnique Study_ID:", dplyr::n_distinct(zoo_db$Study_ID), "\n")

# ES_ID
zoo_db <- zoo_db %>%
  mutate(
    ES_ID = paste0("ZOO_ES_", sprintf("%06d", row_number()))
  )

cat("Unique ES_ID:", dplyr::n_distinct(zoo_db$ES_ID), "\n")
cat("ES_ID unique? ", dplyr::n_distinct(zoo_db$ES_ID) == nrow(zoo_db), "\n")

# Common_id
zoo_db <- zoo_db %>%
  mutate(
    common_key = paste(
      Study_ID,
      Species_clean,
      Group_clean,
      paste0(
        "CO2",
        ifelse(
          is.finite(C_value_CO2_num),
          sprintf("%.2f", C_value_CO2_num),
          "NA"
        )
      ),
      paste0(
        "T",
        ifelse(
          is.finite(C_value_Temp_num),
          sprintf("%.2f", C_value_Temp_num),
          "NA"
        )
      ),
      paste0(
        "NC",
        ifelse(
          is.finite(N_C_num),
          sprintf("%d", as.integer(round(N_C_num))),
          "NA"
        )
      ),
      sep = "__"
    )
  ) %>%
  group_by(common_key) %>%
  mutate(
    Common_id = paste0("ZOO_COMMON_", sprintf("%05d", cur_group_id()))
  ) %>%
  ungroup()

cat("Unique Common_id:", dplyr::n_distinct(zoo_db$Common_id), "\n")

# save the wider cleaned database with IDs
write.csv(
  zoo_db,
  file.path(out_dir, "zoo_db_cleaned_with_ids.csv"),
  row.names = FALSE
)

# model-ready dataset
dat1 <- zoo_db %>%
  mutate(
    Include           = str_squish(as.character(Include)),
    Group             = str_squish(as.character(Group)),
    Treatment         = str_squish(as.character(Treatment)),
    Response_category = str_squish(as.character(Response_category)),
    Broad_subgroup    = str_squish(as.character(Broad_subgroup))
  ) %>%
  filter(Biota_Group == "Zooplankton") %>%
  filter(is.na(Include) | !Include %in% c("No", "out")) %>%
  filter(
    is.na(Group) |
      !str_detect(Group, regex("platy|plathy|elasmo", ignore_case = TRUE))
  ) %>%
  mutate(
    Biota_T    = suppressWarnings(as.numeric(Biota_T)),
    Biota_T_sd = suppressWarnings(as.numeric(Biota_T_sd)),
    N_T        = suppressWarnings(as.numeric(N_T)),
    Biota_C    = suppressWarnings(as.numeric(Biota_C)),
    Biota_C_sd = suppressWarnings(as.numeric(Biota_C_sd)),
    N_C        = suppressWarnings(as.numeric(N_C)),
    Common_id  = str_squish(as.character(Common_id)),
    Study_ID   = str_squish(as.character(Study_ID)),
    ES_ID      = str_squish(as.character(ES_ID)),
    delta_CO2  = suppressWarnings(as.numeric(delta_CO2)),
    delta_T    = suppressWarnings(as.numeric(delta_T))
  ) %>%
  filter(
    is.finite(Biota_C),
    is.finite(Biota_T),
    Biota_C > 0,
    Biota_T > 0,
    is.finite(Biota_C_sd),
    is.finite(Biota_T_sd),
    Biota_C_sd >= 0,
    Biota_T_sd >= 0,
    is.finite(N_C),
    is.finite(N_T),
    N_C > 0,
    N_T > 0,
    !is.na(Common_id), Common_id != "",
    !is.na(Study_ID),  Study_ID  != "",
    !is.na(ES_ID),     ES_ID     != ""
  ) %>%
  mutate(
    Treatment = factor(
      Treatment,
      levels = c("pH", "T", "TpH")
    ),
    Broad_subgroup = case_when(
      str_detect(Broad_subgroup, regex("mero", ignore_case = TRUE)) ~ "Meroplankton",
      str_detect(Broad_subgroup, regex("holo", ignore_case = TRUE)) ~ "Holoplankton",
      TRUE ~ NA_character_
    ),
    Broad_subgroup = factor(
      Broad_subgroup,
      levels = c("Meroplankton", "Holoplankton")
    ),
    Response_category = factor(Response_category)
  ) %>%
  filter(!is.na(Treatment))

write.csv(
  dat1,
  file.path(out_dir, "dat1_model_ready_with_ids.csv"),
  row.names = FALSE
)

cat("\nLoaded shared cleaned dataset from 00B_Dependencies.R\n")
cat("Database path:", db_path, "\n")
cat("Rows in dat1:", nrow(dat1), "\n")
cat("Unique ES_ID:", dplyr::n_distinct(dat1$ES_ID), "\n")
cat("Unique Study_ID:", dplyr::n_distinct(dat1$Study_ID), "\n")
cat("Unique Common_id:", dplyr::n_distinct(dat1$Common_id), "\n")
cat("ES_ID unique? ", dplyr::n_distinct(dat1$ES_ID) == nrow(dat1), "\n\n")

cat("Treatment summary:\n")
print(table(dat1$Treatment, useNA = "ifany"))

cat("\nBroad subgroup summary:\n")
print(table(dat1$Broad_subgroup, useNA = "ifany"))

cat("\nResponse category summary:\n")
print(table(dat1$Response_category, useNA = "ifany"))

cat("\nSaved checking files to:", out_dir, "\n")