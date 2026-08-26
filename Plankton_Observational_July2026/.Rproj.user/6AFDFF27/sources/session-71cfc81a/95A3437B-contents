library(dplyr)
library(readxl)
library(stringr)
library(here)

# ==========================================
# 04_database_check_zooplankton.R
# Database checks for studies, taxa, species, and observation spans
# ==========================================

# ---- Settings ----

script_name <- "04_database_check_zooplankton"

fp_zoo_candidates <- c(
  here::here("Data", "Observational", "Plankton observational.xlsx"),
  here::here("Data", "Plankton observational.xlsx"),
  here::here("Observational", "Plankton observational.xlsx")
)

fp_zoo <- fp_zoo_candidates[file.exists(fp_zoo_candidates)][1]

if (is.na(fp_zoo)) {
  stop(
    paste0(
      "Cannot find Plankton observational.xlsx. Checked:\n",
      paste(fp_zoo_candidates, collapse = "\n")
    )
  )
}

out_dir <- here::here("Outputs", script_name)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

min_span <- 19

apply_span_filter_distribution <- TRUE
apply_span_filter_phenology <- FALSE

# ---- Possible column names ----

possible_study_cols <- c(
  "Study_ID", "Study ID", "StudyID",
  "Study", "Reference", "Citation",
  "Paper", "Author_year", "Author Year"
)

possible_taxa_cols <- c(
  "Scientific_name", "Scientific name", "Scientific.Name",
  "Species", "species",
  "Taxon", "Taxa", "taxon", "taxa",
  "Organism", "Name"
)

# ---- Helper functions ----

get_col <- function(df, possible_cols, label, sheet_name) {
  col <- possible_cols[possible_cols %in% names(df)][1]
  
  if (is.na(col)) {
    stop(
      paste0(
        "No ", label, " column found in the ", sheet_name, " sheet.\n",
        "Check column names using:\n",
        "names(read_xlsx(fp_zoo, sheet = '", sheet_name, "'))"
      )
    )
  }
  
  col
}

clean_text <- function(x) {
  x %>%
    as.character() %>%
    str_squish() %>%
    str_replace_all("_", " ") %>%
    na_if("") %>%
    na_if("NA")
}

extract_genus <- function(x) {
  x <- clean_text(x)
  str_extract(x, "^[A-Z][a-zA-Z-]+")
}

extract_species_binomial <- function(x) {
  x <- clean_text(x)
  
  binomial <- str_extract(
    x,
    "^[A-Z][a-zA-Z-]+\\s+[a-z][a-zA-Z-]+"
  )
  
  bad_species <- str_detect(
    binomial,
    regex("\\b(sp|spp|cf|aff|indet|complex|group)\\b", ignore_case = TRUE)
  )
  
  binomial[bad_species] <- NA_character_
  
  binomial
}

filter_include <- function(df) {
  if ("Include" %in% names(df)) {
    df %>%
      filter(
        is.na(Include) |
          !toupper(str_squish(as.character(Include))) %in% c("NO", "N")
      )
  } else {
    df
  }
}

get_span_col <- function(sheet_name) {
  if (sheet_name == "Phenology") {
    return("Duration")
  }
  
  if (sheet_name == "Distribution") {
    return("Timespan")
  }
  
  stop("Unknown sheet name.")
}

read_database_sheet <- function(sheet_name) {
  df <- read_xlsx(fp_zoo, sheet = sheet_name)
  
  study_col <- get_col(df, possible_study_cols, "study", sheet_name)
  taxa_col  <- get_col(df, possible_taxa_cols, "taxa/species", sheet_name)
  span_col  <- get_span_col(sheet_name)
  
  if (!span_col %in% names(df)) {
    stop(paste0("Column '", span_col, "' not found in the ", sheet_name, " sheet."))
  }
  
  if (!"Rate" %in% names(df)) {
    stop(paste0("Column 'Rate' not found in the ", sheet_name, " sheet."))
  }
  
  df %>%
    filter_include() %>%
    transmute(
      Panel = sheet_name,
      Study = clean_text(.data[[study_col]]),
      Taxon_recorded = clean_text(.data[[taxa_col]]),
      Genus = extract_genus(.data[[taxa_col]]),
      Species_binomial = extract_species_binomial(.data[[taxa_col]]),
      Rate = suppressWarnings(as.numeric(Rate)),
      Span = suppressWarnings(as.numeric(.data[[span_col]]))
    )
}

# ---- Read data ----

phen_data <- read_database_sheet("Phenology")
dist_data <- read_database_sheet("Distribution")

all_data <- bind_rows(phen_data, dist_data)

# ---- Apply validity and span rules ----

all_obs <- all_data %>%
  mutate(
    valid_rate = !is.na(Rate),
    valid_span = !is.na(Span),
    passes_span_filter = case_when(
      Panel == "Distribution" & apply_span_filter_distribution ~ Span >= min_span,
      Panel == "Phenology" & apply_span_filter_phenology ~ Span >= min_span,
      TRUE ~ TRUE
    )
  ) %>%
  filter(
    valid_rate,
    valid_span,
    passes_span_filter
  )

phen_obs <- all_obs %>%
  filter(Panel == "Phenology")

dist_obs <- all_obs %>%
  filter(Panel == "Distribution")

# ---- Unique studies ----

unique_studies <- all_obs %>%
  filter(!is.na(Study)) %>%
  distinct(Study) %>%
  arrange(Study)

# ---- Unique taxa ----

unique_recorded_taxa <- all_obs %>%
  filter(!is.na(Taxon_recorded)) %>%
  distinct(Taxon_recorded) %>%
  arrange(Taxon_recorded)

unique_genera <- all_obs %>%
  filter(!is.na(Genus)) %>%
  distinct(Genus) %>%
  arrange(Genus)

unique_species <- all_obs %>%
  filter(!is.na(Species_binomial)) %>%
  distinct(Species_binomial) %>%
  arrange(Species_binomial)

# ---- Summary tables ----

study_summary <- tibble(
  Category = c(
    "Unique studies in Phenology",
    "Unique studies in Distribution",
    "Unique studies across Phenology + Distribution"
  ),
  Count = c(
    n_distinct(phen_obs$Study, na.rm = TRUE),
    n_distinct(dist_obs$Study, na.rm = TRUE),
    nrow(unique_studies)
  )
)

taxa_summary <- tibble(
  Category = c(
    "Unique recorded taxa in Phenology",
    "Unique recorded taxa in Distribution",
    "Unique recorded taxa across Phenology + Distribution",
    "Unique genera across Phenology + Distribution",
    "Unique species-level taxa across Phenology + Distribution"
  ),
  Count = c(
    n_distinct(phen_obs$Taxon_recorded, na.rm = TRUE),
    n_distinct(dist_obs$Taxon_recorded, na.rm = TRUE),
    nrow(unique_recorded_taxa),
    nrow(unique_genera),
    nrow(unique_species)
  )
)

observation_summary <- all_obs %>%
  count(Panel, name = "n_observations")

span_summary_by_panel <- all_obs %>%
  group_by(Panel) %>%
  summarise(
    n_observations = n(),
    median_span = median(Span, na.rm = TRUE),
    min_span = min(Span, na.rm = TRUE),
    max_span = max(Span, na.rm = TRUE),
    .groups = "drop"
  )

span_summary_overall <- all_obs %>%
  summarise(
    n_observations = n(),
    median_span = median(Span, na.rm = TRUE),
    min_span = min(Span, na.rm = TRUE),
    max_span = max(Span, na.rm = TRUE)
  )

n_obs <- nrow(all_obs)
median_span <- median(all_obs$Span, na.rm = TRUE)
min_span <- min(all_obs$Span, na.rm = TRUE)
max_span <- max(all_obs$Span, na.rm = TRUE)

verification_sentence <- paste0(
  "From these, we compiled ",
  format(n_obs, big.mark = ","),
  " observations (median span ",
  round(median_span),
  " yr, range ",
  round(min_span),
  "–",
  round(max_span),
  " yr)."
)

# ---- Export checks ----

write.csv(
  all_obs,
  file.path(out_dir, "database_check_rows_used.csv"),
  row.names = FALSE
)

write.csv(
  study_summary,
  file.path(out_dir, "study_summary.csv"),
  row.names = FALSE
)

write.csv(
  taxa_summary,
  file.path(out_dir, "taxa_summary.csv"),
  row.names = FALSE
)

write.csv(
  unique_studies,
  file.path(out_dir, "unique_studies.csv"),
  row.names = FALSE
)

write.csv(
  unique_recorded_taxa,
  file.path(out_dir, "unique_recorded_taxa.csv"),
  row.names = FALSE
)

write.csv(
  unique_genera,
  file.path(out_dir, "unique_genera.csv"),
  row.names = FALSE
)

write.csv(
  unique_species,
  file.path(out_dir, "unique_species_level_taxa.csv"),
  row.names = FALSE
)

write.csv(
  observation_summary,
  file.path(out_dir, "observation_summary_by_panel.csv"),
  row.names = FALSE
)

write.csv(
  span_summary_by_panel,
  file.path(out_dir, "span_summary_by_panel.csv"),
  row.names = FALSE
)

write.csv(
  span_summary_overall,
  file.path(out_dir, "span_summary_overall.csv"),
  row.names = FALSE
)

writeLines(
  verification_sentence,
  file.path(out_dir, "verification_sentence.txt")
)

# ---- Print results ----

cat("========================================\n")
cat("DATABASE CHECK SETTINGS\n")
cat("========================================\n")

cat("Input file:\n")
cat(fp_zoo, "\n\n")

cat("Output folder:\n")
cat(out_dir, "\n\n")

cat("Minimum span threshold:", min_span, "years\n")
cat("Distribution span filter applied:", apply_span_filter_distribution, "\n")
cat("Phenology span filter applied:", apply_span_filter_phenology, "\n\n")

cat("========================================\n")
cat("UNIQUE STUDIES\n")
cat("========================================\n")

cat("Unique studies in Phenology:", n_distinct(phen_obs$Study, na.rm = TRUE), "\n")
cat("Unique studies in Distribution:", n_distinct(dist_obs$Study, na.rm = TRUE), "\n")
cat("Unique studies across Phenology + Distribution:", nrow(unique_studies), "\n\n")

cat("========================================\n")
cat("UNIQUE TAXA\n")
cat("========================================\n")

cat("Unique recorded taxa in Phenology:", n_distinct(phen_obs$Taxon_recorded, na.rm = TRUE), "\n")
cat("Unique recorded taxa in Distribution:", n_distinct(dist_obs$Taxon_recorded, na.rm = TRUE), "\n")
cat("Unique recorded taxa across Phenology + Distribution:", nrow(unique_recorded_taxa), "\n\n")

cat("Unique genera across Phenology + Distribution:", nrow(unique_genera), "\n")
cat("Unique species-level taxa across Phenology + Distribution:", nrow(unique_species), "\n\n")

cat("========================================\n")
cat("UNIQUE SPECIES LIST\n")
cat("========================================\n")

cat("Unique species-level taxa, alphabetically arranged:\n")
cat("--------------------------------------------------\n")

if (nrow(unique_species) == 0) {
  cat("No species-level taxa detected.\n")
} else {
  cat(
    paste0(
      seq_len(nrow(unique_species)),
      ". ",
      unique_species$Species_binomial,
      collapse = "\n"
    )
  )
  cat("\n")
}

cat("\n")

cat("========================================\n")
cat("OBSERVATION AND SPAN VERIFICATION\n")
cat("========================================\n")

cat("Observation summary by panel:\n")
print(observation_summary)

cat("\nSpan summary by panel:\n")
print(span_summary_by_panel)

cat("\nOverall span summary:\n")
print(span_summary_overall)

cat("\nVerification sentence:\n")
cat("----------------------\n")
cat(verification_sentence)
cat("\n\n")

cat("========================================\n")
cat("FILES SAVED\n")
cat("========================================\n")

cat(file.path(out_dir, "database_check_rows_used.csv"), "\n")
cat(file.path(out_dir, "study_summary.csv"), "\n")
cat(file.path(out_dir, "taxa_summary.csv"), "\n")
cat(file.path(out_dir, "unique_studies.csv"), "\n")
cat(file.path(out_dir, "unique_recorded_taxa.csv"), "\n")
cat(file.path(out_dir, "unique_genera.csv"), "\n")
cat(file.path(out_dir, "unique_species_level_taxa.csv"), "\n")
cat(file.path(out_dir, "observation_summary_by_panel.csv"), "\n")
cat(file.path(out_dir, "span_summary_by_panel.csv"), "\n")
cat(file.path(out_dir, "span_summary_overall.csv"), "\n")
cat(file.path(out_dir, "verification_sentence.txt"), "\n")