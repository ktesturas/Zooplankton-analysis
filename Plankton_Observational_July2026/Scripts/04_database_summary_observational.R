# ===========================================
# 04_database_summary_observational.R
# Summary of observations used in the
# bottom percentage-consistent panel
# ===========================================

pacman::p_load(dplyr, stringr, openxlsx, here)

# ===========================================
# Source final prepared data from Script 1
# ===========================================

source(
  here::here(
    "Scripts",
    "01_prepare_zooplankton_data.R"
  )
)

# Script 1 already applies:
# Include == "Yes"
# Duration >= 19 years


# ===========================================
# Select observations used in bottom panel
# ===========================================

phen_bottom <- phen_source %>%
  filter(
    Group %in% ZOOP_GROUPS,
    !Subgroup %in% EXCLUDED_SUBGROUPS,
    Consistent %in% c("C", "NC")
  ) %>%
  transmute(
    Panel     = "Phenology",
    Reference = as.character(Reference),
    Species   = as.character(Species),
    Duration  = as.numeric(Duration),
    FirstYear = as.numeric(FirstYear),
    LastYear  = as.numeric(LastYear),
    Consistent = as.character(Consistent)
  )

dist_bottom <- dist_source %>%
  filter(
    Group %in% ZOOP_GROUPS,
    !Subgroup %in% EXCLUDED_SUBGROUPS,
    Consistent %in% c("C", "NC")
  ) %>%
  transmute(
    Panel     = "Distribution",
    Reference = as.character(Reference),
    Species   = as.character(Species),
    Duration  = as.numeric(Duration),
    FirstYear = as.numeric(FirstYear),
    LastYear  = as.numeric(LastYear),
    Consistent = as.character(Consistent)
  )


# ===========================================
# Combine panels and identify taxonomic level
# ===========================================

bottom_data <- bind_rows(
  phen_bottom,
  dist_bottom
) %>%
  mutate(
    Reference = str_squish(Reference),
    Species   = str_squish(Species),
    
    Reference = na_if(Reference, ""),
    Species   = na_if(Species, ""),
    
    Genus = str_extract(
      Species,
      "^[A-Z][a-zA-Z-]+"
    ),
    
    Species_binomial = str_extract(
      Species,
      "^[A-Z][a-zA-Z-]+\\s+[a-z][a-zA-Z-]+"
    ),
    
    Species_binomial = if_else(
      str_detect(
        Species_binomial,
        regex(
          "\\b(sp|spp|cf|aff|indet|complex|group)\\b",
          ignore_case = TRUE
        )
      ),
      NA_character_,
      Species_binomial
    )
  )


# ===========================================
# Summary function
# ===========================================

make_summary <- function(data, panel_name) {
  
  n_consistent <- sum(
    data$Consistent == "C",
    na.rm = TRUE
  )
  
  n_not_consistent <- sum(
    data$Consistent == "NC",
    na.rm = TRUE
  )
  
  denominator <- n_consistent + n_not_consistent
  
  tibble(
    Panel = panel_name,
    
    n_observations = nrow(data),
    
    n_unique_studies = n_distinct(
      data$Reference,
      na.rm = TRUE
    ),
    
    n_unique_recorded_taxa = n_distinct(
      data$Species,
      na.rm = TRUE
    ),
    
    n_unique_genera = n_distinct(
      data$Genus,
      na.rm = TRUE
    ),
    
    n_unique_species_level_taxa = n_distinct(
      data$Species_binomial,
      na.rm = TRUE
    ),
    
    n_consistent = n_consistent,
    
    n_not_consistent = n_not_consistent,
    
    percent_consistent = (
      n_consistent / denominator
    ) * 100,
    
    mean_duration_years = mean(
      data$Duration,
      na.rm = TRUE
    ),
    
    median_duration_years = median(
      data$Duration,
      na.rm = TRUE
    ),
    
    duration_q1_years = quantile(
      data$Duration,
      probs = 0.25,
      na.rm = TRUE,
      names = FALSE
    ),
    
    duration_q3_years = quantile(
      data$Duration,
      probs = 0.75,
      na.rm = TRUE,
      names = FALSE
    ),
    
    minimum_duration_years = min(
      data$Duration,
      na.rm = TRUE
    ),
    
    maximum_duration_years = max(
      data$Duration,
      na.rm = TRUE
    ),
    
    earliest_first_year = min(
      data$FirstYear,
      na.rm = TRUE
    ),
    
    latest_last_year = max(
      data$LastYear,
      na.rm = TRUE
    )
  )
}


# ===========================================
# Create summary table
# ===========================================

database_summary <- bind_rows(
  make_summary(
    bottom_data %>%
      filter(Panel == "Phenology"),
    "Phenology"
  ),
  
  make_summary(
    bottom_data %>%
      filter(Panel == "Distribution"),
    "Distribution"
  ),
  
  make_summary(
    bottom_data,
    "Combined"
  )
) %>%
  mutate(
    across(
      c(
        percent_consistent,
        mean_duration_years,
        median_duration_years,
        duration_q1_years,
        duration_q3_years,
        minimum_duration_years,
        maximum_duration_years
      ),
      ~ round(.x, 1)
    )
  )


# ===========================================
# Export one Excel workbook
# ===========================================

out_dir <- here::here(
  "Outputs",
  "04_database_summary_observational"
)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

output_file <- file.path(
  out_dir,
  "zooplankton_database_summary.xlsx"
)

wb <- createWorkbook()

addWorksheet(
  wb,
  "Database summary"
)

writeData(
  wb,
  sheet = "Database summary",
  x = database_summary
)

header_style <- createStyle(
  textDecoration = "bold",
  fgFill = "#D9EAF7",
  border = "Bottom"
)

addStyle(
  wb,
  sheet = "Database summary",
  style = header_style,
  rows = 1,
  cols = 1:ncol(database_summary),
  gridExpand = TRUE
)

freezePane(
  wb,
  sheet = "Database summary",
  firstRow = TRUE
)

setColWidths(
  wb,
  sheet = "Database summary",
  cols = 1:ncol(database_summary),
  widths = "auto"
)

saveWorkbook(
  wb,
  output_file,
  overwrite = TRUE
)

message(
  "Database summary exported to: ",
  output_file
)

print(database_summary)