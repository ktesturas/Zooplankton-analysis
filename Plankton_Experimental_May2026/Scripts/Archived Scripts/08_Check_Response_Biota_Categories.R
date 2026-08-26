# ============================================================
# 08_Check_Response_Biota_Categories.R
# Purpose:
#   Create a simple Excel workbook for checking the direction
#   of Response_Biota metrics under each lumped Response_category.
#
#   Output:
#   One Excel workbook with 5 sheets:
#     - Growth
#     - Development
#     - Physiology
#     - Reproduction
#     - Survival
#
#   Each sheet contains:
#     Response_Biota | n_records | Direction | Notes
#
#   Direction is intentionally blank so you can manually label:
#     positive_when_higher
#     negative_when_higher
#     check_context
#
# Input:
#   dat1 object from Scripts/00_Dependencies.R
#
# Outputs:
#   Outputs/Tables/08_Check_Response_Biota_Categories/
# ============================================================

suppressPackageStartupMessages(invisible(
  lapply(
    c(
      "dplyr", "tibble", "openxlsx"
    ),
    library,
    character.only = TRUE
  )
))

# ============================================================
# Project dependencies
# ============================================================

source(file.path("Scripts", "00B_Dependencies.R"))

# ============================================================
# Script name and output folder
# ============================================================

script_name <- "08_Check_Response_Biota_Categories"

table_dir <- file.path("Outputs", "Tables", script_name)

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# Check shared cleaned dataset
# ============================================================

if (!exists("dat1")) {
  stop("dat1 not found. Check that Scripts/00_Dependencies.R creates dat1.")
}

required_cols <- c("Response_Biota", "Response_category")

missing_cols <- setdiff(required_cols, names(dat1))

if (length(missing_cols) > 0) {
  stop(
    "The following required columns are missing from dat1: ",
    paste(missing_cols, collapse = ", ")
  )
}

cat("\n==== Checking Response_Biota and Response_category ====\n")
cat("Rows in dat1:", nrow(dat1), "\n")
cat("Unique Response_Biota values:", dplyr::n_distinct(dat1$Response_Biota), "\n")
cat("Unique original Response_category values:", dplyr::n_distinct(dat1$Response_category), "\n\n")

# ============================================================
# Create checked dataset
# ============================================================

dat_check <- dat1 %>%
  mutate(
    Response_Biota = as.character(Response_Biota),
    Response_category_original = as.character(Response_category),
    
    # Lump Behaviour with Physiology
    Response_category_lumped = dplyr::case_when(
      Response_category_original %in% c("Behaviour", "Physiology") ~ "Physiology",
      TRUE ~ Response_category_original
    ),
    
    Response_category_lumped = factor(
      Response_category_lumped,
      levels = c(
        "Growth",
        "Development",
        "Physiology",
        "Reproduction",
        "Survival"
      )
    )
  ) %>%
  filter(!is.na(Response_Biota), !is.na(Response_category_lumped))

# ============================================================
# Create one checking table per lumped response category
# ============================================================

category_levels <- c(
  "Growth",
  "Development",
  "Physiology",
  "Reproduction",
  "Survival"
)

category_tables <- lapply(category_levels, function(cat_name) {
  
  dat_check %>%
    filter(Response_category_lumped == cat_name) %>%
    count(Response_Biota, name = "n_records") %>%
    arrange(Response_Biota) %>%
    mutate(
      Direction = "",
      Notes = ""
    ) %>%
    select(
      Response_Biota,
      n_records,
      Direction,
      Notes
    )
})

names(category_tables) <- category_levels

# ============================================================
# Save Excel workbook with 5 sheets
# ============================================================

output_xlsx <- file.path(
  table_dir,
  "response_biota_direction_check_by_category.xlsx"
)

wb <- openxlsx::createWorkbook()

header_style <- openxlsx::createStyle(
  textDecoration = "bold",
  halign = "center",
  valign = "center",
  border = "Bottom"
)

for (sheet_name in names(category_tables)) {
  
  openxlsx::addWorksheet(wb, sheet_name)
  
  openxlsx::writeData(
    wb,
    sheet = sheet_name,
    x = category_tables[[sheet_name]],
    startRow = 1,
    startCol = 1,
    headerStyle = header_style
  )
  
  openxlsx::freezePane(
    wb,
    sheet = sheet_name,
    firstActiveRow = 2
  )
  
  openxlsx::setColWidths(
    wb,
    sheet = sheet_name,
    cols = 1:ncol(category_tables[[sheet_name]]),
    widths = "auto"
  )
}

openxlsx::saveWorkbook(
  wb,
  output_xlsx,
  overwrite = TRUE
)

# ============================================================
# Optional CSV copies per category
# ============================================================

for (sheet_name in names(category_tables)) {
  
  csv_name <- paste0(
    "response_biota_direction_check_",
    stringr::str_to_lower(sheet_name),
    ".csv"
  )
  
  write.csv(
    category_tables[[sheet_name]],
    file.path(table_dir, csv_name),
    row.names = FALSE
  )
}

# ============================================================
# Console summary
# ============================================================

cat("\n==== Response_Biota direction checking tables ====\n")

for (sheet_name in names(category_tables)) {
  cat(
    sheet_name, ": ",
    nrow(category_tables[[sheet_name]]),
    " unique Response_Biota values\n",
    sep = ""
  )
}

cat("\n==== Saved outputs ====\n")
cat("Tables folder:", table_dir, "\n")
cat("Excel workbook:", output_xlsx, "\n")