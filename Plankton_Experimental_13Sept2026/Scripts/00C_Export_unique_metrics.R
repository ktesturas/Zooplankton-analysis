# ============================================================
# 00C_Export_unique_metrics.R
# Purpose:
#   Use the shared cleaned zooplankton dataset from
#   Scripts/00B_Dependencies.R and export unique Metric values
#   grouped by Response_category and Metric_Type.
#
# Outputs:
#   Outputs/Tables/00C_Export_unique_metrics/
#     - unique_metrics_by_category.xlsx
#
#   Outputs/RDS/00C_Export_unique_metrics/
#     - unique_metrics_all.rds
#     - sheet1_by_response_category.rds
#     - sheet2_by_response_category_and_metric_type.rds
# ============================================================

suppressPackageStartupMessages(invisible(
  lapply(
    c(
      "dplyr", "stringr", "tidyr", "tibble", "openxlsx"
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
# Script name and output folders
# ============================================================

script_name <- "00C_Export_unique_metrics"

table_dir <- file.path("Outputs", "Tables", script_name)
rds_dir   <- file.path("Outputs", "RDS", script_name)

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# Check shared cleaned dataset
# ============================================================

if (!exists("dat1")) {
  stop("dat1 not found. Check Scripts/00B_Dependencies.R.")
}

required_cols <- c("Response_category", "Metric", "Metric_Type")

missing_cols <- setdiff(required_cols, names(dat1))

if (length(missing_cols) > 0) {
  stop(
    "Missing required column(s): ",
    paste(missing_cols, collapse = ", ")
  )
}

cat("\nLoaded dat1 from 00B_Dependencies.R\n")
cat("Rows in dat1:", nrow(dat1), "\n")

# ============================================================
# Clean and extract unique Metric values
# ============================================================

unique_metrics_all <- dat1 %>%
  mutate(
    Response_category = str_squish(as.character(Response_category)),
    Metric            = str_squish(as.character(Metric)),
    Metric_Type       = str_squish(as.character(Metric_Type)),
    Metric_Type       = case_when(
      str_detect(Metric_Type, regex("^pos", ignore_case = TRUE)) ~ "Positive",
      str_detect(Metric_Type, regex("^neg", ignore_case = TRUE)) ~ "Negative",
      str_detect(Metric_Type, regex("^amb", ignore_case = TRUE)) ~ "Ambiguous",
      TRUE ~ Metric_Type
    )
  ) %>%
  filter(
    !is.na(Response_category), Response_category != "",
    !is.na(Metric), Metric != ""
  ) %>%
  distinct(Response_category, Metric, Metric_Type) %>%
  arrange(Response_category, Metric_Type, Metric)

cat("\nUnique Metric values extracted:", nrow(unique_metrics_all), "\n")

# ============================================================
# Helper function: make wide list
# ============================================================

make_wide_list <- function(data, group_col, value_col) {
  
  split_list <- data %>%
    group_by(.data[[group_col]]) %>%
    summarise(
      values = list(sort(unique(.data[[value_col]]))),
      .groups = "drop"
    )
  
  max_n <- max(lengths(split_list$values))
  
  wide_df <- lapply(split_list$values, function(x) {
    length(x) <- max_n
    x
  }) %>%
    as_tibble(.name_repair = "minimal")
  
  names(wide_df) <- split_list[[group_col]]
  
  wide_df
}

# ============================================================
# Sheet 1: Metrics grouped by Response_category
# ============================================================

sheet1_data <- unique_metrics_all %>%
  select(Response_category, Metric) %>%
  distinct(Response_category, Metric)

sheet1_wide <- make_wide_list(
  data      = sheet1_data,
  group_col = "Response_category",
  value_col = "Metric"
)

# ============================================================
# Sheet 2: Metrics grouped by Response_category and Metric_Type
# ============================================================

sheet2_data <- unique_metrics_all %>%
  filter(
    !is.na(Metric_Type), Metric_Type != "",
    Metric_Type %in% c("Positive", "Negative", "Ambiguous")
  ) %>%
  select(Response_category, Metric_Type, Metric) %>%
  distinct(Response_category, Metric_Type, Metric)

response_categories <- sort(unique(sheet2_data$Response_category))
metric_types <- c("Positive", "Negative", "Ambiguous")

sheet2_cols <- list()

for (rc in response_categories) {
  for (mt in metric_types) {
    
    col_name <- paste(rc, mt, sep = "__")
    
    sheet2_cols[[col_name]] <- sheet2_data %>%
      filter(Response_category == rc, Metric_Type == mt) %>%
      arrange(Metric) %>%
      pull(Metric) %>%
      unique()
  }
}

max_n <- max(lengths(sheet2_cols))

sheet2_wide <- lapply(sheet2_cols, function(x) {
  length(x) <- max_n
  x
}) %>%
  as_tibble(.name_repair = "minimal")

names(sheet2_wide) <- names(sheet2_cols)

# ============================================================
# Save RDS objects
# ============================================================

saveRDS(
  unique_metrics_all,
  file.path(rds_dir, "unique_metrics_all.rds")
)

saveRDS(
  sheet1_wide,
  file.path(rds_dir, "sheet1_by_response_category.rds")
)

saveRDS(
  sheet2_wide,
  file.path(rds_dir, "sheet2_by_response_category_and_metric_type.rds")
)

cat("\nSaved RDS objects to:\n")
cat(rds_dir, "\n")

# ============================================================
# Export to Excel
# ============================================================

output_file <- file.path(table_dir, "unique_metrics_by_category.xlsx")

wb <- createWorkbook()

# -----------------------------
# Sheet 1
# -----------------------------

addWorksheet(wb, "By response category")

writeData(
  wb,
  sheet = "By response category",
  x = sheet1_wide,
  startRow = 1,
  startCol = 1
)

# -----------------------------
# Sheet 2
# -----------------------------

addWorksheet(wb, "By category and type")

category_headers <- str_replace(names(sheet2_wide), "__.*$", "")
type_headers     <- str_replace(names(sheet2_wide), "^.*__", "")

writeData(
  wb,
  sheet = "By category and type",
  x = t(category_headers),
  startRow = 1,
  startCol = 1,
  colNames = FALSE
)

writeData(
  wb,
  sheet = "By category and type",
  x = t(type_headers),
  startRow = 2,
  startCol = 1,
  colNames = FALSE
)

writeData(
  wb,
  sheet = "By category and type",
  x = sheet2_wide,
  startRow = 3,
  startCol = 1,
  colNames = FALSE
)

# Merge response-category headers across Positive, Negative, Ambiguous
start_col <- 1

for (rc in response_categories) {
  mergeCells(
    wb,
    sheet = "By category and type",
    cols = start_col:(start_col + 2),
    rows = 1
  )
  
  start_col <- start_col + 3
}

# ============================================================
# Formatting
# ============================================================

header_style <- createStyle(
  textDecoration = "bold",
  halign = "center",
  valign = "center",
  border = "Bottom"
)

subheader_style <- createStyle(
  textDecoration = "bold",
  halign = "center",
  valign = "center"
)

addStyle(
  wb,
  sheet = "By response category",
  style = header_style,
  rows = 1,
  cols = 1:ncol(sheet1_wide),
  gridExpand = TRUE
)

addStyle(
  wb,
  sheet = "By category and type",
  style = header_style,
  rows = 1,
  cols = 1:ncol(sheet2_wide),
  gridExpand = TRUE
)

addStyle(
  wb,
  sheet = "By category and type",
  style = subheader_style,
  rows = 2,
  cols = 1:ncol(sheet2_wide),
  gridExpand = TRUE
)

setColWidths(
  wb,
  sheet = "By response category",
  cols = 1:ncol(sheet1_wide),
  widths = 35
)

setColWidths(
  wb,
  sheet = "By category and type",
  cols = 1:ncol(sheet2_wide),
  widths = 35
)

freezePane(
  wb,
  sheet = "By response category",
  firstRow = TRUE
)

freezePane(
  wb,
  sheet = "By category and type",
  firstActiveRow = 3
)

# ============================================================
# Save workbook
# ============================================================

saveWorkbook(wb, output_file, overwrite = TRUE)

cat("\nSaved Excel file to:\n")
cat(output_file, "\n")

# ============================================================
# Final summary
# ============================================================

cat("\nResponse categories in Sheet 1:\n")
print(names(sheet1_wide))

cat("\nResponse categories and Metric_Type columns in Sheet 2:\n")
print(names(sheet2_wide))