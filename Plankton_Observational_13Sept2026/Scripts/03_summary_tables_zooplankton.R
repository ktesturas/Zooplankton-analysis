# ===========================================
# 03_summary_tables_zooplankton.R
# Exports summary tables for the three-panel zooplankton figure.
# Uses only the current combined zooplankton workbook through Script 1.
# ===========================================

# Load packages
pacman::p_load(
  dplyr,
  tidyr,
  here,
  openxlsx,
  grid,
  gridExtra
)

# Source the current data preparation 
data_script <- here::here("Scripts", "01_prepare_zooplankton_data.R")
source(data_script)

# Script 1 already applies the shared eligibility rules:
# Include == "Yes" and Duration >= 19.
# No separate All_marine_life_data.xlsx file is used here.

# ---- Output folder ----
script_name <- "03_summary_tables_zooplankton"
out_dir <- here::here("Outputs", script_name)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Shared column order from Script 1 ----
column_lookup <- idx_tbl %>%
  transmute(
    x_order = Xnum,
    Xkey = Ykey,
    Block,
    Taxon = sub("^[^|]+\\s\\|\\|\\s", "", as.character(Ykey))
  )

study_label <- function(block) {
  if_else(block == PREV_BLOCK, "Previous work", "This study")
}

# ===========================================
# Rate summaries
# ===========================================

make_rate_summary <- function(data, panel_name, unit_label) {
  data %>%
    filter(is.finite(Rate), !is.na(Xkey)) %>%
    mutate(Study = study_label(as.character(Block))) %>%
    group_by(Study, Block, Taxon, Xkey) %>%
    summarise(
      n = n(),
      mean_value = mean(Rate),
      sd = sd(Rate),
      median_value = median(Rate),
      min_value = min(Rate),
      max_value = max(Rate),
      .groups = "drop"
    ) %>%
    mutate(
      se = if_else(n > 1, sd / sqrt(n), NA_real_),
      tcrit = if_else(n > 1, qt(0.975, df = n - 1), NA_real_),
      ci_low_raw = mean_value - tcrit * se,
      ci_high_raw = mean_value + tcrit * se,
      ci_low = if_else(n >= 3, ci_low_raw, NA_real_),
      ci_high = if_else(n >= 3, ci_high_raw, NA_real_),
      ci_shown = n >= 3
    ) %>%
    right_join(column_lookup, by = c("Xkey", "Block", "Taxon")) %>%
    mutate(
      Panel = panel_name,
      Metric = if_else(panel_name == "Distribution change", "Distribution", "Phenology"),
      Study = coalesce(Study, study_label(as.character(Block))),
      n = replace_na(n, 0L),
      value_type = "Mean",
      unit = unit_label,
      ci_method = "t-based 95% CI around the mean; shown only when n >= 3",
      denominator_note = NA_character_
    ) %>%
    select(
      Panel, Metric, Study, x_order, Block, Taxon, n,
      value_type, mean_value, ci_low, ci_high, ci_shown,
      unit, ci_method, denominator_note,
      sd, se, median_value, min_value, max_value,
      ci_low_raw, ci_high_raw
    ) %>%
    arrange(x_order)
}

distribution_summary <- make_rate_summary(
  dist_all_raw,
  "Distribution change",
  "km decade^-1"
)

phenology_summary <- make_rate_summary(
  phen_all_raw,
  "Phenology change",
  "days decade^-1"
)

# ===========================================
# Percentage-consistent summaries
# ===========================================

wilson_ci <- function(k, n) {
  z <- qnorm(0.975)
  phat <- k / n
  denominator <- 1 + z^2 / n
  centre <- (phat + z^2 / (2 * n)) / denominator
  half <- z * sqrt(phat * (1 - phat) / n + z^2 / (4 * n^2)) / denominator

  tibble(
    mean_value = phat * 100,
    ci_low = pmax(0, (centre - half) * 100),
    ci_high = pmin(100, (centre + half) * 100)
  )
}

prop_source <- bind_rows(
  phen_prop_raw %>% mutate(Metric = "Phenology"),
  dist_prop_raw %>% mutate(Metric = "Distribution")
) %>%
  mutate(
    Block = as.character(Block),

    # The shared figure stores the pooled Zooplankton column in this
    # internal axis block. Match pooled percentage records to that column.
    Block = if_else(
      Block == "Zooplankton" & Taxon == "Zooplankton",
      "Phytoplankton",
      Block
    ),

    Study = study_label(Block),
    Xkey = factor(paste(Block, Taxon, sep = sep), levels = x_levels_master)
  ) %>%
  filter(!is.na(Xkey))

percentage_summary <- prop_source %>%
  count(Metric, Study, Block, Taxon, Xkey, Consistent, name = "n_category") %>%
  group_by(Metric, Study, Block, Taxon, Xkey) %>%
  summarise(
    n = sum(n_category),
    k_consistent = sum(n_category[Consistent == "Consistent"], na.rm = TRUE),
    k_not_consistent = sum(n_category[Consistent == "Not consistent"], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(interval = list(wilson_ci(k_consistent, n))) %>%
  unnest(interval) %>%
  ungroup() %>%
  right_join(
    crossing(Metric = c("Distribution", "Phenology"), column_lookup),
    by = c("Metric", "Xkey", "Block", "Taxon")
  ) %>%
  mutate(
    Panel = "Percentage consistent",
    Study = coalesce(Study, study_label(as.character(Block))),
    n = replace_na(n, 0L),
    k_consistent = replace_na(k_consistent, 0L),
    k_not_consistent = replace_na(k_not_consistent, 0L),
    value_type = "Percentage",
    unit = "%",
    ci_method = "Wilson 95% binomial CI",
    ci_shown = n > 0,
    denominator_note = "C / (C + NC)"
  ) %>%
  select(
    Panel, Metric, Study, x_order, Block, Taxon, n,
    k_consistent, k_not_consistent,
    value_type, mean_value, ci_low, ci_high, ci_shown,
    unit, ci_method, denominator_note
  ) %>%
  arrange(Metric, x_order)

# ===========================================
# Combined tables
# ===========================================

percentage_for_master <- percentage_summary %>%
  mutate(
    sd = NA_real_,
    se = NA_real_,
    median_value = NA_real_,
    min_value = NA_real_,
    max_value = NA_real_,
    ci_low_raw = ci_low,
    ci_high_raw = ci_high
  ) %>%
  select(
    Panel, Metric, Study, x_order, Block, Taxon, n,
    value_type, mean_value, ci_low, ci_high, ci_shown,
    unit, ci_method, denominator_note,
    sd, se, median_value, min_value, max_value,
    ci_low_raw, ci_high_raw
  )

master_all_columns <- bind_rows(
  distribution_summary,
  phenology_summary,
  percentage_for_master
) %>%
  arrange(Panel, Metric, x_order)

master_plotted_only <- master_all_columns %>%
  filter(n > 0)

supervisor_summary <- master_plotted_only %>%
  mutate(
    value_ci = case_when(
      is.na(mean_value) ~ NA_character_,
      is.na(ci_low) | is.na(ci_high) ~ sprintf("%.2f", mean_value),
      TRUE ~ sprintf(
        "%.2f (95%% CI: %.2f to %.2f)",
        mean_value, ci_low, ci_high
      )
    )
  ) %>%
  select(
    Panel, Metric, Study, Block, Taxon, n,
    value_ci, unit, ci_method, denominator_note
  )

# ===========================================
# Export supervisor summary as a PNG table
# ===========================================

supervisor_summary_display <- supervisor_summary %>%
  mutate(
    across(everything(), as.character),
    across(everything(), ~ replace_na(.x, "—"))
  ) %>%
  rename(
    `Panel` = Panel,
    `Metric` = Metric,
    `Study` = Study,
    `Block` = Block,
    `Taxon` = Taxon,
    `n` = n,
    `Value (95% CI)` = value_ci,
    `Unit` = unit,
    `CI method` = ci_method,
    `Denominator` = denominator_note
  )

# Add one row for the column headings
n_table_rows <- nrow(supervisor_summary_display) + 1L

# Dynamically increase the image height according to the number of rows
table_height_in <- max(
  6,
  0.32 * n_table_rows + 0.8
)

supervisor_table_grob <- gridExtra::tableGrob(
  supervisor_summary_display,
  rows = NULL,
  theme = gridExtra::ttheme_minimal(
    base_size = 8.5,
    base_family = "sans",
    padding = grid::unit(c(1.5, 2), "mm"),
    colhead = list(
      fg_params = list(
        fontface = "bold",
        col = "grey10"
      ),
      bg_params = list(
        fill = "grey90",
        col = NA
      )
    ),
    core = list(
      fg_params = list(
        col = "grey10"
      ),
      bg_params = list(
        fill = "white",
        col = NA
      )
    )
  )
)

supervisor_table_title <- grid::textGrob(
  "Zooplankton observational synthesis: supervisor summary",
  gp = grid::gpar(
    fontsize = 15,
    fontface = "bold",
    col = "grey10",
    fontfamily = "sans"
  )
)

supervisor_table_complete <- gridExtra::arrangeGrob(
  supervisor_table_grob,
  top = supervisor_table_title,
  padding = grid::unit(0.25, "in")
)

supervisor_png <- file.path(
  out_dir,
  "supervisor_summary.png"
)

png(
  filename = supervisor_png,
  width = 18,
  height = table_height_in,
  units = "in",
  res = 300,
  bg = "white"
)

grid::grid.newpage()
grid::grid.draw(supervisor_table_complete)

dev.off()

message("Supervisor summary PNG exported to: ", supervisor_png)

# ===========================================
# Export Excel workbook
# ===========================================

wb <- createWorkbook()

addWorksheet(wb, "README")
writeData(
  wb, "README",
  data.frame(
    Notes = c(
      "Summary tables use only the current combined zooplankton workbook.",
      "The separate All_marine_life_data.xlsx dataset is not used.",
      "Script 1 applies Include = Yes and Duration >= 19 to both sheets.",
      "Rate panels require a numerical Rate.",
      "Percentage consistent includes only C and NC classifications.",
      "Percentage denominator = C + NC.",
      "Rate confidence intervals are displayed only when n >= 3."
    )
  )
)

tables <- list(
  "Master all columns" = master_all_columns,
  "Master plotted only" = master_plotted_only,
  "Panel A Distribution" = distribution_summary,
  "Panel B Phenology" = phenology_summary,
  "Panel C Percent" = percentage_summary,
  "Supervisor summary" = supervisor_summary
)

for (sheet_name in names(tables)) {
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, tables[[sheet_name]])
}

for (sheet_name in names(wb)) {
  setColWidths(wb, sheet_name, cols = 1:50, widths = "auto")
  freezePane(wb, sheet_name, firstRow = TRUE)
}

excel_file <- file.path(out_dir, "observational_figure_summary_tables.xlsx")
saveWorkbook(wb, excel_file, overwrite = TRUE)

message("Summary workbook exported to: ", excel_file)
