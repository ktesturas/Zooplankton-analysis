# ============================================================
# 02D_Plot_lnrr_calcifier.R
# Purpose:
#   Plot calcifier vs non-calcifier pooled lnRR estimates
#   using saved outputs from 02C_Model_lnrr_calcifier.R.
#
# Important:
#   This script does NOT run rma.mv().
#   Use this when only editing the plot.
#
# Input:
#   Outputs/RDS/02C_Model_lnrr_calcifier/
#     - dat1_calcifier_model_ready.rds
#     - pooled_estimates_calcifier_one_panel.rds
#     - plot_levels_calcifier.rds
#
# Outputs:
#   Outputs/Figures/02D_Plot_lnrr_calcifier/
#   Outputs/Tables/02D_Plot_lnrr_calcifier/
#   Outputs/RDS/02D_Plot_lnrr_calcifier/
# ============================================================

suppressPackageStartupMessages(invisible(
  lapply(
    c(
      "dplyr", "ggplot2", "readr", "stringr",
      "showtext", "tibble", "grid", "gridExtra",
      "tidyr", "purrr", "scales"
    ),
    library,
    character.only = TRUE
  )
))

# ============================================================
# Fonts and theme
# ============================================================

font_add_google("Inter", "helv")
showtext_auto()
fam <- "helv"

theme_set(theme_minimal(base_size = 18, base_family = fam))
theme_update(text = element_text(family = fam))

# ============================================================
# Script name and output folders
# ============================================================

model_script_name <- "02C_Model_lnrr_calcifier"
plot_script_name  <- "02D_Plot_lnrr_calcifier"

model_rds_dir <- file.path("Outputs", "RDS", model_script_name)

figure_dir <- file.path("Outputs", "Figures", plot_script_name)
table_dir  <- file.path("Outputs", "Tables", plot_script_name)
rds_dir    <- file.path("Outputs", "RDS", plot_script_name)

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# Load saved outputs from 02C
# ============================================================

dat1_calc <- readRDS(
  file.path(model_rds_dir, "dat1_calcifier_model_ready.rds")
)

sum_df_calcifier <- readRDS(
  file.path(model_rds_dir, "pooled_estimates_calcifier_one_panel.rds")
)

plot_levels_calcifier <- readRDS(
  file.path(model_rds_dir, "plot_levels_calcifier.rds")
)

response_levels <- plot_levels_calcifier$response_levels
category_levels <- plot_levels_calcifier$category_levels
treatment_levels <- plot_levels_calcifier$treatment_levels
treatment_label_levels <- plot_levels_calcifier$treatment_label_levels
calcifier_levels <- plot_levels_calcifier$calcifier_levels

cat("\nLoaded saved model outputs from:", model_rds_dir, "\n")
cat("Rows in dat1_calc:", nrow(dat1_calc), "\n")
cat("Rows in sum_df_calcifier:", nrow(sum_df_calcifier), "\n\n")

# ============================================================
# Make sure factors are ordered correctly
# ============================================================

dat1_calc <- dat1_calc %>%
  mutate(
    Response_category = factor(Response_category, levels = response_levels),
    Treatment_label = factor(Treatment_label, levels = treatment_label_levels),
    Calcifier = factor(Calcifier, levels = calcifier_levels)
  )

sum_df_calcifier <- sum_df_calcifier %>%
  mutate(
    Category = factor(Category, levels = category_levels),
    Treatment_label = factor(Treatment_label, levels = treatment_label_levels),
    Calcifier = factor(Calcifier, levels = calcifier_levels)
  )

# ============================================================
# Save pooled means, confidence intervals, percent change, and p-values as a table figure
# ============================================================

pooled_table_plot <- sum_df_calcifier %>%
  mutate(
    Category = as.character(Category),
    Calcifier = as.character(Calcifier),
    Treatment = as.character(Treatment_label),
    
    # lnRR values
    `Mean lnRR` = sprintf("%.3f", mean),
    `95% CI` = sprintf("%.3f to %.3f", ci_low, ci_high),
    
    # p-values for pooled mean lnRR being different from 0
    `p-value` = case_when(
      is.na(pval) ~ NA_character_,
      pval < 0.001 ~ "<0.001",
      TRUE ~ sprintf("%.3f", pval)
    ),
    
    # Percent change values converted from lnRR
    `% change` = sprintf("%.1f%%", (exp(mean) - 1) * 100),
    `% change 95% CI` = sprintf(
      "%.1f%% to %.1f%%",
      (exp(ci_low) - 1) * 100,
      (exp(ci_high) - 1) * 100
    )
  ) %>%
  select(
    Category,
    Calcifier,
    Treatment,
    k,
    `Mean lnRR`,
    `95% CI`,
    `p-value`,
    `% change`,
    `% change 95% CI`
  ) %>%
  arrange(
    factor(Category, levels = category_levels),
    factor(Calcifier, levels = calcifier_levels),
    factor(Treatment, levels = treatment_label_levels)
  )

write.csv(
  pooled_table_plot,
  file.path(table_dir, "pooled_estimates_calcifier_summary_table.csv"),
  row.names = FALSE
)

saveRDS(
  pooled_table_plot,
  file.path(rds_dir, "pooled_estimates_calcifier_summary_table.rds")
)

pooled_table_plot_char <- pooled_table_plot
pooled_table_plot_char[] <- lapply(pooled_table_plot_char, as.character)

tbl_grob <- gridExtra::tableGrob(
  pooled_table_plot_char,
  rows = NULL,
  theme = gridExtra::ttheme_minimal(
    base_size = 12,
    base_family = fam,
    padding = unit(c(4, 4), "mm"),
    colhead = list(
      fg_params = list(fontface = "bold", col = "grey10"),
      bg_params = list(fill = "grey90", col = NA)
    ),
    core = list(
      fg_params = list(col = "grey10"),
      bg_params = list(fill = "white", col = NA)
    )
  )
)

png(
  filename = file.path(figure_dir, "pooled_summary_table_calcifier_lnRR_percent_change_pvalues.png"),
  width = 3800,
  height = 2400,
  res = 300,
  bg = "white"
)

grid.newpage()

grid.text(
  "Pooled mean effect sizes by calcifier status, 95% confidence intervals, p-values, and percent change",
  x = 0.5,
  y = 0.98,
  gp = gpar(
    fontsize = 15,
    fontface = "bold",
    col = "grey10",
    fontfamily = fam
  )
)

grid.draw(
  editGrob(
    tbl_grob,
    vp = viewport(
      x = 0.5,
      y = 0.47,
      width = 0.98,
      height = 0.88
    )
  )
)

dev.off()

# ============================================================
# Plot aesthetics following 02B style
# ============================================================

ymin <- -0.9
ymax <- 0.6

brks <- seq(ymin, ymax, by = 0.3)
brks[abs(brks) < 1e-12] <- 0
brks_no0 <- brks[brks != 0]

category_positions <- tibble(
  Category = factor(category_levels, levels = category_levels),
  category_x = seq_along(category_levels)
)

vbreaks <- seq(0.5, length(category_levels) + 0.5, by = 1)

col_trt_line <- c(
  "pH" = "#8A2D00",
  "T" = "#2C7FB8",
  "T + pH" = "#2CA25F"
)

fill_trt <- c(
  "pH" = "#F4C2A1",
  "T" = "#CFE7F3",
  "T + pH" = "#CFEBD8"
)

# Filled symbols:
# Calcifier = filled circle
# Non-calcifier = filled triangle
calcifier_shapes <- c(
  "Calcifier" = 21,
  "Non-calcifier" = 24
)

# Six streams per category:
# Calcifier pH, Calcifier T, Calcifier T + pH,
# Non-calcifier pH, Non-calcifier T, Non-calcifier T + pH
stream_positions <- tidyr::expand_grid(
  Calcifier = factor(calcifier_levels, levels = calcifier_levels),
  Treatment_label = factor(treatment_label_levels, levels = treatment_label_levels)
) %>%
  mutate(
    stream_id = paste(Calcifier, Treatment_label, sep = "_"),
    offset = case_when(
      Calcifier == "Calcifier"     & Treatment_label == "pH"     ~ -0.30,
      Calcifier == "Calcifier"     & Treatment_label == "T"      ~ -0.18,
      Calcifier == "Calcifier"     & Treatment_label == "T + pH" ~ -0.06,
      Calcifier == "Non-calcifier" & Treatment_label == "pH"     ~  0.06,
      Calcifier == "Non-calcifier" & Treatment_label == "T"      ~  0.18,
      Calcifier == "Non-calcifier" & Treatment_label == "T + pH" ~  0.30,
      TRUE ~ 0
    )
  )

# Overall background
overall_idx  <- which(category_levels == "Overall")
overall_xmin <- overall_idx - 0.5
overall_xmax <- overall_idx + 0.5

summary_fill <- "grey95"

# Separator position between Survival and Overall
response_overall_separator <- which(category_levels == "Survival") + 0.5

# ============================================================
# Raw points for plotting
# ============================================================

raw_all_final <- bind_rows(
  # Raw points for response-category columns
  dat1_calc %>%
    filter(Response_category %in% response_levels, is.finite(yi)) %>%
    transmute(
      Category = factor(as.character(Response_category), levels = category_levels),
      yi_raw = yi,
      Treatment_label = factor(Treatment_label, levels = treatment_label_levels),
      Calcifier = factor(Calcifier, levels = calcifier_levels),
      ES_ID = ES_ID
    ),
  
  # Raw points for the Overall column
  dat1_calc %>%
    filter(is.finite(yi)) %>%
    transmute(
      Category = factor("Overall", levels = category_levels),
      yi_raw = yi,
      Treatment_label = factor(Treatment_label, levels = treatment_label_levels),
      Calcifier = factor(Calcifier, levels = calcifier_levels),
      ES_ID = ES_ID
    )
) %>%
  left_join(category_positions, by = "Category") %>%
  left_join(stream_positions, by = c("Calcifier", "Treatment_label"))

set.seed(123)

raw_all_final <- raw_all_final %>%
  mutate(
    x_pos = category_x + offset,
    x_jitter = x_pos + runif(n(), min = -0.025, max = 0.025)
  )

# ============================================================
# Pooled model summary positions
# ============================================================

sum_df_plot <- sum_df_calcifier %>%
  mutate(
    Category = factor(Category, levels = category_levels),
    Treatment_label = factor(Treatment_label, levels = treatment_label_levels),
    Calcifier = factor(Calcifier, levels = calcifier_levels)
  ) %>%
  left_join(category_positions, by = "Category") %>%
  left_join(stream_positions, by = c("Calcifier", "Treatment_label")) %>%
  mutate(
    x_pos = category_x + offset
  )

# ============================================================
# Final one-panel calcifier plot
# ============================================================

p_calcifier <- ggplot() +
  annotate(
    "rect",
    xmin = overall_xmin,
    xmax = overall_xmax,
    ymin = -Inf,
    ymax = Inf,
    fill = summary_fill,
    colour = NA
  ) +
  geom_hline(yintercept = brks_no0, colour = "grey90", linewidth = 0.3) +
  geom_vline(xintercept = vbreaks, colour = "grey90", linewidth = 0.3) +
  geom_hline(yintercept = 0, colour = "grey50", linewidth = 0.4) +
  
  # Strong separator between Survival and Overall
  geom_vline(
    xintercept = response_overall_separator,
    linetype = "solid",
    colour = "black",
    linewidth = 0.4
  ) +
  
  # Raw points
  # These remain pale, because the pooled estimates should stand out.
  geom_point(
    data = raw_all_final,
    aes(
      x = x_jitter,
      y = yi_raw,
      fill = Treatment_label,
      shape = Calcifier
    ),
    colour = "grey25",
    stroke = 0.10,
    alpha = 0.18,
    size = 2.0,
    show.legend = FALSE
  ) +
  
  # Pooled CI
  # Solid CI bars for both Calcifier and Non-calcifier.
  geom_errorbar(
    data = sum_df_plot,
    aes(
      x = x_pos,
      ymin = ci_low,
      ymax = ci_high,
      colour = Treatment_label
    ),
    width = 0.035,
    linewidth = 0.75
  ) +
  
  # Pooled mean
  # Filled centre symbols with darker outline.
  geom_point(
    data = sum_df_plot,
    aes(
      x = x_pos,
      y = mean,
      colour = Treatment_label,
      fill = Treatment_label,
      shape = Calcifier
    ),
    size = 3.1,
    stroke = 0.85
  ) +
  
  # k labels above plot
  geom_text(
    data = sum_df_plot,
    aes(
      x = x_pos,
      y = Inf,
      label = k,
      colour = Treatment_label
    ),
    vjust = -1.2,
    family = fam,
    size = 7.5,
    show.legend = FALSE
  ) +
  
  scale_colour_manual(
    values = col_trt_line,
    labels = c("pH", "T", "T + pH"),
    name = NULL
  ) +
  scale_fill_manual(
    values = fill_trt,
    labels = c("pH", "T", "T + pH"),
    name = NULL
  ) +
  scale_shape_manual(
    values = calcifier_shapes,
    name = NULL
  ) +
  scale_x_continuous(
    breaks = seq_along(category_levels),
    labels = category_levels,
    limits = c(0.5, length(category_levels) + 0.5),
    expand = expansion(add = 0)
  ) +
  scale_y_continuous(
    limits = c(ymin, ymax),
    breaks = brks,
    labels = label_number(accuracy = 0.1, trim = TRUE),
    expand = expansion(mult = c(0, 0)),
    oob = squish,
    name = "Mean effect size (lnRR)",
    sec.axis = sec_axis(
      ~ (exp(.) - 1) * 100,
      breaks = c(-30, 0, 30),
      labels = c("-30%", "0%", "30%"),
      name = "% change"
    )
  ) +
  xlab(NULL) +
  theme(
    text = element_text(family = fam, colour = "grey10"),
    
    axis.text.y.left = element_text(
      size = 28,
      colour = "grey10",
      margin = margin(r = 8)
    ),
    
    axis.text.y.right = element_text(
      size = 28,
      colour = "grey10",
      margin = margin(l = 8)
    ),
    
    axis.title.y.left = element_text(
      size = 33,
      margin = margin(r = 12)
    ),
    
    axis.title.y.right = element_text(
      size = 33,
      hjust = 0.40,
      margin = margin(l = 12)
    ),
    
    axis.text.x = element_text(
      angle = 30,
      hjust = 1,
      size = 33,
      colour = "grey10",
      vjust = 1.0
    ),
    
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    
    axis.ticks = element_line(
      colour = "grey10",
      linewidth = 0.4
    ),
    
    axis.ticks.length = unit(4, "pt"),
    
    legend.position = "top",
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(0, 0, 0, 0),
    legend.text = element_text(size = 30),
    
    plot.margin = margin(15, 15, 15, 15)
  ) +
  coord_cartesian(clip = "off") +
  geom_rect(
    aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf),
    inherit.aes = FALSE,
    fill = NA,
    colour = "black",
    linewidth = 0.4
  ) +
  guides(
    colour = guide_legend(
      nrow = 1,
      byrow = TRUE,
      override.aes = list(
        shape = 16,
        linewidth = 0.8
      )
    ),
    fill = "none",
    shape = guide_legend(
      nrow = 1,
      byrow = TRUE,
      override.aes = list(
        colour = "grey10",
        fill = "grey80",
        size = 4
      )
    )
  )

# ============================================================
# Save final plot
# ============================================================

ggsave(
  filename = file.path(
    figure_dir,
    "zooplankton_calcifier_noncalcifier_one_panel_clean_StudyID_ESID.png"
  ),
  plot = p_calcifier,
  width = 14.5,
  height = 5.2,
  units = "in",
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = file.path(
    figure_dir,
    "zooplankton_calcifier_noncalcifier_one_panel_clean_StudyID_ESID.pdf"
  ),
  plot = p_calcifier,
  width = 13.5,
  height = 5.0,
  units = "in",
  bg = "white"
)

saveRDS(
  p_calcifier,
  file.path(rds_dir, "p_calcifier_one_panel.rds")
)

# ============================================================
# Console prints
# ============================================================

cat("\n==== Plot saved ====\n")
cat("Figures:", figure_dir, "\n")
cat("Tables :", table_dir, "\n")
cat("RDS    :", rds_dir, "\n")