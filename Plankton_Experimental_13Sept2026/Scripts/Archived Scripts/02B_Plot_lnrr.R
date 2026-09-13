# ============================================================
# 02B_Plot_lnrr.R
# Purpose:
#   Create the final one-panel lnRR plot using saved outputs from
#   02A_Model_lnrr.R.
#
# Important:
#   This script does NOT run rma.mv().
#   Use this when only editing the plot.
#
# Input:
#   Outputs/RDS/02A_Model_lnrr/
#     - dat1_with_lnrr_vi.rds
#     - pooled_estimates_final_one_panel.rds
#     - plot_levels.rds
#
# Outputs:
#   Outputs/Figures/02B_Plot_lnrr/
#   Outputs/RDS/02B_Plot_lnrr/
# ============================================================

suppressPackageStartupMessages(invisible(
  lapply(
    c(
      "dplyr", "ggplot2", "readr", "stringr",
      "showtext", "tibble", "grid", "tidyr",
      "purrr", "scales"
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
# Script name and folders
# ============================================================

model_script_name <- "02A_Model_lnrr"
plot_script_name  <- "02B_Plot_lnrr"

model_rds_dir <- file.path("Outputs", "RDS", model_script_name)

figure_dir <- file.path("Outputs", "Figures", plot_script_name)
rds_dir    <- file.path("Outputs", "RDS", plot_script_name)

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# Load saved model outputs
# ============================================================

dat1_ml <- readRDS(
  file.path(model_rds_dir, "dat1_with_lnrr_vi.rds")
)

sum_df_final <- readRDS(
  file.path(model_rds_dir, "pooled_estimates_final_one_panel.rds")
)

plot_levels <- readRDS(
  file.path(model_rds_dir, "plot_levels.rds")
)

response_levels <- plot_levels$response_levels
treatment_levels <- plot_levels$treatment_levels

# Reorder plotting labels:
# Response categories first, then Holoplankton before Meroplankton, then Overall.
final_label_levels <- c(
  "Growth",
  "Development",
  "Physiology",
  "Reproduction",
  "Survival",
  "Holoplankton",
  "Meroplankton",
  "Overall"
)

subgroup_levels <- c("Holoplankton", "Meroplankton")

cat("\nLoaded saved model outputs from:", model_rds_dir, "\n")
cat("Rows in dat1_ml:", nrow(dat1_ml), "\n")
cat("Rows in sum_df_final:", nrow(sum_df_final), "\n\n")

# ============================================================
# Re-factor model summary labels for new plot order
# ============================================================

sum_df_final <- sum_df_final %>%
  mutate(
    labels = factor(as.character(labels), levels = final_label_levels),
    Treatment = factor(Treatment, levels = treatment_levels)
  )



# ============================================================
# Save pooled means, confidence intervals, percent change, and p-values as a table figure
# ============================================================

pooled_table_plot <- sum_df_final %>%
  mutate(
    Category = as.character(labels),
    Treatment = recode(
      as.character(Treatment),
      "pH"  = "pH",
      "T"   = "T",
      "TpH" = "T + pH"
    ),
    
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
    Treatment,
    k,
    `Mean lnRR`,
    `95% CI`,
    `p-value`,
    `% change`,
    `% change 95% CI`
  ) %>%
  arrange(
    factor(Category, levels = final_label_levels),
    factor(Treatment, levels = c("pH", "Temperature", "Temperature + pH"))
  )

pooled_table_plot[] <- lapply(pooled_table_plot, as.character)

tbl_grob <- gridExtra::tableGrob(
  pooled_table_plot,
  rows = NULL,
  theme = gridExtra::ttheme_minimal(
    base_size = 14,
    base_family = fam,
    padding = unit(c(5, 5), "mm"),
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
  filename = file.path(figure_dir, "pooled_summary_table_lnRR_percent_change_pvalues.png"),
  width = 3400,
  height = 2000,
  res = 300,
  bg = "white"
)

grid.newpage()

grid.text(
  "Pooled mean effect sizes, 95% confidence intervals, p-values, and percent change",
  x = 0.5,
  y = 0.98,
  gp = gpar(
    fontsize = 16,
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
# Raw points for plotting
# ============================================================

raw_all_final <- bind_rows(
  # Raw points for response-category columns
  dat1_ml %>%
    filter(Response_category %in% response_levels, is.finite(yi)) %>%
    transmute(
      labels = factor(Response_category, levels = final_label_levels),
      yi_raw = yi,
      Treatment = factor(Treatment, levels = treatment_levels)
    ),
  
  # Raw points for the two subgroup-overall columns
  dat1_ml %>%
    filter(!is.na(Broad_subgroup), is.finite(yi)) %>%
    transmute(
      labels = factor(as.character(Broad_subgroup), levels = final_label_levels),
      yi_raw = yi,
      Treatment = factor(Treatment, levels = treatment_levels)
    ),
  
  # Raw points for the overall treatment column
  dat1_ml %>%
    filter(is.finite(yi)) %>%
    transmute(
      labels = factor("Overall", levels = final_label_levels),
      yi_raw = yi,
      Treatment = factor(Treatment, levels = treatment_levels)
    )
)

# ============================================================
# Plot aesthetics
# ============================================================

ymin <- -0.9
ymax <- 0.6
brks <- seq(ymin, ymax, by = 0.3)
brks[abs(brks) < 1e-12] <- 0
brks_no0 <- brks[brks != 0]

vbreaks <- seq(0.5, length(final_label_levels) + 0.5, by = 1)

col_trt_line <- c(
  "pH" = "#8A2D00",
  "T" = "#2C7FB8",
  "TpH" = "#2CA25F"
)

fill_trt <- c(
  "pH" = "#F4C2A1",
  "T" = "#CFE7F3",
  "TpH" = "#CFEBD8"
)

treat_lab <- c(
  "pH" = "pH",
  "T" = "T",
  "TpH" = "T + pH"
)

pd <- position_dodge(width = 0.55)

pjd <- position_jitterdodge(
  jitter.width = 0.22,
  jitter.height = 0,
  dodge.width = 0.55
)

# Background panels
subgroup_start_idx <- which(final_label_levels == "Holoplankton")
subgroup_end_idx   <- which(final_label_levels == "Meroplankton")
subgroup_xmin      <- subgroup_start_idx - 0.5
subgroup_xmax      <- subgroup_end_idx + 0.5

overall_idx  <- which(final_label_levels == "Overall")
overall_xmin <- overall_idx - 0.5
overall_xmax <- overall_idx + 0.5

summary_fill <- "grey95"

# Separator positions
response_subgroup_separator <- 5.5
subgroup_overall_separator  <- 7.5

# ============================================================
# Final one-panel plot
# ============================================================

p_final <- ggplot() +
  annotate(
    "rect",
    xmin = subgroup_xmin,
    xmax = subgroup_xmax,
    ymin = -Inf,
    ymax = Inf,
    fill = summary_fill,
    colour = NA
  ) +
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
  
  # Strong separator between Survival and plankton-group summaries
  geom_vline(
    xintercept = response_subgroup_separator,
    linetype = "solid",
    colour = "black",
    linewidth = 0.4
  ) +
  
  # Separator between Holoplankton and Meroplankton
  geom_vline(
    xintercept = 6.5,
    linetype = "dotted",
    colour = "grey70",
    linewidth = 0.4
  ) +
  
  # Separator between plankton-group summaries and Overall
  geom_vline(
    xintercept = subgroup_overall_separator,
    linetype = "dotted",
    colour = "grey70",
    linewidth = 0.4
  ) +
  
  geom_point(
    data = raw_all_final,
    aes(x = labels, y = yi_raw, fill = Treatment),
    position = pjd,
    shape = 21,
    stroke = 0.10,
    alpha = 0.23,
    size = 2.0,
    show.legend = FALSE
  ) +
  geom_errorbar(
    data = sum_df_final,
    aes(x = labels, ymin = ci_low, ymax = ci_high, colour = Treatment),
    width = 0.2,
    linewidth = 0.60,
    position = pd
  ) +
  geom_point(
    data = sum_df_final,
    aes(x = labels, y = mean, colour = Treatment),
    position = pd,
    size = 2.2,
    shape = 16
  ) +
  geom_text(
    data = sum_df_final,
    aes(x = labels, y = Inf, label = k, colour = Treatment),
    position = pd,
    vjust = -1.2,
    family = fam,
    size = 7.5,
    show.legend = FALSE
  ) +
  scale_colour_manual(values = col_trt_line, labels = treat_lab, name = NULL) +
  scale_fill_manual(values = fill_trt, labels = treat_lab, guide = "legend") +
  scale_x_discrete(
    drop = FALSE,
    limits = final_label_levels,
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
    text               = element_text(family = fam, colour = "grey10"),
    axis.text.y.left  = element_text(
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
    axis.text.x        = element_text(angle = 30, hjust = 1, size = 33, colour = "grey10", vjust = 1.0),
    panel.grid.major   = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.ticks         = element_line(colour = "grey10", linewidth = 0.4),
    axis.ticks.length  = unit(4, "pt"),
    legend.position    = "top",
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(0, 0, 0, 0),
    legend.text        = element_text(size = 30),
    plot.margin        = margin(15, 15, 15, 15)
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
    colour = guide_legend(nrow = 1, byrow = TRUE)
  )

# ============================================================
# Save final plot
# ============================================================
# Use ragg for cleaner PNG output and save with proportions closer
# to the RStudio plot panel.

ggsave(
  filename = file.path(
    figure_dir,
    "zooplankton_experimental_treatment_one_panel_clean_StudyID_ESID.png"
  ),
  plot = p_final,
  width = 12.1,
  height = 4.8,
  units = "in",
  dpi = 300,
  bg = "white")

