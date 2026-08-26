# ============================================================
# 07B_Plot_regression_interaction.R
# Purpose:
#   Create interaction plots using saved outputs from
#   07A_Regression_analysis_interaction.R.
#
# Important:
#   This script does NOT run rma.mv().
#   Use this when only editing plots, labels, font sizes, colours,
#   legends, or figure dimensions.
#
# Input:
#   Outputs/RDS/07A_Regression_analysis_interaction/
#     - overall_run.rds
#     - broad_subgroup_runs.rds
#     - response_category_runs.rds
#     - all_interaction_model_summaries_combined.rds
#
# Outputs:
#   Outputs/Figures/07B_Plot_regression_interaction/
#   Outputs/RDS/07B_Plot_regression_interaction/
# ============================================================

suppressPackageStartupMessages(invisible(
  lapply(
    c(
      "dplyr", "ggplot2", "patchwork", "stringr",
      "showtext", "tibble", "grid", "scales"
    ),
    library,
    character.only = TRUE
  )
))

# ============================================================
# Input and output folders
# ============================================================

input_script_name <- "07A_Regression_analysis_interaction"
script_name       <- "07B_Plot_regression_interaction"

input_rds_dir <- file.path("Outputs", "RDS", input_script_name)

rds_dir    <- file.path("Outputs", "RDS", script_name)
figure_dir <- file.path("Outputs", "Figures", script_name)

dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

if (!dir.exists(input_rds_dir)) {
  stop(paste("Input RDS directory not found. Run 07A first:", input_rds_dir))
}

if (!dir.exists(rds_dir)) {
  stop(paste("RDS directory could not be created:", rds_dir))
}

if (!dir.exists(figure_dir)) {
  stop(paste("Figure directory could not be created:", figure_dir))
}

cat("\nWorking directory:\n", getwd(), "\n")
cat("\nInput RDS directory:\n", input_rds_dir, "\n")
cat("Output RDS directory:\n", rds_dir, "\n")
cat("Figure directory:\n", figure_dir, "\n")

# ============================================================
# Check required input files
# ============================================================

overall_run_rds_path <- file.path(input_rds_dir, "overall_run.rds")
subgroup_runs_rds_path <- file.path(input_rds_dir, "broad_subgroup_runs.rds")
response_runs_rds_path <- file.path(input_rds_dir, "response_category_runs.rds")
master_summary_rds_path <- file.path(input_rds_dir, "all_interaction_model_summaries_combined.rds")

if (!file.exists(overall_run_rds_path)) {
  stop(paste("Missing file. Run 07A first:", overall_run_rds_path))
}

if (!file.exists(subgroup_runs_rds_path)) {
  stop(paste("Missing file. Run 07A first:", subgroup_runs_rds_path))
}

if (!file.exists(response_runs_rds_path)) {
  stop(paste("Missing file. Run 07A first:", response_runs_rds_path))
}

if (!file.exists(master_summary_rds_path)) {
  stop(paste("Missing file. Run 07A first:", master_summary_rds_path))
}

# ============================================================
# Fonts and adjustable text sizes
# ============================================================

font_add_google("Inter", "helv")
showtext_auto()
fam <- "helv"

theme_set(theme_minimal(base_size = 18, base_family = fam))
theme_update(text = element_text(family = fam))

# -----------------------------
# EDIT FONT SIZES HERE
# -----------------------------

axis_title_size       <- 93
axis_text_size        <- 61
intext_plot_size      <- 30
panel_title_size      <- axis_title_size

legend_title_size     <- 61
legend_text_size      <- 50
subtitle_size         <- 61

# Notes:
#   axis_title_size  = x-axis and y-axis title size
#   axis_text_size   = x-axis and y-axis number size
#   intext_plot_size = not heavily used here, but retained for consistency
#   panel_title_size = panel labels above each plot
#
# You said you want labels above the plot to be the same size
# as the x and y axis titles, so:
#
#   panel_title_size <- axis_title_size

# ============================================================
# Theme
# ============================================================

theme_meta <- function() {
  theme(
    panel.background = element_rect(
      fill = "white",
      colour = "black",
      linewidth = 0.7
    ),
    panel.grid = element_blank(),
    legend.position = "right",
    
    plot.margin = margin(12, 8, 8, 8),
    plot.background = element_blank(),
    
    axis.ticks = element_line(linewidth = 0.3, colour = "black"),
    axis.ticks.length = unit(-0.15, "lines"),
    
    axis.title.y = element_text(
      size = axis_title_size,
      margin = margin(0, 6, 0, 0),
      face = "bold"
    ),
    axis.title.x = element_text(
      size = axis_title_size,
      margin = margin(6, 0, 0, 0),
      face = "bold"
    ),
    axis.text.y = element_text(
      size = axis_text_size,
      margin = margin(0, 5, 0, 0),
      colour = "#000000"
    ),
    axis.text.x = element_text(
      size = axis_text_size,
      margin = margin(5, 0, 0, 0),
      colour = "#000000"
    ),
    
    legend.title = element_text(
      size = legend_title_size,
      face = "bold"
    ),
    legend.text = element_text(
      size = legend_text_size
    ),
    
    # Panel title above plot
    plot.title = element_text(
      size = panel_title_size,
      face = "bold",
      hjust = 0,
      margin = margin(0, 0, 8, 0)
    ),
    
    plot.subtitle = element_text(
      size = subtitle_size,
      hjust = 0,
      margin = margin(0, 0, 8, 0)
    )
  )
}

# ============================================================
# Load saved model outputs from 07A
# ============================================================

overall_run <- readRDS(overall_run_rds_path)
subgroup_runs <- readRDS(subgroup_runs_rds_path)
response_runs <- readRDS(response_runs_rds_path)
master_summary_tbl <- readRDS(master_summary_rds_path)

cat("\n==== Loaded 07A outputs ====\n")
cat("Overall run loaded:", file.exists(overall_run_rds_path), "\n")
cat("Subgroup runs loaded:", file.exists(subgroup_runs_rds_path), "\n")
cat("Response runs loaded:", file.exists(response_runs_rds_path), "\n")
cat("Master summary loaded:", file.exists(master_summary_rds_path), "\n")

# ============================================================
# Plot interaction function
# ============================================================

plot_metareg_interaction <- function(plot_df,
                                     pred_grid,
                                     summary_row,
                                     title_text = "(a) TpH interaction",
                                     y_lab = "Mean effect size (lnRR)") {
  
  if (nrow(plot_df) == 0 || nrow(pred_grid) == 0) {
    return(NULL)
  }
  
  p_int_lab <- if (is.finite(summary_row$p_interaction)) {
    if (summary_row$p_interaction < 0.0001) {
      "interaction p < 0.0001"
    } else {
      paste0(
        "interaction p = ",
        formatC(summary_row$p_interaction, format = "f", digits = 4)
      )
    }
  } else {
    "interaction p = NA"
  }
  
  b_int_lab <- if (is.finite(summary_row$beta_interaction)) {
    paste0(
      "interaction = ",
      formatC(summary_row$beta_interaction, format = "f", digits = 6)
    )
  } else {
    "interaction = NA"
  }
  
  point_breaks <- quantile(
    plot_df$delta_CO2,
    probs = c(0, 1 / 3, 2 / 3, 1),
    na.rm = TRUE
  )
  
  point_breaks <- unique(point_breaks)
  
  if (length(point_breaks) < 4) {
    plot_df_points <- plot_df %>%
      mutate(CO2_group = "Observed TpH studies")
  } else {
    plot_df_points <- plot_df %>%
      mutate(
        CO2_group = cut(
          delta_CO2,
          breaks = point_breaks,
          include.lowest = TRUE,
          labels = c(
            "Lower observed \u0394CO2",
            "Mid observed \u0394CO2",
            "Higher observed \u0394CO2"
          )
        )
      )
  }
  
  ggplot() +
    geom_point(
      data = plot_df_points,
      aes(x = delta_T, y = yi, colour = CO2_group),
      size = 1.6,
      alpha = 0.28,
      stroke = 0
    ) +
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      colour = "black",
      linewidth = 0.4
    ) +
    geom_line(
      data = pred_grid,
      aes(x = delta_T, y = pred, linetype = CO2_label),
      linewidth = 1.0,
      colour = "black"
    ) +
    scale_colour_manual(
      values = c(
        "Lower observed \u0394CO2" = "#9ecae1",
        "Mid observed \u0394CO2" = "#4292c6",
        "Higher observed \u0394CO2" = "#08519c",
        "Observed TpH studies" = "#4292c6"
      ),
      name = "Observed data"
    ) +
    scale_linetype_manual(
      values = c("solid", "dashed", "dotdash"),
      name = "Predicted lines"
    ) +
    labs(
      title = title_text,
      subtitle = paste0(p_int_lab, "   |   ", b_int_lab),
      x = expression(Delta * "T (" * degree * "C)"),
      y = y_lab
    ) +
    theme_meta()
}

# ============================================================
# Plot runner
# ============================================================

plot_one_analysis <- function(run_obj) {
  
  p <- plot_metareg_interaction(
    plot_df = run_obj$result$plot_df,
    pred_grid = run_obj$result$pred_grid,
    summary_row = run_obj$result$summary,
    title_text = run_obj$title,
    y_lab = "Mean effect size (lnRR)"
  )
  
  if (!is.null(p)) {
    ggsave(
      filename = file.path(figure_dir, paste0(run_obj$id, ".png")),
      plot = p,
      width = 8.8,
      height = 6.8,
      dpi = 600,
      bg = "white"
    )
    
    saveRDS(
      p,
      file.path(rds_dir, paste0(run_obj$id, "_plot.rds"))
    )
  }
  
  p
}

# ============================================================
# Create individual plots
# ============================================================

overall_plot <- plot_one_analysis(overall_run)

subgroup_plots <- lapply(subgroup_runs, plot_one_analysis)
names(subgroup_plots) <- vapply(subgroup_runs, function(x) x$title, character(1))

response_plots <- lapply(response_runs, plot_one_analysis)
names(response_plots) <- vapply(response_runs, function(x) x$title, character(1))

# ============================================================
# Combined panel: overall and broad subgroups
# ============================================================

subgroup_plot_list <- c(
  list(overall_plot),
  subgroup_plots
)

subgroup_plot_list <- subgroup_plot_list[
  !vapply(subgroup_plot_list, is.null, logical(1))
]

if (length(subgroup_plot_list) > 0) {
  
  p_subgroups_combined <- wrap_plots(subgroup_plot_list, nrow = 1)
  
  ggsave(
    filename = file.path(figure_dir, "combined_overall_and_broad_subgroups.png"),
    plot = p_subgroups_combined,
    width = 8 * length(subgroup_plot_list),
    height = 6.8,
    dpi = 600,
    bg = "white"
  )
  
  saveRDS(
    p_subgroups_combined,
    file.path(rds_dir, "combined_overall_and_broad_subgroups.rds")
  )
}

# ============================================================
# Reorder response-category plots
# ============================================================

response_plot_order <- c(
  "Growth",
  "Development",
  "Physiology",
  "Reproduction",
  "Survival"
)

response_runs <- response_runs[
  order(match(
    str_remove(
      vapply(response_runs, function(x) x$title, character(1)),
      "^Response: "
    ),
    response_plot_order
  ))
]

response_runs <- response_runs[
  vapply(response_runs, function(x) {
    str_remove(x$title, "^Response: ") %in% response_plot_order
  }, logical(1))
]

response_plot_list <- lapply(response_runs, function(x) {
  plot_metareg_interaction(
    plot_df = x$result$plot_df,
    pred_grid = x$result$pred_grid,
    summary_row = x$result$summary,
    title_text = str_remove(x$title, "^Response: "),
    y_lab = "Mean effect size (lnRR)"
  )
})

names(response_plot_list) <- vapply(response_runs, function(x) {
  str_remove(x$title, "^Response: ")
}, character(1))

response_plot_list <- response_plot_list[
  !vapply(response_plot_list, is.null, logical(1))
]

# ============================================================
# Combined panel: response categories
# ============================================================

if (length(response_plot_list) > 0) {
  
  ncol_resp <- if (length(response_plot_list) <= 2) {
    length(response_plot_list)
  } else if (length(response_plot_list) <= 4) {
    2
  } else {
    3
  }
  
  p_response_combined <- wrap_plots(response_plot_list, ncol = ncol_resp)
  
  ggsave(
    filename = file.path(figure_dir, "combined_response_categories.png"),
    plot = p_response_combined,
    width = 7.5 * ncol_resp,
    height = 6.5 * ceiling(length(response_plot_list) / ncol_resp),
    dpi = 600,
    bg = "white"
  )
  
  saveRDS(
    p_response_combined,
    file.path(rds_dir, "combined_response_categories.rds")
  )
}

# ============================================================
# Console output
# ============================================================

cat("\n==== Master interaction summary ====\n")
print(master_summary_tbl)

cat("\n==== Saved figures ====\n")
cat("Figures:", figure_dir, "\n")
cat("RDS    :", rds_dir, "\n")

if (exists("p_subgroups_combined")) {
  cat(
    "Combined subgroup figure:",
    file.path(figure_dir, "combined_overall_and_broad_subgroups.png"),
    "\n"
  )
}

if (exists("p_response_combined")) {
  cat(
    "Combined response-category figure:",
    file.path(figure_dir, "combined_response_categories.png"),
    "\n"
  )
}

# ============================================================
# Show plots in console
# ============================================================

if (!is.null(overall_plot)) {
  print(overall_plot)
}

if (length(subgroup_plot_list) > 0 && exists("p_subgroups_combined")) {
  print(p_subgroups_combined)
}

if (length(response_plot_list) > 0 && exists("p_response_combined")) {
  print(p_response_combined)
}

cat("\nDone.\n")