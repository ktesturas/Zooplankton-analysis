# ============================================================
# 05B_Plot_regression_T.R
# Purpose:
#   Create delta_T meta-regression bubble plots using saved
#   outputs from 05A_Regression_analysis_T.R.
#
# Important:
#   This script does NOT run rma.mv().
#   Use this when only editing plots, labels, font sizes, colours,
#   or figure dimensions.
#
# Input:
#   Outputs/RDS/05A_Regression_analysis_T/
#     - overall_deltaT_meta_regression_result.rds
#     - subgroup_deltaT_meta_regression_results.rds
#
# Outputs:
#   Outputs/Figures/05B_Plot_regression_T/
#   Outputs/RDS/05B_Plot_regression_T/
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

input_script_name <- "05A_Regression_analysis_T"
script_name       <- "05B_Plot_regression_T"

input_rds_dir <- file.path("Outputs", "RDS", input_script_name)

rds_dir    <- file.path("Outputs", "RDS", script_name)
figure_dir <- file.path("Outputs", "Figures", script_name)

dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

if (!dir.exists(input_rds_dir)) {
  stop(paste("Input RDS directory not found. Run 05A first:", input_rds_dir))
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

overall_rds_path <- file.path(
  input_rds_dir,
  "overall_deltaT_meta_regression_result.rds"
)

subgroup_rds_path <- file.path(
  input_rds_dir,
  "subgroup_deltaT_meta_regression_results.rds"
)

if (!file.exists(overall_rds_path)) {
  stop(paste("Missing file. Run 05A first:", overall_rds_path))
}

if (!file.exists(subgroup_rds_path)) {
  stop(paste("Missing file. Run 05A first:", subgroup_rds_path))
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

# Notes:
#   axis_title_size  = x-axis and y-axis title size
#   axis_text_size   = x-axis and y-axis number size
#   intext_plot_size = p-value and slope text inside each plot
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
    legend.position = "none",
    
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
    
    # Panel label above each plot
    plot.title = element_text(
      size = panel_title_size,
      face = "plain",
      hjust = 0,
      margin = margin(0, 0, 8, 0)
    )
  )
}

# ============================================================
# Load saved model outputs from 05A
# ============================================================

overall_res <- readRDS(overall_rds_path)
subgroup_results <- readRDS(subgroup_rds_path)

overall_deltaT_summary <- overall_res$summary

cat("\n==== Loaded 05A outputs ====\n")
cat("Overall result loaded:", file.exists(overall_rds_path), "\n")
cat("Subgroup results loaded:", file.exists(subgroup_rds_path), "\n")
cat("Number of subgroup results:", length(subgroup_results), "\n")

# ============================================================
# Plot function
# ============================================================

plot_metareg_deltaT <- function(plot_df,
                                fit_df,
                                ci_poly,
                                summary_row,
                                title_text = "",
                                point_fill = "#8A2D00",
                                line_type = 1,
                                y_limits = NULL,
                                y_breaks = waiver(),
                                x_lab = expression(Delta * "T (" * degree * "C)"),
                                y_lab = "Mean effect size (lnRR)") {
  
  if (nrow(plot_df) == 0 || nrow(fit_df) == 0 || nrow(ci_poly) == 0) {
    return(NULL)
  }
  
  slope_lab <- if (is.finite(summary_row$beta_delta_T)) {
    paste0(
      "slope = ",
      formatC(summary_row$beta_delta_T, format = "f", digits = 4)
    )
  } else {
    "slope = NA"
  }
  
  p_lab <- if (is.finite(summary_row$p_delta_T)) {
    if (summary_row$p_delta_T < 0.0001) {
      "p < 0.0001"
    } else {
      paste0(
        "p = ",
        formatC(summary_row$p_delta_T, format = "f", digits = 4)
      )
    }
  } else {
    "p = NA"
  }
  
  x_min <- min(fit_df$X, na.rm = TRUE)
  x_max <- max(fit_df$X, na.rm = TRUE)
  x_rng <- x_max - x_min
  
  if (is.null(y_limits)) {
    y_min <- min(c(plot_df$yi, ci_poly$y), na.rm = TRUE)
    y_max <- max(c(plot_df$yi, ci_poly$y), na.rm = TRUE)
    y_pad <- 0.08 * (y_max - y_min)
    y_limits <- c(y_min - y_pad, y_max + y_pad)
  }
  
  y_rng <- diff(y_limits)
  
  ggplot(plot_df, aes(x = delta_T, y = yi)) +
    geom_point(
      fill = point_fill,
      colour = "black",
      size = plot_df$BubbleSize,
      shape = 21
    ) +
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      colour = "#000000",
      linewidth = 0.4
    ) +
    geom_polygon(
      data = ci_poly,
      aes(x = x, y = y),
      inherit.aes = FALSE,
      fill = "black",
      alpha = 0.20
    ) +
    geom_line(
      data = fit_df,
      aes(x = X, y = Y),
      inherit.aes = FALSE,
      linewidth = 0.8,
      colour = "black",
      linetype = line_type
    ) +
    annotate(
      "text",
      x = x_max - 0.16 * x_rng,
      y = y_limits[2] - 0.10 * y_rng,
      label = p_lab,
      size = intext_plot_size,
      hjust = 0
    ) +
    annotate(
      "text",
      x = x_max - 0.22 * x_rng,
      y = y_limits[2] - 0.19 * y_rng,
      label = slope_lab,
      size = intext_plot_size,
      hjust = 0
    ) +
    scale_y_continuous(
      limits = y_limits,
      breaks = y_breaks,
      expand = c(0, 0)
    ) +
    labs(
      title = title_text,
      x = x_lab,
      y = y_lab
    ) +
    theme_meta()
}

# ============================================================
# Make overall plot
# ============================================================

p_overall <- plot_metareg_deltaT(
  plot_df = overall_res$plot_df,
  fit_df = overall_res$fit_df,
  ci_poly = overall_res$ci_poly,
  summary_row = overall_deltaT_summary,
  title_text = "(a) Overall",
  point_fill = "#8A2D00",
  line_type = 1,
  x_lab = expression(Delta * "T (" * degree * "C)"),
  y_lab = "Mean effect size (lnRR)"
)

if (!is.null(p_overall)) {
  ggsave(
    filename = file.path(figure_dir, "overall_deltaT_meta_regression_plot.png"),
    plot = p_overall,
    width = 8.5,
    height = 6.8,
    dpi = 600,
    bg = "white"
  )
  
  saveRDS(
    p_overall,
    file.path(rds_dir, "overall_deltaT_meta_regression_plot.rds")
  )
}

# ============================================================
# Make subgroup plots
# ============================================================

subgroup_plots <- list()

for (i in seq_along(subgroup_results)) {
  
  g <- subgroup_results[[i]]$subgroup
  res_g <- subgroup_results[[i]]$result
  
  point_col <- if (g == "Meroplankton") "#1b7c3d" else "#2b6a99"
  line_ty   <- if (g == "Meroplankton") 1 else 2
  panel_lab <- paste0("(", letters[i + 1], ") ", g)
  
  p_g <- plot_metareg_deltaT(
    plot_df = res_g$plot_df,
    fit_df = res_g$fit_df,
    ci_poly = res_g$ci_poly,
    summary_row = res_g$summary,
    title_text = panel_lab,
    point_fill = point_col,
    line_type = line_ty,
    x_lab = expression(Delta * "T (" * degree * "C)"),
    y_lab = "Mean effect size (lnRR)"
  )
  
  subgroup_plots[[g]] <- p_g
  
  if (!is.null(p_g)) {
    
    safe_name <- str_replace_all(g, "[^A-Za-z0-9]+", "_")
    
    ggsave(
      filename = file.path(
        figure_dir,
        paste0("subgroup_", safe_name, "_deltaT_meta_regression_plot.png")
      ),
      plot = p_g,
      width = 8.5,
      height = 6.8,
      dpi = 600,
      bg = "white"
    )
    
    saveRDS(
      p_g,
      file.path(
        rds_dir,
        paste0("subgroup_", safe_name, "_deltaT_meta_regression_plot.rds")
      )
    )
  }
}

# ============================================================
# Combine all plots in one row
# ============================================================

plots_one_row <- c(list(p_overall), subgroup_plots)
plots_one_row <- plots_one_row[!vapply(plots_one_row, is.null, logical(1))]

if (length(plots_one_row) > 0) {
  
  p_combined_row <- wrap_plots(plots_one_row, nrow = 1)
  
  ggsave(
    filename = file.path(figure_dir, "deltaT_meta_regression_all_panels_one_row.png"),
    plot = p_combined_row,
    width = 7 * length(plots_one_row),
    height = 6.8,
    dpi = 600,
    bg = "white"
  )
  
  saveRDS(
    p_combined_row,
    file.path(rds_dir, "deltaT_meta_regression_all_panels_one_row.rds")
  )
}

# ============================================================
# Console output
# ============================================================

cat("\n==== Saved figures ====\n")
cat("Figures:", figure_dir, "\n")
cat("RDS    :", rds_dir, "\n")

if (!is.null(p_overall)) {
  cat(
    "Overall figure:",
    file.path(figure_dir, "overall_deltaT_meta_regression_plot.png"),
    "\n"
  )
}

if (length(subgroup_plots) > 0) {
  cat("Subgroup figures saved for:\n")
  print(names(subgroup_plots))
}

if (exists("p_combined_row")) {
  cat(
    "Combined figure:",
    file.path(figure_dir, "deltaT_meta_regression_all_panels_one_row.png"),
    "\n"
  )
}

# ============================================================
# Show plots in console
# ============================================================

print(p_overall)

if (length(subgroup_plots) > 0) {
  invisible(lapply(subgroup_plots, print))
}

if (exists("p_combined_row")) {
  print(p_combined_row)
}