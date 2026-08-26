# ============================================================
# 04B_Plot_regression_CO2_response_categories.R
# Purpose:
#   Create response-category delta_CO2 meta-regression bubble plots
#   using saved outputs from 04A_Regression_analysis_CO2_response_categories.R.
#
# Important:
#   This script does NOT run rma.mv().
#   Use this when only editing plots, labels, font sizes, colours,
#   or figure dimensions.
#
# Input:
#   Outputs/RDS/04A_Regression_analysis_CO2_response_categories/
#     - response_category_deltaCO2_meta_regression_results.rds
#     - response_category_deltaCO2_meta_regression_summary.rds
#
# Outputs:
#   Outputs/Figures/04B_Plot_regression_CO2_response_categories/
#   Outputs/RDS/04B_Plot_regression_CO2_response_categories/
# ============================================================

suppressPackageStartupMessages(invisible(
  lapply(
    c(
      "dplyr", "ggplot2", "patchwork", "readr", "stringr",
      "showtext", "tibble", "grid", "scales"
    ),
    library,
    character.only = TRUE
  )
))

# ============================================================
# Input and output folders
# ============================================================

input_script_name <- "04A_Regression_analysis_CO2_response_categories"
script_name       <- "04B_Plot_regression_CO2_response_categories"

input_rds_dir <- file.path("Outputs", "RDS", input_script_name)

rds_dir    <- file.path("Outputs", "RDS", script_name)
figure_dir <- file.path("Outputs", "Figures", script_name)

dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

if (!dir.exists(input_rds_dir)) {
  stop(paste("Input RDS directory not found. Run 04A first:", input_rds_dir))
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

response_results_rds_path <- file.path(
  input_rds_dir,
  "response_category_deltaCO2_meta_regression_results.rds"
)

response_summary_rds_path <- file.path(
  input_rds_dir,
  "response_category_deltaCO2_meta_regression_summary.rds"
)

if (!file.exists(response_results_rds_path)) {
  stop(paste("Missing file. Run 04A first:", response_results_rds_path))
}

if (!file.exists(response_summary_rds_path)) {
  stop(paste("Missing file. Run 04A first:", response_summary_rds_path))
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
#combined_title_size   <- 18

# Notes:
#   axis_title_size  = x-axis and y-axis title size
#   axis_text_size   = x-axis and y-axis number size
#   intext_plot_size = p-value and slope text inside each plot
#   panel_title_size = response-category label above each plot
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
    
    plot.margin = margin(10, 8, 8, 8),
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
    
    # Response-category label above the plot
    plot.title = element_text(
      size = panel_title_size,
      face = "bold",
      hjust = 0.5,
      margin = margin(0, 0, 8, 0)
    )
  )
}

# ============================================================
# Load saved model outputs from 04A
# ============================================================

response_results <- readRDS(response_results_rds_path)
response_summary_tbl <- readRDS(response_summary_rds_path)

cat("\n==== Loaded 04A outputs ====\n")
cat("Response results loaded:", file.exists(response_results_rds_path), "\n")
cat("Response summary loaded:", file.exists(response_summary_rds_path), "\n")
cat("Number of response-category results:", length(response_results), "\n")

# ============================================================
# Response-category palette
# ============================================================

response_palette <- c(
  "Behaviour"    = "#A23A00",
  "Development"  = "#A23A00",
  "Growth"       = "#A23A00",
  "Physiology"   = "#A23A00",
  "Reproduction" = "#A23A00",
  "Survival"     = "#A23A00"
)

get_response_colour <- function(x) {
  if (x %in% names(response_palette)) {
    response_palette[[x]]
  } else {
    "#A23A00"
  }
}

# ============================================================
# Plot function
# ============================================================

plot_metareg_deltaCO2 <- function(plot_df,
                                  fit_df,
                                  ci_poly,
                                  summary_row,
                                  title_text = "",
                                  point_fill = "#8A2D00",
                                  line_type = 1,
                                  y_limits = NULL,
                                  y_breaks = waiver(),
                                  x_lab = expression(Delta * "CO"[2] * " (ppm)"),
                                  y_lab = "Effect size (lnRR)") {
  
  if (nrow(plot_df) == 0 || nrow(fit_df) == 0 || nrow(ci_poly) == 0) {
    return(NULL)
  }
  
  slope_lab <- if (is.finite(summary_row$beta_delta_CO2)) {
    paste0(
      "slope = ",
      formatC(summary_row$beta_delta_CO2, format = "f", digits = 4)
    )
  } else {
    "slope = NA"
  }
  
  p_lab <- if (is.finite(summary_row$p_delta_CO2)) {
    if (summary_row$p_delta_CO2 < 0.0001) {
      "p < 0.0001"
    } else {
      paste0(
        "p = ",
        formatC(summary_row$p_delta_CO2, format = "f", digits = 4)
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
  
  ggplot(plot_df, aes(x = delta_CO2, y = yi)) +
    geom_point(
      fill = point_fill,
      colour = "black",
      size = plot_df$BubbleSize,
      shape = 21,
      alpha = 0.60
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
      fill = point_fill,
      alpha = 0.15
    ) +
    geom_line(
      data = fit_df,
      aes(x = X, y = Y),
      inherit.aes = FALSE,
      linewidth = 0.9,
      colour = point_fill,
      linetype = line_type
    ) +
    annotate(
      "text",
      x = x_max - 0.22 * x_rng,
      y = y_limits[2] - 0.10 * y_rng,
      label = p_lab,
      size = intext_plot_size,
      hjust = 0
    ) +
    annotate(
      "text",
      x = x_max - 0.28 * x_rng,
      y = y_limits[2] - 0.20 * y_rng,
      label = slope_lab,
      size = intext_plot_size,
      hjust = 0
    ) +
    scale_y_continuous(
      limits = y_limits,
      breaks = y_breaks,
      expand = c(0.02, 0.02)
    ) +
    labs(
      title = title_text,
      x = x_lab,
      y = y_lab
    ) +
    theme_meta()
}

# ============================================================
# Set response-category plotting order
# ============================================================

response_plot_order <- c(
  "Growth",
  "Development",
  "Physiology",
  "Reproduction",
  "Survival"
)

response_results <- response_results[
  order(match(
    vapply(response_results, function(x) x$response_category, character(1)),
    response_plot_order
  ))
]

response_results <- response_results[
  vapply(response_results, function(x) {
    x$response_category %in% response_plot_order
  }, logical(1))
]

# ============================================================
# Create and save per-category plots
# ============================================================

response_plots <- list()

for (i in seq_along(response_results)) {
  
  rc <- response_results[[i]]$response_category
  res_rc <- response_results[[i]]$result
  
  point_col <- get_response_colour(rc)
  
  p_rc <- plot_metareg_deltaCO2(
    plot_df = res_rc$plot_df,
    fit_df = res_rc$fit_df,
    ci_poly = res_rc$ci_poly,
    summary_row = res_rc$summary,
    title_text = rc,
    point_fill = point_col,
    line_type = 1,
    x_lab = expression(Delta * "CO"[2] * " (ppm)"),
    y_lab = "Effect size (lnRR)"
  )
  
  response_plots[[rc]] <- p_rc
  
  if (!is.null(p_rc)) {
    
    safe_name <- str_replace_all(rc, "[^A-Za-z0-9]+", "_")
    
    ggsave(
      filename = file.path(
        figure_dir,
        paste0("response_", safe_name, "_deltaCO2_plot.png")
      ),
      plot = p_rc,
      width = 6,
      height = 4.5,
      dpi = 600,
      bg = "white"
    )
    
    saveRDS(
      p_rc,
      file.path(
        rds_dir,
        paste0("response_", safe_name, "_deltaCO2_plot.rds")
      )
    )
  }
}

# ============================================================
# Combined response-category panel figure
# ============================================================

response_plots_nonnull <- response_plots[
  !vapply(response_plots, is.null, logical(1))
]

if (length(response_plots_nonnull) > 0) {
  
  p_response_panel <- wrap_plots(response_plots_nonnull, ncol = 2) +
    plot_annotation(
      title = NA,
      theme = theme(
        plot.title = element_text(
          family = fam,
          face = "bold",
          size = combined_title_size,
          hjust = 0.5
        )
      )
    )
  
  ggsave(
    filename = file.path(
      figure_dir,
      "response_category_deltaCO2_meta_regression_panels.png"
    ),
    plot = p_response_panel,
    width = 12,
    height = 12,
    dpi = 600,
    bg = "white"
  )
  
  saveRDS(
    p_response_panel,
    file.path(
      rds_dir,
      "response_category_deltaCO2_meta_regression_panels.rds"
    )
  )
}

# ============================================================
# Console output
# ============================================================

cat("\n==== Response category delta_CO2 meta-regression summary ====\n")
print(response_summary_tbl)

cat("\n==== Saved figures ====\n")
cat("Figures:", figure_dir, "\n")
cat("RDS    :", rds_dir, "\n")

if (length(response_plots_nonnull) > 0) {
  cat("Response-category figures saved for:\n")
  print(names(response_plots_nonnull))
}

if (exists("p_response_panel")) {
  cat(
    "Combined figure:",
    file.path(figure_dir, "response_category_deltaCO2_meta_regression_panels.png"),
    "\n"
  )
}

# ============================================================
# Show plots in console
# ============================================================

if (length(response_plots_nonnull) > 0) {
  invisible(lapply(response_plots_nonnull, print))
}

if (exists("p_response_panel")) {
  print(p_response_panel)
}

