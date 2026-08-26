# ============================================================
# 11B_Plot_lnrr_interaction.R
# Purpose:
#   Create one-panel lnRR plots using saved outputs from

# Notes:
#   The plotted estimates are from the single interaction model:
#
#     yi ~ 0 + pH_stressor + T_stressor + pH_T_interaction
#
#   Therefore:
#     pH  = model-estimated pH effect
#     T   = model-estimated temperature effect
#     TpH = model-predicted T + pH effect
#           = pH + T + pH:T
# ============================================================


# ============================================================
# Load packages
# ============================================================

suppressPackageStartupMessages(invisible(
  lapply(
    c(
      "dplyr", "ggplot2", "readr", "stringr",
      "showtext", "tibble", "grid", "tidyr",
      "purrr", "scales", "gridExtra"
    ),
    library,
    character.only = TRUE
  )
))


# ============================================================
# Fonts and theme
# ============================================================

font_ok <- TRUE

tryCatch(
  {
    font_add_google("Inter", "helv")
  },
  error = function(e) {
    font_ok <<- FALSE
    message("Could not load Google font Inter. Using default sans font instead.")
  }
)

showtext_auto()

fam <- if (font_ok) "helv" else "sans"

theme_set(theme_minimal(base_size = 18, base_family = fam))
theme_update(text = element_text(family = fam))


# ============================================================
# Script name and folders
# ============================================================

model_script_name <- "11A_Model_lnrr_interaction"
plot_script_name  <- "11B_Plot_lnrr_interaction"

model_rds_dir <- file.path("Outputs", "RDS", model_script_name)

figure_dir <- file.path("Outputs", "Figures", plot_script_name)
rds_dir    <- file.path("Outputs", "RDS", plot_script_name)

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)


# ============================================================
# Load saved model-ready data
# ============================================================

dat1_ml_path <- file.path(model_rds_dir, "dat1_with_lnrr_vi.rds")

if (!file.exists(dat1_ml_path)) {
  stop("Missing file: ", dat1_ml_path)
}

dat1_ml <- readRDS(dat1_ml_path)

cat("\nLoaded model-ready data from:", dat1_ml_path, "\n")
cat("Rows in dat1_ml:", nrow(dat1_ml), "\n\n")


# ============================================================
# Plot labels and levels
# ============================================================

response_levels <- c(
  "Growth",
  "Development",
  "Physiology",
  "Reproduction",
  "Survival"
)

subgroup_levels <- c(
  "Holoplankton",
  "Meroplankton"
)

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

treatment_levels <- c("pH", "T", "TpH")

treat_lab <- c(
  "pH"  = "pH",
  "T"   = "T",
  "TpH" = "pH:T"
)


# ============================================================
# Helper functions
# ============================================================

read_required_rds <- function(path) {
  if (!file.exists(path)) {
    stop("Missing file: ", path)
  }
  readRDS(path)
}


get_full_factorial_common_ids <- function(df) {
  df %>%
    mutate(Treatment_chr = as.character(Treatment)) %>%
    group_by(Common_id) %>%
    summarise(
      has_pH  = any(Treatment_chr == "pH",  na.rm = TRUE),
      has_T   = any(Treatment_chr == "T",   na.rm = TRUE),
      has_TpH = any(Treatment_chr == "TpH", na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      is_fully_factorial = has_pH & has_T & has_TpH
    )
}


keep_full_factorial_within_subset <- function(df) {
  factorial_ids <- get_full_factorial_common_ids(df) %>%
    filter(is_fully_factorial) %>%
    select(Common_id)
  
  df %>%
    semi_join(factorial_ids, by = "Common_id")
}


format_p_value <- function(pval) {
  case_when(
    is.na(pval) ~ NA_character_,
    pval < 0.001 ~ "<0.001",
    TRUE ~ sprintf("%.3f", pval)
  )
}


# ============================================================
# Function to create one plot scope
# ============================================================

make_interaction_plot_for_scope <- function(model_scope_to_plot) {
  
  # ----------------------------------------------------------
  # Choose input files
  # ----------------------------------------------------------
  
  if (model_scope_to_plot == "all_data") {
    
    predicted_path <- file.path(
      model_rds_dir,
      "interaction_predicted_effects_all_data.rds"
    )
    
    coefficients_path <- file.path(
      model_rds_dir,
      "interaction_coefficients_all_data.rds"
    )
    
    model_scope_label <- "All data"
    model_scope_file  <- "all_data"
    
  } else if (model_scope_to_plot == "full_factorial_only") {
    
    predicted_path <- file.path(
      model_rds_dir,
      "interaction_predicted_effects_full_factorial_only.rds"
    )
    
    coefficients_path <- file.path(
      model_rds_dir,
      "interaction_coefficients_full_factorial_only.rds"
    )
    
    model_scope_label <- "Full factorial only"
    model_scope_file  <- "full_factorial_only"
    
  } else {
    
    stop("model_scope_to_plot must be either 'all_data' or 'full_factorial_only'.")
  }
  
  
  # ----------------------------------------------------------
  # Load model outputs
  # ----------------------------------------------------------
  
  sum_df_final <- read_required_rds(predicted_path)
  interaction_coefficients <- read_required_rds(coefficients_path)
  
  cat("\n============================================================\n")
  cat("Creating plot for:", model_scope_label, "\n")
  cat("Predicted effects file:", predicted_path, "\n")
  cat("Interaction coefficients file:", coefficients_path, "\n")
  cat("Rows in sum_df_final:", nrow(sum_df_final), "\n")
  cat("============================================================\n\n")
  
  
  # ----------------------------------------------------------
  # Re-factor model summary labels
  # ----------------------------------------------------------
  
  sum_df_final <- sum_df_final %>%
    mutate(
      labels = factor(as.character(labels), levels = final_label_levels),
      Treatment = factor(as.character(Treatment), levels = treatment_levels)
    )
  
  
  # ----------------------------------------------------------
  # Raw points for plotting
  # ----------------------------------------------------------
  
  if (model_scope_to_plot == "all_data") {
    
    raw_response <- dat1_ml %>%
      filter(Response_category %in% response_levels, is.finite(yi)) %>%
      transmute(
        labels = factor(as.character(Response_category), levels = final_label_levels),
        yi_raw = yi,
        Treatment = factor(as.character(Treatment), levels = treatment_levels)
      )
    
    raw_subgroup <- dat1_ml %>%
      filter(!is.na(Broad_subgroup), is.finite(yi)) %>%
      transmute(
        labels = factor(as.character(Broad_subgroup), levels = final_label_levels),
        yi_raw = yi,
        Treatment = factor(as.character(Treatment), levels = treatment_levels)
      )
    
    raw_overall <- dat1_ml %>%
      filter(is.finite(yi)) %>%
      transmute(
        labels = factor("Overall", levels = final_label_levels),
        yi_raw = yi,
        Treatment = factor(as.character(Treatment), levels = treatment_levels)
      )
    
  } else {
    
    raw_response <- purrr::map_dfr(
      response_levels,
      function(rc) {
        dat1_ml %>%
          filter(Response_category == rc, is.finite(yi)) %>%
          keep_full_factorial_within_subset() %>%
          transmute(
            labels = factor(rc, levels = final_label_levels),
            yi_raw = yi,
            Treatment = factor(as.character(Treatment), levels = treatment_levels)
          )
      }
    )
    
    raw_subgroup <- purrr::map_dfr(
      subgroup_levels,
      function(subg) {
        dat1_ml %>%
          filter(Broad_subgroup == subg, is.finite(yi)) %>%
          keep_full_factorial_within_subset() %>%
          transmute(
            labels = factor(subg, levels = final_label_levels),
            yi_raw = yi,
            Treatment = factor(as.character(Treatment), levels = treatment_levels)
          )
      }
    )
    
    raw_overall <- dat1_ml %>%
      filter(is.finite(yi)) %>%
      keep_full_factorial_within_subset() %>%
      transmute(
        labels = factor("Overall", levels = final_label_levels),
        yi_raw = yi,
        Treatment = factor(as.character(Treatment), levels = treatment_levels)
      )
  }
  
  raw_all_final <- bind_rows(
    raw_response,
    raw_subgroup,
    raw_overall
  )
  
  
  # ----------------------------------------------------------
  # Treatment-specific k values for labels above plot
  # ----------------------------------------------------------
  
  raw_k_table <- raw_all_final %>%
    count(labels, Treatment, name = "k_treatment")
  
  sum_df_final <- sum_df_final %>%
    left_join(
      raw_k_table,
      by = c("labels", "Treatment")
    ) %>%
    mutate(
      k_treatment = ifelse(is.na(k_treatment), NA_integer_, k_treatment),
      k_model = k
    )
  
  
  # ----------------------------------------------------------
  # Save model-predicted effects table
  # ----------------------------------------------------------
  
  predicted_table_plot <- sum_df_final %>%
    mutate(
      Category = as.character(labels),
      Treatment_label = recode(
        as.character(Treatment),
        "pH"  = "pH",
        "T"   = "T",
        "TpH" = "pH:T"
      ),
      
      `Mean lnRR` = sprintf("%.3f", mean),
      `95% CI` = sprintf("%.3f to %.3f", ci_low, ci_high),
      
      `p-value` = format_p_value(pval),
      
      `% change` = sprintf("%.1f%%", (exp(mean) - 1) * 100),
      `% change 95% CI` = sprintf(
        "%.1f%% to %.1f%%",
        (exp(ci_low) - 1) * 100,
        (exp(ci_high) - 1) * 100
      )
    ) %>%
    select(
      Category,
      Treatment = Treatment_label,
      `k treatment` = k_treatment,
      `k model` = k_model,
      `Mean lnRR`,
      `95% CI`,
      `p-value`,
      `% change`,
      `% change 95% CI`
    ) %>%
    arrange(
      factor(Category, levels = final_label_levels),
      factor(Treatment, levels = c("pH", "T", "pH:T"))
    )
  
  write.csv(
    predicted_table_plot,
    file.path(
      figure_dir,
      paste0("interaction_model_predicted_effects_table_", model_scope_file, ".csv")
    ),
    row.names = FALSE
  )
  
  predicted_table_plot_for_grob <- predicted_table_plot
  predicted_table_plot_for_grob[] <- lapply(predicted_table_plot_for_grob, as.character)
  
  tbl_grob <- gridExtra::tableGrob(
    predicted_table_plot_for_grob,
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
    filename = file.path(
      figure_dir,
      paste0("interaction_model_predicted_effects_table_", model_scope_file, ".png")
    ),
    width = 3800,
    height = 2200,
    res = 300,
    bg = "white"
  )
  
  grid.newpage()
  
  grid.text(
    paste0(
      "Interaction model-predicted effects, 95% confidence intervals, p-values, and percent change: ",
      model_scope_label
    ),
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
  
  
  # ----------------------------------------------------------
  # Save interaction coefficient table
  # ----------------------------------------------------------
  
  interaction_terms_table <- interaction_coefficients %>%
    filter(term == "pH_T_interaction") %>%
    mutate(
      labels = factor(as.character(labels), levels = final_label_levels),
      Category = as.character(labels),
      `Interaction term` = "pH:T",
      `Mean lnRR` = sprintf("%.3f", estimate),
      `95% CI` = sprintf("%.3f to %.3f", ci_low, ci_high),
      `p-value` = format_p_value(pval),
      `% change` = sprintf("%.1f%%", (exp(estimate) - 1) * 100),
      `% change 95% CI` = sprintf(
        "%.1f%% to %.1f%%",
        (exp(ci_low) - 1) * 100,
        (exp(ci_high) - 1) * 100
      )
    ) %>%
    select(
      Category,
      `Interaction term`,
      k,
      n_studies,
      n_common_id,
      `Mean lnRR`,
      `95% CI`,
      `p-value`,
      `% change`,
      `% change 95% CI`
    ) %>%
    arrange(factor(Category, levels = final_label_levels))
  
  write.csv(
    interaction_terms_table,
    file.path(
      figure_dir,
      paste0("interaction_terms_only_table_", model_scope_file, ".csv")
    ),
    row.names = FALSE
  )
  
  
  # ----------------------------------------------------------
  # Plot aesthetics
  # ----------------------------------------------------------
  
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
  
  pd <- position_dodge(width = 0.55)
  
  pjd <- position_jitterdodge(
    jitter.width = 0.22,
    jitter.height = 0,
    dodge.width = 0.55
  )
  
  subgroup_start_idx <- which(final_label_levels == "Holoplankton")
  subgroup_end_idx   <- which(final_label_levels == "Meroplankton")
  subgroup_xmin      <- subgroup_start_idx - 0.5
  subgroup_xmax      <- subgroup_end_idx + 0.5
  
  overall_idx  <- which(final_label_levels == "Overall")
  overall_xmin <- overall_idx - 0.5
  overall_xmax <- overall_idx + 0.5
  
  summary_fill <- "grey95"
  
  response_subgroup_separator <- 5.5
  subgroup_overall_separator  <- 7.5
  
  
  # ----------------------------------------------------------
  # Final one-panel plot
  # ----------------------------------------------------------
  
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
    
    geom_vline(
      xintercept = response_subgroup_separator,
      linetype = "solid",
      colour = "black",
      linewidth = 0.4
    ) +
    
    geom_vline(
      xintercept = 6.5,
      linetype = "dotted",
      colour = "grey70",
      linewidth = 0.4
    ) +
    
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
      aes(x = labels, y = Inf, label = k_treatment, colour = Treatment),
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
    labs(
      x = NULL
    ) +
    theme(
      text = element_text(family = fam, colour = "grey10"),
      
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      
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
      panel.grid.major   = element_blank(),
      panel.grid.minor   = element_blank(),
      axis.ticks         = element_line(colour = "grey10", linewidth = 0.4),
      axis.ticks.length  = unit(4, "pt"),
      legend.position    = "top",
      legend.margin      = margin(0, 0, 0, 0),
      legend.box.margin  = margin(0, 0, 0, 0),
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
  
  
  # ----------------------------------------------------------
  # Save final plot as PNG only
  # ----------------------------------------------------------
  
  ggsave(
    filename = file.path(
      figure_dir,
      paste0(
        "zooplankton_interaction_model_predicted_treatment_one_panel_",
        model_scope_file,
        ".png"
      )
    ),
    plot = p_final,
    width = 12.1,
    height = 5.3,
    units = "in",
    dpi = 300,
    bg = "white"
  )
  
  saveRDS(
    p_final,
    file.path(
      rds_dir,
      paste0(
        "zooplankton_interaction_model_predicted_treatment_one_panel_",
        model_scope_file,
        ".rds"
      )
    )
  )
  
  saveRDS(
    sum_df_final,
    file.path(
      rds_dir,
      paste0("plotting_summary_interaction_predicted_effects_", model_scope_file, ".rds")
    )
  )
  
  saveRDS(
    raw_all_final,
    file.path(
      rds_dir,
      paste0("plotting_raw_points_interaction_", model_scope_file, ".rds")
    )
  )
  
  cat("\n==== Interaction model predicted effects plotted:", model_scope_label, "====\n")
  print(sum_df_final)
  
  cat("\n==== Interaction terms only:", model_scope_label, "====\n")
  print(interaction_terms_table)
  
  cat("\nSaved PNG plot for:", model_scope_label, "\n")
  
  return(
    list(
      plot = p_final,
      summary = sum_df_final,
      raw = raw_all_final,
      predicted_table = predicted_table_plot,
      interaction_terms = interaction_terms_table
    )
  )
}


# ============================================================
# Produce BOTH plots in one go
# ============================================================

model_scopes_to_plot <- c(
  "all_data",
  "full_factorial_only"
)

plot_outputs <- purrr::map(
  model_scopes_to_plot,
  make_interaction_plot_for_scope
)

names(plot_outputs) <- model_scopes_to_plot

saveRDS(
  plot_outputs,
  file.path(rds_dir, "all_interaction_plot_outputs.rds")
)


# ============================================================
# Console summary
# ============================================================

cat("\n============================================================\n")
cat("DONE: Saved BOTH interaction-model PNG plots without plot titles.\n")
cat("No PDF files were created.\n")
cat("TpH is displayed as pH:T in legends and table summaries.\n")
cat("============================================================\n")

cat("\nMain scientific plots saved as:\n")
cat(
  file.path(
    figure_dir,
    "zooplankton_interaction_model_predicted_treatment_one_panel_all_data.png"
  ),
  "\n"
)
cat(
  file.path(
    figure_dir,
    "zooplankton_interaction_model_predicted_treatment_one_panel_full_factorial_only.png"
  ),
  "\n"
)

cat("\nTable figures saved in:\n")
cat("Figures:", figure_dir, "\n")
cat("RDS    :", rds_dir, "\n")