# ============================================================
# 12D_Plot_lnrr_subgroup.R
# Purpose:
#   Plot Holoplankton vs Meroplankton model-predicted lnRR estimates
#   from the single interaction models fitted in 12C.
#
#   This is the plotting companion to 12C.
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
# Script name and output folders
# ============================================================

model_script_name <- "12C_Model_lnrr_subgroup"
plot_script_name  <- "12D_Plot_lnrr_subgroup"

model_rds_dir <- file.path("Outputs", "RDS", model_script_name)

figure_dir <- file.path("Outputs", "Figures", plot_script_name)
table_dir  <- file.path("Outputs", "Tables", plot_script_name)
rds_dir    <- file.path("Outputs", "RDS", plot_script_name)

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# Helper
# ============================================================

safe_read_rds <- function(path, message_if_missing = NULL) {
  if (!file.exists(path)) {
    if (is.null(message_if_missing)) {
      stop("Could not find: ", path)
    } else {
      stop(message_if_missing, "\nMissing file: ", path)
    }
  }
  
  readRDS(path)
}

safe_text <- function(x) {
  x %>%
    as.character() %>%
    stringr::str_replace_all("[^A-Za-z0-9]+", "_") %>%
    stringr::str_replace_all("^_+|_+$", "")
}

# ============================================================
# Load saved outputs from 12C
# ============================================================

dat1_subgroup <- safe_read_rds(
  file.path(model_rds_dir, "dat1_subgroup_interaction_model_ready.rds"),
  "Please run Scripts/12C_Model_lnrr_subgroup.R first."
)

interaction_predicted_effects_all_data <- safe_read_rds(
  file.path(model_rds_dir, "interaction_predicted_effects_all_data.rds"),
  "Please run Scripts/12C_Model_lnrr_subgroup.R first."
)

interaction_predicted_effects_full_factorial_only <- safe_read_rds(
  file.path(model_rds_dir, "interaction_predicted_effects_full_factorial_only.rds"),
  "Please run Scripts/12C_Model_lnrr_subgroup.R first."
)

interaction_coefficients_all_data <- safe_read_rds(
  file.path(model_rds_dir, "interaction_coefficients_all_data.rds"),
  "Please run Scripts/12C_Model_lnrr_subgroup.R first."
)

interaction_coefficients_full_factorial_only <- safe_read_rds(
  file.path(model_rds_dir, "interaction_coefficients_full_factorial_only.rds"),
  "Please run Scripts/12C_Model_lnrr_subgroup.R first."
)

plot_levels_subgroup_interaction <- safe_read_rds(
  file.path(model_rds_dir, "plot_levels_subgroup_interaction.rds"),
  "Please run Scripts/12C_Model_lnrr_subgroup.R first."
)

response_levels <- plot_levels_subgroup_interaction$response_levels
category_levels <- plot_levels_subgroup_interaction$category_levels
treatment_levels <- plot_levels_subgroup_interaction$treatment_levels
treatment_label_levels <- plot_levels_subgroup_interaction$treatment_label_levels
subgroup_levels <- plot_levels_subgroup_interaction$subgroup_levels
model_scopes <- plot_levels_subgroup_interaction$model_scopes

if (is.null(model_scopes)) {
  model_scopes <- c("all_data", "full_factorial_only")
}

cat("\nLoaded saved model outputs from:", model_rds_dir, "\n")
cat("Rows in dat1_subgroup:", nrow(dat1_subgroup), "\n")
cat("Rows in all-data predicted effects:", nrow(interaction_predicted_effects_all_data), "\n")
cat("Rows in full-factorial predicted effects:", nrow(interaction_predicted_effects_full_factorial_only), "\n\n")

# ============================================================
# Combine model outputs
# ============================================================

interaction_predicted_effects_all_scopes <- bind_rows(
  interaction_predicted_effects_all_data,
  interaction_predicted_effects_full_factorial_only
) %>%
  mutate(
    analysis_scope = factor(as.character(analysis_scope), levels = model_scopes),
    Category = factor(as.character(Category), levels = category_levels),
    Response_category = factor(as.character(Response_category), levels = response_levels),
    Treatment = factor(as.character(Treatment), levels = treatment_levels),
    Treatment_label = factor(as.character(Treatment_label), levels = treatment_label_levels),
    Broad_subgroup = factor(as.character(Broad_subgroup), levels = subgroup_levels)
  )

interaction_coefficients_all_scopes <- bind_rows(
  interaction_coefficients_all_data,
  interaction_coefficients_full_factorial_only
) %>%
  mutate(
    analysis_scope = factor(as.character(analysis_scope), levels = model_scopes),
    Category = factor(as.character(Category), levels = category_levels),
    Response_category = factor(as.character(Response_category), levels = response_levels),
    Broad_subgroup = factor(as.character(Broad_subgroup), levels = subgroup_levels)
  )

# ============================================================
# Make sure raw data factors are ordered correctly
# ============================================================

dat1_subgroup <- dat1_subgroup %>%
  mutate(
    Response_category = factor(as.character(Response_category), levels = response_levels),
    Treatment = factor(as.character(Treatment), levels = treatment_levels),
    Treatment_label = factor(as.character(Treatment_label), levels = treatment_label_levels),
    Broad_subgroup = factor(as.character(Broad_subgroup), levels = subgroup_levels)
  )

# ============================================================
# Full-factorial Common_id helper
# ============================================================

get_full_factorial_common_ids <- function(df) {
  if (nrow(df) == 0) {
    return(character(0))
  }
  
  df %>%
    filter(!is.na(Common_id), !is.na(Treatment)) %>%
    distinct(Common_id, Treatment) %>%
    group_by(Common_id) %>%
    summarise(
      has_pH  = any(as.character(Treatment) == "pH"),
      has_T   = any(as.character(Treatment) == "T"),
      has_TpH = any(as.character(Treatment) == "TpH"),
      .groups = "drop"
    ) %>%
    filter(has_pH, has_T, has_TpH) %>%
    pull(Common_id) %>%
    as.character()
}

# ============================================================
# Save model-predicted effects, interaction coefficients, and
# percent change as table figures
# ============================================================

make_predicted_table <- function(pred_df, scope_name) {
  pred_table <- pred_df %>%
    mutate(
      analysis_scope = as.character(analysis_scope),
      Category = as.character(Category),
      Broad_subgroup = as.character(Broad_subgroup),
      Treatment = as.character(Treatment_label),
      
      `Mean lnRR` = ifelse(
        is.finite(mean),
        sprintf("%.3f", mean),
        NA_character_
      ),
      
      `95% CI` = ifelse(
        is.finite(ci_low) & is.finite(ci_high),
        sprintf("%.3f to %.3f", ci_low, ci_high),
        NA_character_
      ),
      
      `p-value` = case_when(
        is.na(pval) ~ NA_character_,
        pval < 0.001 ~ "<0.001",
        TRUE ~ sprintf("%.3f", pval)
      ),
      
      `% change` = ifelse(
        is.finite(mean),
        sprintf("%.1f%%", (exp(mean) - 1) * 100),
        NA_character_
      ),
      
      `% change 95% CI` = ifelse(
        is.finite(ci_low) & is.finite(ci_high),
        sprintf(
          "%.1f%% to %.1f%%",
          (exp(ci_low) - 1) * 100,
          (exp(ci_high) - 1) * 100
        ),
        NA_character_
      )
    ) %>%
    select(
      Category,
      Broad_subgroup,
      Treatment,
      k,
      n_studies,
      n_common_id,
      n_full_factorial_common_id,
      n_pH,
      n_T,
      n_TpH,
      `Mean lnRR`,
      `95% CI`,
      `p-value`,
      `% change`,
      `% change 95% CI`,
      model_status,
      model_message
    ) %>%
    arrange(
      factor(Category, levels = category_levels),
      factor(Broad_subgroup, levels = subgroup_levels),
      factor(Treatment, levels = treatment_label_levels)
    )
  
  write.csv(
    pred_table,
    file.path(table_dir, paste0("interaction_predicted_effects_table_", scope_name, ".csv")),
    row.names = FALSE
  )
  
  saveRDS(
    pred_table,
    file.path(rds_dir, paste0("interaction_predicted_effects_table_", scope_name, ".rds"))
  )
  
  pred_table_char <- pred_table
  pred_table_char[] <- lapply(pred_table_char, as.character)
  
  tbl_grob <- gridExtra::tableGrob(
    pred_table_char,
    rows = NULL,
    theme = gridExtra::ttheme_minimal(
      base_size = 9,
      base_family = fam,
      padding = unit(c(3, 3), "mm"),
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
  
  table_title <- case_when(
    scope_name == "all_data" ~ "All data: subgroup interaction-model predicted effects, confidence intervals, p-values, and percent change",
    scope_name == "full_factorial_only" ~ "Full factorial only: subgroup interaction-model predicted effects, confidence intervals, p-values, and percent change",
    TRUE ~ paste0(scope_name, ": subgroup interaction-model predicted effects")
  )
  
  png(
    filename = file.path(
      figure_dir,
      paste0("interaction_model_predicted_effects_table_subgroup_", scope_name, ".png")
    ),
    width = 4200,
    height = 1800,
    res = 300,
    bg = "white"
  )
  
  grid.newpage()
  
  grid.text(
    table_title,
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
  
  pred_table
}

# ============================================================
# Plot aesthetics following 12B style
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

subgroup_shapes <- c(
  "Holoplankton" = 21,
  "Meroplankton" = 24
)

stream_positions <- tidyr::expand_grid(
  Broad_subgroup = factor(subgroup_levels, levels = subgroup_levels),
  Treatment_label = factor(treatment_label_levels, levels = treatment_label_levels)
) %>%
  mutate(
    stream_id = paste(Broad_subgroup, Treatment_label, sep = "_"),
    offset = case_when(
      Broad_subgroup == "Holoplankton" & Treatment_label == "pH"     ~ -0.30,
      Broad_subgroup == "Holoplankton" & Treatment_label == "T"      ~ -0.18,
      Broad_subgroup == "Holoplankton" & Treatment_label == "T + pH" ~ -0.06,
      Broad_subgroup == "Meroplankton" & Treatment_label == "pH"     ~  0.06,
      Broad_subgroup == "Meroplankton" & Treatment_label == "T"      ~  0.18,
      Broad_subgroup == "Meroplankton" & Treatment_label == "T + pH" ~  0.30,
      TRUE ~ 0
    )
  )

overall_idx  <- which(category_levels == "Overall")
overall_xmin <- overall_idx - 0.5
overall_xmax <- overall_idx + 0.5

summary_fill <- "grey95"

response_overall_separator <- which(category_levels == "Survival") + 0.5

# ============================================================
# Raw points builder by model scope
# ============================================================

make_raw_points_for_scope <- function(scope_name) {
  analysis_grid <- tidyr::expand_grid(
    Category = category_levels,
    Broad_subgroup = subgroup_levels
  )
  
  raw_df <- purrr::map_dfr(
    seq_len(nrow(analysis_grid)),
    function(i) {
      category_name  <- analysis_grid$Category[i]
      subgroup_name <- analysis_grid$Broad_subgroup[i]
      
      df_base <- dat1_subgroup %>%
        filter(as.character(Broad_subgroup) == subgroup_name)
      
      if (category_name != "Overall") {
        df_base <- df_base %>%
          filter(as.character(Response_category) == category_name)
      }
      
      if (scope_name == "full_factorial_only") {
        full_ids <- get_full_factorial_common_ids(df_base)
        
        df_base <- df_base %>%
          filter(as.character(Common_id) %in% full_ids)
      }
      
      if (nrow(df_base) == 0) {
        return(tibble())
      }
      
      df_base %>%
        filter(is.finite(yi)) %>%
        transmute(
          analysis_scope = scope_name,
          Category = factor(category_name, levels = category_levels),
          yi_raw = yi,
          Treatment = factor(as.character(Treatment), levels = treatment_levels),
          Treatment_label = factor(as.character(Treatment_label), levels = treatment_label_levels),
          Broad_subgroup = factor(as.character(Broad_subgroup), levels = subgroup_levels),
          Study_ID = Study_ID,
          ES_ID = ES_ID,
          Common_id = Common_id
        )
    }
  ) %>%
    left_join(category_positions, by = "Category") %>%
    left_join(stream_positions, by = c("Broad_subgroup", "Treatment_label"))
  
  set.seed(123)
  
  raw_df %>%
    mutate(
      x_pos = category_x + offset,
      x_jitter = x_pos + runif(n(), min = -0.025, max = 0.025)
    )
}

# ============================================================
# Plot builder
# ============================================================

make_subgroup_interaction_plot <- function(scope_name) {
  pred_scope_all <- interaction_predicted_effects_all_scopes %>%
    filter(as.character(analysis_scope) == scope_name)
  
  pred_scope_fitted <- pred_scope_all %>%
    filter(
      model_status == "fitted",
      is.finite(mean),
      is.finite(se),
      is.finite(ci_low),
      is.finite(ci_high)
    )
  
  raw_scope <- make_raw_points_for_scope(scope_name)
  
  raw_k_table <- raw_scope %>%
    group_by(Category, Broad_subgroup, Treatment_label) %>%
    summarise(
      k_treatment = n(),
      .groups = "drop"
    )
  
  sum_df_plot <- pred_scope_fitted %>%
    mutate(
      Category = factor(as.character(Category), levels = category_levels),
      Treatment_label = factor(as.character(Treatment_label), levels = treatment_label_levels),
      Broad_subgroup = factor(as.character(Broad_subgroup), levels = subgroup_levels)
    ) %>%
    left_join(category_positions, by = "Category") %>%
    left_join(stream_positions, by = c("Broad_subgroup", "Treatment_label")) %>%
    left_join(raw_k_table, by = c("Category", "Broad_subgroup", "Treatment_label")) %>%
    mutate(
      x_pos = category_x + offset,
      k_label = ifelse(is.na(k_treatment), "", as.character(k_treatment))
    )
  
  p <- ggplot() +
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
      xintercept = response_overall_separator,
      linetype = "solid",
      colour = "black",
      linewidth = 0.4
    ) +
    
    geom_point(
      data = raw_scope,
      aes(
        x = x_jitter,
        y = yi_raw,
        fill = Treatment_label,
        shape = Broad_subgroup
      ),
      colour = "grey25",
      stroke = 0.10,
      alpha = 0.18,
      size = 2.0,
      show.legend = FALSE
    ) +
    
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
    
    geom_point(
      data = sum_df_plot,
      aes(
        x = x_pos,
        y = mean,
        colour = Treatment_label,
        fill = Treatment_label,
        shape = Broad_subgroup
      ),
      size = 3.1,
      stroke = 0.85
    ) +
    
    geom_text(
      data = sum_df_plot,
      aes(
        x = x_pos,
        y = Inf,
        label = k_label,
        colour = Treatment_label
      ),
      vjust = -1.2,
      family = fam,
      size = 7.0,
      show.legend = FALSE
    ) +
    
    scale_colour_manual(
      values = col_trt_line,
      labels = c("pH", "T", "T:pH"),
      name = NULL
    ) +
    scale_fill_manual(
      values = fill_trt,
      labels = c("pH", "T", "T:pH"),
      name = NULL
    ) +
    scale_shape_manual(
      values = subgroup_shapes,
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
      legend.text = element_text(size = 28),
      
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
  
  # Save plot data
  write.csv(
    raw_scope,
    file.path(table_dir, paste0("raw_points_subgroup_", scope_name, ".csv")),
    row.names = FALSE
  )
  
  write.csv(
    sum_df_plot,
    file.path(table_dir, paste0("plot_summary_subgroup_interaction_", scope_name, ".csv")),
    row.names = FALSE
  )
  
  saveRDS(
    raw_scope,
    file.path(rds_dir, paste0("raw_points_subgroup_", scope_name, ".rds"))
  )
  
  saveRDS(
    sum_df_plot,
    file.path(rds_dir, paste0("plot_summary_subgroup_interaction_", scope_name, ".rds"))
  )
  
  saveRDS(
    p,
    file.path(rds_dir, paste0("p_subgroup_interaction_", scope_name, ".rds"))
  )
  
  # Save figure as PNG only
  ggsave(
    filename = file.path(
      figure_dir,
      paste0("zooplankton_subgroup_interaction_model_predicted_effects_", scope_name, ".png")
    ),
    plot = p,
    width = 12.1,
    height = 5.3,
    units = "in",
    dpi = 300,
    bg = "white"
  )
  
  list(
    plot = p,
    raw_points = raw_scope,
    summary_points = sum_df_plot
  )
}

# ============================================================
# Create tables and plots for each scope
# ============================================================

plot_outputs <- list()
table_outputs <- list()

for (scope_name in model_scopes) {
  cat("\n===============================\n")
  cat("Making plot for:", scope_name, "\n")
  cat("===============================\n")
  
  pred_scope <- interaction_predicted_effects_all_scopes %>%
    filter(as.character(analysis_scope) == scope_name)
  
  table_outputs[[scope_name]] <- make_predicted_table(
    pred_df = pred_scope,
    scope_name = scope_name
  )
  
  plot_outputs[[scope_name]] <- make_subgroup_interaction_plot(
    scope_name = scope_name
  )
}

# ============================================================
# Save combined plotting objects
# ============================================================

saveRDS(
  plot_outputs,
  file.path(rds_dir, "plot_outputs_subgroup_interaction_all_scopes.rds")
)

saveRDS(
  table_outputs,
  file.path(rds_dir, "table_outputs_subgroup_interaction_all_scopes.rds")
)

write.csv(
  interaction_predicted_effects_all_scopes,
  file.path(table_dir, "interaction_predicted_effects_all_scopes_for_plotting.csv"),
  row.names = FALSE
)

write.csv(
  interaction_coefficients_all_scopes,
  file.path(table_dir, "interaction_coefficients_all_scopes_for_plotting.csv"),
  row.names = FALSE
)

saveRDS(
  interaction_predicted_effects_all_scopes,
  file.path(rds_dir, "interaction_predicted_effects_all_scopes_for_plotting.rds")
)

saveRDS(
  interaction_coefficients_all_scopes,
  file.path(rds_dir, "interaction_coefficients_all_scopes_for_plotting.rds")
)

# ============================================================
# Console prints
# ============================================================

cat("\n==== Model status counts ====\n")
print(table(
  interaction_predicted_effects_all_scopes$analysis_scope,
  interaction_predicted_effects_all_scopes$model_status,
  useNA = "ifany"
))

cat("\n==== Plot files saved ====\n")
cat("Figures:", figure_dir, "\n")
cat("Tables :", table_dir, "\n")
cat("RDS    :", rds_dir, "\n")

cat("\nMain figure outputs:\n")
cat(
  "- ",
  file.path(
    figure_dir,
    "zooplankton_subgroup_interaction_model_predicted_effects_all_data.png"
  ),
  "\n",
  sep = ""
)
cat(
  "- ",
  file.path(
    figure_dir,
    "zooplankton_subgroup_interaction_model_predicted_effects_full_factorial_only.png"
  ),
  "\n",
  sep = ""
)

cat("\nDone: 12D subgroup interaction plots completed. PNG only, no plot titles, no PDF outputs.\n")
