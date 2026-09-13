# ============================================================
# 07_Regression_analysis_interaction.R
# Purpose:
#   Use the shared cleaned zooplankton dataset from
#   Scripts/00_Dependencies.R, calculate lnRR and sampling variance,
#   fit TpH-only multilevel interaction meta-regressions with
#   delta_T * delta_CO2, and create interaction plots.
#
# Input:
#   dat1 object from Scripts/00_Dependencies.R
#
# Outputs:
#   Outputs/Tables/07_Regression_analysis_interaction/
#   Outputs/RDS/07_Regression_analysis_interaction/
#   Outputs/Figures/07_Regression_analysis_interaction/
# ============================================================

suppressPackageStartupMessages(invisible(
  lapply(
    c(
      "dplyr", "ggplot2", "patchwork", "readr", "metafor", "stringr",
      "showtext", "tibble", "grid", "Matrix", "matrixcalc",
      "tidyr", "purrr", "scales"
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
# Fonts and theme
# ============================================================

font_add_google("Inter", "helv")
showtext_auto()
fam <- "helv"

theme_set(theme_minimal(base_size = 18, base_family = fam))
theme_update(text = element_text(family = fam))

theme_meta <- function() {
  theme(
    panel.background = element_rect(fill = "white", colour = "black", linewidth = 0.7),
    panel.grid = element_blank(),
    legend.position = "right",
    plot.margin = margin(8, 8, 8, 8),
    plot.background = element_blank(),
    axis.ticks = element_line(linewidth = 0.3, colour = "black"),
    axis.ticks.length = unit(-0.15, "lines"),
    axis.title.y = element_text(size = 13, margin = margin(0, 6, 0, 0), face = "bold"),
    axis.title.x = element_text(size = 13, margin = margin(6, 0, 0, 0), face = "bold"),
    axis.text.y = element_text(size = 11, margin = margin(0, 5, 0, 0), colour = "#000000"),
    axis.text.x = element_text(size = 11, margin = margin(5, 0, 0, 0), colour = "#000000"),
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 11)
  )
}

# ============================================================
# Script name and output folders
# ============================================================

script_name <- "07_Regression_analysis_interaction"

table_dir  <- file.path("Outputs", "Tables", script_name)
rds_dir    <- file.path("Outputs", "RDS", script_name)
figure_dir <- file.path("Outputs", "Figures", script_name)

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

if (!dir.exists(table_dir)) stop(paste("Table directory could not be created:", table_dir))
if (!dir.exists(rds_dir)) stop(paste("RDS directory could not be created:", rds_dir))
if (!dir.exists(figure_dir)) stop(paste("Figure directory could not be created:", figure_dir))

cat("\nWorking directory:\n", getwd(), "\n")
cat("\nTable directory:\n", table_dir, "\n")
cat("RDS directory:\n", rds_dir, "\n")
cat("Figure directory:\n", figure_dir, "\n")

# ============================================================
# Check shared cleaned dataset
# ============================================================

if (!exists("dat1")) {
  stop("dat1 not found. Check that Scripts/00_Dependencies.R creates dat1.")
}

if (!"delta_T" %in% names(dat1)) {
  stop("delta_T not found in dat1. Check that Scripts/00_Dependencies.R creates delta_T.")
}

if (!"delta_CO2" %in% names(dat1)) {
  stop("delta_CO2 not found in dat1. Check that Scripts/00_Dependencies.R creates delta_CO2.")
}

if (!"Response_category" %in% names(dat1)) {
  stop("Response_category not found in dat1. Check that Scripts/00_Dependencies.R creates Response_category.")
}

if (!"Broad_subgroup" %in% names(dat1)) {
  stop("Broad_subgroup not found in dat1. Check that Scripts/00_Dependencies.R creates Broad_subgroup.")
}

cat("\n==== Shared cleaned dataset ====\n")
cat("Rows in dat1:", nrow(dat1), "\n")
cat("Unique ES_ID:", n_distinct(dat1$ES_ID), "\n")
cat("Unique Study_ID:", n_distinct(dat1$Study_ID), "\n")
cat("Unique Common_id:", n_distinct(dat1$Common_id), "\n")
cat("ES_ID unique? ", n_distinct(dat1$ES_ID) == nrow(dat1), "\n\n")

# ============================================================
# Prepare shared dataset for this script
# ============================================================

dat1 <- dat1 %>%
  mutate(
    Treatment = factor(Treatment, levels = c("pH", "T", "TpH")),
    Broad_subgroup = str_squish(as.character(Broad_subgroup)),
    Broad_subgroup = case_when(
      str_detect(Broad_subgroup, regex("mero", ignore_case = TRUE)) ~ "Meroplankton",
      str_detect(Broad_subgroup, regex("holo", ignore_case = TRUE)) ~ "Holoplankton",
      TRUE ~ Broad_subgroup
    ),
    Broad_subgroup = factor(Broad_subgroup, levels = c("Meroplankton", "Holoplankton")),
    Response_category = str_squish(as.character(Response_category)),
    delta_T = suppressWarnings(as.numeric(delta_T)),
    delta_CO2 = suppressWarnings(as.numeric(delta_CO2))
  ) %>%
  filter(
    is.finite(delta_T), delta_T > 0, delta_T <= 10,
    is.finite(delta_CO2), delta_CO2 > 0, delta_CO2 < 5000
  )

# ============================================================
# Compute lnRR and vi
# ============================================================

escalc_out <- metafor::escalc(
  measure = "ROM",
  m1i  = Biota_T,
  sd1i = Biota_T_sd,
  n1i  = N_T,
  m2i  = Biota_C,
  sd2i = Biota_C_sd,
  n2i  = N_C,
  data = dat1
)

dat1_ml <- dat1 %>%
  mutate(
    yi = escalc_out$yi,
    vi = escalc_out$vi
  ) %>%
  filter(is.finite(yi), is.finite(vi), vi > 0)

saveRDS(
  dat1_ml,
  file.path(rds_dir, "dat1_with_lnrr_vi.rds")
)

write_csv(
  dat1_ml,
  file.path(table_dir, "dat1_with_lnrr_vi.csv")
)

# ============================================================
# Keep only TpH treatment
# ============================================================

dat_interaction <- dat1_ml %>%
  filter(Treatment == "TpH") %>%
  filter(
    is.finite(delta_T), delta_T > 0, delta_T <= 10,
    is.finite(delta_CO2), delta_CO2 > 0, delta_CO2 < 5000,
    is.finite(yi), is.finite(vi), vi > 0
  )

# ============================================================
# Trim extreme vi values
# ============================================================

q_low  <- quantile(dat_interaction$vi, 0.01, na.rm = TRUE)
q_high <- quantile(dat_interaction$vi, 0.99, na.rm = TRUE)

dat_interaction <- dat_interaction %>%
  filter(vi >= q_low, vi <= q_high)

saveRDS(
  dat_interaction,
  file.path(rds_dir, "dat_interaction_TpH_only_trimmed.rds")
)

write_csv(
  dat_interaction,
  file.path(table_dir, "dat_interaction_TpH_only_trimmed.csv")
)

cat("\n==== TpH interaction meta-regression data ====\n")
cat("Rows available:", nrow(dat_interaction), "\n")
cat("Unique studies:", dplyr::n_distinct(dat_interaction$Study_ID), "\n")
cat("Unique effect sizes:", dplyr::n_distinct(dat_interaction$ES_ID), "\n")
cat("Unique Common_id:", dplyr::n_distinct(dat_interaction$Common_id), "\n")
cat("delta_T range:", min(dat_interaction$delta_T, na.rm = TRUE), "to", max(dat_interaction$delta_T, na.rm = TRUE), "\n")
cat("delta_CO2 range:", min(dat_interaction$delta_CO2, na.rm = TRUE), "to", max(dat_interaction$delta_CO2, na.rm = TRUE), "\n")
cat("vi range:", min(dat_interaction$vi, na.rm = TRUE), "to", max(dat_interaction$vi, na.rm = TRUE), "\n")

# ============================================================
# VCV matrix builder
# ============================================================

make_V_matrix <- function(df) {
  split_list <- split(df, df$Common_id)
  
  block_list <- lapply(split_list, function(x) {
    n <- nrow(x)
    
    shared_cov <- (x$Biota_C_sd[1]^2) / (x$N_C[1] * x$Biota_C[1]^2)
    
    v <- matrix(shared_cov, nrow = n, ncol = n)
    diag(v) <- x$vi
    
    v
  })
  
  V <- as.matrix(Matrix::bdiag(block_list))
  
  if (!isSymmetric(V)) {
    V <- (V + t(V)) / 2
  }
  
  if (!matrixcalc::is.positive.definite(V)) {
    V <- as.matrix(Matrix::nearPD(V, ensureSymmetry = TRUE)$mat)
  }
  
  V
}

# ============================================================
# Fit interaction meta-regression
# ============================================================

fit_metareg_interaction <- function(df) {
  df <- df %>%
    arrange(Study_ID, ES_ID) %>%
    mutate(
      Study_ID = factor(Study_ID),
      ES_ID = factor(ES_ID),
      Common_id = factor(Common_id)
    )
  
  if (nrow(df) < 4) return(NULL)
  if (length(unique(df$delta_T)) < 2) return(NULL)
  if (length(unique(df$delta_CO2)) < 2) return(NULL)
  if (dplyr::n_distinct(df$Study_ID) < 2) return(NULL)
  
  V_matrix <- make_V_matrix(df)
  
  fit1 <- tryCatch(
    rma.mv(
      yi = yi,
      V = V_matrix,
      mods = ~ delta_T * delta_CO2,
      random = ~ 1 | Study_ID/ES_ID,
      data = df,
      method = "REML",
      sparse = TRUE,
      control = list(
        optimizer = "nlminb",
        rel.tol = 1e-8,
        iter.max = 1000,
        eval.max = 2000
      )
    ),
    error = function(e) NULL
  )
  
  if (!is.null(fit1)) {
    return(list(fit = fit1, V_matrix = V_matrix, data = df, optimizer = "nlminb"))
  }
  
  fit2 <- tryCatch(
    rma.mv(
      yi = yi,
      V = V_matrix,
      mods = ~ delta_T * delta_CO2,
      random = ~ 1 | Study_ID/ES_ID,
      data = df,
      method = "REML",
      sparse = TRUE,
      control = list(
        optimizer = "optim",
        optmethod = "BFGS",
        reltol = 1e-8,
        maxit = 1000
      )
    ),
    error = function(e) NULL
  )
  
  if (!is.null(fit2)) {
    return(list(fit = fit2, V_matrix = V_matrix, data = df, optimizer = "optim_BFGS"))
  }
  
  fit3 <- tryCatch(
    rma.mv(
      yi = yi,
      V = V_matrix,
      mods = ~ delta_T * delta_CO2,
      random = ~ 1 | Study_ID/ES_ID,
      data = df,
      method = "REML",
      sparse = TRUE,
      control = list(
        optimizer = "optim",
        optmethod = "Nelder-Mead",
        reltol = 1e-8,
        maxit = 2000
      )
    ),
    error = function(e) NULL
  )
  
  if (!is.null(fit3)) {
    return(list(fit = fit3, V_matrix = V_matrix, data = df, optimizer = "optim_NelderMead"))
  }
  
  NULL
}

# ============================================================
# Prediction grid
# ============================================================

make_interaction_prediction_grid <- function(fit, df, n_points = 100) {
  co2_vals <- quantile(
    df$delta_CO2,
    probs = c(0.25, 0.50, 0.75),
    na.rm = TRUE
  )
  
  co2_vals <- as.numeric(co2_vals)
  
  co2_tbl <- tibble(
    delta_CO2 = co2_vals,
    CO2_label = paste0("\u0394CO2 = ", round(co2_vals, 0), " ppm")
  )
  
  delta_T_seq <- seq(
    min(df$delta_T, na.rm = TRUE),
    max(df$delta_T, na.rm = TRUE),
    length.out = n_points
  )
  
  pred_grid <- tidyr::crossing(
    delta_T = delta_T_seq,
    co2_tbl
  ) %>%
    mutate(
      interaction_term = delta_T * delta_CO2
    )
  
  pred_mat <- cbind(
    pred_grid$delta_T,
    pred_grid$delta_CO2,
    pred_grid$interaction_term
  )
  
  preds <- predict(fit, newmods = pred_mat)
  
  pred_grid %>%
    mutate(
      pred = as.numeric(preds$pred),
      se = as.numeric(preds$se),
      ci_lb = pred - qnorm(0.975) * se,
      ci_ub = pred + qnorm(0.975) * se
    )
}

# ============================================================
# Summarise model
# ============================================================

summarise_metareg_interaction <- function(df, analysis_id, save_v = FALSE) {
  k <- nrow(df)
  
  if (
    k < 4 ||
    length(unique(df$delta_T)) < 2 ||
    length(unique(df$delta_CO2)) < 2 ||
    dplyr::n_distinct(df$Study_ID) < 2
  ) {
    return(list(
      summary = tibble(
        analysis_id = analysis_id,
        optimizer = NA_character_,
        k = k,
        studies = dplyr::n_distinct(df$Study_ID),
        beta_0 = NA_real_,
        se_0 = NA_real_,
        z_0 = NA_real_,
        p_0 = NA_real_,
        beta_delta_T = NA_real_,
        se_delta_T = NA_real_,
        z_delta_T = NA_real_,
        p_delta_T = NA_real_,
        beta_delta_CO2 = NA_real_,
        se_delta_CO2 = NA_real_,
        z_delta_CO2 = NA_real_,
        p_delta_CO2 = NA_real_,
        beta_interaction = NA_real_,
        se_interaction = NA_real_,
        z_interaction = NA_real_,
        p_interaction = NA_real_,
        QM = NA_real_,
        QM_p = NA_real_,
        tau2_1 = NA_real_,
        tau2_2 = NA_real_
      ),
      fit = NULL,
      pred_grid = tibble(),
      plot_df = tibble()
    ))
  }
  
  fit_obj <- fit_metareg_interaction(df)
  
  if (is.null(fit_obj)) {
    return(list(
      summary = tibble(
        analysis_id = analysis_id,
        optimizer = NA_character_,
        k = k,
        studies = dplyr::n_distinct(df$Study_ID),
        beta_0 = NA_real_,
        se_0 = NA_real_,
        z_0 = NA_real_,
        p_0 = NA_real_,
        beta_delta_T = NA_real_,
        se_delta_T = NA_real_,
        z_delta_T = NA_real_,
        p_delta_T = NA_real_,
        beta_delta_CO2 = NA_real_,
        se_delta_CO2 = NA_real_,
        z_delta_CO2 = NA_real_,
        p_delta_CO2 = NA_real_,
        beta_interaction = NA_real_,
        se_interaction = NA_real_,
        z_interaction = NA_real_,
        p_interaction = NA_real_,
        QM = NA_real_,
        QM_p = NA_real_,
        tau2_1 = NA_real_,
        tau2_2 = NA_real_
      ),
      fit = NULL,
      pred_grid = tibble(),
      plot_df = tibble()
    ))
  }
  
  fit <- fit_obj$fit
  
  coef_tab <- as.data.frame(coef(summary(fit)))
  coef_tab$term <- rownames(coef_tab)
  rownames(coef_tab) <- NULL
  
  intercept_row   <- coef_tab %>% filter(term == "intrcpt")
  deltaT_row      <- coef_tab %>% filter(term == "delta_T")
  deltaCO2_row    <- coef_tab %>% filter(term == "delta_CO2")
  interaction_row <- coef_tab %>% filter(term == "delta_T:delta_CO2")
  
  if (save_v) {
    write.csv(
      fit_obj$V_matrix,
      file.path(table_dir, paste0("V_matrix_", analysis_id, ".csv")),
      row.names = TRUE
    )
    
    saveRDS(
      fit_obj$V_matrix,
      file.path(rds_dir, paste0("V_matrix_", analysis_id, ".rds"))
    )
  }
  
  pred_grid <- make_interaction_prediction_grid(fit, df, n_points = 100)
  
  tau2_vals <- fit$sigma2
  tau2_1 <- if (length(tau2_vals) >= 1) tau2_vals[1] else NA_real_
  tau2_2 <- if (length(tau2_vals) >= 2) tau2_vals[2] else NA_real_
  
  summary_tbl <- tibble(
    analysis_id = analysis_id,
    optimizer = fit_obj$optimizer,
    k = fit$k,
    studies = dplyr::n_distinct(df$Study_ID),
    beta_0 = if (nrow(intercept_row) == 1) intercept_row$estimate else NA_real_,
    se_0 = if (nrow(intercept_row) == 1) intercept_row$se else NA_real_,
    z_0 = if (nrow(intercept_row) == 1) intercept_row$zval else NA_real_,
    p_0 = if (nrow(intercept_row) == 1) intercept_row$pval else NA_real_,
    beta_delta_T = if (nrow(deltaT_row) == 1) deltaT_row$estimate else NA_real_,
    se_delta_T = if (nrow(deltaT_row) == 1) deltaT_row$se else NA_real_,
    z_delta_T = if (nrow(deltaT_row) == 1) deltaT_row$zval else NA_real_,
    p_delta_T = if (nrow(deltaT_row) == 1) deltaT_row$pval else NA_real_,
    beta_delta_CO2 = if (nrow(deltaCO2_row) == 1) deltaCO2_row$estimate else NA_real_,
    se_delta_CO2 = if (nrow(deltaCO2_row) == 1) deltaCO2_row$se else NA_real_,
    z_delta_CO2 = if (nrow(deltaCO2_row) == 1) deltaCO2_row$zval else NA_real_,
    p_delta_CO2 = if (nrow(deltaCO2_row) == 1) deltaCO2_row$pval else NA_real_,
    beta_interaction = if (nrow(interaction_row) == 1) interaction_row$estimate else NA_real_,
    se_interaction = if (nrow(interaction_row) == 1) interaction_row$se else NA_real_,
    z_interaction = if (nrow(interaction_row) == 1) interaction_row$zval else NA_real_,
    p_interaction = if (nrow(interaction_row) == 1) interaction_row$pval else NA_real_,
    QM = if (!is.null(fit$QM)) as.numeric(fit$QM) else NA_real_,
    QM_p = if (!is.null(fit$QMp)) as.numeric(fit$QMp) else NA_real_,
    tau2_1 = tau2_1,
    tau2_2 = tau2_2
  )
  
  list(
    summary = summary_tbl,
    fit = fit,
    pred_grid = pred_grid,
    plot_df = df
  )
}

# ============================================================
# Plot interaction
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
      paste0("interaction p = ", formatC(summary_row$p_interaction, format = "f", digits = 4))
    }
  } else {
    "interaction p = NA"
  }
  
  b_int_lab <- if (is.finite(summary_row$beta_interaction)) {
    paste0("interaction = ", formatC(summary_row$beta_interaction, format = "f", digits = 6))
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
# Safe runner
# ============================================================

run_one_analysis <- function(df, analysis_id, panel_title, save_v = FALSE) {
  cat("\n-----------------------------\n")
  cat("Running:", analysis_id, "\n")
  cat("Rows:", nrow(df), "\n")
  cat("Studies:", dplyr::n_distinct(df$Study_ID), "\n")
  cat("-----------------------------\n")
  
  res <- summarise_metareg_interaction(
    df = df,
    analysis_id = analysis_id,
    save_v = save_v
  )
  
  write_csv(
    res$summary,
    file.path(table_dir, paste0(analysis_id, "_summary.csv"))
  )
  
  write_csv(
    res$pred_grid,
    file.path(table_dir, paste0(analysis_id, "_prediction_grid.csv"))
  )
  
  write_csv(
    res$plot_df,
    file.path(table_dir, paste0(analysis_id, "_plot_data.csv"))
  )
  
  saveRDS(
    res,
    file.path(rds_dir, paste0(analysis_id, "_result.rds"))
  )
  
  if (!is.null(res$fit)) {
    saveRDS(
      res$fit,
      file.path(rds_dir, paste0(analysis_id, "_model.rds"))
    )
    
    sink(file.path(table_dir, paste0(analysis_id, "_model_summary.txt")))
    print(summary(res$fit))
    sink()
  }
  
  p <- plot_metareg_interaction(
    plot_df = res$plot_df,
    pred_grid = res$pred_grid,
    summary_row = res$summary,
    title_text = panel_title
  )
  
  if (!is.null(p)) {
    ggsave(
      filename = file.path(figure_dir, paste0(analysis_id, ".png")),
      plot = p,
      width = 8.8,
      height = 6.8,
      dpi = 600,
      bg = "white"
    )
    
    saveRDS(
      p,
      file.path(rds_dir, paste0(analysis_id, "_plot.rds"))
    )
  }
  
  list(
    id = analysis_id,
    title = panel_title,
    result = res,
    plot = p
  )
}

# ============================================================
# 1. Overall interaction model
# ============================================================

overall_run <- run_one_analysis(
  df = dat_interaction,
  analysis_id = "overall_TpH_deltaT_x_deltaCO2",
  panel_title = "(a) Overall TpH interaction",
  save_v = TRUE
)

# ============================================================
# 2. Broad subgroup runs
# ============================================================

subgroup_levels <- dat_interaction %>%
  filter(!is.na(Broad_subgroup)) %>%
  distinct(Broad_subgroup) %>%
  pull(Broad_subgroup) %>%
  as.character()

subgroup_runs <- lapply(seq_along(subgroup_levels), function(i) {
  g <- subgroup_levels[i]
  
  df_g <- dat_interaction %>%
    filter(Broad_subgroup == g)
  
  run_one_analysis(
    df = df_g,
    analysis_id = paste0("broad_subgroup_", str_replace_all(g, "[^A-Za-z0-9]+", "_")),
    panel_title = paste0("(", letters[i + 1], ") ", g),
    save_v = FALSE
  )
})

subgroup_summary_tbl <- bind_rows(lapply(subgroup_runs, function(x) {
  x$result$summary %>%
    mutate(Broad_subgroup = x$title)
}))

write_csv(
  subgroup_summary_tbl,
  file.path(table_dir, "broad_subgroup_summaries_combined.csv")
)

saveRDS(
  subgroup_runs,
  file.path(rds_dir, "broad_subgroup_runs.rds")
)

# ============================================================
# 3. Response category runs
# ============================================================

response_levels <- dat_interaction %>%
  filter(!is.na(Response_category), Response_category != "") %>%
  distinct(Response_category) %>%
  arrange(Response_category) %>%
  pull(Response_category)

response_runs <- lapply(seq_along(response_levels), function(i) {
  rc <- response_levels[i]
  
  df_rc <- dat_interaction %>%
    filter(Response_category == rc)
  
  run_one_analysis(
    df = df_rc,
    analysis_id = paste0("response_category_", str_replace_all(rc, "[^A-Za-z0-9]+", "_")),
    panel_title = paste0("Response: ", rc),
    save_v = FALSE
  )
})

response_summary_tbl <- bind_rows(lapply(response_runs, function(x) {
  x$result$summary %>%
    mutate(Response_category_panel = x$title)
}))

write_csv(
  response_summary_tbl,
  file.path(table_dir, "response_category_summaries_combined.csv")
)

saveRDS(
  response_runs,
  file.path(rds_dir, "response_category_runs.rds")
)

# ============================================================
# Combined panel: overall and broad subgroups
# ============================================================

subgroup_plots <- c(
  list(overall_run$plot),
  lapply(subgroup_runs, function(x) x$plot)
)

subgroup_plots <- subgroup_plots[
  !vapply(subgroup_plots, is.null, logical(1))
]

if (length(subgroup_plots) > 0) {
  p_subgroups_combined <- wrap_plots(subgroup_plots, nrow = 1)
  
  ggsave(
    filename = file.path(figure_dir, "combined_overall_and_broad_subgroups.png"),
    plot = p_subgroups_combined,
    width = 8 * length(subgroup_plots),
    height = 6.8,
    dpi = 600,
    bg = "white"
  )
  
  saveRDS(
    p_subgroups_combined,
    file.path(rds_dir, "combined_overall_and_broad_subgroups.rds")
  )
  
  print(p_subgroups_combined)
}

# ============================================================
# Combined panel: response categories
# ============================================================

response_plot_list <- lapply(response_runs, function(x) x$plot)

response_plot_list <- response_plot_list[
  !vapply(response_plot_list, is.null, logical(1))
]

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
  
  print(p_response_combined)
}

# ============================================================
# Master summary table
# ============================================================

master_summary_tbl <- bind_rows(
  overall_run$result$summary %>%
    mutate(
      panel_type = "Overall",
      panel_name = "Overall TpH interaction"
    ),
  bind_rows(lapply(subgroup_runs, function(x) {
    x$result$summary %>%
      mutate(
        panel_type = "Broad_subgroup",
        panel_name = x$title
      )
  })),
  bind_rows(lapply(response_runs, function(x) {
    x$result$summary %>%
      mutate(
        panel_type = "Response_category",
        panel_name = x$title
      )
  }))
)

write_csv(
  master_summary_tbl,
  file.path(table_dir, "all_interaction_model_summaries_combined.csv")
)

saveRDS(
  master_summary_tbl,
  file.path(rds_dir, "all_interaction_model_summaries_combined.rds")
)

# ============================================================
# Console output
# ============================================================

cat("\n==== Overall interaction summary ====\n")
print(overall_run$result$summary)

if (nrow(subgroup_summary_tbl) > 0) {
  cat("\n==== Broad subgroup interaction summaries ====\n")
  print(subgroup_summary_tbl)
}

if (nrow(response_summary_tbl) > 0) {
  cat("\n==== Response category interaction summaries ====\n")
  print(response_summary_tbl)
}

cat("\n==== Saved outputs ====\n")
cat("Tables :", table_dir, "\n")
cat("RDS    :", rds_dir, "\n")
cat("Figures:", figure_dir, "\n")

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

cat("\nDone.\n")