# ============================================================
# 05_Regression_analysis_T.R
# Purpose:
#   Use the shared cleaned zooplankton dataset from
#   Scripts/00_Dependencies.R, calculate lnRR and sampling variance,
#   fit T-only multilevel meta-regressions with delta_T,
#   and create bubble plots.
#
# Input:
#   dat1 object from Scripts/00_Dependencies.R
#
# Outputs:
#   Outputs/Tables/05_Regression_analysis_T/
#   Outputs/RDS/05_Regression_analysis_T/
#   Outputs/Figures/05_Regression_analysis_T/
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
    legend.position = "none",
    plot.margin = margin(8, 8, 8, 8),
    plot.background = element_blank(),
    axis.ticks = element_line(linewidth = 0.3, colour = "black"),
    axis.ticks.length = unit(-0.15, "lines"),
    axis.title.y = element_text(size = 13, margin = margin(0, 6, 0, 0), face = "bold"),
    axis.title.x = element_text(size = 13, margin = margin(6, 0, 0, 0), face = "bold"),
    axis.text.y = element_text(size = 11, margin = margin(0, 5, 0, 0), colour = "#000000"),
    axis.text.x = element_text(size = 11, margin = margin(5, 0, 0, 0), colour = "#000000")
  )
}

# ============================================================
# Script name and output folders
# ============================================================

script_name <- "05_Regression_analysis_T"

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
    Broad_subgroup = factor(Broad_subgroup, levels = c("Meroplankton", "Holoplankton")),
    Response_category = factor(Response_category),
    delta_T = suppressWarnings(as.numeric(delta_T))
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
# Keep only T treatment for delta_T meta-regression
# ============================================================

dat_deltaT <- dat1_ml %>%
  filter(Treatment == "T") %>%
  filter(is.finite(delta_T), delta_T > 0, delta_T <= 10) %>%
  filter(is.finite(yi), is.finite(vi), vi > 0)

# ============================================================
# Trim extreme vi values for stability
# ============================================================

q_low  <- quantile(dat_deltaT$vi, 0.01, na.rm = TRUE)
q_high <- quantile(dat_deltaT$vi, 0.99, na.rm = TRUE)

dat_deltaT <- dat_deltaT %>%
  filter(vi >= q_low, vi <= q_high)

saveRDS(
  dat_deltaT,
  file.path(rds_dir, "dat_deltaT_T_only_trimmed.rds")
)

write_csv(
  dat_deltaT,
  file.path(table_dir, "dat_deltaT_T_only_trimmed.csv")
)

cat("\n==== T-only delta_T meta-regression data ====\n")
cat("Rows available:", nrow(dat_deltaT), "\n")
cat("Unique studies:", dplyr::n_distinct(dat_deltaT$Study_ID), "\n")
cat("Unique effect sizes:", dplyr::n_distinct(dat_deltaT$ES_ID), "\n")
cat("Unique Common_id:", dplyr::n_distinct(dat_deltaT$Common_id), "\n")
cat("delta_T range:", min(dat_deltaT$delta_T, na.rm = TRUE), "to", max(dat_deltaT$delta_T, na.rm = TRUE), "\n")
cat("vi range:", min(dat_deltaT$vi, na.rm = TRUE), "to", max(dat_deltaT$vi, na.rm = TRUE), "\n")

cat("\nSummary of vi:\n")
print(summary(dat_deltaT$vi))

cat("\nSummary of delta_T:\n")
print(summary(dat_deltaT$delta_T))

cat("\nStudy size summary:\n")
print(summary(table(dat_deltaT$Study_ID)))

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
# Fit T-only delta_T meta-regression
# ============================================================

fit_metareg_deltaT <- function(df) {
  df <- df %>%
    arrange(Study_ID, ES_ID) %>%
    mutate(
      Study_ID = factor(Study_ID),
      ES_ID = factor(ES_ID),
      Common_id = factor(Common_id)
    )
  
  if (nrow(df) < 2) return(NULL)
  if (length(unique(df$delta_T)) < 2) return(NULL)
  
  V_matrix <- make_V_matrix(df)
  
  fit1 <- tryCatch(
    rma.mv(
      yi = yi,
      V = V_matrix,
      mods = ~ delta_T,
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
      mods = ~ delta_T,
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
      mods = ~ delta_T,
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
# Plot object helpers
# ============================================================

make_bubble_size <- function(df) {
  if (nrow(df) == 0) {
    return(df %>% mutate(BubbleSize = numeric()))
  }
  
  Weight <- 1 / sqrt(df$vi)
  Breaks <- quantile(Weight, c(0.2, 0.4, 0.6, 0.8, 1), na.rm = TRUE)
  
  BubbleSize <- rep(1, length(Weight))
  BubbleSize[Weight < Breaks[1]] <- 1
  
  for (i in seq(2, length(Breaks))) {
    idx <- which(Weight < Breaks[i] & Weight > Breaks[i - 1])
    BubbleSize[idx] <- 0.5 + (i - 1) * 0.3
  }
  
  BubbleSize <- BubbleSize * 2.5
  
  df %>%
    mutate(
      Weight = Weight,
      BubbleSize = BubbleSize
    )
}

trim_plot_yi <- function(df, lower_q = 0.02, upper_q = 0.98) {
  if (nrow(df) < 3) return(df)
  
  q_low_yi  <- quantile(df$yi, lower_q, na.rm = TRUE)
  q_high_yi <- quantile(df$yi, upper_q, na.rm = TRUE)
  
  df %>%
    filter(yi > q_low_yi, yi < q_high_yi)
}

make_prediction_objects <- function(fit, df, xvar = "delta_T", n_points = 100) {
  x_seq <- seq(
    min(df[[xvar]], na.rm = TRUE),
    max(df[[xvar]], na.rm = TRUE),
    length.out = n_points
  )
  
  fit_df <- tibble(
    X = x_seq,
    Y = NA_real_,
    Upper = NA_real_,
    Lower = NA_real_
  )
  
  for (i in seq_along(x_seq)) {
    pred_i <- predict(fit, newmods = matrix(x_seq[i], nrow = 1))
    
    fit_df$Y[i]     <- as.numeric(pred_i$pred)
    fit_df$Upper[i] <- as.numeric(pred_i$pred + qnorm(0.975) * pred_i$se)
    fit_df$Lower[i] <- as.numeric(pred_i$pred - qnorm(0.975) * pred_i$se)
  }
  
  ci_poly <- tibble(
    x = c(fit_df$X, rev(fit_df$X), fit_df$X[1]),
    y = c(fit_df$Upper, rev(fit_df$Lower), fit_df$Upper[1])
  )
  
  list(fit_df = fit_df, ci_poly = ci_poly)
}

# ============================================================
# Summarise T-only delta_T meta-regression
# ============================================================

summarise_metareg_deltaT <- function(df, analysis_id, save_v = FALSE) {
  k <- nrow(df)
  
  if (k < 2 || length(unique(df$delta_T)) < 2) {
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
        QM = NA_real_,
        QM_p = NA_real_,
        tau2_1 = NA_real_,
        tau2_2 = NA_real_
      ),
      fit = NULL,
      fit_df = tibble(),
      ci_poly = tibble(),
      plot_df = tibble()
    ))
  }
  
  fit_obj <- fit_metareg_deltaT(df)
  
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
        QM = NA_real_,
        QM_p = NA_real_,
        tau2_1 = NA_real_,
        tau2_2 = NA_real_
      ),
      fit = NULL,
      fit_df = tibble(),
      ci_poly = tibble(),
      plot_df = tibble()
    ))
  }
  
  fit <- fit_obj$fit
  
  coef_tab <- as.data.frame(coef(summary(fit)))
  coef_tab$term <- rownames(coef_tab)
  rownames(coef_tab) <- NULL
  
  intercept_row <- coef_tab %>% filter(term == "intrcpt")
  slope_row     <- coef_tab %>% filter(term == "delta_T")
  
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
  
  plot_df <- df %>%
    make_bubble_size() %>%
    trim_plot_yi()
  
  pred_objs <- make_prediction_objects(
    fit,
    df,
    xvar = "delta_T",
    n_points = 100
  )
  
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
    beta_delta_T = if (nrow(slope_row) == 1) slope_row$estimate else NA_real_,
    se_delta_T = if (nrow(slope_row) == 1) slope_row$se else NA_real_,
    z_delta_T = if (nrow(slope_row) == 1) slope_row$zval else NA_real_,
    p_delta_T = if (nrow(slope_row) == 1) slope_row$pval else NA_real_,
    QM = if (!is.null(fit$QM)) as.numeric(fit$QM) else NA_real_,
    QM_p = if (!is.null(fit$QMp)) as.numeric(fit$QMp) else NA_real_,
    tau2_1 = tau2_1,
    tau2_2 = tau2_2
  )
  
  list(
    summary = summary_tbl,
    fit = fit,
    fit_df = pred_objs$fit_df,
    ci_poly = pred_objs$ci_poly,
    plot_df = plot_df
  )
}

# ============================================================
# Plot function
# ============================================================

plot_metareg_deltaT <- function(plot_df,
                                fit_df,
                                ci_poly,
                                summary_row,
                                title_text = "",
                                point_fill = "#2b6a99",
                                line_type = 1,
                                y_limits = NULL,
                                y_breaks = waiver(),
                                x_lab = expression(Delta * "T (" * degree * "C)"),
                                y_lab = "Mean effect size (lnRR)") {
  if (nrow(plot_df) == 0 || nrow(fit_df) == 0 || nrow(ci_poly) == 0) {
    return(NULL)
  }
  
  slope_lab <- if (is.finite(summary_row$beta_delta_T)) {
    paste0("slope = ", formatC(summary_row$beta_delta_T, format = "f", digits = 4))
  } else {
    "slope = NA"
  }
  
  p_lab <- if (is.finite(summary_row$p_delta_T)) {
    if (summary_row$p_delta_T < 0.0001) {
      "p < 0.0001"
    } else {
      paste0("p = ", formatC(summary_row$p_delta_T, format = "f", digits = 4))
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
      size = 5,
      hjust = 0
    ) +
    annotate(
      "text",
      x = x_max - 0.22 * x_rng,
      y = y_limits[2] - 0.19 * y_rng,
      label = slope_lab,
      size = 5,
      hjust = 0
    ) +
    annotate(
      "text",
      x = x_min + 0.01 * x_rng,
      y = y_limits[2] - 0.10 * y_rng,
      label = title_text,
      size = 5,
      hjust = 0
    ) +
    scale_y_continuous(
      limits = y_limits,
      breaks = y_breaks,
      expand = c(0, 0)
    ) +
    labs(
      x = x_lab,
      y = y_lab
    ) +
    theme_meta()
}

# ============================================================
# Run overall model
# ============================================================

overall_res <- summarise_metareg_deltaT(
  df = dat_deltaT,
  analysis_id = "overall_T_deltaT",
  save_v = TRUE
)

overall_deltaT_summary <- overall_res$summary
overall_deltaT_fit     <- overall_res$fit

write_csv(
  overall_deltaT_summary,
  file.path(table_dir, "overall_deltaT_meta_regression_summary.csv")
)

saveRDS(
  overall_res,
  file.path(rds_dir, "overall_deltaT_meta_regression_result.rds")
)

if (!is.null(overall_deltaT_fit)) {
  saveRDS(
    overall_deltaT_fit,
    file.path(rds_dir, "overall_deltaT_meta_regression_model.rds")
  )
  
  sink(file.path(table_dir, "overall_deltaT_model_summary.txt"))
  print(summary(overall_deltaT_fit))
  sink()
}

write_csv(
  overall_res$fit_df,
  file.path(table_dir, "overall_deltaT_prediction_grid.csv")
)

write_csv(
  overall_res$ci_poly,
  file.path(table_dir, "overall_deltaT_ci_polygon.csv")
)

write_csv(
  overall_res$plot_df,
  file.path(table_dir, "overall_deltaT_plot_data.csv")
)

# ============================================================
# Make overall plot
# ============================================================

p_overall <- plot_metareg_deltaT(
  plot_df = overall_res$plot_df,
  fit_df = overall_res$fit_df,
  ci_poly = overall_res$ci_poly,
  summary_row = overall_deltaT_summary,
  title_text = "(a) Overall",
  point_fill = "#2b6a99",
  line_type = 1,
  x_lab = expression(Delta * "T (" * degree * "C)"),
  y_lab = "Mean effect size (lnRR)"
)

if (!is.null(p_overall)) {
  ggsave(
    filename = file.path(figure_dir, "overall_deltaT_meta_regression_plot.png"),
    plot = p_overall,
    width = 8.5,
    height = 6.5,
    dpi = 600,
    bg = "white"
  )
  
  saveRDS(
    p_overall,
    file.path(rds_dir, "overall_deltaT_meta_regression_plot.rds")
  )
}

# ============================================================
# Optional subgroup models and plots
# ============================================================

subgroup_levels <- dat_deltaT %>%
  filter(!is.na(Broad_subgroup)) %>%
  distinct(Broad_subgroup) %>%
  pull(Broad_subgroup) %>%
  as.character()

subgroup_results <- lapply(subgroup_levels, function(g) {
  df_g <- dat_deltaT %>%
    filter(Broad_subgroup == g)
  
  res_g <- summarise_metareg_deltaT(
    df = df_g,
    analysis_id = paste0("subgroup_", g),
    save_v = FALSE
  )
  
  list(
    subgroup = g,
    data = df_g,
    result = res_g
  )
})

subgroup_summary_tbl <- bind_rows(lapply(subgroup_results, function(x) {
  x$result$summary %>%
    mutate(Broad_subgroup = x$subgroup)
}))

write_csv(
  subgroup_summary_tbl,
  file.path(table_dir, "subgroup_deltaT_meta_regression_summary.csv")
)

saveRDS(
  subgroup_results,
  file.path(rds_dir, "subgroup_deltaT_meta_regression_results.rds")
)

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
      height = 6.5,
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
    
    if (!is.null(res_g$fit)) {
      saveRDS(
        res_g$fit,
        file.path(
          rds_dir,
          paste0("subgroup_", safe_name, "_deltaT_meta_regression_model.rds")
        )
      )
      
      sink(file.path(table_dir, paste0("subgroup_", safe_name, "_model_summary.txt")))
      print(summary(res_g$fit))
      sink()
    }
    
    write_csv(
      res_g$fit_df,
      file.path(table_dir, paste0("subgroup_", safe_name, "_prediction_grid.csv"))
    )
    
    write_csv(
      res_g$ci_poly,
      file.path(table_dir, paste0("subgroup_", safe_name, "_ci_polygon.csv"))
    )
    
    write_csv(
      res_g$plot_df,
      file.path(table_dir, paste0("subgroup_", safe_name, "_plot_data.csv"))
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
    height = 6.5,
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

cat("\n==== Overall delta_T meta-regression summary ====\n")
print(overall_deltaT_summary)

if (!is.null(overall_deltaT_fit)) {
  cat("\n==== Overall model output ====\n")
  print(summary(overall_deltaT_fit))
}

if (nrow(subgroup_summary_tbl) > 0) {
  cat("\n==== Subgroup delta_T meta-regression summary ====\n")
  print(subgroup_summary_tbl)
}

cat("\n==== Quick interpretation ====\n")

if (is.finite(overall_deltaT_summary$p_delta_T)) {
  if (
    overall_deltaT_summary$p_delta_T < 0.05 &&
    overall_deltaT_summary$beta_delta_T < 0
  ) {
    cat("Larger temperature increases are associated with more negative lnRR.\n")
  } else if (
    overall_deltaT_summary$p_delta_T < 0.05 &&
    overall_deltaT_summary$beta_delta_T > 0
  ) {
    cat("Larger temperature increases are associated with more positive lnRR.\n")
  } else {
    cat("No slope detected for delta_T.\n")
  }
} else {
  cat("Model did not return an interpretable slope estimate.\n")
}

cat("\n==== Saved outputs ====\n")
cat("Tables :", table_dir, "\n")
cat("RDS    :", rds_dir, "\n")
cat("Figures:", figure_dir, "\n")

if (!is.null(p_overall)) {
  cat(
    "Overall figure:",
    file.path(figure_dir, "overall_deltaT_meta_regression_plot.png"),
    "\n"
  )
}

if (length(plots_one_row) > 0) {
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