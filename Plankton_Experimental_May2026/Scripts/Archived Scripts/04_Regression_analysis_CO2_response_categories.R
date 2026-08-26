# ============================================================
# 04_Regression_analysis_CO2_response_categories.R
# Purpose:
#   Use the shared cleaned zooplankton dataset from
#   Scripts/00_Dependencies.R, calculate lnRR and sampling variance,
#   fit pH-only multilevel meta-regressions with delta_CO2
#   by response category, and create response-category panels.
#
# Input:
#   dat1 object from Scripts/00_Dependencies.R
#
# Outputs:
#   Outputs/Tables/04_Regression_analysis_CO2_response_categories/
#   Outputs/RDS/04_Regression_analysis_CO2_response_categories/
#   Outputs/Figures/04_Regression_analysis_CO2_response_categories/
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
    axis.text.x = element_text(size = 11, margin = margin(5, 0, 0, 0), colour = "#000000"),
    strip.text = element_text(face = "bold", size = 12)
  )
}

# ============================================================
# Script name and output folders
# ============================================================

script_name <- "04_Regression_analysis_CO2_response_categories"

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

if (!"delta_CO2" %in% names(dat1)) {
  stop("delta_CO2 not found in dat1. Check that Scripts/00_Dependencies.R creates delta_CO2.")
}

if (!"Response_category" %in% names(dat1)) {
  stop("Response_category not found in dat1. Check that Scripts/00_Dependencies.R creates Response_category.")
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
    Response_category = str_squish(as.character(Response_category)),
    Response_category = factor(Response_category),
    delta_CO2 = suppressWarnings(as.numeric(delta_CO2))
  ) %>%
  filter(
    !is.na(Response_category),
    Response_category != ""
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
# Keep only pH treatment for delta_CO2 meta-regression
# ============================================================

dat_deltaCO2 <- dat1_ml %>%
  filter(Treatment == "pH") %>%
  filter(is.finite(delta_CO2), delta_CO2 > 0) %>%
  filter(is.finite(yi), is.finite(vi), vi > 0)

# ============================================================
# Trim extreme vi values for stability
# ============================================================

q_low  <- quantile(dat_deltaCO2$vi, 0.01, na.rm = TRUE)
q_high <- quantile(dat_deltaCO2$vi, 0.99, na.rm = TRUE)

dat_deltaCO2 <- dat_deltaCO2 %>%
  filter(vi >= q_low, vi <= q_high)

saveRDS(
  dat_deltaCO2,
  file.path(rds_dir, "dat_deltaCO2_pH_only_trimmed.rds")
)

write_csv(
  dat_deltaCO2,
  file.path(table_dir, "dat_deltaCO2_pH_only_trimmed.csv")
)

cat("\n==== pH-only delta_CO2 response-category data ====\n")
cat("Rows available:", nrow(dat_deltaCO2), "\n")
cat("Unique studies:", dplyr::n_distinct(dat_deltaCO2$Study_ID), "\n")
cat("Unique effect sizes:", dplyr::n_distinct(dat_deltaCO2$ES_ID), "\n")
cat("Unique Common_id:", dplyr::n_distinct(dat_deltaCO2$Common_id), "\n")
cat("Response categories:", paste(sort(unique(as.character(dat_deltaCO2$Response_category))), collapse = ", "), "\n")
cat("delta_CO2 range:", min(dat_deltaCO2$delta_CO2, na.rm = TRUE), "to", max(dat_deltaCO2$delta_CO2, na.rm = TRUE), "\n")
cat("vi range:", min(dat_deltaCO2$vi, na.rm = TRUE), "to", max(dat_deltaCO2$vi, na.rm = TRUE), "\n")

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
# Fit pH-only delta_CO2 meta-regression
# ============================================================

fit_metareg_deltaCO2 <- function(df) {
  df <- df %>%
    arrange(Study_ID, ES_ID) %>%
    mutate(
      Study_ID = factor(Study_ID),
      ES_ID = factor(ES_ID),
      Common_id = factor(Common_id)
    )
  
  if (nrow(df) < 2) return(NULL)
  if (length(unique(df$delta_CO2)) < 2) return(NULL)
  if (dplyr::n_distinct(df$Study_ID) < 2) return(NULL)
  
  V_matrix <- make_V_matrix(df)
  
  fit1 <- tryCatch(
    rma.mv(
      yi = yi,
      V = V_matrix,
      mods = ~ delta_CO2,
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
      mods = ~ delta_CO2,
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
      mods = ~ delta_CO2,
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
# Plot-data helpers
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

make_prediction_objects <- function(fit, df, xvar = "delta_CO2", n_points = 100) {
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
# Summary/model wrapper
# ============================================================

summarise_metareg_deltaCO2 <- function(df, analysis_id, save_v = FALSE) {
  k <- nrow(df)
  
  if (k < 2 || length(unique(df$delta_CO2)) < 2 || dplyr::n_distinct(df$Study_ID) < 2) {
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
        beta_delta_CO2 = NA_real_,
        se_delta_CO2 = NA_real_,
        z_delta_CO2 = NA_real_,
        p_delta_CO2 = NA_real_,
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
  
  fit_obj <- fit_metareg_deltaCO2(df)
  
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
        beta_delta_CO2 = NA_real_,
        se_delta_CO2 = NA_real_,
        z_delta_CO2 = NA_real_,
        p_delta_CO2 = NA_real_,
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
  slope_row     <- coef_tab %>% filter(term == "delta_CO2")
  
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
    xvar = "delta_CO2",
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
    beta_delta_CO2 = if (nrow(slope_row) == 1) slope_row$estimate else NA_real_,
    se_delta_CO2 = if (nrow(slope_row) == 1) slope_row$se else NA_real_,
    z_delta_CO2 = if (nrow(slope_row) == 1) slope_row$zval else NA_real_,
    p_delta_CO2 = if (nrow(slope_row) == 1) slope_row$pval else NA_real_,
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

plot_metareg_deltaCO2 <- function(plot_df,
                                  fit_df,
                                  ci_poly,
                                  summary_row,
                                  title_text = "",
                                  point_fill = "#8A2D00",
                                  line_type = 1,
                                  x_lab = expression(Delta * "CO"[2] * " (ppm)"),
                                  y_lab = "Effect size (lnRR)") {
  if (nrow(plot_df) == 0 || nrow(fit_df) == 0 || nrow(ci_poly) == 0) {
    return(NULL)
  }
  
  slope_lab <- if (is.finite(summary_row$beta_delta_CO2)) {
    paste0("slope = ", formatC(summary_row$beta_delta_CO2, format = "f", digits = 4))
  } else {
    "slope = NA"
  }
  
  p_lab <- if (is.finite(summary_row$p_delta_CO2)) {
    if (summary_row$p_delta_CO2 < 0.0001) {
      "p < 0.0001"
    } else {
      paste0("p = ", formatC(summary_row$p_delta_CO2, format = "f", digits = 4))
    }
  } else {
    "p = NA"
  }
  
  x_min <- min(fit_df$X, na.rm = TRUE)
  x_max <- max(fit_df$X, na.rm = TRUE)
  x_rng <- x_max - x_min
  
  y_min <- min(c(plot_df$yi, ci_poly$y), na.rm = TRUE)
  y_max <- max(c(plot_df$yi, ci_poly$y), na.rm = TRUE)
  y_pad <- 0.08 * (y_max - y_min)
  y_limits <- c(y_min - y_pad, y_max + y_pad)
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
      size = 3.8,
      hjust = 0
    ) +
    annotate(
      "text",
      x = x_max - 0.28 * x_rng,
      y = y_limits[2] - 0.20 * y_rng,
      label = slope_lab,
      size = 3.8,
      hjust = 0
    ) +
    labs(
      title = title_text,
      x = x_lab,
      y = y_lab
    ) +
    scale_y_continuous(expand = c(0.02, 0.02)) +
    theme_meta() +
    theme(
      plot.title = element_text(face = "bold", size = 13, hjust = 0.5)
    )
}

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
# Run per response category
# ============================================================

response_levels <- dat_deltaCO2 %>%
  filter(!is.na(Response_category), Response_category != "") %>%
  distinct(Response_category) %>%
  pull(Response_category) %>%
  as.character()

response_results <- lapply(response_levels, function(rc) {
  df_rc <- dat_deltaCO2 %>%
    filter(Response_category == rc)
  
  res_rc <- summarise_metareg_deltaCO2(
    df = df_rc,
    analysis_id = paste0("response_", make.names(rc)),
    save_v = FALSE
  )
  
  list(
    response_category = rc,
    data = df_rc,
    result = res_rc
  )
})

response_summary_tbl <- bind_rows(lapply(response_results, function(x) {
  x$result$summary %>%
    mutate(Response_category = x$response_category)
}))

write_csv(
  response_summary_tbl,
  file.path(table_dir, "response_category_deltaCO2_meta_regression_summary.csv")
)

saveRDS(
  response_results,
  file.path(rds_dir, "response_category_deltaCO2_meta_regression_results.rds")
)

# ============================================================
# Save per-category model outputs and plots
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
    x_lab = expression(Delta * "CO"[2]),
    y_lab = "Effect size (lnRR)"
  )
  
  response_plots[[rc]] <- p_rc
  
  safe_name <- str_replace_all(rc, "[^A-Za-z0-9]+", "_")
  
  write_csv(
    res_rc$summary,
    file.path(table_dir, paste0("response_", safe_name, "_summary.csv"))
  )
  
  write_csv(
    res_rc$fit_df,
    file.path(table_dir, paste0("response_", safe_name, "_prediction_grid.csv"))
  )
  
  write_csv(
    res_rc$ci_poly,
    file.path(table_dir, paste0("response_", safe_name, "_ci_polygon.csv"))
  )
  
  write_csv(
    res_rc$plot_df,
    file.path(table_dir, paste0("response_", safe_name, "_plot_data.csv"))
  )
  
  saveRDS(
    res_rc,
    file.path(rds_dir, paste0("response_", safe_name, "_result.rds"))
  )
  
  if (!is.null(res_rc$fit)) {
    saveRDS(
      res_rc$fit,
      file.path(rds_dir, paste0("response_", safe_name, "_model.rds"))
    )
    
    sink(file.path(table_dir, paste0("response_", safe_name, "_model_summary.txt")))
    print(summary(res_rc$fit))
    sink()
  }
  
  if (!is.null(p_rc)) {
    ggsave(
      filename = file.path(figure_dir, paste0("response_", safe_name, "_deltaCO2_plot.png")),
      plot = p_rc,
      width = 6,
      height = 4.3,
      dpi = 600,
      bg = "white"
    )
    
    saveRDS(
      p_rc,
      file.path(rds_dir, paste0("response_", safe_name, "_deltaCO2_plot.rds"))
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
      title = "pH-only meta-regression with delta_CO2 by response category",
      theme = theme(
        plot.title = element_text(
          family = fam,
          face = "bold",
          size = 18,
          hjust = 0.5
        )
      )
    )
  
  ggsave(
    filename = file.path(figure_dir, "response_category_deltaCO2_meta_regression_panels.png"),
    plot = p_response_panel,
    width = 12,
    height = 9,
    dpi = 600,
    bg = "white"
  )
  
  saveRDS(
    p_response_panel,
    file.path(rds_dir, "response_category_deltaCO2_meta_regression_panels.rds")
  )
  
  print(p_response_panel)
}

# ============================================================
# Console output
# ============================================================

cat("\n==== Response category delta_CO2 meta-regression summary ====\n")
print(response_summary_tbl)

cat("\n==== Saved outputs ====\n")
cat("Tables :", table_dir, "\n")
cat("RDS    :", rds_dir, "\n")
cat("Figures:", figure_dir, "\n")

if (exists("p_response_panel")) {
  cat(
    "Combined figure:",
    file.path(figure_dir, "response_category_deltaCO2_meta_regression_panels.png"),
    "\n"
  )
}