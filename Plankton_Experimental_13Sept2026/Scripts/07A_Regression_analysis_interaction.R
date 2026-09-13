# ============================================================
# 07A_Regression_analysis_interaction.R
# Purpose:
#   Use the shared cleaned zooplankton dataset from
#   Scripts/00B_Dependencies.R, calculate lnRR and sampling variance,
#   and fit TpH-only multilevel interaction meta-regressions with
#   delta_T * delta_CO2.
#
# Important:
#   This script runs the slow rma.mv() models.
#   Run this only when the data/model inputs change.
#
# Input:
#   dat1 object from Scripts/00B_Dependencies.R
#
# Outputs:
#   Outputs/Tables/07A_Regression_analysis_interaction/
#   Outputs/RDS/07A_Regression_analysis_interaction/
# ============================================================

suppressPackageStartupMessages(invisible(
  lapply(
    c(
      "dplyr", "readr", "metafor", "stringr",
      "tibble", "Matrix", "matrixcalc",
      "tidyr", "purrr"
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
# Script name and output folders
# ============================================================

script_name <- "07A_Regression_analysis_interaction"

table_dir <- file.path("Outputs", "Tables", script_name)
rds_dir   <- file.path("Outputs", "RDS", script_name)

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)

if (!dir.exists(table_dir)) stop(paste("Table directory could not be created:", table_dir))
if (!dir.exists(rds_dir)) stop(paste("RDS directory could not be created:", rds_dir))

cat("\nWorking directory:\n", getwd(), "\n")
cat("\nTable directory:\n", table_dir, "\n")
cat("RDS directory:\n", rds_dir, "\n")

# ============================================================
# Check shared cleaned dataset
# ============================================================

if (!exists("dat1")) {
  stop("dat1 not found. Check that Scripts/00B_Dependencies.R creates dat1.")
}

if (!"delta_T" %in% names(dat1)) {
  stop("delta_T not found in dat1. Check that Scripts/00B_Dependencies.R creates delta_T.")
}

if (!"delta_CO2" %in% names(dat1)) {
  stop("delta_CO2 not found in dat1. Check that Scripts/00B_Dependencies.R creates delta_CO2.")
}

if (!"Response_category" %in% names(dat1)) {
  stop("Response_category not found in dat1. Check that Scripts/00B_Dependencies.R creates Response_category.")
}

if (!"Broad_subgroup" %in% names(dat1)) {
  stop("Broad_subgroup not found in dat1. Check that Scripts/00B_Dependencies.R creates Broad_subgroup.")
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
# Summarise interaction model
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
  
  list(
    id = analysis_id,
    title = panel_title,
    result = res
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

saveRDS(
  overall_run,
  file.path(rds_dir, "overall_run.rds")
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

saveRDS(
  subgroup_summary_tbl,
  file.path(rds_dir, "broad_subgroup_summaries_combined.rds")
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

saveRDS(
  response_summary_tbl,
  file.path(rds_dir, "response_category_summaries_combined.rds")
)

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
cat("Tables:", table_dir, "\n")
cat("RDS   :", rds_dir, "\n")

cat("\nDone.\n")