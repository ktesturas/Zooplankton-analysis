# ============================================================
# 11A_Model_lnrr_interaction.R
# Purpose:
#   Use the shared cleaned zooplankton dataset from
#   Scripts/00B_Dependencies.R, calculate lnRR and sampling variance,
#   adjust lnRR direction using Metric_Type,
#   and fit single multilevel meta-analytic models with:
#
#       pH main effect + Temperature main effect + pH:Temperature interaction
#
#   Main model:
#     Uses all final filtered model-ready data.
#
#   Sensitivity model:
#     Uses only full factorial Common_id sets, where pH, T, and TpH
#     are all present within the same Common_id.

# ============================================================


# ============================================================
# Load packages
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

script_name <- "11A_Model_lnrr_interaction"

table_dir <- file.path("Outputs", "Tables", script_name)
rds_dir   <- file.path("Outputs", "RDS", script_name)

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)


# ============================================================
# Check shared cleaned dataset
# ============================================================

if (!exists("dat1")) {
  stop("dat1 not found. Check that Scripts/00B_Dependencies.R creates dat1.")
}

if (!"Metric_Type" %in% names(dat1)) {
  stop("Metric_Type column not found in dat1. Check that the column exists in the database and is retained by 00B_Dependencies.R.")
}

cat("\n==== Shared cleaned dataset ====\n")
cat("Rows in dat1:", nrow(dat1), "\n")
cat("Unique ES_ID:", dplyr::n_distinct(dat1$ES_ID), "\n")
cat("Unique Study_ID:", dplyr::n_distinct(dat1$Study_ID), "\n")
cat("Unique Common_id:", dplyr::n_distinct(dat1$Common_id), "\n")
cat("ES_ID unique? ", dplyr::n_distinct(dat1$ES_ID) == nrow(dat1), "\n\n")

cat("Metric_Type values in dat1 before cleaning:\n")
print(table(dat1$Metric_Type, useNA = "ifany"))


# ============================================================
# Clean and check Metric_Type
# ============================================================

dat1 <- dat1 %>%
  mutate(
    Metric_Type_clean = stringr::str_to_lower(
      stringr::str_squish(as.character(Metric_Type))
    ),
    Metric_Type_clean = case_when(
      is.na(Metric_Type_clean) | Metric_Type_clean == "" ~ "ambiguous",
      Metric_Type_clean %in% c("positive", "negative", "ambiguous") ~ Metric_Type_clean,
      TRUE ~ Metric_Type_clean
    )
  )

invalid_metric_types <- dat1 %>%
  filter(!Metric_Type_clean %in% c("positive", "negative", "ambiguous")) %>%
  distinct(Metric_Type, Metric_Type_clean)

if (nrow(invalid_metric_types) > 0) {
  print(invalid_metric_types)
  stop(
    "Invalid Metric_Type values found. Allowed values are: positive, negative, ambiguous. ",
    "Please fix these values in the database or extend the recoding rules."
  )
}

metric_type_summary <- dat1 %>%
  count(Metric_Type_clean, name = "n") %>%
  mutate(
    pct = round(100 * n / sum(n), 2)
  )

cat("\nMetric_Type summary after cleaning:\n")
print(metric_type_summary)

write.csv(
  metric_type_summary,
  file.path(table_dir, "metric_type_summary.csv"),
  row.names = FALSE
)

saveRDS(
  metric_type_summary,
  file.path(rds_dir, "metric_type_summary.rds")
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
    yi_raw = escalc_out$yi,
    vi = escalc_out$vi,
    
    metric_direction_multiplier = case_when(
      Metric_Type_clean == "negative" ~ -1,
      Metric_Type_clean %in% c("positive", "ambiguous") ~ 1,
      TRUE ~ NA_real_
    ),
    
    yi = yi_raw * metric_direction_multiplier,
    
    # Lump Behaviour into Physiology
    Response_category = as.character(Response_category),
    Response_category = dplyr::case_when(
      Response_category %in% c("Behaviour", "Physiology") ~ "Physiology",
      TRUE ~ Response_category
    ),
    Response_category = factor(
      Response_category,
      levels = c(
        "Growth",
        "Development",
        "Physiology",
        "Reproduction",
        "Survival"
      )
    ),
    Treatment = factor(Treatment, levels = c("pH", "T", "TpH")),
    Broad_subgroup = factor(
      Broad_subgroup,
      levels = c("Meroplankton", "Holoplankton")
    )
  ) %>%
  filter(is.finite(yi), is.finite(vi), vi > 0)

metric_type_used_summary <- dat1_ml %>%
  count(Metric_Type_clean, metric_direction_multiplier, name = "n") %>%
  mutate(
    pct = round(100 * n / sum(n), 2)
  )

cat("\nMetric_Type summary used in model-ready data:\n")
print(metric_type_used_summary)

write.csv(
  metric_type_used_summary,
  file.path(table_dir, "metric_type_used_in_model_summary.csv"),
  row.names = FALSE
)

saveRDS(
  metric_type_used_summary,
  file.path(rds_dir, "metric_type_used_in_model_summary.rds")
)

saveRDS(
  dat1_ml,
  file.path(rds_dir, "dat1_with_lnrr_vi.rds")
)

write.csv(
  dat1_ml,
  file.path(table_dir, "dat1_with_lnrr_vi.csv"),
  row.names = FALSE
)


# ============================================================
# Labels
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

plot_levels <- list(
  response_levels = response_levels,
  subgroup_levels = subgroup_levels,
  final_label_levels = final_label_levels,
  treatment_levels = treatment_levels
)

saveRDS(
  plot_levels,
  file.path(rds_dir, "plot_levels_interaction.rds")
)


# ============================================================
# Add stressor variables
# ============================================================
# Coding:
#
#   Treatment   pH_stressor   T_stressor   pH_T_interaction
#   pH          1             0            0
#   T           0             1            0
#   TpH         1             1            1
#
# Model:
#   yi = beta_pH + beta_T + beta_pH_T_interaction
#
# For TpH rows:
#   predicted TpH = beta_pH + beta_T + beta_pH_T_interaction
# ============================================================

add_stressor_variables <- function(df) {
  df %>%
    mutate(
      Treatment_chr = as.character(Treatment),
      pH_stressor = ifelse(Treatment_chr %in% c("pH", "TpH"), 1, 0),
      T_stressor  = ifelse(Treatment_chr %in% c("T", "TpH"), 1, 0),
      pH_T_interaction = pH_stressor * T_stressor
    )
}


# ============================================================
# Full factorial checker
# ============================================================

get_full_factorial_common_ids <- function(df) {
  df %>%
    mutate(Treatment_chr = as.character(Treatment)) %>%
    group_by(Common_id) %>%
    summarise(
      has_pH  = any(Treatment_chr == "pH",  na.rm = TRUE),
      has_T   = any(Treatment_chr == "T",   na.rm = TRUE),
      has_TpH = any(Treatment_chr == "TpH", na.rm = TRUE),
      treatment_list = paste(sort(unique(Treatment_chr)), collapse = ", "),
      n_effect_sizes = n(),
      .groups = "drop"
    ) %>%
    mutate(
      is_fully_factorial = has_pH & has_T & has_TpH
    )
}

factorial_check_common_all <- get_full_factorial_common_ids(dat1_ml)

cat("\nFull factorial check at Common_id level in dat1_ml:\n")
print(table(factorial_check_common_all$is_fully_factorial, useNA = "ifany"))

write.csv(
  factorial_check_common_all,
  file.path(table_dir, "full_factorial_check_common_id_all_data.csv"),
  row.names = FALSE
)

saveRDS(
  factorial_check_common_all,
  file.path(rds_dir, "full_factorial_check_common_id_all_data.rds")
)


# ============================================================
# VCV matrix builder
# ============================================================
# Important:
#   df must be arranged by Common_id before creating V_matrix,
#   because V_matrix is block diagonal by Common_id.
# ============================================================

make_V_matrix <- function(df) {
  df <- df %>%
    arrange(Common_id, Study_ID, ES_ID)
  
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
# Interaction model-fitting function
# ============================================================

fit_interaction_model <- function(df,
                                  analysis_id,
                                  label,
                                  data_scope,
                                  save_v = FALSE) {
  
  label <- as.character(label)
  data_scope <- as.character(data_scope)
  
  k <- nrow(df)
  
  if (k == 0) {
    return(list(
      predictions = tibble(),
      coefficients = tibble(),
      model_stats = tibble(
        analysis_id = analysis_id,
        labels = label,
        data_scope = data_scope,
        model_status = "no_data",
        k = 0L,
        n_studies = 0L,
        n_common_id = 0L
      ),
      fit = NULL
    ))
  }
  
  treatment_present <- df %>%
    mutate(Treatment_chr = as.character(Treatment)) %>%
    summarise(
      has_pH  = any(Treatment_chr == "pH",  na.rm = TRUE),
      has_T   = any(Treatment_chr == "T",   na.rm = TRUE),
      has_TpH = any(Treatment_chr == "TpH", na.rm = TRUE),
      .groups = "drop"
    )
  
  if (!all(treatment_present$has_pH,
           treatment_present$has_T,
           treatment_present$has_TpH)) {
    
    missing_treatments <- c(
      if (!treatment_present$has_pH) "pH" else NULL,
      if (!treatment_present$has_T) "T" else NULL,
      if (!treatment_present$has_TpH) "TpH" else NULL
    )
    
    return(list(
      predictions = tibble(),
      coefficients = tibble(),
      model_stats = tibble(
        analysis_id = analysis_id,
        labels = label,
        data_scope = data_scope,
        model_status = paste0(
          "not_fitted_missing_treatments_",
          paste(missing_treatments, collapse = "_")
        ),
        k = as.integer(k),
        n_studies = dplyr::n_distinct(df$Study_ID),
        n_common_id = dplyr::n_distinct(df$Common_id)
      ),
      fit = NULL
    ))
  }
  
  df_fit <- df %>%
    add_stressor_variables() %>%
    arrange(Common_id, Study_ID, ES_ID) %>%
    mutate(
      Study_ID = factor(Study_ID),
      ES_ID = factor(ES_ID),
      Common_id = factor(Common_id),
      Broad_subgroup = factor(Broad_subgroup),
      Response_category = factor(Response_category)
    )
  
  V_matrix <- make_V_matrix(df_fit)
  
  if (save_v) {
    write.csv(
      V_matrix,
      file.path(table_dir, paste0("V_matrix_", data_scope, "_", analysis_id, ".csv")),
      row.names = TRUE
    )
    
    saveRDS(
      V_matrix,
      file.path(rds_dir, paste0("V_matrix_", data_scope, "_", analysis_id, ".rds"))
    )
  }
  
  fit <- tryCatch(
    metafor::rma.mv(
      yi = yi,
      V = V_matrix,
      mods = ~ 0 + pH_stressor + T_stressor + pH_T_interaction,
      random = ~ 1 | Study_ID/ES_ID,
      data = df_fit
    ),
    error = function(e) e
  )
  
  if (inherits(fit, "error")) {
    return(list(
      predictions = tibble(),
      coefficients = tibble(),
      model_stats = tibble(
        analysis_id = analysis_id,
        labels = label,
        data_scope = data_scope,
        model_status = paste0("model_error: ", fit$message),
        k = as.integer(k),
        n_studies = dplyr::n_distinct(df$Study_ID),
        n_common_id = dplyr::n_distinct(df$Common_id)
      ),
      fit = NULL
    ))
  }
  
  # ----------------------------------------------------------
  # Coefficients:
  #   pH_stressor       = pH main effect
  #   T_stressor        = temperature main effect
  #   pH_T_interaction  = extra combined effect beyond pH + T
  # ----------------------------------------------------------
  
  coef_table <- as.data.frame(coef(summary(fit))) %>%
    tibble::rownames_to_column("term") %>%
    as_tibble() %>%
    transmute(
      analysis_id = analysis_id,
      labels = label,
      data_scope = data_scope,
      term = term,
      estimate = estimate,
      se = se,
      ci_low = ci.lb,
      ci_high = ci.ub,
      zval = zval,
      pval = pval,
      k = as.integer(fit$k),
      n_studies = dplyr::n_distinct(df$Study_ID),
      n_common_id = dplyr::n_distinct(df$Common_id)
    )
  
  # ----------------------------------------------------------
  # Predicted effects for plotting in the same style as before:
  #
  #   pH       = beta_pH
  #   T        = beta_T
  #   TpH      = beta_pH + beta_T + beta_pH_T_interaction
  #
  # This means the green point is model-predicted T + pH,
  # not the separate observed pooled TpH mean.
  # ----------------------------------------------------------
  
  newmods <- rbind(
    c(1, 0, 0),  # pH
    c(0, 1, 0),  # T
    c(1, 1, 1)   # predicted TpH
  )
  
  pred <- predict(fit, newmods = newmods)
  
  prediction_table <- tibble(
    analysis_id = analysis_id,
    labels = label,
    data_scope = data_scope,
    Treatment = factor(c("pH", "T", "TpH"), levels = treatment_levels),
    model_term = c(
      "pH main effect",
      "T main effect",
      "predicted T + pH"
    ),
    mean = as.numeric(pred$pred),
    se = as.numeric(pred$se),
    ci_low = as.numeric(pred$ci.lb),
    ci_high = as.numeric(pred$ci.ub)
  ) %>%
    mutate(
      zval = mean / se,
      pval = 2 * pnorm(abs(zval), lower.tail = FALSE),
      k = as.integer(fit$k),
      n_studies = dplyr::n_distinct(df$Study_ID),
      n_common_id = dplyr::n_distinct(df$Common_id),
      labels = factor(labels, levels = final_label_levels)
    )
  
  sigma2_1 <- ifelse(length(fit$sigma2) >= 1, fit$sigma2[1], NA_real_)
  sigma2_2 <- ifelse(length(fit$sigma2) >= 2, fit$sigma2[2], NA_real_)
  
  model_stats <- tibble(
    analysis_id = analysis_id,
    labels = label,
    data_scope = data_scope,
    model_status = "fitted",
    model = "yi ~ 0 + pH_stressor + T_stressor + pH_T_interaction + random(~ 1 | Study_ID/ES_ID)",
    k = as.integer(fit$k),
    n_studies = dplyr::n_distinct(df$Study_ID),
    n_common_id = dplyr::n_distinct(df$Common_id),
    sigma2_1 = sigma2_1,
    sigma2_2 = sigma2_2,
    QE = fit$QE,
    QE_df = fit$k - fit$p,
    QEp = fit$QEp,
    logLik = as.numeric(logLik(fit)),
    AIC = AIC(fit),
    BIC = BIC(fit)
  )
  
  list(
    predictions = prediction_table,
    coefficients = coef_table,
    model_stats = model_stats,
    fit = fit
  )
}


# ============================================================
# Helper function to fit all plot panels
# ============================================================

fit_all_panels <- function(df, data_scope) {
  
  out <- list()
  
  # Overall
  out[["overall"]] <- fit_interaction_model(
    df = df,
    analysis_id = "overall",
    label = "Overall",
    data_scope = data_scope,
    save_v = FALSE
  )
  
  # Response categories
  for (rc in response_levels) {
    out[[paste0("response_", rc)]] <- fit_interaction_model(
      df = df %>% filter(Response_category == rc),
      analysis_id = paste0("response_", rc),
      label = rc,
      data_scope = data_scope,
      save_v = FALSE
    )
  }
  
  # Broad subgroups
  for (subg in subgroup_levels) {
    out[[paste0("subgroup_", subg)]] <- fit_interaction_model(
      df = df %>% filter(Broad_subgroup == subg),
      analysis_id = paste0("subgroup_", subg),
      label = subg,
      data_scope = data_scope,
      save_v = FALSE
    )
  }
  
  out
}


# ============================================================
# Main Model 1: all final filtered data
# ============================================================

cat("\n============================================================\n")
cat("Fitting Model 1: all-data interaction models\n")
cat("============================================================\n")

model_all_data <- fit_all_panels(
  df = dat1_ml,
  data_scope = "all_data"
)

predictions_all_data <- purrr::map_dfr(model_all_data, "predictions") %>%
  mutate(
    Treatment = factor(Treatment, levels = treatment_levels),
    labels = factor(labels, levels = final_label_levels)
  )

coefficients_all_data <- purrr::map_dfr(model_all_data, "coefficients") %>%
  mutate(
    labels = factor(labels, levels = final_label_levels)
  )

model_stats_all_data <- purrr::map_dfr(model_all_data, "model_stats") %>%
  mutate(
    labels = factor(labels, levels = final_label_levels)
  )

write.csv(
  predictions_all_data,
  file.path(table_dir, "interaction_predicted_effects_all_data.csv"),
  row.names = FALSE
)

write.csv(
  coefficients_all_data,
  file.path(table_dir, "interaction_coefficients_all_data.csv"),
  row.names = FALSE
)

write.csv(
  model_stats_all_data,
  file.path(table_dir, "interaction_model_stats_all_data.csv"),
  row.names = FALSE
)

saveRDS(
  predictions_all_data,
  file.path(rds_dir, "interaction_predicted_effects_all_data.rds")
)

saveRDS(
  coefficients_all_data,
  file.path(rds_dir, "interaction_coefficients_all_data.rds")
)

saveRDS(
  model_stats_all_data,
  file.path(rds_dir, "interaction_model_stats_all_data.rds")
)

saveRDS(
  model_all_data,
  file.path(rds_dir, "interaction_model_objects_all_data.rds")
)


# ============================================================
# Sensitivity Model 2: full factorial Common_id only
# ============================================================
# For each analysis panel, full factorial Common_id sets are defined
# within that analysis subset.
#
# Example:
#   For Growth, a Common_id is retained only if Growth rows contain
#   pH, T, and TpH.
# ============================================================

keep_full_factorial_within_subset <- function(df) {
  factorial_ids <- get_full_factorial_common_ids(df) %>%
    filter(is_fully_factorial) %>%
    select(Common_id)
  
  df %>%
    semi_join(factorial_ids, by = "Common_id")
}


fit_all_panels_full_factorial <- function(df, data_scope) {
  
  out <- list()
  
  # Overall
  out[["overall"]] <- fit_interaction_model(
    df = keep_full_factorial_within_subset(df),
    analysis_id = "overall",
    label = "Overall",
    data_scope = data_scope,
    save_v = FALSE
  )
  
  # Response categories
  for (rc in response_levels) {
    df_sub <- df %>%
      filter(Response_category == rc) %>%
      keep_full_factorial_within_subset()
    
    out[[paste0("response_", rc)]] <- fit_interaction_model(
      df = df_sub,
      analysis_id = paste0("response_", rc),
      label = rc,
      data_scope = data_scope,
      save_v = FALSE
    )
  }
  
  # Broad subgroups
  for (subg in subgroup_levels) {
    df_sub <- df %>%
      filter(Broad_subgroup == subg) %>%
      keep_full_factorial_within_subset()
    
    out[[paste0("subgroup_", subg)]] <- fit_interaction_model(
      df = df_sub,
      analysis_id = paste0("subgroup_", subg),
      label = subg,
      data_scope = data_scope,
      save_v = FALSE
    )
  }
  
  out
}

cat("\n============================================================\n")
cat("Fitting Model 2: full-factorial-only sensitivity models\n")
cat("============================================================\n")

model_full_factorial <- fit_all_panels_full_factorial(
  df = dat1_ml,
  data_scope = "full_factorial_only"
)

predictions_full_factorial <- purrr::map_dfr(model_full_factorial, "predictions") %>%
  mutate(
    Treatment = factor(Treatment, levels = treatment_levels),
    labels = factor(labels, levels = final_label_levels)
  )

coefficients_full_factorial <- purrr::map_dfr(model_full_factorial, "coefficients") %>%
  mutate(
    labels = factor(labels, levels = final_label_levels)
  )

model_stats_full_factorial <- purrr::map_dfr(model_full_factorial, "model_stats") %>%
  mutate(
    labels = factor(labels, levels = final_label_levels)
  )

write.csv(
  predictions_full_factorial,
  file.path(table_dir, "interaction_predicted_effects_full_factorial_only.csv"),
  row.names = FALSE
)

write.csv(
  coefficients_full_factorial,
  file.path(table_dir, "interaction_coefficients_full_factorial_only.csv"),
  row.names = FALSE
)

write.csv(
  model_stats_full_factorial,
  file.path(table_dir, "interaction_model_stats_full_factorial_only.csv"),
  row.names = FALSE
)

saveRDS(
  predictions_full_factorial,
  file.path(rds_dir, "interaction_predicted_effects_full_factorial_only.rds")
)

saveRDS(
  coefficients_full_factorial,
  file.path(rds_dir, "interaction_coefficients_full_factorial_only.rds")
)

saveRDS(
  model_stats_full_factorial,
  file.path(rds_dir, "interaction_model_stats_full_factorial_only.rds")
)

saveRDS(
  model_full_factorial,
  file.path(rds_dir, "interaction_model_objects_full_factorial_only.rds")
)


# ============================================================
# Combined outputs
# ============================================================

interaction_predicted_effects_combined <- bind_rows(
  predictions_all_data,
  predictions_full_factorial
) %>%
  mutate(
    Treatment = factor(Treatment, levels = treatment_levels),
    labels = factor(labels, levels = final_label_levels)
  )

interaction_coefficients_combined <- bind_rows(
  coefficients_all_data,
  coefficients_full_factorial
) %>%
  mutate(
    labels = factor(labels, levels = final_label_levels)
  )

interaction_model_stats_combined <- bind_rows(
  model_stats_all_data,
  model_stats_full_factorial
) %>%
  mutate(
    labels = factor(labels, levels = final_label_levels)
  )

write.csv(
  interaction_predicted_effects_combined,
  file.path(table_dir, "interaction_predicted_effects_combined.csv"),
  row.names = FALSE
)

write.csv(
  interaction_coefficients_combined,
  file.path(table_dir, "interaction_coefficients_combined.csv"),
  row.names = FALSE
)

write.csv(
  interaction_model_stats_combined,
  file.path(table_dir, "interaction_model_stats_combined.csv"),
  row.names = FALSE
)

saveRDS(
  interaction_predicted_effects_combined,
  file.path(rds_dir, "interaction_predicted_effects_combined.rds")
)

saveRDS(
  interaction_coefficients_combined,
  file.path(rds_dir, "interaction_coefficients_combined.rds")
)

saveRDS(
  interaction_model_stats_combined,
  file.path(rds_dir, "interaction_model_stats_combined.rds")
)


# ============================================================
# Pull out interaction terms only
# ============================================================

interaction_terms_only <- interaction_coefficients_combined %>%
  filter(term == "pH_T_interaction") %>%
  arrange(data_scope, labels)

write.csv(
  interaction_terms_only,
  file.path(table_dir, "interaction_terms_only.csv"),
  row.names = FALSE
)

saveRDS(
  interaction_terms_only,
  file.path(rds_dir, "interaction_terms_only.rds")
)


# ============================================================
# Console prints
# ============================================================

cat("\n==== Model 1: all-data predicted effects ====\n")
print(predictions_all_data)

cat("\n==== Model 1: all-data interaction coefficients ====\n")
print(
  coefficients_all_data %>%
    filter(term == "pH_T_interaction")
)

cat("\n==== Model 1: all-data model stats ====\n")
print(model_stats_all_data)

cat("\n==== Model 2: full-factorial-only predicted effects ====\n")
print(predictions_full_factorial)

cat("\n==== Model 2: full-factorial-only interaction coefficients ====\n")
print(
  coefficients_full_factorial %>%
    filter(term == "pH_T_interaction")
)

cat("\n==== Model 2: full-factorial-only model stats ====\n")
print(model_stats_full_factorial)

cat("\n==== Saved interaction model outputs ====\n")
cat("Tables:", table_dir, "\n")
cat("RDS   :", rds_dir, "\n")