# ============================================================
# 02A_Model_lnrr_calcifier.R
# Purpose:
#   Use the shared cleaned dataset from 00B_Dependencies.R,
#   calculate lnRR and sampling variance directly, and fit single
#   interaction models for Calcifier vs Non-calcifier groups.
#
#   This is the calcifier equivalent of 11A:
#     Model 1: all final model-ready data
#     Model 2: full-factorial Common_id only
#
#   Interaction model:
#     yi ~ 0 + pH_stressor + T_stressor + pH_T_interaction
#
#   Predicted treatment effects:
#     pH      = pH_stressor
#     T       = T_stressor
#     T + pH  = pH_stressor + T_stressor + pH_T_interaction
#
# Important:
#   Common_id is used through the V matrix to account for shared
#   controls/non-independence.
#
#   This script does not depend on outputs from 02A_Model_lnrr.R.
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
# Script name and output folders
# ============================================================

script_name <- "02A_Model_lnrr_calcifier"

table_dir <- file.path("Outputs", "Tables", script_name)
rds_dir   <- file.path("Outputs", "RDS", script_name)

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# Load cleaned data and calculate lnRR and sampling variance
# ============================================================

source(file.path("Scripts", "00B_Dependencies.R"))

if (!exists("dat1")) {
  stop("dat1 not found. Check Scripts/00B_Dependencies.R.")
}

required_cols <- c(
  "Biota_C", "Biota_T",
  "Biota_C_sd", "Biota_T_sd",
  "N_C", "N_T",
  "Study_ID", "ES_ID", "Common_id",
  "Treatment", "Response_category", "Metric_Type", "Calcifier"
)

missing_cols <- setdiff(required_cols, names(dat1))

if (length(missing_cols) > 0) {
  stop(
    "Missing required column(s): ",
    paste(missing_cols, collapse = ", ")
  )
}

dat1_ml <- dat1 %>%
  mutate(
    Biota_C    = suppressWarnings(as.numeric(Biota_C)),
    Biota_T    = suppressWarnings(as.numeric(Biota_T)),
    Biota_C_sd = suppressWarnings(as.numeric(Biota_C_sd)),
    Biota_T_sd = suppressWarnings(as.numeric(Biota_T_sd)),
    N_C        = suppressWarnings(as.numeric(N_C)),
    N_T        = suppressWarnings(as.numeric(N_T)),
    Metric_Type = stringr::str_squish(as.character(Metric_Type)),
    Metric_Type = case_when(
      stringr::str_detect(Metric_Type, regex("^pos", ignore_case = TRUE)) ~ "Positive",
      stringr::str_detect(Metric_Type, regex("^neg", ignore_case = TRUE)) ~ "Negative",
      stringr::str_detect(Metric_Type, regex("^amb", ignore_case = TRUE)) ~ "Ambiguous",
      TRUE ~ NA_character_
    ),
    direction_multiplier = case_when(
      Metric_Type == "Positive" ~  1,
      Metric_Type == "Negative" ~ -1,
      Metric_Type == "Ambiguous" ~ 1,
      TRUE ~ NA_real_
    )
  ) %>%
  filter(
    is.finite(Biota_C),
    is.finite(Biota_T),
    Biota_C > 0,
    Biota_T > 0,
    is.finite(Biota_C_sd),
    is.finite(Biota_T_sd),
    Biota_C_sd >= 0,
    Biota_T_sd >= 0,
    is.finite(N_C),
    is.finite(N_T),
    N_C > 0,
    N_T > 0
  ) %>%
  mutate(
    yi_raw = log(Biota_T / Biota_C),
    yi = yi_raw * direction_multiplier,
    vi =
      (Biota_T_sd^2) / (N_T * Biota_T^2) +
      (Biota_C_sd^2) / (N_C * Biota_C^2)
  ) %>%
  filter(
    is.finite(yi),
    is.finite(vi),
    vi > 0
  )

cat("\n==== Created model-ready data from 00B_Dependencies.R ====\n")
cat("Rows in dat1_ml:", nrow(dat1_ml), "\n")
cat("Unique ES_ID:", dplyr::n_distinct(dat1_ml$ES_ID), "\n")
cat("Unique Study_ID:", dplyr::n_distinct(dat1_ml$Study_ID), "\n")
cat("Unique Common_id:", dplyr::n_distinct(dat1_ml$Common_id), "\n\n")

cat("Metric direction handling:\n")
cat("Positive  : raw lnRR retained\n")
cat("Negative  : raw lnRR multiplied by -1\n")
cat("Ambiguous : raw lnRR retained\n\n")
print(table(dat1_ml$Metric_Type, useNA = "ifany"))

cat("Raw Calcifier values:\n")
print(table(dat1_ml$Calcifier, useNA = "ifany"))

# ============================================================
# Labels and levels
# ============================================================

response_levels <- c(
  "Growth",
  "Development",
  "Physiology",
  "Reproduction",
  "Survival"
)

category_levels <- c(
  "Growth",
  "Development",
  "Physiology",
  "Reproduction",
  "Survival",
  "Overall"
)

treatment_levels <- c("pH", "T", "TpH")

treatment_label_levels <- c(
  "pH",
  "T",
  "T + pH"
)

calcifier_levels <- c(
  "Calcifier",
  "Non-calcifier"
)

model_scopes <- c(
  "all_data",
  "full_factorial_only"
)

plot_levels_calcifier_interaction <- list(
  response_levels = response_levels,
  category_levels = category_levels,
  treatment_levels = treatment_levels,
  treatment_label_levels = treatment_label_levels,
  calcifier_levels = calcifier_levels,
  model_scopes = model_scopes
)

saveRDS(
  plot_levels_calcifier_interaction,
  file.path(rds_dir, "plot_levels_calcifier_interaction.rds")
)

# ============================================================
# Helper functions
# ============================================================

safe_text <- function(x) {
  x %>%
    as.character() %>%
    stringr::str_replace_all("[^A-Za-z0-9]+", "_") %>%
    stringr::str_replace_all("^_+|_+$", "")
}

calc_p_value <- function(estimate, se) {
  ifelse(
    is.finite(estimate) & is.finite(se) & se > 0,
    2 * stats::pnorm(abs(estimate / se), lower.tail = FALSE),
    NA_real_
  )
}

# ============================================================
# Clean Calcifier, Treatment, and add stressor coding
# ============================================================

dat1_calc <- dat1_ml %>%
  mutate(
    Calcifier_clean = stringr::str_squish(as.character(Calcifier)),
    Calcifier_clean_lower = stringr::str_to_lower(Calcifier_clean),
    Calcifier_clean = case_when(
      Calcifier_clean_lower %in% c("calcifier", "calcifiers") ~ "Calcifier",
      Calcifier_clean_lower %in% c(
        "non-calcifier",
        "non calcifier",
        "noncalcifier",
        "non-calcifiers",
        "non calcifiers",
        "noncalcifiers"
      ) ~ "Non-calcifier",
      TRUE ~ NA_character_
    ),
    Calcifier = factor(Calcifier_clean, levels = calcifier_levels),
    
    Treatment = factor(as.character(Treatment), levels = treatment_levels),
    Treatment_label = dplyr::recode(
      as.character(Treatment),
      "pH"  = "pH",
      "T"   = "T",
      "TpH" = "T + pH"
    ),
    Treatment_label = factor(Treatment_label, levels = treatment_label_levels),
    
    Response_category = as.character(Response_category),
    Response_category = case_when(
      Response_category %in% c("Behaviour", "Physiology") ~ "Physiology",
      TRUE ~ Response_category
    ),
    Response_category = factor(Response_category, levels = response_levels),
    
    Treatment_chr = as.character(Treatment),
    pH_stressor = ifelse(Treatment_chr %in% c("pH", "TpH"), 1, 0),
    T_stressor  = ifelse(Treatment_chr %in% c("T", "TpH"), 1, 0),
    pH_T_interaction = pH_stressor * T_stressor
  ) %>%
  filter(
    !is.na(Calcifier),
    !is.na(Treatment),
    !is.na(Treatment_label),
    is.finite(yi),
    is.finite(vi),
    vi > 0
  )

cat("\n==== Cleaned Calcifier counts ====\n")
print(table(dat1_calc$Calcifier, useNA = "ifany"))

cat("\n==== Treatment counts by Calcifier ====\n")
print(table(dat1_calc$Calcifier, dat1_calc$Treatment_label, useNA = "ifany"))

write.csv(
  dat1_calc,
  file.path(table_dir, "dat1_calcifier_interaction_model_ready.csv"),
  row.names = FALSE
)

saveRDS(
  dat1_calc,
  file.path(rds_dir, "dat1_calcifier_interaction_model_ready.rds")
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

make_common_id_design_flags <- function(df, category_name, calcifier_name) {
  if (nrow(df) == 0) {
    return(tibble(
      Category = category_name,
      Response_category = ifelse(category_name == "Overall", NA_character_, category_name),
      Calcifier = calcifier_name,
      Common_id = NA_character_,
      has_pH = NA,
      has_T = NA,
      has_TpH = NA,
      full_factorial_common_id = NA
    ))
  }
  
  df %>%
    filter(!is.na(Common_id), !is.na(Treatment)) %>%
    distinct(Common_id, Treatment) %>%
    group_by(Common_id) %>%
    summarise(
      has_pH  = any(as.character(Treatment) == "pH"),
      has_T   = any(as.character(Treatment) == "T"),
      has_TpH = any(as.character(Treatment) == "TpH"),
      full_factorial_common_id = has_pH & has_T & has_TpH,
      .groups = "drop"
    ) %>%
    mutate(
      Category = category_name,
      Response_category = ifelse(category_name == "Overall", NA_character_, category_name),
      Calcifier = calcifier_name
    ) %>%
    select(
      Category,
      Response_category,
      Calcifier,
      Common_id,
      has_pH,
      has_T,
      has_TpH,
      full_factorial_common_id
    )
}

# ============================================================
# VCV matrix builder
# ============================================================

make_V_matrix <- function(df) {
  split_list <- split(df, df$Common_id)
  
  block_list <- lapply(split_list, function(x) {
    n <- nrow(x)
    
    shared_cov <- suppressWarnings(
      (x$Biota_C_sd[1]^2) / (x$N_C[1] * x$Biota_C[1]^2)
    )
    
    if (!is.finite(shared_cov)) {
      shared_cov <- 0
    }
    
    direction <- x$direction_multiplier
    
    v <- outer(direction, direction) * shared_cov
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
# Empty output templates
# ============================================================

make_empty_predicted_effects <- function(
    df,
    analysis_scope,
    analysis_id,
    category_name,
    calcifier_name,
    model_status,
    model_message = NA_character_
) {
  tibble(
    analysis_scope = analysis_scope,
    analysis_id = analysis_id,
    Category = category_name,
    Response_category = ifelse(category_name == "Overall", NA_character_, category_name),
    Calcifier = calcifier_name,
    Treatment = treatment_levels,
    Treatment_label = treatment_label_levels,
    mean = NA_real_,
    se = NA_real_,
    ci_low = NA_real_,
    ci_high = NA_real_,
    zval = NA_real_,
    pval = NA_real_,
    k = nrow(df),
    n_studies = dplyr::n_distinct(df$Study_ID),
    n_common_id = dplyr::n_distinct(df$Common_id),
    n_full_factorial_common_id = length(get_full_factorial_common_ids(df)),
    n_pH = sum(as.character(df$Treatment) == "pH", na.rm = TRUE),
    n_T = sum(as.character(df$Treatment) == "T", na.rm = TRUE),
    n_TpH = sum(as.character(df$Treatment) == "TpH", na.rm = TRUE),
    model_status = model_status,
    model_message = model_message
  )
}

make_empty_coefficients <- function(
    df,
    analysis_scope,
    analysis_id,
    category_name,
    calcifier_name,
    model_status,
    model_message = NA_character_
) {
  tibble(
    analysis_scope = analysis_scope,
    analysis_id = analysis_id,
    Category = category_name,
    Response_category = ifelse(category_name == "Overall", NA_character_, category_name),
    Calcifier = calcifier_name,
    term = c("pH_stressor", "T_stressor", "pH_T_interaction"),
    estimate = NA_real_,
    se = NA_real_,
    ci_low = NA_real_,
    ci_high = NA_real_,
    zval = NA_real_,
    pval = NA_real_,
    k = nrow(df),
    n_studies = dplyr::n_distinct(df$Study_ID),
    n_common_id = dplyr::n_distinct(df$Common_id),
    n_full_factorial_common_id = length(get_full_factorial_common_ids(df)),
    n_pH = sum(as.character(df$Treatment) == "pH", na.rm = TRUE),
    n_T = sum(as.character(df$Treatment) == "T", na.rm = TRUE),
    n_TpH = sum(as.character(df$Treatment) == "TpH", na.rm = TRUE),
    model_status = model_status,
    model_message = model_message
  )
}

make_model_stats <- function(
    df,
    analysis_scope,
    analysis_id,
    category_name,
    calcifier_name,
    model_status,
    model_message = NA_character_,
    fit = NULL
) {
  if (is.null(fit)) {
    return(tibble(
      analysis_scope = analysis_scope,
      analysis_id = analysis_id,
      Category = category_name,
      Response_category = ifelse(category_name == "Overall", NA_character_, category_name),
      Calcifier = calcifier_name,
      model = "yi ~ 0 + pH_stressor + T_stressor + pH_T_interaction + random(~1 | Study_ID/ES_ID)",
      model_status = model_status,
      model_message = model_message,
      k = nrow(df),
      n_studies = dplyr::n_distinct(df$Study_ID),
      n_common_id = dplyr::n_distinct(df$Common_id),
      n_full_factorial_common_id = length(get_full_factorial_common_ids(df)),
      n_pH = sum(as.character(df$Treatment) == "pH", na.rm = TRUE),
      n_T = sum(as.character(df$Treatment) == "T", na.rm = TRUE),
      n_TpH = sum(as.character(df$Treatment) == "TpH", na.rm = TRUE),
      rank_X = NA_integer_,
      sigma2 = NA_character_,
      QE = NA_real_,
      QE_df = NA_real_,
      QEp = NA_real_,
      logLik = NA_real_,
      AIC = NA_real_,
      BIC = NA_real_
    ))
  }
  
  X <- model.matrix(
    ~ 0 + pH_stressor + T_stressor + pH_T_interaction,
    data = df
  )
  
  sigma2_val <- if (!is.null(fit$sigma2) && length(fit$sigma2) > 0) {
    paste(round(fit$sigma2, 6), collapse = "; ")
  } else {
    NA_character_
  }
  
  tibble(
    analysis_scope = analysis_scope,
    analysis_id = analysis_id,
    Category = category_name,
    Response_category = ifelse(category_name == "Overall", NA_character_, category_name),
    Calcifier = calcifier_name,
    model = "yi ~ 0 + pH_stressor + T_stressor + pH_T_interaction + random(~1 | Study_ID/ES_ID)",
    model_status = model_status,
    model_message = model_message,
    k = as.integer(fit$k),
    n_studies = dplyr::n_distinct(df$Study_ID),
    n_common_id = dplyr::n_distinct(df$Common_id),
    n_full_factorial_common_id = length(get_full_factorial_common_ids(df)),
    n_pH = sum(as.character(df$Treatment) == "pH", na.rm = TRUE),
    n_T = sum(as.character(df$Treatment) == "T", na.rm = TRUE),
    n_TpH = sum(as.character(df$Treatment) == "TpH", na.rm = TRUE),
    rank_X = qr(X)$rank,
    sigma2 = sigma2_val,
    QE = fit$QE,
    QE_df = fit$k - fit$p,
    QEp = fit$QEp,
    logLik = as.numeric(logLik(fit)),
    AIC = AIC(fit),
    BIC = BIC(fit)
  )
}

# ============================================================
# Main interaction model-fitting function
# ============================================================

fit_interaction_model <- function(
    df,
    analysis_scope,
    analysis_id,
    category_name,
    calcifier_name,
    save_v = FALSE
) {
  if (nrow(df) == 0) {
    return(list(
      predicted_effects = make_empty_predicted_effects(
        df,
        analysis_scope,
        analysis_id,
        category_name,
        calcifier_name,
        model_status = "no_data"
      ),
      coefficients = make_empty_coefficients(
        df,
        analysis_scope,
        analysis_id,
        category_name,
        calcifier_name,
        model_status = "no_data"
      ),
      model_stats = make_model_stats(
        df,
        analysis_scope,
        analysis_id,
        category_name,
        calcifier_name,
        model_status = "no_data"
      )
    ))
  }
  
  df_fit <- df %>%
    arrange(Common_id, Study_ID, ES_ID) %>%
    mutate(
      Study_ID = factor(Study_ID),
      ES_ID = factor(ES_ID),
      Common_id = factor(Common_id),
      Calcifier = factor(Calcifier, levels = calcifier_levels),
      Treatment = factor(Treatment, levels = treatment_levels),
      Treatment_label = factor(Treatment_label, levels = treatment_label_levels),
      Response_category = factor(Response_category, levels = response_levels),
      pH_stressor = as.numeric(pH_stressor),
      T_stressor = as.numeric(T_stressor),
      pH_T_interaction = as.numeric(pH_T_interaction)
    )
  
  X <- model.matrix(
    ~ 0 + pH_stressor + T_stressor + pH_T_interaction,
    data = df_fit
  )
  
  rank_X <- qr(X)$rank
  treatments_present <- sort(unique(as.character(df_fit$Treatment)))
  
  if (nrow(df_fit) < 3 || rank_X < 3) {
    msg <- paste0(
      "Insufficient treatment structure for full pH + T + pH:T model. ",
      "Treatments present: ",
      paste(treatments_present, collapse = ", "),
      "; rank_X = ",
      rank_X,
      "."
    )
    
    return(list(
      predicted_effects = make_empty_predicted_effects(
        df_fit,
        analysis_scope,
        analysis_id,
        category_name,
        calcifier_name,
        model_status = "insufficient_treatment_structure",
        model_message = msg
      ),
      coefficients = make_empty_coefficients(
        df_fit,
        analysis_scope,
        analysis_id,
        category_name,
        calcifier_name,
        model_status = "insufficient_treatment_structure",
        model_message = msg
      ),
      model_stats = make_model_stats(
        df_fit,
        analysis_scope,
        analysis_id,
        category_name,
        calcifier_name,
        model_status = "insufficient_treatment_structure",
        model_message = msg
      )
    ))
  }
  
  V_matrix <- make_V_matrix(df_fit)
  
  if (save_v) {
    write.csv(
      V_matrix,
      file.path(table_dir, paste0("V_matrix_", analysis_id, "_", analysis_scope, ".csv")),
      row.names = TRUE
    )
    
    saveRDS(
      V_matrix,
      file.path(rds_dir, paste0("V_matrix_", analysis_id, "_", analysis_scope, ".rds"))
    )
  }
  
  fit_try <- tryCatch(
    metafor::rma.mv(
      yi = yi,
      V = V_matrix,
      mods = ~ 0 + pH_stressor + T_stressor + pH_T_interaction,
      random = ~ 1 | Study_ID/ES_ID,
      data = df_fit
    ),
    error = function(e) e
  )
  
  if (inherits(fit_try, "error")) {
    msg <- paste0("Model error: ", fit_try$message)
    
    return(list(
      predicted_effects = make_empty_predicted_effects(
        df_fit,
        analysis_scope,
        analysis_id,
        category_name,
        calcifier_name,
        model_status = "model_error",
        model_message = msg
      ),
      coefficients = make_empty_coefficients(
        df_fit,
        analysis_scope,
        analysis_id,
        category_name,
        calcifier_name,
        model_status = "model_error",
        model_message = msg
      ),
      model_stats = make_model_stats(
        df_fit,
        analysis_scope,
        analysis_id,
        category_name,
        calcifier_name,
        model_status = "model_error",
        model_message = msg
      )
    ))
  }
  
  fit <- fit_try
  
  # ------------------------------------------------------------
  # Model-predicted treatment effects
  # ------------------------------------------------------------
  
  newmods <- rbind(
    c(1, 0, 0),  # pH
    c(0, 1, 0),  # T
    c(1, 1, 1)   # predicted T + pH
  )
  
  colnames(newmods) <- c(
    "pH_stressor",
    "T_stressor",
    "pH_T_interaction"
  )
  
  pred_try <- tryCatch(
    predict(fit, newmods = newmods),
    error = function(e) e
  )
  
  if (inherits(pred_try, "error")) {
    msg <- paste0("Prediction error: ", pred_try$message)
    
    return(list(
      predicted_effects = make_empty_predicted_effects(
        df_fit,
        analysis_scope,
        analysis_id,
        category_name,
        calcifier_name,
        model_status = "prediction_error",
        model_message = msg
      ),
      coefficients = make_empty_coefficients(
        df_fit,
        analysis_scope,
        analysis_id,
        category_name,
        calcifier_name,
        model_status = "prediction_error",
        model_message = msg
      ),
      model_stats = make_model_stats(
        df_fit,
        analysis_scope,
        analysis_id,
        category_name,
        calcifier_name,
        model_status = "prediction_error",
        model_message = msg,
        fit = fit
      )
    ))
  }
  
  pred <- pred_try
  
  predicted_effects <- tibble(
    analysis_scope = analysis_scope,
    analysis_id = analysis_id,
    Category = category_name,
    Response_category = ifelse(category_name == "Overall", NA_character_, category_name),
    Calcifier = calcifier_name,
    Treatment = treatment_levels,
    Treatment_label = treatment_label_levels,
    mean = as.numeric(pred$pred),
    se = as.numeric(pred$se),
    ci_low = as.numeric(pred$ci.lb),
    ci_high = as.numeric(pred$ci.ub)
  ) %>%
    mutate(
      zval = ifelse(is.finite(se) & se > 0, mean / se, NA_real_),
      pval = calc_p_value(mean, se),
      k = as.integer(fit$k),
      n_studies = dplyr::n_distinct(df_fit$Study_ID),
      n_common_id = dplyr::n_distinct(df_fit$Common_id),
      n_full_factorial_common_id = length(get_full_factorial_common_ids(df_fit)),
      n_pH = sum(as.character(df_fit$Treatment) == "pH", na.rm = TRUE),
      n_T = sum(as.character(df_fit$Treatment) == "T", na.rm = TRUE),
      n_TpH = sum(as.character(df_fit$Treatment) == "TpH", na.rm = TRUE),
      model_status = "fitted",
      model_message = NA_character_
    )
  
  # ------------------------------------------------------------
  # Model coefficients
  # ------------------------------------------------------------
  
  coef_df <- as.data.frame(coef(summary(fit))) %>%
    tibble::rownames_to_column("term") %>%
    as_tibble() %>%
    transmute(
      analysis_scope = analysis_scope,
      analysis_id = analysis_id,
      Category = category_name,
      Response_category = ifelse(category_name == "Overall", NA_character_, category_name),
      Calcifier = calcifier_name,
      term = term,
      estimate = estimate,
      se = se,
      ci_low = ci.lb,
      ci_high = ci.ub,
      zval = zval,
      pval = pval,
      k = as.integer(fit$k),
      n_studies = dplyr::n_distinct(df_fit$Study_ID),
      n_common_id = dplyr::n_distinct(df_fit$Common_id),
      n_full_factorial_common_id = length(get_full_factorial_common_ids(df_fit)),
      n_pH = sum(as.character(df_fit$Treatment) == "pH", na.rm = TRUE),
      n_T = sum(as.character(df_fit$Treatment) == "T", na.rm = TRUE),
      n_TpH = sum(as.character(df_fit$Treatment) == "TpH", na.rm = TRUE),
      model_status = "fitted",
      model_message = NA_character_
    )
  
  model_stats <- make_model_stats(
    df_fit,
    analysis_scope,
    analysis_id,
    category_name,
    calcifier_name,
    model_status = "fitted",
    model_message = NA_character_,
    fit = fit
  )
  
  list(
    predicted_effects = predicted_effects,
    coefficients = coef_df,
    model_stats = model_stats
  )
}

# ============================================================
# Analysis grid
# ============================================================

analysis_grid <- tidyr::expand_grid(
  Category = category_levels,
  Calcifier = calcifier_levels
)

get_analysis_subset <- function(dat, category_name, calcifier_name) {
  x <- dat %>%
    filter(as.character(Calcifier) == calcifier_name)
  
  if (category_name != "Overall") {
    x <- x %>%
      filter(as.character(Response_category) == category_name)
  }
  
  x
}

# ============================================================
# Common_id design summary
# ============================================================

common_id_design_flags <- purrr::map_dfr(
  seq_len(nrow(analysis_grid)),
  function(i) {
    category_name  <- analysis_grid$Category[i]
    calcifier_name <- analysis_grid$Calcifier[i]
    
    df_base <- get_analysis_subset(
      dat1_calc,
      category_name,
      calcifier_name
    )
    
    make_common_id_design_flags(
      df_base,
      category_name,
      calcifier_name
    )
  }
)

common_id_design_summary <- common_id_design_flags %>%
  filter(!is.na(Common_id)) %>%
  group_by(Category, Response_category, Calcifier) %>%
  summarise(
    n_common_id_total = n_distinct(Common_id),
    n_full_factorial_common_id = sum(full_factorial_common_id, na.rm = TRUE),
    n_incomplete_common_id = n_common_id_total - n_full_factorial_common_id,
    .groups = "drop"
  )

write.csv(
  common_id_design_flags,
  file.path(table_dir, "common_id_design_flags_calcifier.csv"),
  row.names = FALSE
)

write.csv(
  common_id_design_summary,
  file.path(table_dir, "common_id_design_summary_calcifier.csv"),
  row.names = FALSE
)

saveRDS(
  common_id_design_flags,
  file.path(rds_dir, "common_id_design_flags_calcifier.rds")
)

saveRDS(
  common_id_design_summary,
  file.path(rds_dir, "common_id_design_summary_calcifier.rds")
)

cat("\n==== Common_id design summary ====\n")
print(common_id_design_summary)

# ============================================================
# Fit interaction models
# ============================================================

all_fit_results <- list()
counter <- 1

for (scope_name in model_scopes) {
  cat("\n===============================\n")
  cat("Running model scope:", scope_name, "\n")
  cat("===============================\n")
  
  for (i in seq_len(nrow(analysis_grid))) {
    category_name  <- analysis_grid$Category[i]
    calcifier_name <- analysis_grid$Calcifier[i]
    
    cat("\n--- Category:", category_name, "| Calcifier:", calcifier_name, "---\n")
    
    df_base <- get_analysis_subset(
      dat1_calc,
      category_name,
      calcifier_name
    )
    
    if (scope_name == "full_factorial_only") {
      full_ids <- get_full_factorial_common_ids(df_base)
      
      df_scope <- df_base %>%
        filter(as.character(Common_id) %in% full_ids)
    } else {
      df_scope <- df_base
    }
    
    cat("Rows used:", nrow(df_scope), "\n")
    cat("Studies:", dplyr::n_distinct(df_scope$Study_ID), "\n")
    cat("Common_id:", dplyr::n_distinct(df_scope$Common_id), "\n")
    cat("Treatment counts:\n")
    print(table(df_scope$Treatment_label, useNA = "ifany"))
    
    analysis_id <- paste0(
      safe_text(category_name),
      "_",
      safe_text(calcifier_name)
    )
    
    all_fit_results[[counter]] <- fit_interaction_model(
      df = df_scope,
      analysis_scope = scope_name,
      analysis_id = analysis_id,
      category_name = category_name,
      calcifier_name = calcifier_name,
      save_v = FALSE
    )
    
    counter <- counter + 1
  }
}

# ============================================================
# Combine model outputs
# ============================================================

interaction_predicted_effects_all_scopes <- purrr::map_dfr(
  all_fit_results,
  "predicted_effects"
) %>%
  mutate(
    analysis_scope = factor(analysis_scope, levels = model_scopes),
    Category = factor(Category, levels = category_levels),
    Response_category = factor(Response_category, levels = response_levels),
    Calcifier = factor(Calcifier, levels = calcifier_levels),
    Treatment = factor(Treatment, levels = treatment_levels),
    Treatment_label = factor(Treatment_label, levels = treatment_label_levels)
  )

interaction_coefficients_all_scopes <- purrr::map_dfr(
  all_fit_results,
  "coefficients"
) %>%
  mutate(
    analysis_scope = factor(analysis_scope, levels = model_scopes),
    Category = factor(Category, levels = category_levels),
    Response_category = factor(Response_category, levels = response_levels),
    Calcifier = factor(Calcifier, levels = calcifier_levels)
  )

interaction_model_stats_all_scopes <- purrr::map_dfr(
  all_fit_results,
  "model_stats"
) %>%
  mutate(
    analysis_scope = factor(analysis_scope, levels = model_scopes),
    Category = factor(Category, levels = category_levels),
    Response_category = factor(Response_category, levels = response_levels),
    Calcifier = factor(Calcifier, levels = calcifier_levels)
  )

# ============================================================
# Save outputs by model scope
# ============================================================

for (scope_name in model_scopes) {
  pred_scope <- interaction_predicted_effects_all_scopes %>%
    filter(as.character(analysis_scope) == scope_name)
  
  coef_scope <- interaction_coefficients_all_scopes %>%
    filter(as.character(analysis_scope) == scope_name)
  
  stats_scope <- interaction_model_stats_all_scopes %>%
    filter(as.character(analysis_scope) == scope_name)
  
  write.csv(
    pred_scope,
    file.path(table_dir, paste0("interaction_predicted_effects_", scope_name, ".csv")),
    row.names = FALSE
  )
  
  write.csv(
    coef_scope,
    file.path(table_dir, paste0("interaction_coefficients_", scope_name, ".csv")),
    row.names = FALSE
  )
  
  write.csv(
    stats_scope,
    file.path(table_dir, paste0("interaction_model_stats_", scope_name, ".csv")),
    row.names = FALSE
  )
  
  saveRDS(
    pred_scope,
    file.path(rds_dir, paste0("interaction_predicted_effects_", scope_name, ".rds"))
  )
  
  saveRDS(
    coef_scope,
    file.path(rds_dir, paste0("interaction_coefficients_", scope_name, ".rds"))
  )
  
  saveRDS(
    stats_scope,
    file.path(rds_dir, paste0("interaction_model_stats_", scope_name, ".rds"))
  )
}

# ============================================================
# Save combined outputs
# ============================================================

write.csv(
  interaction_predicted_effects_all_scopes,
  file.path(table_dir, "interaction_predicted_effects_all_scopes.csv"),
  row.names = FALSE
)

write.csv(
  interaction_coefficients_all_scopes,
  file.path(table_dir, "interaction_coefficients_all_scopes.csv"),
  row.names = FALSE
)

write.csv(
  interaction_model_stats_all_scopes,
  file.path(table_dir, "interaction_model_stats_all_scopes.csv"),
  row.names = FALSE
)

saveRDS(
  interaction_predicted_effects_all_scopes,
  file.path(rds_dir, "interaction_predicted_effects_all_scopes.rds")
)

saveRDS(
  interaction_coefficients_all_scopes,
  file.path(rds_dir, "interaction_coefficients_all_scopes.rds")
)

saveRDS(
  interaction_model_stats_all_scopes,
  file.path(rds_dir, "interaction_model_stats_all_scopes.rds")
)

# ============================================================
# Interaction terms only
# ============================================================

interaction_terms_only <- interaction_coefficients_all_scopes %>%
  filter(term == "pH_T_interaction") %>%
  mutate(
    interaction_interpretation = case_when(
      is.na(estimate) ~ NA_character_,
      estimate < 0 ~ "Combined effect is more negative than additive expectation",
      estimate > 0 ~ "Combined effect is more positive than additive expectation",
      estimate == 0 ~ "Combined effect matches additive expectation"
    )
  )

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
# Compact summary table with percent change
# ============================================================

interaction_predicted_effects_summary_table <- interaction_predicted_effects_all_scopes %>%
  mutate(
    Category = as.character(Category),
    Response_category = as.character(Response_category),
    Calcifier = as.character(Calcifier),
    Treatment = as.character(Treatment_label),
    `Mean lnRR` = ifelse(is.finite(mean), sprintf("%.3f", mean), NA_character_),
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
    analysis_scope,
    Category,
    Calcifier,
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
    factor(analysis_scope, levels = model_scopes),
    factor(Category, levels = category_levels),
    factor(Calcifier, levels = calcifier_levels),
    factor(Treatment, levels = treatment_label_levels)
  )

write.csv(
  interaction_predicted_effects_summary_table,
  file.path(table_dir, "interaction_predicted_effects_summary_table_all_scopes.csv"),
  row.names = FALSE
)

saveRDS(
  interaction_predicted_effects_summary_table,
  file.path(rds_dir, "interaction_predicted_effects_summary_table_all_scopes.rds")
)

# ============================================================
# Console prints
# ============================================================

cat("\n==== Interaction predicted effects: all scopes ====\n")
print(interaction_predicted_effects_all_scopes)

cat("\n==== Interaction coefficients: all scopes ====\n")
print(interaction_coefficients_all_scopes)

cat("\n==== Interaction terms only ====\n")
print(interaction_terms_only)

cat("\n==== Model stats: all scopes ====\n")
print(interaction_model_stats_all_scopes)

cat("\n==== Model status counts ====\n")
print(table(
  interaction_model_stats_all_scopes$analysis_scope,
  interaction_model_stats_all_scopes$model_status,
  useNA = "ifany"
))

cat("\n==== Saved model outputs ====\n")
cat("Tables:", table_dir, "\n")
cat("RDS   :", rds_dir, "\n")

cat("\nDone: 02A calcifier interaction models completed.\n")
