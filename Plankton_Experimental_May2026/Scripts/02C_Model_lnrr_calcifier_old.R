# ============================================================
# 02C_Model_lnrr_calcifier.R
# Purpose:
#   Use model-ready lnRR data from 02A_Model_lnrr.R and calculate
#   pooled mean lnRRs for Calcifier vs Non-calcifier groups.
#
#   Outputs are saved for plotting in 02D_Plot_lnrr_calcifier.R.
#
# Important:
#   This script runs rma.mv() models.
#   Run this only when data/model inputs change.
#
# Input:
#   Outputs/RDS/02A_Model_lnrr/dat1_with_lnrr_vi.rds
#
# Outputs:
#   Outputs/Tables/02C_Model_lnrr_calcifier/
#   Outputs/RDS/02C_Model_lnrr_calcifier/
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

script_name <- "02C_Model_lnrr_calcifier"

table_dir <- file.path("Outputs", "Tables", script_name)
rds_dir   <- file.path("Outputs", "RDS", script_name)

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# Load model-ready data from 02A
# ============================================================

input_rds <- file.path(
  "Outputs", "RDS", "02A_Model_lnrr",
  "dat1_with_lnrr_vi.rds"
)

if (!file.exists(input_rds)) {
  stop(
    "Could not find: ", input_rds, "\n",
    "Please run Scripts/02A_Model_lnrr.R first."
  )
}

dat1_ml <- readRDS(input_rds)

if (!"Calcifier" %in% names(dat1_ml)) {
  stop(
    "Column 'Calcifier' not found in dat1_ml. ",
    "Please check that the Calcifier column is retained in 00B_Dependencies.R and 02A_Model_lnrr.R."
  )
}

cat("\n==== Loaded model-ready data from 02A ====\n")
cat("Rows in dat1_ml:", nrow(dat1_ml), "\n")
cat("Unique ES_ID:", n_distinct(dat1_ml$ES_ID), "\n")
cat("Unique Study_ID:", n_distinct(dat1_ml$Study_ID), "\n")
cat("Unique Common_id:", n_distinct(dat1_ml$Common_id), "\n\n")

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

plot_levels_calcifier <- list(
  response_levels = response_levels,
  category_levels = category_levels,
  treatment_levels = treatment_levels,
  treatment_label_levels = treatment_label_levels,
  calcifier_levels = calcifier_levels
)

saveRDS(
  plot_levels_calcifier,
  file.path(rds_dir, "plot_levels_calcifier.rds")
)

# ============================================================
# Clean Calcifier and Treatment columns
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
    Treatment_label = recode(
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
    Response_category = factor(Response_category, levels = response_levels)
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
  file.path(table_dir, "dat1_calcifier_model_ready.csv"),
  row.names = FALSE
)

saveRDS(
  dat1_calc,
  file.path(rds_dir, "dat1_calcifier_model_ready.rds")
)

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
# Main model-fitting function
# ============================================================

fit_pooled_model <- function(df, analysis_id, save_v = FALSE) {
  k <- nrow(df)
  
  if (k == 0) {
    return(tibble(
      analysis_id = analysis_id,
      model = NA_character_,
      mean = NA_real_,
      se = NA_real_,
      ci_low = NA_real_,
      ci_high = NA_real_,
      zval = NA_real_,
      pval = NA_real_,
      k = 0L,
      n_studies = 0L,
      n_common_id = 0L,
      sigma2 = NA_real_,
      QE = NA_real_,
      QE_df = NA_real_,
      QEp = NA_real_,
      logLik = NA_real_,
      AIC = NA_real_,
      BIC = NA_real_
    ))
  }
  
  if (k == 1) {
    return(tibble(
      analysis_id = analysis_id,
      model = "single_effect",
      mean = df$yi[1],
      se = sqrt(df$vi[1]),
      ci_low = df$yi[1] - 1.96 * sqrt(df$vi[1]),
      ci_high = df$yi[1] + 1.96 * sqrt(df$vi[1]),
      zval = NA_real_,
      pval = NA_real_,
      k = 1L,
      n_studies = dplyr::n_distinct(df$Study_ID),
      n_common_id = dplyr::n_distinct(df$Common_id),
      sigma2 = NA_real_,
      QE = NA_real_,
      QE_df = NA_real_,
      QEp = NA_real_,
      logLik = NA_real_,
      AIC = NA_real_,
      BIC = NA_real_
    ))
  }
  
  df_fit <- df %>%
    arrange(Study_ID, ES_ID) %>%
    mutate(
      Study_ID = factor(Study_ID),
      ES_ID = factor(ES_ID),
      Common_id = factor(Common_id),
      Calcifier = factor(Calcifier),
      Treatment = factor(Treatment),
      Treatment_label = factor(Treatment_label),
      Response_category = factor(Response_category)
    )
  
  V_matrix <- make_V_matrix(df_fit)
  
  if (save_v) {
    write.csv(
      V_matrix,
      file.path(table_dir, paste0("V_matrix_", analysis_id, ".csv")),
      row.names = TRUE
    )
    
    saveRDS(
      V_matrix,
      file.path(rds_dir, paste0("V_matrix_", analysis_id, ".rds"))
    )
  }
  
  fit <- metafor::rma.mv(
    yi = yi,
    V = V_matrix,
    data = df_fit,
    random = ~ 1 | Study_ID/ES_ID
  )
  
  coefs <- as.data.frame(coef(summary(fit)))
  
  sigma2_val <- if (!is.null(fit$sigma2) && length(fit$sigma2) > 0) {
    fit$sigma2[1]
  } else {
    NA_real_
  }
  
  tibble(
    analysis_id = analysis_id,
    model = "(~ 1 | Study_ID/ES_ID)",
    mean = coefs$estimate[1],
    se = coefs$se[1],
    ci_low = coefs$ci.lb[1],
    ci_high = coefs$ci.ub[1],
    zval = coefs$zval[1],
    pval = coefs$pval[1],
    k = as.integer(fit$k),
    n_studies = dplyr::n_distinct(df$Study_ID),
    n_common_id = dplyr::n_distinct(df$Common_id),
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
# 1. Treatment × response category × Calcifier pooled means
# ============================================================

calcifier_response_summary <- dat1_calc %>%
  filter(Response_category %in% response_levels) %>%
  group_by(Response_category, Calcifier, Treatment, Treatment_label) %>%
  group_split() %>%
  purrr::map_dfr(function(df_sub) {
    rc    <- as.character(df_sub$Response_category[1])
    calc  <- as.character(df_sub$Calcifier[1])
    trt   <- as.character(df_sub$Treatment[1])
    trt_l <- as.character(df_sub$Treatment_label[1])
    
    fit_pooled_model(
      df_sub,
      paste0("response_", rc, "_", calc, "_", trt),
      save_v = FALSE
    ) %>%
      mutate(
        Category = rc,
        Response_category = rc,
        Calcifier = calc,
        Treatment = trt,
        Treatment_label = trt_l
      )
  }) %>%
  mutate(
    Category = factor(Category, levels = category_levels),
    Response_category = factor(Response_category, levels = response_levels),
    Calcifier = factor(Calcifier, levels = calcifier_levels),
    Treatment = factor(Treatment, levels = treatment_levels),
    Treatment_label = factor(Treatment_label, levels = treatment_label_levels)
  )

write.csv(
  calcifier_response_summary,
  file.path(table_dir, "model_summary_treatment_response_calcifier.csv"),
  row.names = FALSE
)

saveRDS(
  calcifier_response_summary,
  file.path(rds_dir, "model_summary_treatment_response_calcifier.rds")
)

# ============================================================
# 2. Treatment × Overall × Calcifier pooled means
# ============================================================

calcifier_overall_summary <- dat1_calc %>%
  group_by(Calcifier, Treatment, Treatment_label) %>%
  group_split() %>%
  purrr::map_dfr(function(df_sub) {
    calc  <- as.character(df_sub$Calcifier[1])
    trt   <- as.character(df_sub$Treatment[1])
    trt_l <- as.character(df_sub$Treatment_label[1])
    
    fit_pooled_model(
      df_sub,
      paste0("overall_", calc, "_", trt),
      save_v = FALSE
    ) %>%
      mutate(
        Category = "Overall",
        Response_category = NA_character_,
        Calcifier = calc,
        Treatment = trt,
        Treatment_label = trt_l
      )
  }) %>%
  mutate(
    Category = factor(Category, levels = category_levels),
    Calcifier = factor(Calcifier, levels = calcifier_levels),
    Treatment = factor(Treatment, levels = treatment_levels),
    Treatment_label = factor(Treatment_label, levels = treatment_label_levels)
  )

write.csv(
  calcifier_overall_summary,
  file.path(table_dir, "model_summary_treatment_overall_calcifier.csv"),
  row.names = FALSE
)

saveRDS(
  calcifier_overall_summary,
  file.path(rds_dir, "model_summary_treatment_overall_calcifier.rds")
)

# ============================================================
# 3. Final pooled estimates for calcifier plot
# ============================================================

sum_df_calcifier <- bind_rows(
  calcifier_response_summary,
  calcifier_overall_summary
) %>%
  filter(is.finite(mean), is.finite(se), k > 0) %>%
  mutate(
    Category = factor(Category, levels = category_levels),
    Calcifier = factor(Calcifier, levels = calcifier_levels),
    Treatment = factor(Treatment, levels = treatment_levels),
    Treatment_label = factor(Treatment_label, levels = treatment_label_levels)
  )

write.csv(
  sum_df_calcifier,
  file.path(table_dir, "pooled_estimates_calcifier_one_panel.csv"),
  row.names = FALSE
)

saveRDS(
  sum_df_calcifier,
  file.path(rds_dir, "pooled_estimates_calcifier_one_panel.rds")
)

# ============================================================
# Compact summary table with percent change
# ============================================================

summary_table_calcifier <- sum_df_calcifier %>%
  mutate(
    Category = as.character(Category),
    Calcifier = as.character(Calcifier),
    Treatment = as.character(Treatment_label),
    `Mean lnRR` = sprintf("%.3f", mean),
    `95% CI` = sprintf("%.3f to %.3f", ci_low, ci_high),
    `p-value` = case_when(
      is.na(pval) ~ NA_character_,
      pval < 0.001 ~ "<0.001",
      TRUE ~ sprintf("%.3f", pval)
    ),
    `% change` = sprintf("%.1f%%", (exp(mean) - 1) * 100),
    `% change 95% CI` = sprintf(
      "%.1f%% to %.1f%%",
      (exp(ci_low) - 1) * 100,
      (exp(ci_high) - 1) * 100
    )
  ) %>%
  select(
    Category,
    Calcifier,
    Treatment,
    k,
    `Mean lnRR`,
    `95% CI`,
    `p-value`,
    `% change`,
    `% change 95% CI`
  ) %>%
  arrange(
    factor(Category, levels = category_levels),
    factor(Calcifier, levels = calcifier_levels),
    factor(Treatment, levels = treatment_label_levels)
  )

write.csv(
  summary_table_calcifier,
  file.path(table_dir, "pooled_estimates_calcifier_summary_table.csv"),
  row.names = FALSE
)

saveRDS(
  summary_table_calcifier,
  file.path(rds_dir, "pooled_estimates_calcifier_summary_table.rds")
)

# ============================================================
# Console prints
# ============================================================

cat("\n==== Calcifier response-category summaries ====\n")
print(calcifier_response_summary)

cat("\n==== Calcifier overall summaries ====\n")
print(calcifier_overall_summary)

cat("\n==== Final calcifier one-panel pooled estimates ====\n")
print(sum_df_calcifier)

cat("\n==== Saved model outputs ====\n")
cat("Tables:", table_dir, "\n")
cat("RDS   :", rds_dir, "\n")