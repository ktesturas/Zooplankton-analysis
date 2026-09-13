# ============================================================
# 02_Panel_Plots_Model_lnrr.R
# Purpose:
#   Use the shared cleaned zooplankton dataset from
#   Scripts/00_Dependencies.R, calculate lnRR and sampling variance,
#   fit fixed multilevel meta-analytic models, and create one final plot.
#
# Key updates:
#   1. Behaviour is lumped with Physiology.
#   2. The final plot has one panel only.
#   3. The rightmost "Overall" column from the old bottom panel
#      is pasted to the right side of the top panel as:
#      Meroplankton and Holoplankton.
#
# Input:
#   dat1 object from Scripts/00_Dependencies.R
#
# Outputs:
#   Outputs/Tables/02_Panel_Plots_Model_lnrr/
#   Outputs/RDS/02_Panel_Plots_Model_lnrr/
#   Outputs/Figures/02_Panel_Plots_Model_lnrr/
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

# ============================================================
# Script name and output folders
# ============================================================

script_name <- "02_Panel_Plots_Model_lnrr_new"

table_dir  <- file.path("Outputs", "Tables", script_name)
rds_dir    <- file.path("Outputs", "RDS", script_name)
figure_dir <- file.path("Outputs", "Figures", script_name)

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# Check shared cleaned dataset
# ============================================================

if (!exists("dat1")) {
  stop("dat1 not found. Check that Scripts/00B_Dependencies.R creates dat1.")
}

cat("\n==== Shared cleaned dataset ====\n")
cat("Rows in dat1:", nrow(dat1), "\n")
cat("Unique ES_ID:", n_distinct(dat1$ES_ID), "\n")
cat("Unique Study_ID:", n_distinct(dat1$Study_ID), "\n")
cat("Unique Common_id:", n_distinct(dat1$Common_id), "\n")
cat("ES_ID unique? ", n_distinct(dat1$ES_ID) == nrow(dat1), "\n\n")

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
    vi = escalc_out$vi,
    
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

final_label_levels <- c(
  "Growth",
  "Development",
  "Physiology",
  "Reproduction",
  "Survival",
  "Overall",
  "Meroplankton",
  "Holoplankton"
)

treatment_levels <- c("pH", "T", "TpH")
subgroup_levels  <- c("Meroplankton", "Holoplankton")

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
# Fixed multilevel model: Study_ID/ES_ID
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
      Broad_subgroup = factor(Broad_subgroup),
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
# 1. Overall treatment pooled means
#    This becomes the "Overall" column in the final plot.
# ============================================================

overall_treatment_summary <- purrr::map_dfr(
  treatment_levels,
  function(trt) {
    df_sub <- dat1_ml %>% filter(Treatment == trt)
    
    fit_pooled_model(
      df_sub,
      paste0("overall_", trt),
      save_v = TRUE
    ) %>%
      mutate(
        Treatment = trt,
        labels = "Overall"
      )
  }
) %>%
  mutate(
    Treatment = factor(Treatment, levels = treatment_levels),
    labels = factor(labels, levels = final_label_levels)
  )

write.csv(
  overall_treatment_summary,
  file.path(table_dir, "model_summary_overall_treatments.csv"),
  row.names = FALSE
)

saveRDS(
  overall_treatment_summary,
  file.path(rds_dir, "model_summary_overall_treatments.rds")
)

# ============================================================
# 2. Overall subgroup pooled means
#    These are across treatments and saved only as reference.
# ============================================================

overall_subgroup_summary <- purrr::map_dfr(
  subgroup_levels,
  function(subg) {
    df_sub <- dat1_ml %>% filter(Broad_subgroup == subg)
    
    fit_pooled_model(
      df_sub,
      paste0("overall_", subg),
      save_v = FALSE
    ) %>%
      mutate(Broad_subgroup = subg)
  }
) %>%
  mutate(
    Broad_subgroup = factor(Broad_subgroup, levels = subgroup_levels)
  )

write.csv(
  overall_subgroup_summary,
  file.path(table_dir, "model_summary_overall_subgroups.csv"),
  row.names = FALSE
)

saveRDS(
  overall_subgroup_summary,
  file.path(rds_dir, "model_summary_overall_subgroups.rds")
)

# ============================================================
# 3. Treatment × response-category pooled means
#    Behaviour has already been merged into Physiology.
# ============================================================

response_summary <- dat1_ml %>%
  filter(Response_category %in% response_levels) %>%
  group_by(Treatment, Response_category) %>%
  group_split() %>%
  purrr::map_dfr(function(df_sub) {
    trt <- as.character(df_sub$Treatment[1])
    rc  <- as.character(df_sub$Response_category[1])
    
    fit_pooled_model(
      df_sub,
      paste0("response_", trt, "_", rc),
      save_v = FALSE
    ) %>%
      mutate(
        Treatment = trt,
        labels = rc
      )
  }) %>%
  mutate(
    Treatment = factor(Treatment, levels = treatment_levels),
    labels = factor(labels, levels = final_label_levels)
  )

write.csv(
  response_summary,
  file.path(table_dir, "model_summary_treatment_response_category.csv"),
  row.names = FALSE
)

saveRDS(
  response_summary,
  file.path(rds_dir, "model_summary_treatment_response_category.rds")
)

# ============================================================
# 4. Treatment × subgroup overall pooled means
#    This is the rightmost "Overall" column from the old bottom panel.
#    These are pasted to the right side of the final one-panel plot.
# ============================================================

treatment_subgroup_overall_summary <- dat1_ml %>%
  filter(!is.na(Broad_subgroup)) %>%
  group_by(Treatment, Broad_subgroup) %>%
  group_split() %>%
  purrr::map_dfr(function(df_sub) {
    trt  <- as.character(df_sub$Treatment[1])
    subg <- as.character(df_sub$Broad_subgroup[1])
    
    fit_pooled_model(
      df_sub,
      paste0("subgroup_overall_", trt, "_", subg),
      save_v = FALSE
    ) %>%
      mutate(
        Treatment = trt,
        Broad_subgroup = subg,
        
        # Key step:
        # The old bottom-panel Overall estimates are relabelled as
        # Meroplankton and Holoplankton so they appear on the right.
        labels = subg
      )
  }) %>%
  mutate(
    Treatment = factor(Treatment, levels = treatment_levels),
    Broad_subgroup = factor(Broad_subgroup, levels = subgroup_levels),
    labels = factor(labels, levels = final_label_levels)
  )

write.csv(
  treatment_subgroup_overall_summary,
  file.path(table_dir, "model_summary_treatment_subgroup_overall.csv"),
  row.names = FALSE
)

saveRDS(
  treatment_subgroup_overall_summary,
  file.path(rds_dir, "model_summary_treatment_subgroup_overall.rds")
)

# ============================================================
# 5. Final pooled estimates for the one-panel plot
# ============================================================

sum_df_final <- bind_rows(
  response_summary,
  overall_treatment_summary,
  treatment_subgroup_overall_summary
) %>%
  filter(is.finite(mean), is.finite(se), k > 0) %>%
  mutate(
    Treatment = factor(Treatment, levels = treatment_levels),
    labels = factor(labels, levels = final_label_levels)
  )

write.csv(
  sum_df_final,
  file.path(table_dir, "pooled_estimates_final_one_panel.csv"),
  row.names = FALSE
)

saveRDS(
  sum_df_final,
  file.path(rds_dir, "pooled_estimates_final_one_panel.rds")
)

# ============================================================
# Combined model summary table
# ============================================================

all_model_summaries <- bind_rows(
  overall_treatment_summary,
  overall_subgroup_summary,
  response_summary,
  treatment_subgroup_overall_summary
)

write.csv(
  all_model_summaries,
  file.path(table_dir, "model_summary_all_pooled_means.csv"),
  row.names = FALSE
)

saveRDS(
  all_model_summaries,
  file.path(rds_dir, "model_summary_all_pooled_means.rds")
)

# ============================================================
# Raw points for plotting
# ============================================================

raw_all_final <- bind_rows(
  # Raw points for response-category columns
  dat1_ml %>%
    filter(Response_category %in% response_levels, is.finite(yi)) %>%
    transmute(
      labels = factor(Response_category, levels = final_label_levels),
      yi_raw = yi,
      Treatment = factor(Treatment, levels = treatment_levels)
    ),
  
  # Raw points for the overall treatment column
  dat1_ml %>%
    filter(is.finite(yi)) %>%
    transmute(
      labels = factor("Overall", levels = final_label_levels),
      yi_raw = yi,
      Treatment = factor(Treatment, levels = treatment_levels)
    ),
  
  # Raw points for the two subgroup-overall columns
  dat1_ml %>%
    filter(!is.na(Broad_subgroup), is.finite(yi)) %>%
    transmute(
      labels = factor(as.character(Broad_subgroup), levels = final_label_levels),
      yi_raw = yi,
      Treatment = factor(Treatment, levels = treatment_levels)
    )
)

# ============================================================
# Plot aesthetics
# ============================================================

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

treat_lab <- c(
  "pH" = "pH",
  "T" = "T",
  "TpH" = "T × pH"
)

pd <- position_dodge(width = 0.55)

pjd <- position_jitterdodge(
  jitter.width = 0.22,
  jitter.height = 0,
  dodge.width = 0.55
)

overall_idx  <- which(final_label_levels == "Overall")
overall_xmin <- overall_idx - 0.5
overall_xmax <- overall_idx + 0.5

subgroup_start_idx <- which(final_label_levels == "Meroplankton")
subgroup_end_idx   <- which(final_label_levels == "Holoplankton")
subgroup_xmin      <- subgroup_start_idx - 0.5
subgroup_xmax      <- subgroup_end_idx + 0.5

overall_fill  <- "grey95"
subgroup_fill <- "grey98"

# ============================================================
# Final one-panel plot
# ============================================================

p_final <- ggplot() +
  annotate(
    "rect",
    xmin = overall_xmin,
    xmax = overall_xmax,
    ymin = -Inf,
    ymax = Inf,
    fill = overall_fill,
    colour = NA
  ) +
  annotate(
    "rect",
    xmin = subgroup_xmin,
    xmax = subgroup_xmax,
    ymin = -Inf,
    ymax = Inf,
    fill = subgroup_fill,
    colour = NA
  ) +
  geom_hline(yintercept = brks_no0, colour = "grey90", linewidth = 0.3) +
  geom_vline(xintercept = vbreaks, colour = "grey90", linewidth = 0.3) +
  geom_hline(yintercept = 0, colour = "grey50", linewidth = 0.6) +
  
  # Separator between response categories and Overall
  geom_vline(
    xintercept = 5.5,
    linetype = "dotted",
    colour = "grey70",
    linewidth = 0.5
  ) +
  
  # Separator between Overall and plankton-group overall summaries
  geom_vline(
    xintercept = 6.5,
    linetype = "dotted",
    colour = "grey70",
    linewidth = 0.5
  ) +
  
  geom_point(
    data = raw_all_final,
    aes(x = labels, y = yi_raw, fill = Treatment),
    position = pjd,
    shape = 21,
    stroke = 0,
    alpha = 0.10,
    size = 2.4,
    show.legend = FALSE
  ) +
  geom_errorbar(
    data = sum_df_final,
    aes(x = labels, ymin = ci_low, ymax = ci_high, colour = Treatment),
    width = 0.2,
    linewidth = 0.9,
    position = pd
  ) +
  geom_point(
    data = sum_df_final,
    aes(x = labels, y = mean, colour = Treatment),
    position = pd,
    size = 3.0,
    shape = 16
  ) +
  geom_text(
    data = sum_df_final,
    aes(x = labels, y = Inf, label = k, colour = Treatment),
    position = pd,
    vjust = -1.2,
    family = fam,
    size = 4.1,
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
  xlab(NULL) +
  theme(
    text               = element_text(family = fam, colour = "grey10"),
    axis.text.y        = element_text(size = 16, colour = "grey10"),
    axis.title.y.left  = element_text(size = 18),
    axis.title.y.right = element_text(size = 18, hjust = 0.40, margin = margin(l = 10)),
    axis.text.x        = element_text(angle = 30, hjust = 1, size = 18, colour = "grey10"),
    panel.grid.major   = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.ticks         = element_line(colour = "grey10", linewidth = 0.6),
    axis.ticks.length  = unit(4, "pt"),
    legend.position    = "top",
    legend.margin      = margin(0, 0, 6, 0),
    plot.margin        = margin(20, 20, 20, 15)
  ) +
  coord_cartesian(clip = "off") +
  geom_rect(
    aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf),
    inherit.aes = FALSE,
    fill = NA,
    colour = "black",
    linewidth = 0.6
  ) +
  guides(
    colour = guide_legend(nrow = 1, byrow = TRUE)
  )

print(p_final)

# ============================================================
# Save final plot
# ============================================================

ggsave(
  filename = file.path(
    figure_dir,
    "zooplankton_experimental_treatment_one_panel_clean_StudyID_ESID.png"
  ),
  plot = p_final,
  width = 11,
  height = 4.8,
  dpi = 300,
  bg = "white",
  device = "png"
)

saveRDS(
  p_final,
  file.path(rds_dir, "zooplankton_experimental_treatment_one_panel_clean_StudyID_ESID.rds")
)

# ============================================================
# Console prints
# ============================================================

cat("\n==== Overall treatment model summaries ====\n")
print(overall_treatment_summary)

cat("\n==== Overall subgroup model summaries ====\n")
print(overall_subgroup_summary)

cat("\n==== Treatment × response-category summaries ====\n")
print(response_summary)

cat("\n==== Treatment × subgroup overall summaries ====\n")
print(treatment_subgroup_overall_summary)

cat("\n==== Final one-panel pooled estimates ====\n")
print(sum_df_final)

cat("\n==== Saved outputs ====\n")
cat("Tables :", table_dir, "\n")
cat("RDS    :", rds_dir, "\n")
cat("Figures:", figure_dir, "\n")
cat(
  "Final figure:",
  file.path(figure_dir, "zooplankton_experimental_treatment_one_panel_clean_StudyID_ESID.png"),
  "\n"
)

