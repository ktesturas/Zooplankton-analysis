# ============================================================
# 02_Panel_Plots_Model_lnrr.R
# Purpose:
#   Use the shared cleaned zooplankton dataset from
#   Scripts/00_Dependencies.R, calculate lnRR and sampling variance,
#   fit fixed multilevel meta-analytic models, and create panel plots.
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

source(file.path("Scripts", "00_Dependencies.R"))

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

script_name <- "02_Panel_Plots_Model_lnrr"

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
  stop("dat1 not found. Check that Scripts/00_Dependencies.R creates dat1.")
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
    vi = escalc_out$vi
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

label_levels <- c(
  "Growth",
  "Behaviour",
  "Development",
  "Physiology",
  "Reproduction",
  "Survival",
  "Overall"
)

cats <- label_levels[label_levels != "Overall"]

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
# ============================================================

overall_treatment_summary <- purrr::map_dfr(
  levels(dat1_ml$Treatment),
  function(trt) {
    df_sub <- dat1_ml %>% filter(Treatment == trt)
    fit_pooled_model(df_sub, paste0("overall_", trt), save_v = TRUE)
  }
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
# ============================================================

overall_subgroup_summary <- purrr::map_dfr(
  levels(dat1_ml$Broad_subgroup),
  function(subg) {
    if (is.na(subg)) return(NULL)
    
    df_sub <- dat1_ml %>% filter(Broad_subgroup == subg)
    
    fit_pooled_model(
      df_sub,
      paste0("overall_", subg),
      save_v = FALSE
    )
  }
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
# ============================================================

panel1_summary <- dat1_ml %>%
  filter(Response_category %in% cats) %>%
  group_by(Treatment, Response_category) %>%
  group_split() %>%
  purrr::map_dfr(function(df_sub) {
    trt <- as.character(df_sub$Treatment[1])
    rc  <- as.character(df_sub$Response_category[1])
    
    fit_pooled_model(
      df_sub,
      paste0("panel1_", trt, "_", rc),
      save_v = FALSE
    ) %>%
      mutate(
        Treatment = trt,
        labels = rc
      )
  }) %>%
  mutate(
    Treatment = factor(Treatment, levels = c("pH", "T", "TpH")),
    labels = factor(labels, levels = label_levels)
  )

panel1_overall <- overall_treatment_summary %>%
  mutate(
    Treatment = sub("^overall_", "", analysis_id),
    labels = "Overall"
  ) %>%
  mutate(
    Treatment = factor(Treatment, levels = c("pH", "T", "TpH")),
    labels = factor(labels, levels = label_levels)
  )

sum_df_trt <- bind_rows(panel1_summary, panel1_overall) %>%
  filter(is.finite(mean), is.finite(se), k > 0)

write.csv(
  sum_df_trt,
  file.path(table_dir, "pooled_estimates_panel1_treatment.csv"),
  row.names = FALSE
)

saveRDS(
  sum_df_trt,
  file.path(rds_dir, "pooled_estimates_panel1_treatment.rds")
)

# ============================================================
# 4. Treatment × subgroup pooled means
# ============================================================

subgroup_overall_summary <- dat1_ml %>%
  filter(!is.na(Broad_subgroup)) %>%
  group_by(Treatment, Broad_subgroup) %>%
  group_split() %>%
  purrr::map_dfr(function(df_sub) {
    trt  <- as.character(df_sub$Treatment[1])
    subg <- as.character(df_sub$Broad_subgroup[1])
    
    fit_pooled_model(
      df_sub,
      paste0("overall_", trt, "_", subg),
      save_v = FALSE
    ) %>%
      mutate(
        Treatment = trt,
        Broad_subgroup = subg
      )
  }) %>%
  mutate(
    Treatment = factor(Treatment, levels = c("pH", "T", "TpH")),
    Broad_subgroup = factor(Broad_subgroup, levels = c("Meroplankton", "Holoplankton"))
  )

write.csv(
  subgroup_overall_summary,
  file.path(table_dir, "model_summary_treatment_subgroup.csv"),
  row.names = FALSE
)

saveRDS(
  subgroup_overall_summary,
  file.path(rds_dir, "model_summary_treatment_subgroup.rds")
)

# ============================================================
# 5. Treatment × subgroup × response-category pooled means
# ============================================================

panel2_summary <- dat1_ml %>%
  filter(!is.na(Broad_subgroup), Response_category %in% cats) %>%
  group_by(Treatment, Broad_subgroup, Response_category) %>%
  group_split() %>%
  purrr::map_dfr(function(df_sub) {
    trt  <- as.character(df_sub$Treatment[1])
    subg <- as.character(df_sub$Broad_subgroup[1])
    rc   <- as.character(df_sub$Response_category[1])
    
    fit_pooled_model(
      df_sub,
      paste0("panel2_", trt, "_", subg, "_", rc),
      save_v = FALSE
    ) %>%
      mutate(
        Treatment = trt,
        Broad_subgroup = subg,
        labels = rc
      )
  }) %>%
  mutate(
    Treatment = factor(Treatment, levels = c("pH", "T", "TpH")),
    Broad_subgroup = factor(Broad_subgroup, levels = c("Meroplankton", "Holoplankton")),
    labels = factor(labels, levels = label_levels)
  )

panel2_overall <- subgroup_overall_summary %>%
  mutate(labels = "Overall") %>%
  mutate(
    labels = factor(labels, levels = label_levels)
  )

sum_df_trt_sub <- bind_rows(panel2_summary, panel2_overall) %>%
  filter(is.finite(mean), is.finite(se), k > 0)

write.csv(
  sum_df_trt_sub,
  file.path(table_dir, "pooled_estimates_panel2_treatment_subgroup.csv"),
  row.names = FALSE
)

saveRDS(
  sum_df_trt_sub,
  file.path(rds_dir, "pooled_estimates_panel2_treatment_subgroup.rds")
)

# ============================================================
# Combined model summary table
# ============================================================

all_model_summaries <- bind_rows(
  overall_treatment_summary,
  overall_subgroup_summary,
  panel1_summary,
  panel1_overall,
  subgroup_overall_summary,
  panel2_summary,
  panel2_overall
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

raw_all_trt <- bind_rows(
  dat1_ml %>%
    filter(Response_category %in% cats, is.finite(yi)) %>%
    transmute(
      labels    = factor(Response_category, levels = label_levels),
      yi_raw    = yi,
      Treatment = factor(Treatment, levels = c("pH", "T", "TpH"))
    ),
  dat1_ml %>%
    filter(is.finite(yi)) %>%
    transmute(
      labels    = factor("Overall", levels = label_levels),
      yi_raw    = yi,
      Treatment = factor(Treatment, levels = c("pH", "T", "TpH"))
    )
)

raw_all_trt_sub <- bind_rows(
  dat1_ml %>%
    filter(!is.na(Broad_subgroup), Response_category %in% cats, is.finite(yi)) %>%
    transmute(
      labels = factor(Response_category, levels = label_levels),
      yi_raw = yi,
      Treatment = factor(Treatment, levels = c("pH", "T", "TpH")),
      Broad_subgroup = factor(Broad_subgroup, levels = c("Meroplankton", "Holoplankton"))
    ),
  dat1_ml %>%
    filter(!is.na(Broad_subgroup), is.finite(yi)) %>%
    transmute(
      labels = factor("Overall", levels = label_levels),
      yi_raw = yi,
      Treatment = factor(Treatment, levels = c("pH", "T", "TpH")),
      Broad_subgroup = factor(Broad_subgroup, levels = c("Meroplankton", "Holoplankton"))
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
vbreaks  <- seq(0.5, length(label_levels) + 0.5, by = 1)

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

pd  <- position_dodge(width = 0.55)
pjd <- position_jitterdodge(
  jitter.width = 0.22,
  jitter.height = 0,
  dodge.width = 0.55
)

pd2  <- position_dodge(width = 0.75)
pjd2 <- position_jitterdodge(
  jitter.width = 0.22,
  jitter.height = 0,
  dodge.width = 0.75
)

overall_idx  <- which(label_levels == "Overall")
overall_xmin <- overall_idx - 0.5
overall_xmax <- overall_idx + 0.5
overall_fill <- "grey95"

# ============================================================
# Panel 1: Treatment × response category
# ============================================================

p_base <- ggplot() +
  annotate(
    "rect",
    xmin = overall_xmin,
    xmax = overall_xmax,
    ymin = -Inf,
    ymax = Inf,
    fill = overall_fill,
    colour = NA
  ) +
  geom_hline(yintercept = brks_no0, colour = "grey90", linewidth = 0.3) +
  geom_vline(xintercept = vbreaks, colour = "grey90", linewidth = 0.3) +
  geom_hline(yintercept = 0, colour = "grey50", linewidth = 0.6) +
  geom_vline(xintercept = 6.5, linetype = "dotted", colour = "grey70", linewidth = 0.5) +
  geom_point(
    data = raw_all_trt,
    aes(x = labels, y = yi_raw, fill = Treatment),
    position = pjd,
    shape = 21,
    stroke = 0,
    alpha = 0.10,
    size = 2.4,
    show.legend = FALSE
  ) +
  geom_errorbar(
    data = sum_df_trt,
    aes(x = labels, ymin = ci_low, ymax = ci_high, colour = Treatment),
    width = 0.2,
    linewidth = 0.9,
    position = pd
  ) +
  geom_point(
    data = sum_df_trt,
    aes(x = labels, y = mean, colour = Treatment),
    position = pd,
    size = 3.0,
    shape = 16
  ) +
  geom_text(
    data = sum_df_trt,
    aes(x = labels, y = Inf, label = k, colour = Treatment),
    position = pd,
    vjust = -1.2,
    family = fam,
    size = 4.1,
    show.legend = FALSE
  ) +
  scale_colour_manual(values = col_trt_line, labels = treat_lab, name = NULL) +
  scale_fill_manual(values = fill_trt, labels = treat_lab, guide = "legend") +
  scale_x_discrete(drop = FALSE, limits = label_levels, expand = expansion(add = 0)) +
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
  theme(
    text               = element_text(family = fam, colour = "grey10"),
    axis.text.y        = element_text(size = 16, colour = "grey10"),
    axis.title.y.left  = element_text(size = 18),
    axis.title.y.right = element_text(size = 18, hjust = 0.40, margin = margin(l = 10)),
    panel.grid.major   = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.ticks         = element_line(colour = "grey10", linewidth = 0.6),
    axis.ticks.length  = unit(4, "pt"),
    legend.position    = "top",
    legend.margin      = margin(0, 0, 6, 0),
    plot.margin        = margin(20, 20, 8, 15)
  ) +
  coord_cartesian(clip = "off") +
  geom_rect(
    aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf),
    inherit.aes = FALSE,
    fill = NA,
    colour = "black",
    linewidth = 0.6
  )

p1 <- p_base +
  theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none"
  )

# ============================================================
# Panel 2: Treatment × subgroup × response category
# ============================================================

p2 <- ggplot() +
  annotate(
    "rect",
    xmin = overall_xmin,
    xmax = overall_xmax,
    ymin = -Inf,
    ymax = Inf,
    fill = overall_fill,
    colour = NA
  ) +
  geom_hline(yintercept = brks_no0, colour = "grey90", linewidth = 0.3) +
  geom_vline(xintercept = vbreaks, colour = "grey90", linewidth = 0.3) +
  geom_hline(yintercept = 0, colour = "grey50", linewidth = 0.6) +
  geom_vline(xintercept = 6.5, linetype = "dotted", colour = "grey70", linewidth = 0.5) +
  geom_point(
    data = raw_all_trt_sub,
    aes(
      x = labels,
      y = yi_raw,
      fill = Treatment,
      group = interaction(Treatment, Broad_subgroup)
    ),
    position = pjd2,
    shape = 21,
    stroke = 0,
    alpha = 0.10,
    size = 2.4,
    show.legend = FALSE
  ) +
  geom_errorbar(
    data = sum_df_trt_sub,
    aes(
      x = labels,
      ymin = ci_low,
      ymax = ci_high,
      colour = Treatment,
      group = interaction(Treatment, Broad_subgroup)
    ),
    width = 0.2,
    linewidth = 0.9,
    position = pd2
  ) +
  geom_point(
    data = sum_df_trt_sub,
    aes(
      x = labels,
      y = mean,
      colour = Treatment,
      shape = Broad_subgroup,
      group = interaction(Treatment, Broad_subgroup)
    ),
    position = pd2,
    size = 3.0
  ) +
  geom_text(
    data = sum_df_trt_sub,
    aes(
      x = labels,
      y = Inf,
      label = k,
      colour = Treatment,
      group = interaction(Treatment, Broad_subgroup)
    ),
    position = pd2,
    vjust = -1.2,
    family = fam,
    size = 4.1,
    show.legend = FALSE
  ) +
  scale_colour_manual(values = col_trt_line, labels = treat_lab, name = NULL) +
  scale_fill_manual(values = fill_trt, labels = treat_lab, guide = "legend") +
  scale_shape_manual(
    values = c("Meroplankton" = 17, "Holoplankton" = 4),
    drop = FALSE,
    name = NULL
  ) +
  scale_x_discrete(drop = FALSE, limits = label_levels, expand = expansion(add = 0)) +
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
    legend.box         = "horizontal",
    plot.margin        = margin(8, 20, 20, 15)
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
    colour = guide_legend(nrow = 1, byrow = TRUE),
    shape  = guide_legend(nrow = 1, byrow = TRUE)
  )

# ============================================================
# Stack panels
# ============================================================

gap <- plot_spacer()

p_stack <- p1 / gap / p2 +
  plot_layout(heights = c(1, 0.12, 1), guides = "collect") +
  plot_annotation(
    theme = theme(
      legend.position = "top",
      legend.box = "horizontal",
      text = element_text(family = fam)
    )
  )

print(p_stack)

# ============================================================
# Save final plot
# ============================================================

# DO MANUAL EXPORTING
# ON THE PLOT PANEL, CLICK EXPORT. BUT BEFORE THAT, MAKE SURE THE WHOLE PLOT PANELS IS STRETCHED ACROSS ENTIRE SCREEN
# THIS IS BECAUSE THE EXPORT FUNCTION WOULD JUST EXPORT AS TO HOW THE PLOT CURRENTLY LOOKS LIKE ON SCREEN

ggsave(
  filename = file.path(
    figure_dir,
    "zooplankton_experimental_treatment_panels_clean_StudyID_ESID.png"
  ),
  plot = p_stack,
  width = 10,
  height = 6,
  dpi = 300,
  bg = "white",
  device = "png"
)

saveRDS(
  p_stack,
  file.path(rds_dir, "zooplankton_experimental_treatment_panels_clean_StudyID_ESID.rds")
)

# ============================================================
# Console prints
# ============================================================

cat("\n==== Overall treatment model summaries ====\n")
print(overall_treatment_summary)

cat("\n==== Overall subgroup model summaries ====\n")
print(overall_subgroup_summary)

cat("\n==== First rows of pooled estimates for panel 1 ====\n")
print(head(sum_df_trt, 10))

cat("\n==== First rows of pooled estimates for panel 2 ====\n")
print(head(sum_df_trt_sub, 10))

cat("\n==== Saved outputs ====\n")
cat("Tables :", table_dir, "\n")
cat("RDS    :", rds_dir, "\n")
cat("Figures:", figure_dir, "\n")
cat(
  "Final figure:",
  file.path(figure_dir, "zooplankton_experimental_treatment_panels_clean_StudyID_ESID.png"),
  "\n"
)