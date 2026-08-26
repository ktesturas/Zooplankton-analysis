# ===========================================
# 03_summary_tables_zooplankton.R
# Summary tables for the three-panel observational figure
# ===========================================

# ---- Packages ----
packages <- c(
  "dplyr", "readxl", "stringr", "tidyr", "ggplot2",
  "scales", "forcats", "here", "openxlsx", "gridExtra", "grid"
)

to_install <- packages[!packages %in% rownames(installed.packages())]
if (length(to_install) > 0) {
  install.packages(to_install)
}

library(dplyr)
library(readxl)
library(stringr)
library(tidyr)
library(ggplot2)
library(scales)
library(forcats)
library(here)
library(openxlsx)
library(gridExtra)
library(grid)

# ---- File paths ----
script_name <- "03_summary_tables_zooplankton"

find_input_file <- function(paths, label) {
  existing <- paths[file.exists(paths)]

  if (length(existing) > 0) {
    return(existing[1])
  }

  stop(
    paste0(
      "Cannot find ", label, ". Checked:\n",
      paste(paths, collapse = "\n")
    ),
    call. = FALSE
  )
}

fp_zoo <- find_input_file(
  c(
    here::here("Data", "Observational", "Plankton observational.xlsx"),
    here::here("Data", "Plankton observational.xlsx")
  ),
  "Plankton observational.xlsx"
)

fp_gen <- find_input_file(
  c(
    here::here("Data", "All marine life", "All_marine_life_data.xlsx"),
    here::here("Data", "All_marine_life", "All_marine_life_data.xlsx"),
    here::here("Data", "All_marine_life_data.xlsx")
  ),
  "All_marine_life_data.xlsx"
)

# ---- Output folder ----
out_dir <- here::here("Outputs", script_name)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ===========================================
# Helper functions
# ===========================================

pluralise_custom <- function(x) {
  out <- str_squish(as.character(x))
  xl  <- str_to_lower(out)
  
  out[grepl("\\bfish\\s*larva(e|es)?\\b", xl, perl = TRUE)] <- "Fish larvae"
  out[xl == "other"] <- "Other zooplankton"
  
  keep_end_re <- "(fish|larvae|zooplankton|phytoplankton|plankton|algae|foraminifera)\\s*$"
  need_s <- !grepl(keep_end_re, xl) & !stringr::str_ends(out, "s")
  out[need_s] <- paste0(out[need_s], "s")
  
  out
}

fix_taxon <- function(x) {
  out <- str_squish(as.character(x))
  out <- str_replace_all(out, "\\s*\\.(?=\\s|$)", "")
  out <- str_replace_all(out, regex("\\bFish\\s*larvaes\\b", ignore_case = TRUE), "Fish larvae")
  out <- str_replace_all(out, regex("\\bBenthic\\s+crustaceas\\b", ignore_case = TRUE), "Benthic crustaceans")
  out <- str_replace_all(out, regex("\\bBenthic\\s+algaes\\b", ignore_case = TRUE), "Benthic algae")
  out <- str_replace_all(
    out,
    regex("\\bBenthic\\s+invertebrates\\s*\\(other\\)?s?\\)?\\b", ignore_case = TRUE),
    "Other benthic invertebrates"
  )
  out <- str_replace_all(out, c(
    "\\bChateognaths\\b"           = "Chaetognaths",
    "\\bPolychaete\\b"             = "Polychaetes",
    "\\bDecapod\\b"                = "Decapods",
    "\\bCalanoid\\b"               = "Calanoids",
    "\\bNon[ -]?calanoid(s)?\\b"   = "Non-calanoids",
    "\\bEchinoderm\\b"             = "Echinoderms",
    "\\bBenthic invert\\.?\\b"     = "Benthic invertebrates",
    "\\bBenthic crustacea\\b"      = "Benthic crustaceans",
    "\\bBenthic cnidaria(ns)?\\b"  = "Benthic cnidarians",
    "\\bZooplankotn\\b"            = "Zooplankton",
    "\\bPhytoplanktoN\\b"          = "Phytoplankton"
  ))
  out
}

normalise_taxon <- function(x) {
  x |> fix_taxon() |> pluralise_custom() |> fix_taxon()
}

pretty_taxon <- function(x) {
  stringr::str_replace_all(x, c(
    "\\bPolychates\\b"     = "Polychaetes",
    "\\bPhytoplanktons\\b" = "Phytoplankton"
  ))
}

recode_leftblock_subgroup <- function(group, subgroup) {
  g_is_zoo <- grepl("^zooplankton", group, ignore.case = TRUE)
  out <- subgroup
  
  out[g_is_zoo & grepl("^pteropod", subgroup, ignore.case = TRUE)] <- "Molluscs"
  
  out[g_is_zoo & (
    grepl("^mysid", subgroup, ignore.case = TRUE) |
      grepl("^cirrip", subgroup, ignore.case = TRUE) |
      grepl("^amphipod", subgroup, ignore.case = TRUE)
  )] <- "Other crustaceans"
  
  out
}

leftblock_exclude_flag <- function(group, subgroup) {
  g_is_zoo <- grepl("^zooplankton", group, ignore.case = TRUE)
  
  bad <- grepl("^rotifer", subgroup, ignore.case = TRUE) |
    grepl("^dinoflagellate", subgroup, ignore.case = TRUE)
  
  g_is_zoo & bad
}

zoo_kind <- function(g) {
  case_when(
    grepl("Zooplankton\\s*-\\s*holo", g, ignore.case = TRUE) ~ "Holoplankton",
    grepl("Zooplankton\\s*-\\s*mero", g, ignore.case = TRUE) ~ "Meroplankton",
    TRUE ~ NA_character_
  )
}

norm_consistency <- function(x) {
  ux <- toupper(str_squish(as.character(x)))
  ux[is.na(x)] <- NA_character_
  
  dplyr::case_when(
    ux %in% c("Y", "YES", "C", "CONSISTENT") ~ "Consistent",
    ux %in% c("N", "NC", "NO", "NO CHANGE") ~ "Not consistent",
    ux %in% c("EQUIVOCAL", "NEUTRAL", "0", "-") ~ "Unsure",
    ux %in% c("", "NA") ~ NA_character_,
    TRUE ~ "Unsure"
  )
}

canonical_taxon <- function(x, for_metric = c("phenology", "distribution")) {
  for_metric <- match.arg(for_metric)
  
  s <- str_squish(as.character(x))
  s <- str_replace_all(s, "\\s*\\.(?=\\s|$)", "")
  s <- str_replace_all(s, "^Amphipod$", "Amphipods")
  s <- str_replace_all(s, "^Decapod$", "Decapods")
  s <- str_replace_all(s, "^Cladoceran$", "Cladocerans")
  s <- str_replace_all(s, "^Rotifer$", "Rotifers")
  s <- str_replace_all(s, "^Polychaete$", "Polychaetes")
  s <- str_replace_all(s, "^Ctenophores?$", "Ctenophores")
  s <- str_replace_all(s, "^Cha?etognaths?$", "Chaetognaths")
  s <- str_replace_all(s, "^Non-?calanoids?$", "Non-calanoids")
  s <- str_replace_all(s, "^Calanoid[s]?$", "Calanoids")
  s <- str_replace_all(s, "^Cnidarian[s]?$", "Cnidarians")
  s <- str_replace_all(s, "^Bivalve[s]?$", "Bivalves")
  s <- str_replace_all(s, "^Gastropod[s]?$", "Gastropods")
  s <- str_replace_all(s, "^Foraminifera$", "Foraminifera")
  s <- str_replace_all(s, "^Euphausiids?$", "Euphausiids")
  s <- str_replace_all(s, "^Benthic invertebrate[s]?$", "Benthic invertebrates")
  s <- str_replace_all(s, "^Benthic crustacea$", "Benthic crustaceans")
  s <- str_replace_all(s, "^Benthic cnidarian[s]?$", "Benthic cnidarians")
  s <- str_replace_all(s, "^Benthic mollusc[s]?$", "Benthic molluscs")
  s <- str_replace_all(s, "^Benthic algae$", "Benthic algae")
  s <- str_replace_all(s, "^Appendicularians?$", "Appendicularians")
  s <- str_replace_all(s, "^Seabird[s]?$", "Seabirds")
  s <- str_replace_all(s, "^Mammal[s]?$", "Mammals")
  s <- str_replace_all(s, "^Reptile[s]?$", "Reptiles")
  s <- str_replace_all(s, "^Fish larvae$", "Fish larvae")
  s <- str_replace_all(s, "^Zooplankton$", "Zooplankton")
  s <- str_replace_all(s, "^Phytoplankton$", "Phytoplankton")
  s <- str_replace_all(s, "^Bony fish$", "Bony fish")
  s <- str_replace_all(s, "^Non-bony fish$", "Non-bony fish")
  s <- str_replace_all(s, "^Fish$", "Fish")
  
  if (for_metric == "distribution") {
    s <- ifelse(s == "Bony fish", "Fish", s)
  }
  
  s
}

summarise_mean_ci <- function(df) {
  df %>%
    summarise(
      n = sum(is.finite(Rate)),
      mean_value = mean(Rate, na.rm = TRUE),
      sd = sd(Rate, na.rm = TRUE),
      median_value = median(Rate, na.rm = TRUE),
      min_value = min(Rate, na.rm = TRUE),
      max_value = max(Rate, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      se = ifelse(n > 1, sd / sqrt(n), NA_real_),
      tcrit = ifelse(n > 1, qt(0.975, df = n - 1), NA_real_),
      ci_low_raw = ifelse(n > 1, mean_value - tcrit * se, NA_real_),
      ci_high_raw = ifelse(n > 1, mean_value + tcrit * se, NA_real_),
      ci_low = ifelse(n >= 3, ci_low_raw, NA_real_),
      ci_high = ifelse(n >= 3, ci_high_raw, NA_real_),
      ci_shown = n >= 3
    )
}

wilson_ci <- function(k, n, conf = 0.95) {
  z <- qnorm(1 - (1 - conf) / 2)
  phat <- ifelse(n > 0, k / n, NA_real_)
  
  denom <- 1 + (z^2 / n)
  centre <- (phat + (z^2 / (2 * n))) / denom
  half <- (z * sqrt((phat * (1 - phat) / n) + (z^2 / (4 * n^2)))) / denom
  
  tibble(
    proportion = phat,
    percentage = phat * 100,
    ci_low = pmax(0, (centre - half) * 100),
    ci_high = pmin(100, (centre + half) * 100)
  )
}

# ===========================================
# Distribution data
# ===========================================

dist_raw0 <- read_xlsx(fp_zoo, sheet = "Distribution") %>%
  mutate(
    Group       = str_squish(as.character(Group)),
    Subgroup    = str_squish(as.character(Subgroup)),
    Subgroup2   = recode_leftblock_subgroup(Group, Subgroup),
    ExcludeLeft = leftblock_exclude_flag(Group, Subgroup),
    Rate        = suppressWarnings(as.numeric(Rate)),
    Timespan    = suppressWarnings(as.numeric(Timespan)),
    ZooKind     = zoo_kind(Group)
  ) %>%
  filter(!is.na(Rate), !is.na(Subgroup), Subgroup != "NA")

dist_zoo <- dist_raw0 %>%
  filter(grepl("^Zooplankton", Group), !ExcludeLeft) %>%
  { if ("Timespan" %in% names(.)) filter(., Timespan >= 19) else . } %>%
  transmute(
    Taxon = normalise_taxon(ifelse(str_to_lower(Subgroup2) == "other", "Other zooplankton", Subgroup2)),
    Rate,
    Block = "Zooplankton"
  )

dist_zoo_kind <- dist_raw0 %>%
  filter(grepl("^Zooplankton", Group), !ExcludeLeft, !is.na(ZooKind)) %>%
  { if ("Timespan" %in% names(.)) filter(., Timespan >= 19) else . } %>%
  transmute(
    Taxon = ZooKind,
    Rate,
    Block = "Zooplankton"
  )

dist_pool <- dist_raw0 %>%
  filter(grepl("^Zooplankton", Group), !ExcludeLeft) %>%
  { if ("Timespan" %in% names(.)) filter(., Timespan >= 19) else . } %>%
  transmute(
    Taxon = "Zooplankton",
    Rate,
    Block = "Phytoplankton"
  )

gen_raw_dist <- read_xlsx(fp_gen, sheet = "Distribution") %>%
  mutate(
    Taxa = str_squish(as.character(Taxa)),
    Rate = suppressWarnings(as.numeric(Rate))
  ) %>%
  filter(!is.na(Taxa), Taxa != "NA", !is.na(Rate)) %>%
  filter(Database != "Kris")

dist_gen <- gen_raw_dist %>%
  mutate(
    merged = case_when(
      str_detect(Taxa, regex("\\bZooplankton\\b.*\\b(mero|holo)\\b", ignore_case = TRUE)) ~ "Zooplankton",
      TRUE ~ Taxa
    )
  ) %>%
  transmute(
    Taxon = normalise_taxon(merged),
    Rate,
    Block = "General marine taxa"
  )

dist_all <- gen_raw_dist %>%
  transmute(
    Taxon = "All marine taxa",
    Rate,
    Block = "General marine taxa"
  )

dist_all_raw <- bind_rows(
  dist_zoo,
  dist_zoo_kind,
  dist_pool,
  dist_gen,
  dist_all
) %>%
  mutate(
    Panel = "Distribution change",
    Study = case_when(
      Block %in% c("Zooplankton", "Phytoplankton") ~ "This study",
      Block == "General marine taxa" ~ "Poloczanska et al. (2016)",
      TRUE ~ NA_character_
    )
  )

# ===========================================
# Phenology data
# ===========================================

zp_raw <- read_xlsx(fp_zoo, sheet = "Phenology") %>%
  filter(is.na(Include) | Include != "No") %>%
  mutate(
    Rate        = suppressWarnings(as.numeric(Rate)),
    Group       = str_squish(as.character(Group)),
    Subgroup    = str_squish(as.character(Subgroup)),
    Subgroup2   = recode_leftblock_subgroup(Group, Subgroup),
    ExcludeLeft = leftblock_exclude_flag(Group, Subgroup),
    ZooKind     = zoo_kind(Group)
  ) %>%
  filter(!is.na(Rate), !is.na(Subgroup), Subgroup != "NA")

phen_zoo <- zp_raw %>%
  filter(grepl("^Zooplankton", Group), !ExcludeLeft) %>%
  transmute(
    Taxon = Subgroup2 %>%
      (\(x) ifelse(str_to_lower(x) == "other", "Other zooplankton", x))() %>%
      pluralise_custom() %>%
      fix_taxon(),
    Rate,
    Block = "Zooplankton"
  )

phen_zoo_kind <- zp_raw %>%
  filter(grepl("^Zooplankton", Group), !ExcludeLeft, !is.na(ZooKind)) %>%
  transmute(
    Taxon = ZooKind,
    Rate,
    Block = "Zooplankton"
  )

phen_pool <- zp_raw %>%
  filter(grepl("^Zooplankton", Group), !ExcludeLeft) %>%
  transmute(
    Taxon = "Zooplankton",
    Rate,
    Block = "Phytoplankton"
  )

ipcc_raw <- read_xlsx(fp_gen, sheet = "Phenology") %>%
  filter(is.na(Include) | Include != "N") %>%
  filter(Database != "Kris") %>%
  mutate(
    Rate = suppressWarnings(as.numeric(Shift)),
    Group = str_squish(as.character(Group))
  ) %>%
  filter(!is.na(Rate), !is.na(Group), Group != "NA")

phen_gen <- ipcc_raw %>%
  mutate(
    Group = case_when(
      str_detect(Group, regex("Mero", ignore_case = TRUE)) ~ "Zooplankton",
      str_detect(Group, regex("Zooplankton", ignore_case = TRUE)) ~ "Zooplankton",
      TRUE ~ Group
    )
  ) %>%
  transmute(
    Taxon = pluralise_custom(Group) %>% fix_taxon(),
    Rate,
    Block = "General marine taxa"
  )

phen_all <- ipcc_raw %>%
  transmute(
    Taxon = "All marine taxa",
    Rate,
    Block = "General marine taxa"
  )

phen_all_raw <- bind_rows(
  phen_gen,
  phen_all,
  phen_zoo,
  phen_zoo_kind,
  phen_pool
) %>%
  mutate(
    Panel = "Phenology change",
    Study = case_when(
      Block %in% c("Zooplankton", "Phytoplankton") ~ "This study",
      Block == "General marine taxa" ~ "Cooley et al. (2022)",
      TRUE ~ NA_character_
    )
  )

# ===========================================
# Master column order
# ===========================================

sep <- " || "

sum_dist <- dist_all_raw %>%
  group_by(Block, Taxon) %>%
  summarise(mean_rate = mean(Rate, na.rm = TRUE), .groups = "drop")

sum_phen <- phen_all_raw %>%
  group_by(Block, Taxon) %>%
  summarise(mean_rate = mean(Rate, na.rm = TRUE), .groups = "drop")

z_order_phen <- sum_phen %>%
  filter(Block == "Zooplankton") %>%
  arrange(mean_rate, Taxon) %>%
  pull(Taxon)

p_order_phen <- sum_phen %>%
  filter(Block == "Phytoplankton") %>%
  arrange(mean_rate, Taxon) %>%
  pull(Taxon)

g_order_phen <- sum_phen %>%
  filter(Block == "General marine taxa") %>%
  arrange(mean_rate, Taxon) %>%
  pull(Taxon)

append_block <- function(block_name, base_order, sum_other) {
  present <- base_order
  
  candidates <- sum_other %>%
    filter(Block == block_name) %>%
    pull(Taxon)
  
  only_other <- setdiff(candidates, present)
  
  if (!length(only_other)) return(base_order)
  
  add <- sum_other %>%
    filter(Block == block_name, Taxon %in% only_other) %>%
    arrange(mean_rate, Taxon) %>%
    pull(Taxon)
  
  c(base_order, add)
}

z_order <- append_block("Zooplankton", z_order_phen, sum_dist)
p_order <- append_block("Phytoplankton", p_order_phen, sum_dist)

g_order_tmp <- append_block("General marine taxa", g_order_phen, sum_dist)

g_order <- g_order_tmp %>%
  { c("Zooplankton", setdiff(., c("Zooplankton", "All marine taxa")), "All marine taxa") } %>%
  unique()

move_items_before <- function(vec, items, before_item = "Zooplankton") {
  items <- items[items %in% vec]
  
  if (!length(items)) return(vec)
  
  vec_wo <- vec[!vec %in% items]
  pos <- match(before_item, vec_wo)
  
  if (is.na(pos)) {
    c(vec_wo, items)
  } else {
    c(vec_wo[seq_len(pos - 1)], items, vec_wo[pos:length(vec_wo)])
  }
}

z_order <- move_items_before(z_order, c("Holoplankton", "Meroplankton"), "Zooplankton")

z_levels <- paste("Zooplankton", z_order, sep = sep)
p_levels <- paste("Phytoplankton", p_order, sep = sep)
g_levels <- paste("General marine taxa", g_order, sep = sep)

x_levels_master <- c(z_levels, p_levels, g_levels)

column_lookup <- tibble(
  x_order = seq_along(x_levels_master),
  Xkey = factor(x_levels_master, levels = x_levels_master),
  Block = sub("\\s\\|\\|\\s.*$", "", as.character(Xkey)),
  Taxon = sub("^[^|]+\\s\\|\\|\\s", "", as.character(Xkey))
)

# ===========================================
# Mean-rate summary tables
# ===========================================

make_mean_summary <- function(raw_df, panel_name, unit_label, ci_method_label) {
  raw_df %>%
    mutate(
      Xkey = factor(paste(Block, Taxon, sep = sep), levels = x_levels_master)
    ) %>%
    filter(!is.na(Xkey)) %>%
    group_by(Panel, Study, Block, Taxon, Xkey) %>%
    summarise_mean_ci() %>%
    ungroup() %>%
    right_join(
      column_lookup,
      by = c("Xkey", "Block", "Taxon")
    ) %>%
    mutate(
      Panel = panel_name,
      Study = case_when(
        panel_name == "Distribution change" & Block %in% c("Zooplankton", "Phytoplankton") ~ "This study",
        panel_name == "Distribution change" & Block == "General marine taxa" ~ "Poloczanska et al. (2016)",
        panel_name == "Phenology change" & Block %in% c("Zooplankton", "Phytoplankton") ~ "This study",
        panel_name == "Phenology change" & Block == "General marine taxa" ~ "Cooley et al. (2022)",
        TRUE ~ Study
      ),
      n = replace_na(n, 0L),
      unit = unit_label,
      ci_method = ci_method_label,
      value_type = "Mean",
      denominator_note = NA_character_
    ) %>%
    select(
      Panel, Study, x_order, Block, Taxon, n,
      value_type, mean_value, ci_low, ci_high, ci_shown,
      unit, ci_method,
      sd, se, median_value, min_value, max_value,
      ci_low_raw, ci_high_raw,
      denominator_note
    ) %>%
    arrange(x_order)
}

distribution_summary <- make_mean_summary(
  raw_df = dist_all_raw,
  panel_name = "Distribution change",
  unit_label = "km decade^-1",
  ci_method_label = "t-based 95% CI around the mean; shown only when n >= 3"
)

phenology_summary <- make_mean_summary(
  raw_df = phen_all_raw,
  panel_name = "Phenology change",
  unit_label = "days decade^-1",
  ci_method_label = "t-based 95% CI around the mean; shown only when n >= 3"
)

# ===========================================
# Percentage-consistent data
# ===========================================

# This study phenology consistency
phen_raw_cons <- read_xlsx(fp_zoo, sheet = "Phenology") %>%
  filter(is.na(Include) | Include != "No") %>%
  filter(Group != "Phytoplankton") %>%
  mutate(
    Group       = str_squish(as.character(Group)),
    Subgroup    = na_if(str_squish(as.character(Subgroup)), "NA"),
    Subgroup2   = recode_leftblock_subgroup(Group, Subgroup),
    ExcludeLeft = leftblock_exclude_flag(Group, Subgroup),
    Consistent  = case_when(
      toupper(as.character(Consistent)) == "C" ~ "Consistent",
      toupper(as.character(Consistent)) == "NC" ~ "Not consistent",
      as.character(Consistent) %in% c("-", "No change", "no change", "Neutral", "0") ~ "Unsure",
      TRUE ~ as.character(Consistent)
    ),
    ZooKind = zoo_kind(Group)
  ) %>%
  filter(!is.na(Subgroup2), !ExcludeLeft) %>%
  filter(Consistent %in% c("Consistent", "Not consistent", "Unsure"))

phen_zoo_sub_cons <- phen_raw_cons %>%
  filter(grepl("^Zooplankton", Group)) %>%
  transmute(
    Taxon = canonical_taxon(Subgroup2, "phenology"),
    Consistent,
    Block = "Zooplankton",
    Study = "This study"
  )

phen_zoo_coll_cons <- phen_zoo_sub_cons %>%
  mutate(Taxon = "Zooplankton")

phen_zoo_kind_cons <- phen_raw_cons %>%
  filter(!is.na(ZooKind)) %>%
  transmute(
    Taxon = ZooKind,
    Consistent,
    Block = "Zooplankton",
    Study = "This study"
  )

# Cooley phenology consistency: zooplankton only
ipcc_raw_cons <- read_xlsx(fp_gen, sheet = "Phenology") %>%
  filter(is.na(Include) | Include != "N") %>%
  filter(Database != "Kris") %>%
  mutate(
    Group = str_squish(as.character(Group)),
    Consistent = toupper(str_squish(as.character(Consistent)))
  ) %>%
  filter(!is.na(Group), Group != "NA")

ipcc_prop_z <- ipcc_raw_cons %>%
  mutate(
    ConsBucket = case_when(
      Consistent %in% c("Y", "YES", "C") ~ "Consistent",
      Consistent %in% c("N", "NC", "NO") ~ "Not consistent",
      Consistent %in% c("", "NA") ~ NA_character_,
      TRUE ~ "Unsure"
    ),
    Taxon = case_when(
      Group %in% c("Zooplankton", "Meroplankton") ~ "Zooplankton",
      str_detect(Group, regex("Zooplankton|Mero", ignore_case = TRUE)) ~ "Zooplankton",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(ConsBucket), Taxon == "Zooplankton") %>%
  transmute(
    Taxon,
    Consistent = ConsBucket,
    Block = "General marine taxa",
    Study = "Cooley et al. (2022)"
  )

phen_prop_raw <- bind_rows(
  phen_zoo_sub_cons,
  phen_zoo_coll_cons,
  phen_zoo_kind_cons,
  ipcc_prop_z
) %>%
  mutate(
    Metric = "Phenology",
    Consistent = factor(Consistent, levels = c("Consistent", "Unsure", "Not consistent"))
  )

# This study distribution consistency
dist_raw_cons <- read_xlsx(fp_zoo, sheet = "Distribution") %>%
  mutate(
    Group       = str_squish(as.character(Group)),
    Subgroup    = na_if(str_squish(as.character(Subgroup)), "NA"),
    Subgroup2   = recode_leftblock_subgroup(Group, Subgroup),
    ExcludeLeft = leftblock_exclude_flag(Group, Subgroup),
    Consistent  = if ("Consistent" %in% names(.)) norm_consistency(Consistent) else NA_character_,
    ZooKind     = zoo_kind(Group)
  ) %>%
  { if ("Include" %in% names(.)) filter(., is.na(Include) | Include != "No") else . } %>%
  filter(!is.na(Subgroup2), !ExcludeLeft) %>%
  filter(Consistent %in% c("Consistent", "Not consistent", "Unsure"))

dist_zoo_sub_cons <- dist_raw_cons %>%
  filter(grepl("^Zooplankton", Group)) %>%
  transmute(
    Taxon = canonical_taxon(Subgroup2, "distribution"),
    Consistent,
    Block = "Zooplankton",
    Study = "This study"
  )

dist_zoo_coll_cons <- dist_zoo_sub_cons %>%
  mutate(Taxon = "Zooplankton")

dist_zoo_kind_cons <- dist_raw_cons %>%
  filter(!is.na(ZooKind)) %>%
  transmute(
    Taxon = ZooKind,
    Consistent,
    Block = "Zooplankton",
    Study = "This study"
  )

# Poloczanska distribution consistency: zooplankton only
gen_raw_cons <- read_xlsx(fp_gen, sheet = "Distribution") %>%
  filter(Database != "Kris")

cand_cols <- c("Consistent", "Consistency", "Consistent_with_CC", "Consistence", "Consistent?")
cons_cols_found <- cand_cols[cand_cols %in% names(gen_raw_cons)]
cons_col <- if (length(cons_cols_found) > 0) cons_cols_found[1] else NA_character_

gen_prop_z <- gen_raw_cons %>%
  mutate(
    Taxa = str_squish(as.character(Taxa)),
    Taxon = case_when(
      Taxa %in% c("Zooplankton", "Zooplankton - holo", "Zooplankton - mero") ~ "Zooplankton",
      str_detect(Taxa, regex("Zooplankton", ignore_case = TRUE)) ~ "Zooplankton",
      TRUE ~ NA_character_
    ),
    Consistent = if (!is.na(cons_col)) norm_consistency(.data[[cons_col]]) else NA_character_
  ) %>%
  filter(!is.na(Taxon), Taxon == "Zooplankton") %>%
  filter(Consistent %in% c("Consistent", "Not consistent", "Unsure")) %>%
  transmute(
    Taxon,
    Consistent,
    Block = "General marine taxa",
    Study = "Poloczanska et al. (2016)"
  )

dist_prop_raw <- bind_rows(
  dist_zoo_sub_cons,
  dist_zoo_coll_cons,
  dist_zoo_kind_cons,
  gen_prop_z
) %>%
  mutate(
    Metric = "Distribution",
    Consistent = factor(Consistent, levels = c("Consistent", "Unsure", "Not consistent"))
  )

prop_raw_all <- bind_rows(
  phen_prop_raw,
  dist_prop_raw
) %>%
  mutate(
    Xkey = factor(paste(Block, Taxon, sep = sep), levels = x_levels_master)
  ) %>%
  filter(!is.na(Xkey))

# ===========================================
# Percentage-consistent summary table
# ===========================================

make_percentage_summary <- function(raw_df, metric_name) {
  counted <- raw_df %>%
    filter(Metric == metric_name) %>%
    count(Metric, Study, Block, Taxon, Xkey, Consistent, name = "n_cat") %>%
    group_by(Metric, Study, Block, Taxon, Xkey) %>%
    summarise(
      n = sum(n_cat),
      k_consistent = sum(n_cat[Consistent == "Consistent"], na.rm = TRUE),
      k_unsure = sum(n_cat[Consistent == "Unsure"], na.rm = TRUE),
      k_not_consistent = sum(n_cat[Consistent == "Not consistent"], na.rm = TRUE),
      .groups = "drop"
    )
  
  counted %>%
    rowwise() %>%
    mutate(
      wilson = list(if (n > 0) wilson_ci(k_consistent, n) else tibble(
        proportion = NA_real_,
        percentage = NA_real_,
        ci_low = NA_real_,
        ci_high = NA_real_
      ))
    ) %>%
    unnest(wilson) %>%
    ungroup() %>%
    right_join(
      column_lookup,
      by = c("Xkey", "Block", "Taxon")
    ) %>%
    mutate(
      Panel = "Percentage consistent",
      Metric = metric_name,
      Study = case_when(
        Metric == "Distribution" & Block %in% c("Zooplankton", "Phytoplankton") ~ "This study",
        Metric == "Distribution" & Block == "General marine taxa" ~ "Poloczanska et al. (2016)",
        Metric == "Phenology" & Block %in% c("Zooplankton", "Phytoplankton") ~ "This study",
        Metric == "Phenology" & Block == "General marine taxa" ~ "Cooley et al. (2022)",
        TRUE ~ Study
      ),
      n = replace_na(n, 0L),
      k_consistent = replace_na(k_consistent, 0L),
      k_unsure = replace_na(k_unsure, 0L),
      k_not_consistent = replace_na(k_not_consistent, 0L),
      value_type = paste0(Metric, " percentage consistent"),
      mean_value = percentage,
      unit = "%",
      ci_method = "Wilson 95% CI for binomial proportion",
      ci_shown = n > 0,
      denominator_note = "consistent / (consistent + not consistent + unsure) x 100"
    ) %>%
    select(
      Panel, Metric, Study, x_order, Block, Taxon,
      n, k_consistent, k_not_consistent, k_unsure,
      value_type, mean_value, ci_low, ci_high, ci_shown,
      unit, ci_method, denominator_note
    ) %>%
    arrange(x_order)
}

percentage_distribution_summary <- make_percentage_summary(prop_raw_all, "Distribution")
percentage_phenology_summary <- make_percentage_summary(prop_raw_all, "Phenology")

percentage_summary <- bind_rows(
  percentage_distribution_summary,
  percentage_phenology_summary
) %>%
  arrange(x_order, Metric)

# ===========================================
# Master summary tables
# ===========================================

distribution_summary_for_master <- distribution_summary %>%
  mutate(Metric = "Distribution") %>%
  select(
    Panel, Metric, Study, x_order, Block, Taxon, n,
    value_type, mean_value, ci_low, ci_high, ci_shown,
    unit, ci_method, denominator_note,
    sd, se, median_value, min_value, max_value,
    ci_low_raw, ci_high_raw
  )

phenology_summary_for_master <- phenology_summary %>%
  mutate(Metric = "Phenology") %>%
  select(
    Panel, Metric, Study, x_order, Block, Taxon, n,
    value_type, mean_value, ci_low, ci_high, ci_shown,
    unit, ci_method, denominator_note,
    sd, se, median_value, min_value, max_value,
    ci_low_raw, ci_high_raw
  )

percentage_summary_for_master <- percentage_summary %>%
  mutate(
    sd = NA_real_,
    se = NA_real_,
    median_value = NA_real_,
    min_value = NA_real_,
    max_value = NA_real_,
    ci_low_raw = ci_low,
    ci_high_raw = ci_high
  ) %>%
  select(
    Panel, Metric, Study, x_order, Block, Taxon, n,
    value_type, mean_value, ci_low, ci_high, ci_shown,
    unit, ci_method, denominator_note,
    sd, se, median_value, min_value, max_value,
    ci_low_raw, ci_high_raw
  )

master_all_columns <- bind_rows(
  distribution_summary_for_master,
  phenology_summary_for_master,
  percentage_summary_for_master
) %>%
  arrange(Panel, Metric, x_order)

master_plotted_only <- master_all_columns %>%
  filter(n > 0)

# ===========================================
# Short summary table
# ===========================================

supervisor_summary <- master_plotted_only %>%
  mutate(
    value_ci = case_when(
      is.na(mean_value) ~ NA_character_,
      is.na(ci_low) | is.na(ci_high) ~ sprintf("%.2f", mean_value),
      TRUE ~ sprintf("%.2f (95%% CI: %.2f to %.2f)", mean_value, ci_low, ci_high)
    )
  ) %>%
  select(
    Panel, Metric, Study, Block, Taxon, n,
    value_ci, unit, ci_method, denominator_note
  ) %>%
  arrange(Panel, Metric, Study, Block, Taxon)

# ===========================================
# Export CSV files
# ===========================================

write.csv(
  master_all_columns,
  file.path(out_dir, "master_summary_all_columns.csv"),
  row.names = FALSE
)

write.csv(
  master_plotted_only,
  file.path(out_dir, "master_summary_plotted_only.csv"),
  row.names = FALSE
)

write.csv(
  distribution_summary,
  file.path(out_dir, "panel_A_distribution_summary.csv"),
  row.names = FALSE
)

write.csv(
  phenology_summary,
  file.path(out_dir, "panel_B_phenology_summary.csv"),
  row.names = FALSE
)

write.csv(
  percentage_summary,
  file.path(out_dir, "panel_C_percentage_consistent_summary.csv"),
  row.names = FALSE
)

write.csv(
  supervisor_summary,
  file.path(out_dir, "supervisor_summary_easy_reading.csv"),
  row.names = FALSE
)

# ===========================================
# Export Excel workbook
# ===========================================

wb <- createWorkbook()

addWorksheet(wb, "README")
writeData(
  wb,
  "README",
  data.frame(
    Notes = c(
      "This workbook summarises all values behind the three-panel observational figure.",
      "Panel A = Distribution change, using t-based 95% CI around the mean.",
      "Panel B = Phenology change, using t-based 95% CI around the mean.",
      "Panel C = Percentage consistent, using Wilson 95% CI for binomial proportions.",
      "For Panels A and B, CIs are only shown when n >= 3.",
      "For Panel C, denominator = consistent + not consistent + unsure.",
      "Poloczanska is labelled as Poloczanska et al. (2016), based on database review.",
      "master_summary_all_columns includes blank/non-plotted columns to match the full figure axis.",
      "master_summary_plotted_only includes only rows where n > 0."
    )
  )
)

addWorksheet(wb, "Master all columns")
writeData(wb, "Master all columns", master_all_columns)

addWorksheet(wb, "Master plotted only")
writeData(wb, "Master plotted only", master_plotted_only)

addWorksheet(wb, "Panel A Distribution")
writeData(wb, "Panel A Distribution", distribution_summary)

addWorksheet(wb, "Panel B Phenology")
writeData(wb, "Panel B Phenology", phenology_summary)

addWorksheet(wb, "Panel C Percent")
writeData(wb, "Panel C Percent", percentage_summary)

addWorksheet(wb, "Supervisor summary")
writeData(wb, "Supervisor summary", supervisor_summary)

sheets <- names(wb)

for (s in sheets) {
  setColWidths(wb, s, cols = 1:50, widths = "auto")
  freezePane(wb, s, firstRow = TRUE)
}

saveWorkbook(
  wb,
  file.path(out_dir, "observational_figure_summary_tables.xlsx"),
  overwrite = TRUE
)

# ===========================================
# Export PDF table
# ===========================================

round_numeric_for_print <- function(df) {
  df %>%
    mutate(
      across(
        where(is.numeric),
        ~ ifelse(is.na(.x), NA, round(.x, 2))
      )
    )
}

make_print_table <- function(df, title, max_rows_per_page = 28) {
  df_print <- df %>%
    round_numeric_for_print() %>%
    mutate(across(everything(), as.character)) %>%
    mutate(across(everything(), ~ replace_na(.x, "")))
  
  if (nrow(df_print) == 0) {
    df_print <- data.frame(Message = "No rows to display")
  }
  
  chunks <- split(df_print, ceiling(seq_len(nrow(df_print)) / max_rows_per_page))
  
  grobs <- lapply(seq_along(chunks), function(i) {
    chunk <- chunks[[i]]
    
    table_grob <- tableGrob(
      chunk,
      rows = NULL,
      theme = ttheme_minimal(
        base_size = 7,
        core = list(fg_params = list(hjust = 0, x = 0.02)),
        colhead = list(fg_params = list(fontface = "bold", hjust = 0, x = 0.02))
      )
    )
    
    title_grob <- textGrob(
      paste0(title, " — page ", i, " of ", length(chunks)),
      gp = gpar(fontsize = 14, fontface = "bold"),
      x = 0,
      hjust = 0
    )
    
    arrangeGrob(
      title_grob,
      table_grob,
      ncol = 1,
      heights = c(0.08, 0.92)
    )
  })
  
  grobs
}

# Selected columns for PDF
distribution_print <- distribution_summary %>%
  select(Study, x_order, Block, Taxon, n, mean_value, ci_low, ci_high, ci_shown, unit)

phenology_print <- phenology_summary %>%
  select(Study, x_order, Block, Taxon, n, mean_value, ci_low, ci_high, ci_shown, unit)

percentage_print <- percentage_summary %>%
  select(Metric, Study, x_order, Block, Taxon, n,
         k_consistent, k_not_consistent, k_unsure,
         mean_value, ci_low, ci_high, unit)

supervisor_print <- supervisor_summary

pdf_file <- file.path(out_dir, "observational_figure_summary_tables_PRINTED.pdf")

pdf(pdf_file, width = 16.5, height = 11.7)

all_grobs <- c(
  make_print_table(supervisor_print, "Supervisor summary", max_rows_per_page = 24),
  make_print_table(distribution_print, "Panel A — Distribution change", max_rows_per_page = 30),
  make_print_table(phenology_print, "Panel B — Phenology change", max_rows_per_page = 30),
  make_print_table(percentage_print, "Panel C — Percentage consistent", max_rows_per_page = 30)
)

for (g in all_grobs) {
  grid.newpage()
  grid.draw(g)
}

dev.off()

# ===========================================
# Print output locations
# ===========================================

cat("\n===========================================\n")
cat("SUMMARY TABLE EXPORT COMPLETE\n")
cat("===========================================\n")

cat("\nOutput folder:\n")
cat(out_dir, "\n")

cat("\nExcel workbook:\n")
cat(file.path(out_dir, "observational_figure_summary_tables.xlsx"), "\n")

cat("\nPrinted PDF table:\n")
cat(pdf_file, "\n")

cat("\nMain CSVs:\n")
cat(file.path(out_dir, "master_summary_all_columns.csv"), "\n")
cat(file.path(out_dir, "master_summary_plotted_only.csv"), "\n")
cat(file.path(out_dir, "supervisor_summary_easy_reading.csv"), "\n")

cat("\nQuick check: Poloczanska et al. (2016) zooplankton distribution summary:\n")
print(
  distribution_summary %>%
    filter(Study == "Poloczanska et al. (2016)", Taxon == "Zooplankton") %>%
    select(Study, Block, Taxon, n, mean_value, ci_low, ci_high, unit, ci_method)
)

cat("\nQuick check: This study zooplankton total distribution summary:\n")
print(
  distribution_summary %>%
    filter(Study == "This study", Block == "Phytoplankton", Taxon == "Zooplankton") %>%
    select(Study, Block, Taxon, n, mean_value, ci_low, ci_high, unit, ci_method)
)