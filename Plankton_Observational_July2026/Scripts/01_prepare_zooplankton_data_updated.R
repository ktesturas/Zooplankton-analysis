# ===========================================
# 01_prepare_zooplankton_data.R
# Prepares all data objects needed for the zooplankton observational plot.
# ===========================================

# ---- Packages ----
library(dplyr)
library(readxl)
library(stringr)
library(tidyr)
library(here)

# ---- Project folder system ----
dir_data    <- here::here("Data")
dir_outputs <- here::here("Outputs")
dir_docs    <- here::here("Docs")
dir_scripts <- here::here("Scripts")

# Script 1 does not export anything, but these checks keep the project
# folder structure standard for scripts that source this file later.
invisible(lapply(
  c(dir_data, dir_outputs, dir_docs, dir_scripts),
  dir.create,
  recursive = TRUE,
  showWarnings = FALSE
))

# ---- Input file paths ----

find_existing_file <- function(candidates, label) {
  hit <- candidates[file.exists(candidates)][1]

  if (length(hit) == 0 || is.na(hit)) {
    stop(
      paste0(
        "Missing ", label, ". Expected one of:\n",
        paste0("- ", candidates, collapse = "\n"),
        "\n\nMove/copy the needed Excel file into Data/ and rerun the script."
      ),
      call. = FALSE
    )
  }

  hit
}

fp_zoo <- find_existing_file(
  c(
    here::here("Data", "Observational", "Plankton observational.xlsx"),
    here::here("Data", "Plankton observational.xlsx")
  ),
  "zooplankton observational workbook"
)

fp_gen <- find_existing_file(
  c(
    here::here("Data", "All marine life", "All_marine_life_data.xlsx"),
    here::here("Data", "All_marine_life", "All_marine_life_data.xlsx"),
    here::here("Data", "All_marine_life_data.xlsx")
  ),
  "all-marine-life / previous-work workbook"
)

# ---- Shared data constants ----
PREV_BLOCK <- "Previous work"
sep <- " || "

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

lab_x <- function(x) {
  taxon <- sub("^[^|]+\\s\\|\\|\\s", "", as.character(x))
  pretty_taxon(taxon)
}

pretty_taxon_prop <- function(lbl) {
  taxon <- sub("^[^|]+\\s\\|\\|\\s", "", as.character(lbl))
  pretty_taxon(taxon)
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
  
  case_when(
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
  s <- str_replace_all(s, "^Appendicularians?$", "Appendicularians")
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

# ===========================================
# Distribution rate data
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
    Taxon = normalise_taxon(
      ifelse(str_to_lower(Subgroup2) == "other", "Other zooplankton", Subgroup2)
    ),
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
  filter(!is.na(Taxa), Taxa != "NA", !is.na(Rate), Database != "Kris")

dist_gen <- gen_raw_dist %>%
  mutate(
    Taxon = case_when(
      Taxa %in% c("Zooplankton", "Zooplankton - holo", "Zooplankton - mero") ~ "Zooplankton",
      str_detect(Taxa, regex("Zooplankton", ignore_case = TRUE)) ~ "Zooplankton",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(Taxon == "Zooplankton") %>%
  transmute(
    Taxon,
    Rate,
    Block = PREV_BLOCK
  )

dist_all_raw <- bind_rows(dist_zoo, dist_zoo_kind, dist_pool, dist_gen) %>%
  mutate(
    sign = ifelse(Rate >= 0, "Positive", "Negative"),
    Panel = "Distribution"
  )

# ===========================================
# Phenology rate data
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
  filter(is.na(Include) | Include != "N", Database != "Kris") %>%
  mutate(
    Rate = suppressWarnings(as.numeric(Shift)),
    Group = str_squish(as.character(Group))
  ) %>%
  filter(!is.na(Rate), !is.na(Group), Group != "NA")

phen_gen <- ipcc_raw %>%
  mutate(
    Taxon = case_when(
      Group %in% c("Zooplankton", "Meroplankton") ~ "Zooplankton",
      str_detect(Group, regex("Zooplankton|Mero", ignore_case = TRUE)) ~ "Zooplankton",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(Taxon == "Zooplankton") %>%
  transmute(
    Taxon,
    Rate,
    Block = PREV_BLOCK
  )

phen_all_raw <- bind_rows(phen_gen, phen_zoo, phen_zoo_kind, phen_pool) %>%
  mutate(
    sign = ifelse(Rate >= 0, "Positive", "Negative"),
    Panel = "Phenology"
  )

# ===========================================
# Master x-axis levels and shared layout objects
# ===========================================

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

append_block <- function(block_name, base_order, sum_other) {
  candidates <- sum_other %>%
    filter(Block == block_name) %>%
    pull(Taxon)
  
  only_other <- setdiff(candidates, base_order)
  
  if (!length(only_other)) {
    return(base_order)
  }
  
  add <- sum_other %>%
    filter(Block == block_name, Taxon %in% only_other) %>%
    arrange(mean_rate, Taxon) %>%
    pull(Taxon)
  
  c(base_order, add)
}

move_items_before <- function(vec, items, before_item = "Zooplankton") {
  items <- items[items %in% vec]
  
  if (!length(items)) {
    return(vec)
  }
  
  vec_wo <- vec[!vec %in% items]
  pos <- match(before_item, vec_wo)
  
  if (is.na(pos)) {
    c(vec_wo, items)
  } else {
    c(vec_wo[seq_len(pos - 1)], items, vec_wo[pos:length(vec_wo)])
  }
}

z_order <- append_block("Zooplankton", z_order_phen, sum_dist)
p_order <- append_block("Phytoplankton", p_order_phen, sum_dist)

z_order <- move_items_before(
  z_order,
  c("Holoplankton", "Meroplankton"),
  "Zooplankton"
)

z_levels <- paste("Zooplankton", z_order, sep = sep)
p_levels <- paste("Phytoplankton", p_order, sep = sep)
g_levels <- paste(PREV_BLOCK, "Zooplankton", sep = sep)

x_levels_master <- c(z_levels, p_levels, g_levels)

dist_all_raw <- dist_all_raw %>%
  mutate(
    Xkey = factor(paste(Block, Taxon, sep = sep), levels = x_levels_master)
  )

phen_all_raw <- phen_all_raw %>%
  mutate(
    Xkey = factor(paste(Block, Taxon, sep = sep), levels = x_levels_master)
  )

combined_raw <- bind_rows(dist_all_raw, phen_all_raw) %>%
  tidyr::complete(
    Panel,
    Xkey,
    fill = list(Rate = NA_real_, sign = NA_character_)
  ) %>%
  mutate(
    Block = factor(sub("\\s\\|\\|\\s.*$", "", as.character(Xkey))),
    Taxon = factor(sub("^[^|]+\\s\\|\\|\\s", "", as.character(Xkey)))
  )

idx_tbl <- tibble(
  Ykey  = factor(x_levels_master, levels = x_levels_master),
  Xnum  = seq_along(x_levels_master),
  Block = sub("\\s\\|\\|\\s.*$", "", as.character(Ykey))
)

is_block <- function(prefix) {
  grepl(paste0("^", prefix, "\\s\\|\\|\\s"), x_levels_master)
}

n_zoop  <- sum(is_block("Zooplankton"))
n_phyto <- sum(is_block("Phytoplankton"))

xmin_grey <- n_zoop + n_phyto + 0.5
xmax_grey <- length(x_levels_master) + 0.5

x_left_center  <- (0.5 + xmin_grey) / 2
x_right_center <- (xmin_grey + xmax_grey) / 2

vbreaks_x <- seq(0.5, length(x_levels_master) + 0.5, by = 1)

get_pos <- function(block, taxon) {
  match(paste(block, taxon, sep = sep), x_levels_master)
}

pos_holo     <- get_pos("Zooplankton", "Holoplankton")
pos_mero     <- get_pos("Zooplankton", "Meroplankton")
pos_zoop_col <- get_pos("Zooplankton", "Zooplankton")

sep_positions <- c(
  if (!is.na(pos_holo)) pos_holo - 0.5 else NA_real_,
  if (!is.na(pos_mero)) pos_mero - 0.5 else NA_real_,
  if (!is.na(pos_zoop_col)) pos_zoop_col - 0.5 else NA_real_,
  n_zoop + 0.5,
  n_zoop + n_phyto + 0.5
)

sep_positions <- sort(unique(sep_positions[is.finite(sep_positions)]))

x_n_right <- max(idx_tbl$Xnum, na.rm = TRUE) + 0.5

# ===========================================
# Rate summaries: t-based 95% CI
# ===========================================

sum_panel <- combined_raw %>%
  filter(is.finite(Rate)) %>%
  group_by(Panel, Xkey) %>%
  summarise(
    mean_rate = mean(Rate, na.rm = TRUE),
    sd_rate   = sd(Rate, na.rm = TRUE),
    n         = sum(!is.na(Rate)),
    .groups   = "drop"
  ) %>%
  mutate(
    se      = ifelse(n > 1, sd_rate / sqrt(n), NA_real_),
    tcrit   = ifelse(n > 1, qt(0.975, df = n - 1), NA_real_),
    ci_low  = ifelse(n >= 3, mean_rate - tcrit * se, NA_real_),
    ci_high = ifelse(n >= 3, mean_rate + tcrit * se, NA_real_)
  )

n_df <- sum_panel %>%
  select(Panel, Xkey, n)

zoop_means <- combined_raw %>%
  filter(is.finite(Rate), Block == "Zooplankton") %>%
  group_by(Panel) %>%
  summarise(
    zoop_mean = mean(Rate, na.rm = TRUE),
    .groups = "drop"
  )

# ===========================================
# Percentage-consistent data
# ===========================================

phen_raw <- read_xlsx(fp_zoo, sheet = "Phenology") %>%
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
  filter(
    !is.na(Subgroup2),
    !ExcludeLeft,
    Consistent %in% c("Consistent", "Not consistent", "Unsure")
  )

phen_zoo_sub <- phen_raw %>%
  filter(grepl("^Zooplankton", Group)) %>%
  transmute(
    Taxon = canonical_taxon(Subgroup2, "phenology"),
    Consistent,
    Block = "Zooplankton"
  )

phen_zoo_coll <- phen_zoo_sub %>%
  mutate(Taxon = "Zooplankton")

phen_zoo_kind_prop <- phen_raw %>%
  filter(!is.na(ZooKind)) %>%
  transmute(
    Taxon = ZooKind,
    Consistent,
    Block = "Zooplankton"
  )

ipcc_raw2 <- read_xlsx(fp_gen, sheet = "Phenology") %>%
  filter(is.na(Include) | Include != "N", Database != "Kris") %>%
  mutate(
    Group = str_squish(as.character(Group)),
    Consistent = toupper(str_squish(as.character(Consistent)))
  ) %>%
  filter(!is.na(Group), Group != "NA")

ipcc_prop_z <- ipcc_raw2 %>%
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
    Block = PREV_BLOCK
  )

phen_prop_raw <- bind_rows(
  phen_zoo_sub,
  phen_zoo_coll,
  phen_zoo_kind_prop,
  ipcc_prop_z
) %>%
  mutate(
    Consistent = factor(
      Consistent,
      levels = c("Consistent", "Unsure", "Not consistent")
    )
  )

phen_prop <- phen_prop_raw %>%
  count(Block, Taxon, Consistent, name = "n_cat") %>%
  group_by(Block, Taxon) %>%
  complete(Consistent, fill = list(n_cat = 0)) %>%
  mutate(
    n = sum(n_cat),
    prop = if_else(n > 0, n_cat / n, 0)
  ) %>%
  ungroup()

dist_raw_p <- read_xlsx(fp_zoo, sheet = "Distribution") %>%
  mutate(
    Group       = str_squish(as.character(Group)),
    Subgroup    = na_if(str_squish(as.character(Subgroup)), "NA"),
    Subgroup2   = recode_leftblock_subgroup(Group, Subgroup),
    ExcludeLeft = leftblock_exclude_flag(Group, Subgroup),
    Consistent  = if ("Consistent" %in% names(.)) norm_consistency(Consistent) else NA_character_,
    ZooKind     = zoo_kind(Group)
  ) %>%
  { if ("Include" %in% names(.)) filter(., is.na(Include) | Include != "No") else . } %>%
  filter(
    !is.na(Subgroup2),
    !ExcludeLeft,
    Consistent %in% c("Consistent", "Not consistent", "Unsure")
  )

dist_zoo_sub <- dist_raw_p %>%
  filter(grepl("^Zooplankton", Group)) %>%
  transmute(
    Taxon = canonical_taxon(Subgroup2, "distribution"),
    Consistent,
    Block = "Zooplankton"
  )

dist_zoo_coll <- dist_zoo_sub %>%
  mutate(Taxon = "Zooplankton")

dist_zoo_kind_prop <- dist_raw_p %>%
  filter(!is.na(ZooKind)) %>%
  transmute(
    Taxon = ZooKind,
    Consistent,
    Block = "Zooplankton"
  )

gen_raw <- read_xlsx(fp_gen, sheet = "Distribution") %>%
  filter(Database != "Kris")

cons_col_candidates <- c(
  "Consistent",
  "Consistency",
  "Consistent_with_CC",
  "Consistence",
  "Consistent?"
)

cons_col <- cons_col_candidates[cons_col_candidates %in% names(gen_raw)][1]

if (length(cons_col) == 0 || is.na(cons_col)) {
  cons_col <- NA_character_
}

gen_prop_z <- gen_raw %>%
  mutate(
    Taxa = str_squish(as.character(Taxa)),
    Taxon = case_when(
      Taxa %in% c("Zooplankton", "Zooplankton - holo", "Zooplankton - mero") ~ "Zooplankton",
      str_detect(Taxa, regex("Zooplankton", ignore_case = TRUE)) ~ "Zooplankton",
      TRUE ~ NA_character_
    ),
    Consistent = if (!is.na(cons_col)) {
      norm_consistency(.data[[cons_col]])
    } else {
      NA_character_
    }
  ) %>%
  filter(
    !is.na(Taxon),
    Taxon == "Zooplankton",
    Consistent %in% c("Consistent", "Not consistent", "Unsure")
  ) %>%
  transmute(
    Taxon,
    Consistent,
    Block = PREV_BLOCK
  )

dist_prop_raw <- bind_rows(
  dist_zoo_sub,
  dist_zoo_coll,
  dist_zoo_kind_prop,
  gen_prop_z
) %>%
  mutate(
    Consistent = factor(
      Consistent,
      levels = c("Consistent", "Unsure", "Not consistent")
    )
  )

dist_prop <- dist_prop_raw %>%
  count(Block, Taxon, Consistent, name = "n_cat") %>%
  group_by(Block, Taxon) %>%
  complete(Consistent, fill = list(n_cat = 0)) %>%
  mutate(
    n = sum(n_cat),
    prop = if_else(n > 0, n_cat / n, 0)
  ) %>%
  ungroup()

# ===========================================
# Percentage-consistent summaries: Wilson 95% binomial CI
# ===========================================

x_levels_prop <- x_levels_master
idx_tbl_prop  <- idx_tbl

phen_prop_all <- phen_prop %>%
  mutate(
    Ykey = factor(paste(Block, Taxon, sep = sep), levels = x_levels_prop)
  )

dist_prop_all <- dist_prop %>%
  mutate(
    Ykey = factor(paste(Block, Taxon, sep = sep), levels = x_levels_prop)
  )

phen_prop_plot <- phen_prop_all %>%
  filter(!is.na(Ykey))

dist_prop_plot <- dist_prop_all %>%
  filter(!is.na(Ykey))

summarise_wilson_prop <- function(df) {
  df %>%
    filter(Consistent == "Consistent", n > 0) %>%
    transmute(
      Block,
      Taxon,
      Ykey,
      n,
      k = n_cat,
      phat = prop
    ) %>%
    mutate(
      z      = qnorm(0.975),
      denom  = 1 + (z^2 / n),
      centre = (phat + (z^2 / (2 * n))) / denom,
      half   = (z * sqrt((phat * (1 - phat) / n) + (z^2 / (4 * n^2)))) / denom,
      p      = phat * 100,
      lo     = pmax(0, (centre - half) * 100),
      hi     = pmin(100, (centre + half) * 100)
    ) %>%
    select(Block, Taxon, Ykey, n, p, lo, hi)
}

se_phen_all <- summarise_wilson_prop(phen_prop_all) %>%
  mutate(Metric = "Phenology")

se_dist_all <- summarise_wilson_prop(dist_prop_all) %>%
  mutate(Metric = "Distribution")

se_phen <- se_phen_all %>%
  filter(!is.na(Ykey))

se_dist <- se_dist_all %>%
  filter(!is.na(Ykey))

get_p_zoop <- function(se_df) {
  out <- se_df %>%
    filter(Block == "Zooplankton", Taxon == "Zooplankton") %>%
    summarise(p = dplyr::first(p)) %>%
    pull(p)
  
  if (length(out) == 0) {
    NA_real_
  } else {
    out
  }
}

p_zoop_phen <- get_p_zoop(se_phen_all)
p_zoop_dist <- get_p_zoop(se_dist_all)

n_phen <- phen_prop_plot %>%
  distinct(Ykey, n) %>%
  rename(n_phen = n)

n_dist <- dist_prop_plot %>%
  distinct(Ykey, n) %>%
  rename(n_dist = n)

offset <- 0.13

se_dist2 <- se_dist %>%
  left_join(idx_tbl_prop, by = "Ykey") %>%
  mutate(xplot = Xnum - offset)

se_phen2 <- se_phen %>%
  left_join(idx_tbl_prop, by = "Ykey") %>%
  mutate(xplot = Xnum + offset)

n_phen2 <- n_phen %>%
  left_join(idx_tbl_prop, by = "Ykey")

n_dist2 <- n_dist %>%
  left_join(idx_tbl_prop, by = "Ykey")

# Echo collapsed Zooplankton stats into the Phytoplankton || Zooplankton placeholder column
phyto_zoop_key <- paste("Phytoplankton", "Zooplankton", sep = sep)

phyto_anchor <- idx_tbl_prop %>%
  filter(as.character(Ykey) == phyto_zoop_key) %>%
  pull(Xnum) %>%
  { if (length(.) == 0) NA_real_ else .[1] }

zoop_overall_phen <- se_phen_all %>%
  filter(Block == "Zooplankton", Taxon == "Zooplankton")

zoop_overall_dist <- se_dist_all %>%
  filter(Block == "Zooplankton", Taxon == "Zooplankton")

se_dist_phytoEcho <- if (nrow(zoop_overall_dist) && is.finite(phyto_anchor)) {
  zoop_overall_dist %>%
    mutate(xplot = phyto_anchor - offset)
} else {
  NULL
}

se_phen_phytoEcho <- if (nrow(zoop_overall_phen) && is.finite(phyto_anchor)) {
  zoop_overall_phen %>%
    mutate(xplot = phyto_anchor + offset)
} else {
  NULL
}

n_zoop_phen_val <- phen_prop_all %>%
  filter(Block == "Zooplankton", Taxon == "Zooplankton") %>%
  distinct(n) %>%
  pull(n) %>%
  { if (length(.) == 0) NA_integer_ else .[1] }

n_zoop_dist_val <- dist_prop_all %>%
  filter(Block == "Zooplankton", Taxon == "Zooplankton") %>%
  distinct(n) %>%
  pull(n) %>%
  { if (length(.) == 0) NA_integer_ else .[1] }

n_phen_phytoEcho <- if (is.finite(phyto_anchor) && !is.na(n_zoop_phen_val)) {
  tibble(Xnum = phyto_anchor, n_phen = n_zoop_phen_val)
} else {
  NULL
}

n_dist_phytoEcho <- if (is.finite(phyto_anchor) && !is.na(n_zoop_dist_val)) {
  tibble(Xnum = phyto_anchor, n_dist = n_zoop_dist_val)
} else {
  NULL
}

message("Data preparation complete. Plot design/export should happen in Script 2.")

