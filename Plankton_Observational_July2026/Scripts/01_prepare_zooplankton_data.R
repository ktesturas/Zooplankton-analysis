# ===========================================
# 01_prepare_zooplankton_data.R
# Prepares all data objects used by Script 2 and Script 3.
# ===========================================

# Packages
pacman::p_load(dplyr, readxl, tidyr, here)

# Project folder system
dir_data <- here::here("Data")
dir_outputs <- here::here("Outputs")
dir_docs <- here::here("Docs")
dir_scripts <- here::here("Scripts")

# Input workbook
fp_zoo <- file.path(
  dir_data, "Zooplankton_observational_08September2026.xlsx")

# Shared data constants
PREV_BLOCK    <- "Previous work"
PREV_DATABASE <- "Poloczanska et al. 2013"
sep           <- " || "

ZOOP_GROUPS <- c("Zooplankton - holo", "Zooplankton - mero")
EXCLUDED_SUBGROUPS <- c("Dinoflagellates", "Rotifers")

# Read the worksheets
dist_source <- read_xlsx(
  fp_zoo,
  sheet = "Distribution",
  na = "NA")

phen_source <- read_xlsx(
  fp_zoo,
  sheet = "Phenology",
  na = "NA")

# Columns required by the analysis
required_distribution <- c(
  "Database", "Group", "Subgroup", "Duration",
  "Rate", "Consistent", "Include")

required_phenology <- c(
  "Database", "Group", "Subgroup", "Duration",
  "Rate", "Consistent", "Include")

# Filtering the dataset
dist_source <- dist_source %>%
  mutate(
    Include = as.character(Include),
    Duration = as.numeric(Duration)
  ) %>%
  filter(Include == "Yes", Duration >= 19)

phen_source <- phen_source %>%
  mutate(
    Include = as.character(Include),
    Duration = as.numeric(Duration)
  ) %>%
  filter(Include == "Yes", Duration >= 19)

# ===========================================
# Helper functions
# ===========================================
is_previous_work <- function(x) {
  x == PREV_DATABASE
}

# Organizing the taxa
group_subgroup <- function(x) {
  case_when(
    x %in% c("Amphipods", "Cirripedes", "Mysids") ~ "Other crustaceans",
    x == "Cladoceran" ~ "Cladocerans",
    x == "Pteropods" ~ "Molluscs",
    x == "Hydromedusae" ~ "Cnidarians",
    x %in% c("Phoronids", "Platyhelminthes", "Other zooplankton") ~ "Other metazoans",
    TRUE ~ x
    )
  }

zoo_kind <- function(x) {
  case_when(
    x == "Zooplankton - holo" ~ "Holoplankton",
    x == "Zooplankton - mero" ~ "Meroplankton",
    TRUE ~ NA_character_
  )
}

classify_consistency <- function(x) {
  case_when(
    x == "C" ~ "Consistent",
    x == "NC" ~ "Not consistent",
    x == "NA" ~ NA_character_,
    TRUE ~ NA_character_
  )
}

lab_x <- function(x) sub("^[^|]+\\s\\|\\|\\s", "", as.character(x))
pretty_taxon_prop <- lab_x

# ===========================================
# Distribution rate data
# ===========================================
# Clean distribution data
dist_raw0 <- dist_source %>%
  mutate(
    Group       = as.character(Group),
    Subgroup    = as.character(Subgroup),
    Subgroup2   = group_subgroup(Subgroup),
    ExcludeLeft = Subgroup %in% EXCLUDED_SUBGROUPS,
    Rate        = as.numeric(Rate),
    ZooKind     = zoo_kind(Group)
  ) %>%
  filter(!is.na(Rate), !is.na(Subgroup))

# Data by zooplankton taxonomic grouops
dist_zoo <- dist_raw0 %>%
  filter(Group %in% ZOOP_GROUPS, !ExcludeLeft) %>%
  transmute(
    Taxon = Subgroup2,
    Rate,
    Block = "Zooplankton"
  )

# Data by mero/holo
dist_zoo_kind <- dist_raw0 %>%
  filter(Group %in% ZOOP_GROUPS, !ExcludeLeft, !is.na(ZooKind)) %>%
  transmute(
    Taxon = ZooKind,
    Rate,
    Block = "Zooplankton"
  )

# Data by overall zoo
dist_pool <- dist_raw0 %>%
  filter(Group %in% ZOOP_GROUPS, !ExcludeLeft) %>%
  transmute(
    Taxon = "Zooplankton",
    Rate,
    Block = "Phytoplankton"
  )

# Select previous work observations
gen_raw_dist <- dist_source %>%
  mutate(
    Group = as.character(Group),
    Rate = as.numeric(Rate)
  ) %>%
  filter(is_previous_work(Database), !is.na(Group), !is.na(Rate))

# Data by previous work
dist_gen <- gen_raw_dist %>%
  mutate(
    Taxon = if_else(Group %in% ZOOP_GROUPS, "Zooplankton", NA_character_)
  ) %>%
  filter(Taxon == "Zooplankton") %>%
  transmute(
    Taxon,
    Rate,
    Block = PREV_BLOCK
  )

# Distribution data columns: taxa, holo/meron, overall, previous
dist_all_raw <- bind_rows(dist_zoo, dist_zoo_kind, dist_pool, dist_gen) %>%
  mutate(
    sign = ifelse(Rate >= 0, "Positive", "Negative"),
    Panel = "Distribution"
  )

# ===========================================
# Phenology rate data
# ===========================================
# Clean phenology data
zp_raw <- phen_source %>%
  mutate(
    Rate        = as.numeric(Rate),
    Group       = as.character(Group),
    Subgroup    = as.character(Subgroup),
    Subgroup2   = group_subgroup(Subgroup),
    ExcludeLeft = Subgroup %in% EXCLUDED_SUBGROUPS,
    ZooKind     = zoo_kind(Group)
  ) %>%
  filter(!is.na(Rate), !is.na(Subgroup))

# Data by zooplankton taxonomic grouops
phen_zoo <- zp_raw %>%
  filter(Group %in% ZOOP_GROUPS, !ExcludeLeft) %>%
  transmute(
    Taxon = Subgroup2,
    Rate,
    Block = "Zooplankton"
  )

# Data by holo/mero
phen_zoo_kind <- zp_raw %>%
  filter(Group %in% ZOOP_GROUPS, !ExcludeLeft, !is.na(ZooKind)) %>%
  transmute(
    Taxon = ZooKind,
    Rate,
    Block = "Zooplankton"
  )

# Data by holo/mero
phen_pool <- zp_raw %>%
  filter(Group %in% ZOOP_GROUPS, !ExcludeLeft) %>%
  transmute(
    Taxon = "Zooplankton",
    Rate,
    Block = "Phytoplankton"
  )

# Data by overall zoo
ipcc_raw <- phen_source %>%
  mutate(
    Rate = as.numeric(Rate),
    Group = as.character(Group)
  ) %>%
  filter(is_previous_work(Database), !is.na(Rate), !is.na(Group))

# Data by previous work
phen_gen <- ipcc_raw %>%
  mutate(
    Taxon = if_else(Group %in% ZOOP_GROUPS, "Zooplankton", NA_character_)
  ) %>%
  filter(Taxon == "Zooplankton") %>%
  transmute(
    Taxon,
    Rate,
    Block = PREV_BLOCK
  )

# Phenology data columns: taxa, holo/meron, overall, previous
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
  filter(Block == "Phytoplankton") %>% # just to separate the overall zooplankton data for the plot
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

# Keep pooled zooplankton categories at the end
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
    mean_rate = mean(Rate),
    sd_rate   = sd(Rate),
    n         = n(),
    .groups   = "drop"
  ) %>%
  mutate(
    se = if_else(
      n > 1,
      sd_rate / sqrt(n),
      NA_real_
    ),
    
    # pmax prevents qt() from receiving df = 0
    tcrit_raw = qt(
      0.975,
      df = pmax(n - 1, 1)
    ),
    
    tcrit = if_else(
      n > 1,
      tcrit_raw,
      NA_real_
    ),
    
    ci_low = if_else(
      n >= 3,
      mean_rate - tcrit * se,
      NA_real_
    ),
    
    ci_high = if_else(
      n >= 3,
      mean_rate + tcrit * se,
      NA_real_
    )
  ) %>%
  select(-tcrit_raw)
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

phen_raw <- phen_source %>%
  mutate(
    Group       = as.character(Group),
    Subgroup    = as.character(Subgroup),
    Subgroup2   = group_subgroup(Subgroup),
    ExcludeLeft = Subgroup %in% EXCLUDED_SUBGROUPS,
    Consistent  = classify_consistency(as.character(Consistent)),
    ZooKind = zoo_kind(Group)
  ) %>%
  filter(
    !is.na(Subgroup2),
    !ExcludeLeft,
    Consistent %in% c("Consistent", "Not consistent")
  )

phen_zoo_sub <- phen_raw %>%
  filter(Group %in% ZOOP_GROUPS) %>%
  transmute(
    Taxon = Subgroup2,
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

ipcc_raw2 <- phen_source %>%
  mutate(
    Group = as.character(Group),
    Consistent = as.character(Consistent)
  ) %>%
  filter(is_previous_work(Database), !is.na(Group))

ipcc_prop_z <- ipcc_raw2 %>%
  mutate(
    ConsBucket = classify_consistency(Consistent),
    Taxon = if_else(Group %in% ZOOP_GROUPS, "Zooplankton", NA_character_)
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
      levels = c("Consistent", "Not consistent")
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

dist_raw_p <- dist_source %>%
  mutate(
    Group       = as.character(Group),
    Subgroup    = as.character(Subgroup),
    Subgroup2   = group_subgroup(Subgroup),
    ExcludeLeft = Subgroup %in% EXCLUDED_SUBGROUPS,
    Consistent  = classify_consistency(as.character(Consistent)),
    ZooKind     = zoo_kind(Group)
  ) %>%
  filter(
    !is.na(Subgroup2),
    !ExcludeLeft,
    Consistent %in% c("Consistent", "Not consistent")
  )

dist_zoo_sub <- dist_raw_p %>%
  filter(Group %in% ZOOP_GROUPS) %>%
  transmute(
    Taxon = Subgroup2,
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

gen_raw <- dist_source %>%
  filter(is_previous_work(Database))


gen_prop_z <- gen_raw %>%
  mutate(
    Group = as.character(Group),
    Taxon = if_else(Group %in% ZOOP_GROUPS, "Zooplankton", NA_character_),
    Consistent = classify_consistency(as.character(Consistent))
  ) %>%
  filter(
    !is.na(Taxon),
    Taxon == "Zooplankton",
    Consistent %in% c("Consistent", "Not consistent")
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
      levels = c("Consistent", "Not consistent")
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