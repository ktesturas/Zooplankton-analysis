# ===========================================
# 02_plot_and_export_zooplankton.R
# Designs and exports the three-panel zooplankton observational plot.
#
# Exports TWO versions:
#   1. With grid lines
#   2. Without grid lines
#
# Background sections:
#   - Individual zooplankton groups = white
#   - Pooled groups                = light grey
#   - Previous work                = darker grey
# ===========================================

# Load Packages
pacman::p_load(dplyr, ggplot2, patchwork, scales, showtext, here, grid)


# Project paths
script_name <- "02_plot_and_export_zooplankton"
scripts_dir <- here::here("Scripts")
figures_dir <- here::here("Outputs", script_name)
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

# Source prepared data from Script 1
data_script <- here::here("Scripts","01_prepare_zooplankton_data.R")
source(data_script)

# ===========================================
# Font setup
# ===========================================

#font_ok <- TRUE

#tryCatch(
#  {
#    font_add_google("Inter", "helv")
#  },
#  error = function(e) {
#    font_ok <<- FALSE
#    message(
#      "Could not load Google font Inter. Using default sans font instead."
#    )
#  }
#)

#showtext_auto()

#fam <- if (font_ok) "helv" else "sans"

sysfonts::font_add(
  family = "helv",
  regular = "/System/Library/Fonts/Helvetica.ttc"
)

sysfonts::font_add(
  family = "header_bold",
  regular = "/System/Library/Fonts/Supplemental/Arial Bold.ttf"
)

showtext::showtext_auto()

fam <- "helv"
# ===========================================
# Plot styling constants
# ===========================================

# ---- Text size controls ----

BASE_SIZE <- 10

Y_TICK_TEXT_SIZE  <- 24
X_TICK_TEXT_SIZE  <- 24

Y_AXIS_TITLE_SIZE <- 28

HEADER_TEXT_SIZE  <- 10
N_LABEL_SIZE      <- 7
BOTTOM_LABEL_SIZE <- 7


# ---- Header label controls ----

HEADER_Y <- 0.63

PREV_LABEL_TOP    <- "Previous"
PREV_LABEL_BOTTOM <- "work"

PREV_LABEL_TOP_Y    <- 0.84
PREV_LABEL_BOTTOM_Y <- 0.48


# ---- Point and line styling ----

MEAN_PT_SIZE <- 1.7
RAW_PT_SIZE  <- 0.8
RAW_PT_ALPHA <- 0.35
RAW_STROKE   <- 0.25

ERR_LW          <- 0.35
LWD_ZERO        <- 0.30
LWD_DASH        <- 0.30
LWD_BOT_MIDLINE <- 0.30

GRID_LW   <- 0.20
BORDER_LW <- 0.20

MINOR_SEP_LW <- 0.25
MAJOR_SEP_LW <- 0.25


# ---- Layout controls ----

GAP_HEIGHT <- 0.15


# ===========================================
# Colours
# ===========================================

# Raw observations
col_map <- c(
  "Negative" = "#FDB462",
  "Positive" = "#009E73"
)

# Grid
grid_col <- "grey90"

# NEW:
# Background colours for the three sections
#
# Individual groups = no rectangle added = white
# Pooled groups     = light grey
# Previous work     = darker grey

pooled_fill   <- "#F2F2F2"
previous_fill <- "#D9D9D9"

# Other line colours
border_col <- "grey10"
zero_col   <- "grey35"

minor_sep_col <- "grey35"
major_sep_col <- "grey10"

# Mean colours
mean_col_dist <- "#2B6CB0"
mean_col_phen <- "black"

# ===========================================
# Axis limits and breaks
# ===========================================

dist_ymin <- -100
dist_ymax <- 200

phen_ymin <- -20
phen_ymax <- 10

brks_y_dist <- seq(-100, 200, 50)
brks_y_phen <- seq(-20, 10, 5)
brks_y_prop <- seq(0, 100, 20)


# ===========================================
# Separator positions
# ===========================================

major_sep_positions <- c(
  n_zoop + n_phyto + 0.5
)

minor_sep_positions <- setdiff(
  sep_positions,
  major_sep_positions
)


# ===========================================
# Background section positions
# ===========================================

# There are three pooled categories in this study:
#   Holoplankton
#   Meroplankton
#   Zooplankton
#
# The major separator marks the boundary between
# "This study" and "Previous work".

N_POOLED_GROUPS <- 3

previous_xmin <- major_sep_positions[1]

previous_xmax_master <- length(x_levels_master) + 0.5
previous_xmax_prop   <- length(x_levels_prop) + 0.5

pooled_xmax <- previous_xmin
pooled_xmin <- pooled_xmax - N_POOLED_GROUPS


# ===========================================
# One-time top header
# ===========================================

p_header <- ggplot() +
  
  annotate(
    "text",
    x = x_left_center,
    y = HEADER_Y,
    label = "This study",
    hjust = 0.5,
    vjust = 0.5,
    family = "header_bold",
    fontface = "plain",
    size = HEADER_TEXT_SIZE,
    colour = border_col
  ) +
  
  annotate(
    "text",
    x = x_right_center,
    y = PREV_LABEL_TOP_Y,
    label = PREV_LABEL_TOP,
    hjust = 0.5,
    vjust = 0.5,
    family = "header_bold",
    fontface = "plain",
    size = HEADER_TEXT_SIZE,
    colour = border_col
  ) +
  
  annotate(
    "text",
    x = x_right_center,
    y = PREV_LABEL_BOTTOM_Y,
    label = PREV_LABEL_BOTTOM,
    hjust = 0.5,
    vjust = 0.5,
    family = "header_bold",
    fontface = "plain",
    size = HEADER_TEXT_SIZE,
    colour = border_col
  ) +
  
  scale_x_continuous(
    limits = c(
      0.5,
      length(x_levels_master) + 0.5
    ),
    expand = expansion(mult = c(0, 0))
  ) +
  
  scale_y_continuous(
    limits = c(0, 1),
    expand = expansion(mult = c(0, 0))
  ) +
  
  coord_cartesian(clip = "off") +
  
  theme_void() +
  
  theme(
    plot.margin = margin(8, 14, -2, 25)
  )


# ===========================================
# Rate panel plotting function
# ===========================================

make_rate_panel <- function(
    panel_name,
    ymin,
    ymax,
    y_breaks,
    mean_col,
    y_lab,
    top_margin = 6,
    show_grid = TRUE
) {
  
  panel_mean <- zoop_means$zoop_mean[
    zoop_means$Panel == panel_name
  ]
  
  
  # -----------------------------------------
  # Raw observations
  # -----------------------------------------
  
  raw2 <- combined_raw %>%
    dplyr::filter(
      Panel == panel_name,
      is.finite(Rate)
    ) %>%
    dplyr::mutate(
      Rate_plot = pmin(
        pmax(Rate, ymin),
        ymax
      )
    )
  
  
  # -----------------------------------------
  # Summary data
  # -----------------------------------------
  
  sum2 <- sum_panel %>%
    dplyr::filter(
      Panel == panel_name
    ) %>%
    dplyr::mutate(
      mean_plot = pmin(
        pmax(mean_rate, ymin),
        ymax
      ),
      ci_low_plot = pmax(
        ci_low,
        ymin
      ),
      ci_high_plot = pmin(
        ci_high,
        ymax
      )
    )
  
  
  n2 <- n_df %>%
    dplyr::filter(
      Panel == panel_name
    )
  
  
  # -----------------------------------------
  # Plot
  # -----------------------------------------
  
  p <- ggplot() +
    
    
    # =======================================
  # BACKGROUND SECTION 1
  # Individual groups
  #
  # No rectangle is required:
  # panel remains WHITE
  # =======================================
  
  
  # =======================================
  # BACKGROUND SECTION 2
  # Pooled groups
  # =======================================
  
  annotate(
    "rect",
    xmin = pooled_xmin,
    xmax = pooled_xmax,
    ymin = -Inf,
    ymax = Inf,
    fill = pooled_fill,
    colour = NA
  ) +
    
    
    # =======================================
  # BACKGROUND SECTION 3
  # Previous work
  # =======================================
  
  annotate(
    "rect",
    xmin = previous_xmin,
    xmax = previous_xmax_master,
    ymin = -Inf,
    ymax = Inf,
    fill = previous_fill,
    colour = NA
  )
  
  
  # -----------------------------------------
  # Optional grid lines
  # -----------------------------------------
  
  if (show_grid) {
    
    p <- p +
      
      geom_hline(
        yintercept = y_breaks,
        colour = grid_col,
        linewidth = GRID_LW
      ) +
      
      geom_vline(
        xintercept = vbreaks_x,
        colour = grid_col,
        linewidth = GRID_LW
      )
  }
  
  
  # -----------------------------------------
  # Main plot layers
  # -----------------------------------------
  
  p <- p +
    
    # Zero line
    geom_hline(
      yintercept = 0,
      colour = zero_col,
      linetype = "solid",
      linewidth = LWD_ZERO
    ) +
    
    # Overall zooplankton mean
    geom_hline(
      yintercept = panel_mean,
      colour = mean_col,
      linetype = "dashed",
      linewidth = LWD_DASH
    )
  
  
  # -----------------------------------------
  # Minor dashed separators
  # -----------------------------------------
  
  if (length(minor_sep_positions) > 0) {
    
    p <- p +
      
      geom_vline(
        xintercept = minor_sep_positions,
        colour = minor_sep_col,
        linewidth = MINOR_SEP_LW,
        linetype = "longdash"
      )
  }
  
  
  # -----------------------------------------
  # Major separator:
  # This study vs previous work
  # -----------------------------------------
  
  if (length(major_sep_positions) > 0) {
    
    p <- p +
      
      geom_vline(
        xintercept = major_sep_positions,
        colour = major_sep_col,
        linewidth = MAJOR_SEP_LW,
        linetype = "solid"
      )
  }
  
  
  # -----------------------------------------
  # Observations + means
  # -----------------------------------------
  
  p <- p +
    
    geom_point(
      data = raw2,
      aes(
        x = Xkey,
        y = Rate_plot,
        colour = sign
      ),
      position = position_jitter(
        width = 0.17,
        height = 0
      ),
      fill = if (panel_name == "Distribution") "#D6E6F3" else "grey80",
      colour = if (panel_name == "Distribution") "#7FAED3" else "grey60",
      alpha = RAW_PT_ALPHA,
      size = RAW_PT_SIZE,
      stroke = RAW_STROKE,
      shape = 21,
      show.legend = FALSE
    ) +
    
    geom_errorbar(
      data = sum2 %>%
        dplyr::filter(
          !is.na(ci_low_plot),
          !is.na(ci_high_plot)
        ),
      aes(
        x = Xkey,
        ymin = ci_low_plot,
        ymax = ci_high_plot
      ),
      width = 0.15,
      linewidth = ERR_LW,
      colour = mean_col
    ) +
    
    geom_point(
      data = sum2,
      aes(
        x = Xkey,
        y = mean_plot
      ),
      shape = 16,
      colour = mean_col,
      size = MEAN_PT_SIZE
    ) +
    
    
    # Sample sizes
    geom_text(
      data = n2,
      aes(
        x = Xkey,
        y = Inf,
        label = n
      ),
      vjust = -0.90,
      family = fam,
      size = N_LABEL_SIZE,
      colour = mean_col,
      inherit.aes = FALSE
    ) +
    # ---------------------------------------
  # Scales
  # ---------------------------------------
  
  scale_colour_manual(
    values = col_map
  ) +
    
    scale_x_discrete(
      drop = FALSE,
      limits = x_levels_master,
      labels = lab_x
    ) +
    
    scale_y_continuous(
      limits = c(ymin, ymax),
      breaks = y_breaks,
      expand = expansion(
        mult = c(0, 0)
      ),
      oob = scales::squish
    ) +
    
    
    # ---------------------------------------
  # Labels + theme
  # ---------------------------------------
  
  labs(
    y = y_lab,
    x = NULL
  ) +
    
    theme_minimal(
      base_size = BASE_SIZE,
      base_family = fam
    ) +
    
    theme(
      text = element_text(
        family = fam,
        colour = border_col
      ),
      
      axis.text.x = element_blank(),
      
      axis.ticks.x = element_blank(),
      
      axis.title.x = element_blank(),
      
      axis.text.y = element_text(
        size = Y_TICK_TEXT_SIZE,
        colour = border_col,
        margin = margin(r = 5)
      ),
      
      axis.title.y = element_text(
        size = Y_AXIS_TITLE_SIZE,
        colour = border_col,
        margin = margin(r = 8)
      ),
      
      axis.ticks = element_line(
        colour = border_col,
        linewidth = BORDER_LW
      ),
      
      axis.ticks.length = unit(
        4,
        "pt"
      ),
      
      # Important:
      # ggplot's automatic grid remains off,
      # because we manually control grid lines above
      panel.grid = element_blank(),
      
      panel.border = element_rect(
        colour = border_col,
        fill = NA,
        linewidth = BORDER_LW
      ),
      
      plot.margin = margin(
        top_margin,
        14,
        4,
        10
      ),
      
      legend.position = "none"
    ) +
    
    coord_cartesian(
      xlim = c(0.5, length(x_levels_master) + 0.5),
      expand = FALSE,
      clip = "off"
    )
  
  
  return(p)
}


# ===========================================
# Bottom panel function
# ===========================================

make_bottom_panel <- function(show_grid = TRUE) {
  
  
  # -----------------------------------------
  # Background colours
  # -----------------------------------------
  
  p <- ggplot() +
    
    # Pooled groups = light grey
    annotate(
      "rect",
      xmin = pooled_xmin,
      xmax = pooled_xmax,
      ymin = -Inf,
      ymax = Inf,
      fill = pooled_fill,
      colour = NA
    ) +
    
    # Previous work = darker grey
    annotate(
      "rect",
      xmin = previous_xmin,
      xmax = previous_xmax_prop,
      ymin = -Inf,
      ymax = Inf,
      fill = previous_fill,
      colour = NA
    )
  
  
  # -----------------------------------------
  # Optional grid
  # -----------------------------------------
  
  if (show_grid) {
    
    p <- p +
      
      geom_hline(
        yintercept = brks_y_prop,
        colour = grid_col,
        linewidth = GRID_LW
      ) +
      
      geom_vline(
        xintercept = vbreaks_x,
        colour = grid_col,
        linewidth = GRID_LW
      )
  }
  
  
  # -----------------------------------------
  # Main reference lines
  # -----------------------------------------
  
  p <- p +
    
    geom_hline(
      yintercept = 50,
      colour = zero_col,
      linewidth = LWD_BOT_MIDLINE
    )
  
  
  if (is.finite(p_zoop_dist)) {
    
    p <- p +
      
      geom_hline(
        yintercept = p_zoop_dist,
        linetype = "dashed",
        colour = mean_col_dist,
        linewidth = LWD_DASH
      )
  }
  
  
  if (is.finite(p_zoop_phen)) {
    
    p <- p +
      
      geom_hline(
        yintercept = p_zoop_phen,
        linetype = "dashed",
        colour = mean_col_phen,
        linewidth = LWD_DASH
      )
  }
  
  
  # -----------------------------------------
  # Vertical separators
  # -----------------------------------------
  
  if (length(minor_sep_positions) > 0) {
    
    p <- p +
      
      geom_vline(
        xintercept = minor_sep_positions,
        colour = minor_sep_col,
        linewidth = MINOR_SEP_LW,
        linetype = "longdash"
      )
  }
  
  
  if (length(major_sep_positions) > 0) {
    
    p <- p +
      
      geom_vline(
        xintercept = major_sep_positions,
        colour = major_sep_col,
        linewidth = MAJOR_SEP_LW,
        linetype = "solid"
      )
  }
  
  
  # =========================================
  # Main bottom-panel data
  # =========================================
  
  p <- p +
    
    # Phenology
    geom_errorbar(
      data = se_phen2,
      aes(
        x = xplot,
        ymin = lo,
        ymax = hi
      ),
      width = 0.15,
      linewidth = ERR_LW,
      colour = mean_col_phen
    ) +
    
    geom_point(
      data = se_phen2,
      aes(
        x = xplot,
        y = p
      ),
      size = MEAN_PT_SIZE,
      shape = 16,
      colour = mean_col_phen
    ) +
    
    
    # Distribution
    geom_errorbar(
      data = se_dist2,
      aes(
        x = xplot,
        ymin = lo,
        ymax = hi
      ),
      width = 0.15,
      linewidth = ERR_LW,
      colour = mean_col_dist
    ) +
    
    geom_point(
      data = se_dist2,
      aes(
        x = xplot,
        y = p
      ),
      size = MEAN_PT_SIZE,
      shape = 16,
      colour = mean_col_dist
    )
  
  
  # =========================================
  # Optional phytoplankton echo data
  # =========================================
  
  if (!is.null(se_phen_phytoEcho)) {
    
    p <- p +
      
      geom_errorbar(
        data = se_phen_phytoEcho,
        aes(
          x = xplot,
          ymin = lo,
          ymax = hi
        ),
        width = 0.15,
        linewidth = ERR_LW,
        colour = mean_col_phen
      ) +
      
      geom_point(
        data = se_phen_phytoEcho,
        aes(
          x = xplot,
          y = p
        ),
        size = MEAN_PT_SIZE,
        shape = 16,
        colour = mean_col_phen
      )
  }
  
  
  if (!is.null(se_dist_phytoEcho)) {
    
    p <- p +
      
      geom_errorbar(
        data = se_dist_phytoEcho,
        aes(
          x = xplot,
          ymin = lo,
          ymax = hi
        ),
        width = 0.15,
        linewidth = ERR_LW,
        colour = mean_col_dist
      ) +
      
      geom_point(
        data = se_dist_phytoEcho,
        aes(
          x = xplot,
          y = p
        ),
        size = MEAN_PT_SIZE,
        shape = 16,
        colour = mean_col_dist
      )
  }
  
  
  # =========================================
  # Axes
  # =========================================
  
  p <- p +
    
    scale_x_continuous(
      breaks = idx_tbl_prop$Xnum,
      labels = pretty_taxon_prop(
        as.character(idx_tbl_prop$Ykey)
      ),
      limits = c(
        0.5,
        length(x_levels_prop) + 0.5
      ),
      expand = expansion(
        mult = c(0, 0)
      )
    ) +
    
    scale_y_continuous(
      limits = c(0, 100),
      breaks = brks_y_prop,
      labels = label_number(
        accuracy = 1
      ),
      expand = expansion(
        mult = c(0, 0)
      )
    ) +
    
    labs(
      y = "Percentage consistent (%)",
      x = NULL
    ) +
    
    
    # =======================================
  # Theme
  # =======================================
  
  theme_minimal(
    base_size = BASE_SIZE,
    base_family = fam
  ) +
    
    theme(
      text = element_text(
        family = fam,
        colour = border_col
      ),
      
      axis.text.x = element_text(
        size = X_TICK_TEXT_SIZE,
        angle = 35,
        hjust = 1,
        vjust = 1,
        colour = border_col
      ),
      
      axis.text.y = element_text(
        size = Y_TICK_TEXT_SIZE,
        colour = border_col,
        margin = margin(r = 5)
      ),
      
      axis.title.y = element_text(
        size = Y_AXIS_TITLE_SIZE,
        colour = border_col,
        margin = margin(r = 8)
      ),
      
      # Automatic ggplot grid off;
      # manually controlled by show_grid
      panel.grid = element_blank(),
      
      axis.ticks = element_line(
        colour = border_col,
        linewidth = BORDER_LW
      ),
      
      axis.ticks.length = unit(
        4,
        "pt"
      ),
      
      plot.margin = margin(
        12,
        18,
        5,
        0
      ),
      
      legend.position = "none",
      
      panel.border = element_rect(
        colour = border_col,
        fill = NA,
        linewidth = BORDER_LW
      )
    ) +
    
    coord_cartesian(
      clip = "off"
    )
  
  
  # =========================================
  # Sample-size labels
  # =========================================
  
  p <- p +
    
    geom_text(
      data = n_phen2,
      aes(
        x = Xnum,
        y = Inf,
        label = n_phen
      ),
      vjust = -1.55,
      hjust = 0.5,
      family = fam,
      size = BOTTOM_LABEL_SIZE,
      colour = mean_col_phen
    ) +
    
    geom_text(
      data = n_dist2,
      aes(
        x = Xnum,
        y = Inf,
        label = n_dist
      ),
      vjust = -3.25,
      hjust = 0.5,
      family = fam,
      size = BOTTOM_LABEL_SIZE,
      colour = mean_col_dist
    )
  
  
  if (!is.null(n_phen_phytoEcho)) {
    
    p <- p +
      
      geom_text(
        data = n_phen_phytoEcho,
        aes(
          x = Xnum,
          y = Inf,
          label = n_phen
        ),
        vjust = -1.55,
        hjust = 0.5,
        family = fam,
        size = BOTTOM_LABEL_SIZE,
        colour = mean_col_phen
      )
  }
  
  
  if (!is.null(n_dist_phytoEcho)) {
    
    p <- p +
      
      geom_text(
        data = n_dist_phytoEcho,
        aes(
          x = Xnum,
          y = Inf,
          label = n_dist
        ),
        vjust = -3.25,
        hjust = 0.5,
        family = fam,
        size = BOTTOM_LABEL_SIZE,
        colour = mean_col_dist
      )
  }
  
  
  # =========================================
  # Distribution / Phenology labels
  # =========================================
  
  p <- p +
    
    annotate(
      "text",
      x = 0.5,
      y = Inf,
      label = "Distribution",
      vjust = -3.25,
      hjust = 1,
      family = "header_bold",
      fontface = "plain",
      size = BOTTOM_LABEL_SIZE,
      colour = mean_col_dist
    ) +
    
    annotate(
      "text",
      x = 0.5,
      y = Inf,
      label = "Phenology",
      vjust = -1.55,
      hjust = 1,
      family = "header_bold",
      fontface = "plain",
      size = BOTTOM_LABEL_SIZE,
      colour = mean_col_phen
    )
  
  
  return(p)
}


# ===========================================
# Function to build complete figure
# ===========================================

make_full_plot <- function(show_grid = TRUE) {
  
  # Annotation placed outside the right edge.
  # Positions use panel-relative coordinates.
  direction_annotation <- function(upper_label, lower_label) {
    
    grid::grobTree(
      grid::segmentsGrob(
        x0 = grid::unit(1, "npc") + grid::unit(21, "pt"),
        x1 = grid::unit(1, "npc") + grid::unit(21, "pt"),
        y0 = grid::unit(0.30, "npc"),
        y1 = grid::unit(0.70, "npc"),
        arrow = grid::arrow(
          ends = "both",
          type = "closed",
          length = grid::unit(2.5, "pt")
        ),
        gp = grid::gpar(
          col = border_col,
          fill = border_col,
          lwd = 0.4
        )
      ),
      
      grid::textGrob(
        label = upper_label,
        x = grid::unit(1, "npc") + grid::unit(21, "pt"),
        y = grid::unit(0.77, "npc"),
        gp = grid::gpar(
          fontfamily = fam,
          fontsize = X_TICK_TEXT_SIZE,
          col = border_col
        )
      ),
      
      grid::textGrob(
        label = lower_label,
        x = grid::unit(1, "npc") + grid::unit(21, "pt"),
        y = grid::unit(0.23, "npc"),
        gp = grid::gpar(
          fontfamily = fam,
          fontsize = X_TICK_TEXT_SIZE,
          col = border_col
        )
      )
    )
  }
  
  
  
  # Panel a: Percentage consistent
  p_consistent <- make_bottom_panel(
    show_grid = show_grid
  ) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      plot.margin = margin(12, 40, 5, 20)
    )
  
  # Panel b: Distribution
  p_distribution <- make_rate_panel(
    panel_name = "Distribution",
    ymin = dist_ymin,
    ymax = dist_ymax,
    y_breaks = brks_y_dist,
    mean_col = mean_col_dist,
    y_lab = expression(
      "Distribution change (km dec"^{-1}*")"
    ),
    top_margin = 4,
    show_grid = show_grid
  ) +
    annotation_custom(
      grob = direction_annotation(
        upper_label = "Poleward",
        lower_label = "Equatorward"
      ),
      xmin = -Inf,
      xmax = Inf,
      ymin = -Inf,
      ymax = Inf
    ) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      plot.margin = margin(4, 40, 4, 20)
    )
  
  # Panel c: Phenology
  # Positive values = later; negative values = earlier.
  p_phenology <- make_rate_panel(
    panel_name = "Phenology",
    ymin = phen_ymin,
    ymax = phen_ymax,
    y_breaks = brks_y_phen,
    mean_col = mean_col_phen,
    y_lab = expression(
      "Phenology change (days dec"^{-1}*")"
    ),
    top_margin = 8,
    show_grid = show_grid
  ) +
    annotation_custom(
      grob = direction_annotation(
        upper_label = "Later",
        lower_label = "Earlier"
      ),
      xmin = -Inf,
      xmax = Inf,
      ymin = -Inf,
      ymax = Inf
    ) +
    theme(
      axis.text.x = element_text(
        size = X_TICK_TEXT_SIZE,
        angle = 35,
        hjust = 1,
        vjust = 1,
        colour = border_col
      ),
      axis.ticks.x = element_line(
        colour = border_col,
        linewidth = BORDER_LW
      ),
      plot.margin = margin(8, 40, 4, 20)
    )
  
  # Match the header margin to the panels.
  header <- p_header +
    theme(
      plot.margin = margin(8, 40, -2, 20)
    )
  
  gap <- plot_spacer() +
    theme(
      plot.margin = margin(0, 0, 0, 0)
    )
  
  # Bold panel labels at the top-left
  panel_tag_theme <- theme(
    plot.tag = element_text(
      family = "header_bold",
      fontface = "plain",
      size = 34,
      colour = border_col,
      hjust = 3.75,
      vjust = -1.75
    ),
    plot.tag.position = c(0, 1)
  )
  
  p_consistent <- p_consistent +
    labs(tag = "a") +
    panel_tag_theme
  
  p_distribution <- p_distribution +
    labs(tag = "b") +
    panel_tag_theme
  
  p_phenology <- p_phenology +
    labs(tag = "c") +
    panel_tag_theme
  
  p_final <- header /
    p_consistent /
    gap /
    p_distribution /
    gap /
    p_phenology +
    plot_layout(
      heights = c(
        0.22,
        1.08,
        GAP_HEIGHT,
        1,
        0.20,
        1
      )
    )
  
  return(p_final)
}
# ===========================================
# Generate BOTH versions
# ===========================================

# Version 1:
# Background shading + grid lines

p_final_grid <- make_full_plot(
  show_grid = TRUE
)


# Version 2:
# Background shading + NO grid lines

p_final_no_grid <- make_full_plot(
  show_grid = FALSE
)


# ===========================================
# Preview figures
# ===========================================

print(p_final_grid)

print(p_final_no_grid)


# ===========================================
# Export
# ===========================================


# -------------------------------------------
# Version 1: WITH GRID
# -------------------------------------------

png_grid_out <- file.path(
  figures_dir,
  "zooplankton_observational_three_panel_GRID.png"
)


ggsave(
  filename = png_grid_out,
  plot = p_final_grid,
  width = 6.4,
  height = 5.6,
  units = "in",
  dpi = 400,
  bg = "white"
)



# -------------------------------------------
# Version 2: NO GRID
# -------------------------------------------

png_no_grid_out <- file.path(
  figures_dir,
  "zooplankton_observational_three_panel_NO_GRID.png"
)

ggsave(
  filename = png_no_grid_out,
  plot = p_final_no_grid,
  width = 6.4,
  height = 5.6,
  units = "in",
  dpi = 400,
  bg = "white"
)


# ===========================================
# Done
# ===========================================
 
message(
  "Saved WITH-GRID PNG to: ",
  png_grid_out
)

message(
  "Saved NO-GRID PNG to: ",
  png_no_grid_out
)

