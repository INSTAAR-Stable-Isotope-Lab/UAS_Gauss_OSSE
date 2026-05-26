#!/usr/bin/env Rscript
###############################################################################
# make_gaussian_plume_stability_comparison.R
#
# Standalone generator for the Pasquill stability comparison plume figure used
# in the full-analysis figure set. The plume is the analytical disk-integrated
# Gaussian plume over the small 100 m lake.
###############################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

find_project_root <- function() {
  wd <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  candidates <- c(
    wd,
    dirname(wd),
    normalizePath(file.path(wd, ".."), winslash = "/", mustWork = FALSE),
    normalizePath(file.path(wd, "../.."), winslash = "/", mustWork = FALSE)
  )
  for (candidate in unique(candidates)) {
    if (file.exists(file.path(candidate, "R", "utils.R"))) {
      return(normalizePath(candidate, winslash = "/", mustWork = TRUE))
    }
  }
  stop("Cannot find project root. Run from the repository root or scripts/full_analysis/.")
}

proj_root <- find_project_root()
options(uas_gauss_osse_root = proj_root)

suppressMessages(source(file.path(proj_root, "R", "utils.R")))
suppressMessages(source(file.path(proj_root, "R", "02b_gaussian_forward.R")))
suppressMessages(source(file.path(proj_root, "R", "06_visualization.R")))


# -----------------------------------------------------------------------------
# 1. Figure settings
# -----------------------------------------------------------------------------
LAKE_NAME <- "small"
DESIGN_NAME <- "design_24az_5h_30min"
FLUX_MG_M2_D <- 50
U_REF_M_S <- 2
WD_FROM_DEG <- 225
CLASSES <- c("A", "B", "C", "D", "E", "F")

X_RANGE_M <- c(1, 500)
Z_RANGE_M <- c(0.1, 100)
GRID_NX <- 250L
GRID_NZ <- 180L

CONC_MIN_PPB <- 1e-4
RECEPTOR_DISTANCE_FROM_EDGE_M <- 10

lk <- LAKE_CONFIGS[[LAKE_NAME]]
design_cfg <- SAMPLING_CONFIGS[[DESIGN_NAME]]

LAKE_RADIUS_M <- lk$diameter_m / 2
RECEPTOR_X_M <- LAKE_RADIUS_M + RECEPTOR_DISTANCE_FROM_EDGE_M
RECEPTOR_Z_MAX_M <- Z_RANGE_M[2]

FIG_DIR <- file.path(proj_root, "output")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

FIG_PNG <- file.path(FIG_DIR, "gaussian_plume_stability_comparison.png")


# -----------------------------------------------------------------------------
# 2. Build plume grids
# -----------------------------------------------------------------------------
message("Building Gaussian plume stability comparison...")
message(sprintf("  Flux:      %g mg CH4 m-2 d-1", FLUX_MG_M2_D))
message(sprintf("  Wind:      %g m/s", U_REF_M_S))
message(sprintf("  Downwind:  %g-%g m", X_RANGE_M[1], X_RANGE_M[2]))
message(sprintf("  Output:    %s", normalizePath(FIG_PNG, winslash = "/", mustWork = FALSE)))

grids <- lapply(CLASSES, function(pc) {
  g <- gaussian_plume_grid(
    source_radius = LAKE_RADIUS_M,
    flux_mg_m2_d = FLUX_MG_M2_D,
    U_ref = U_REF_M_S,
    wd_deg = WD_FROM_DEG,
    pasquill_class = pc,
    plane = "xz",
    x_range = X_RANGE_M,
    z_range = Z_RANGE_M,
    nx = GRID_NX,
    ny = GRID_NZ
  )
  g[, pasquill_class := pc]
  g
})

dt <- rbindlist(grids)
dt[conc_ppb < CONC_MIN_PPB, conc_ppb := NA_real_]
dt[, class_label := sprintf(
  "Class %s (%s)",
  pasquill_class,
  PASQUILL_CLASSES[[pasquill_class]]$description
), by = pasquill_class]
dt[, class_label := factor(
  class_label,
  levels = sapply(CLASSES, function(pc) {
    sprintf("Class %s (%s)", pc, PASQUILL_CLASSES[[pc]]$description)
  })
)]

lake_edge_marks <- unique(dt[, .(class_label)])
lake_edge_marks[, `:=`(
  x = LAKE_RADIUS_M,
  z = 0
)]

receptor_line <- unique(dt[, .(class_label)])
receptor_line[, `:=`(
  x = RECEPTOR_X_M,
  z0 = 0,
  z1 = RECEPTOR_Z_MAX_M
)]


# -----------------------------------------------------------------------------
# 3. Plot and save
# -----------------------------------------------------------------------------
p <- ggplot(dt, aes(x = x_down, y = z, fill = conc_ppb)) +
  geom_tile() +
  geom_point(
    data = lake_edge_marks,
    aes(x = x, y = z),
    inherit.aes = FALSE,
    shape = 25,
    size = 2.5,
    stroke = 0.35,
    colour = "black",
    fill = "black"
  ) +
  geom_segment(
    data = receptor_line,
    aes(x = x, xend = x, y = z0, yend = z1),
    inherit.aes = FALSE,
    colour = "grey15",
    linewidth = 0.75,
    linetype = "longdash"
  ) +
  facet_wrap(~class_label, ncol = 3L, scales = "fixed") +
  scale_fill_viridis_c(
    option = "C",
    na.value = "white",
    trans = "log10",
    name = expression(Delta*CH[4]~"(ppb)")
  ) +
  coord_cartesian(
    xlim = c(0, X_RANGE_M[2]),
    ylim = c(0, Z_RANGE_M[2]),
    expand = FALSE,
    clip = "off"
  ) +
  labs(
    title = "Gaussian plume: Pasquill stability class comparison (A-F)",
    subtitle = sprintf(
      "F = %g mg m-2 d-1, U = %g m/s, R = %g m, centerline x-z",
      FLUX_MG_M2_D, U_REF_M_S, LAKE_RADIUS_M
    ),
    x = "Downwind distance from lake center (m)",
    y = "Height (m)"
  ) +
  theme_osse(base_size = 10) +
  theme(
    legend.position = "right",
    strip.text = element_text(size = 9),
    plot.margin = margin(8, 8, 8, 8)
  )

ggsave_osse(FIG_PNG, p, width = 14, height = 8, dpi = 220)

message("Saved: ", normalizePath(FIG_PNG, winslash = "/", mustWork = FALSE))

