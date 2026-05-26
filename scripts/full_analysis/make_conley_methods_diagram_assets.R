#!/usr/bin/env Rscript
###############################################################################
# make_conley_methods_diagram_assets.R
#
# Conceptual, transparent-background method-diagram assets for explaining the
# Conley receptor geometry and control-volume interpretation. These figures are
# intentionally schematic rather than simulation outputs.
###############################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(grid)
})

proj_root <- if (file.exists("R/utils.R")) {
  normalizePath(".")
} else if (file.exists("../R/utils.R")) {
  normalizePath("..")
} else {
  stop("Cannot find project root. Run from the repository root or scripts/full_analysis/.")
}
options(uas_gauss_osse_root = proj_root)

suppressMessages(source(file.path(proj_root, "R", "utils.R")))
suppressMessages(source(file.path(proj_root, "R", "06_visualization.R")))


# -----------------------------------------------------------------------------
# 1. Settings
# -----------------------------------------------------------------------------
OUTPUT_DIR <- proj_path("output")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

AZIMUTHAL_RESOLUTION_PNG <- file.path(
  OUTPUT_DIR,
  "conley_methods_azimuthal_resolution.png"
)
ISOMETRIC_RETRIEVAL_PNG <- file.path(
  OUTPUT_DIR,
  "conley_methods_isometric_retrieval.png"
)
LOW_RES_CONTROL_SURFACE_PNG <- file.path(
  OUTPUT_DIR,
  "conley_methods_low_resolution_control_surface.png"
)

COL_LAKE <- "#3F9FB5"
COL_LAKE_DARK <- "#1F6E7C"
COL_RING <- "#2F3B46"
COL_RECEPTOR <- "#253342"
COL_RECEPTOR_FAINT <- "#5E7487"
COL_ARROW <- "#1F2933"
COL_FLUX <- "#238B45"
COL_SURFACE <- "#8FB6D8"
COL_SURFACE_EDGE <- "#3B596F"
COL_NORMAL <- "#D05B21"
COL_WIND <- "#58636F"

LAKE_RADIUS_M <- 50
CONTROL_RADIUS_M <- 65


# -----------------------------------------------------------------------------
# 2. Geometry helpers
# -----------------------------------------------------------------------------
transparent_theme <- function() {
  theme_void() +
    theme(
      legend.position = "none",
      plot.background = element_rect(fill = "transparent", colour = NA),
      panel.background = element_rect(fill = "transparent", colour = NA),
      plot.margin = margin(6, 6, 6, 6)
    )
}

save_asset <- function(path, plot, width, height, dpi = 300) {
  ggsave_osse(
    path,
    plot,
    width = width,
    height = height,
    dpi = dpi,
    bg = "transparent"
  )
  message(sprintf("Saved asset: %s", normalizePath(path, mustWork = FALSE)))
}

circle_path <- function(center_x = 0, center_y = 0, radius = 1, n = 240L) {
  theta <- seq(0, 2 * pi, length.out = n)
  data.table(
    x = center_x + radius * cos(theta),
    y = center_y + radius * sin(theta)
  )
}

azimuth_path <- function(radius,
                         az_start_deg,
                         az_end_deg,
                         n = 120L,
                         center_x = 0,
                         center_y = 0) {
  az <- seq(az_start_deg, az_end_deg, length.out = n)
  data.table(
    azimuth_deg = az,
    x = center_x + radius * sin(az * pi / 180),
    y = center_y + radius * cos(az * pi / 180)
  )
}

iso_project <- function(x, y, z, z_scale = 1.35) {
  data.table(
    x2 = 0.86 * (x - y),
    y2 = 0.46 * (x + y) + z * z_scale
  )
}

iso_ring <- function(radius, z, n = 240L) {
  az <- seq(0, 360, length.out = n)
  xy <- data.table(
    azimuth_deg = az,
    x = radius * sin(az * pi / 180),
    y = radius * cos(az * pi / 180),
    z = z
  )
  cbind(xy, iso_project(xy$x, xy$y, xy$z))
}

iso_ring_points <- function(radius, z, n_azimuth, ring_name) {
  az <- seq(0, 360 - 360 / n_azimuth, length.out = n_azimuth)
  xy <- data.table(
    azimuth_deg = az,
    x = radius * sin(az * pi / 180),
    y = radius * cos(az * pi / 180),
    z = z,
    ring_name = ring_name
  )
  cbind(xy, iso_project(xy$x, xy$y, xy$z))
}

iso_arc <- function(radius, z, az_start_deg, az_end_deg, n = 140L) {
  az <- seq(az_start_deg, az_end_deg, length.out = n)
  xy <- data.table(
    azimuth_deg = az,
    x = radius * sin(az * pi / 180),
    y = radius * cos(az * pi / 180),
    z = z
  )
  cbind(xy, iso_project(xy$x, xy$y, xy$z))
}

tail_segment <- function(path_dt, x_col = "x", y_col = "y", n_tail = 5L) {
  n <- nrow(path_dt)
  i0 <- max(1L, n - n_tail)
  data.table(
    x = path_dt[[x_col]][i0],
    y = path_dt[[y_col]][i0],
    xend = path_dt[[x_col]][n],
    yend = path_dt[[y_col]][n]
  )
}

uas_symbol_layers <- function(center_x, center_y, size = 10,
                              colour = COL_ARROW) {
  arm <- size * 0.62
  rotor <- size * 0.27
  body <- circle_path(center_x, center_y, size * 0.17, n = 40L)
  rotors <- rbindlist(lapply(
    list(c(-arm, -arm), c(-arm, arm), c(arm, -arm), c(arm, arm)),
    function(offset) {
      dt <- circle_path(center_x + offset[1], center_y + offset[2], rotor, n = 40L)
      dt[, rotor_id := paste(offset, collapse = "_")]
      dt
    }
  ))
  arms <- data.table(
    x = c(center_x - arm, center_x - arm, center_x + arm, center_x + arm),
    y = c(center_y - arm, center_y + arm, center_y - arm, center_y + arm),
    xend = center_x,
    yend = center_y
  )

  list(
    geom_segment(
      data = arms,
      aes(x = x, y = y, xend = xend, yend = yend),
      inherit.aes = FALSE,
      colour = colour,
      linewidth = 0.8,
      lineend = "round"
    ),
    geom_path(
      data = rotors,
      aes(x = x, y = y, group = rotor_id),
      inherit.aes = FALSE,
      colour = colour,
      linewidth = 0.75
    ),
    geom_polygon(
      data = body,
      aes(x = x, y = y),
      inherit.aes = FALSE,
      fill = colour,
      colour = colour,
      linewidth = 0.3
    )
  )
}


# -----------------------------------------------------------------------------
# 3. Asset 1: top-down azimuthal receptor resolution
# -----------------------------------------------------------------------------
build_azimuthal_resolution_plot <- function() {
  lake <- circle_path(radius = LAKE_RADIUS_M)
  control_ring <- circle_path(radius = CONTROL_RADIUS_M)

  raw_increments <- seq(28, 7, length.out = 28)
  increments <- raw_increments * 360 / sum(raw_increments)
  azimuths <- (215 + c(0, cumsum(increments[-length(increments)]))) %% 360
  receptors <- data.table(
    receptor_order = seq_along(azimuths),
    azimuth_deg = azimuths
  )
  receptors[, `:=`(
    x = CONTROL_RADIUS_M * sin(azimuth_deg * pi / 180),
    y = CONTROL_RADIUS_M * cos(azimuth_deg * pi / 180),
    receptor_colour = grDevices::colorRampPalette(c("#AFC3D2", COL_RECEPTOR))(length(azimuth_deg))
  )]

  resolution_arrow <- azimuth_path(
    radius = CONTROL_RADIUS_M + 20,
    az_start_deg = 255,
    az_end_deg = 375,
    n = 110L
  )
  arrow_tail <- tail_segment(resolution_arrow)

  ggplot() +
    geom_polygon(
      data = lake,
      aes(x = x, y = y),
      fill = COL_LAKE,
      alpha = 0.64,
      colour = COL_LAKE_DARK,
      linewidth = 0.9
    ) +
    geom_path(
      data = control_ring,
      aes(x = x, y = y),
      colour = COL_RING,
      linewidth = 0.9,
      alpha = 0.95
    ) +
    geom_point(
      data = receptors,
      aes(x = x, y = y, fill = receptor_colour),
      shape = 21,
      size = 3.9,
      colour = "white",
      stroke = 0.55
    ) +
    geom_path(
      data = resolution_arrow,
      aes(x = x, y = y),
      colour = COL_ARROW,
      linewidth = 1.25,
      lineend = "round"
    ) +
    geom_segment(
      data = arrow_tail,
      aes(x = x, y = y, xend = xend, yend = yend),
      colour = COL_ARROW,
      linewidth = 1.25,
      arrow = arrow(length = unit(0.18, "inches"), type = "closed"),
      inherit.aes = FALSE
    ) +
    scale_fill_identity() +
    coord_equal(xlim = c(-100, 100), ylim = c(-100, 100), expand = FALSE) +
    transparent_theme()
}


# -----------------------------------------------------------------------------
# 4. Asset 2: isometric Conley retrieval profile
# -----------------------------------------------------------------------------
build_isometric_retrieval_plot <- function() {
  lake <- iso_ring(LAKE_RADIUS_M, z = 0)
  lower_ring <- iso_ring(CONTROL_RADIUS_M, z = 18)
  middle_ring <- iso_ring(CONTROL_RADIUS_M, z = 42)
  upper_ring <- iso_ring(CONTROL_RADIUS_M, z = 66)
  top_ring <- iso_ring(CONTROL_RADIUS_M, z = 90)
  context_rings <- rbindlist(list(
    copy(lower_ring)[, ring_alpha := 0.22],
    copy(upper_ring)[, ring_alpha := 0.22],
    copy(top_ring)[, ring_alpha := 0.18]
  ))

  receptors <- rbindlist(list(
    iso_ring_points(CONTROL_RADIUS_M, 18, 18, "context"),
    iso_ring_points(CONTROL_RADIUS_M, 42, 32, "active"),
    iso_ring_points(CONTROL_RADIUS_M, 66, 18, "context"),
    iso_ring_points(CONTROL_RADIUS_M, 90, 18, "context")
  ))
  receptors[, `:=`(
    point_alpha = fifelse(ring_name == "active", 0.98, 0.28),
    point_size = fifelse(ring_name == "active", 2.3, 1.7),
    point_colour = fifelse(ring_name == "active", COL_RECEPTOR, COL_RECEPTOR_FAINT)
  )]

  active_arc <- iso_arc(CONTROL_RADIUS_M, z = 42, az_start_deg = 205, az_end_deg = 330)
  active_tail <- tail_segment(active_arc, x_col = "x2", y_col = "y2")
  flux_arrow_pts <- cbind(
    data.table(name = c("start", "end"), x = c(0, 0), y = c(0, 0), z = c(0, 56)),
    iso_project(c(0, 0), c(0, 0), c(0, 56))
  )
  flux_arrow <- data.table(
    x = flux_arrow_pts$x2[1],
    y = flux_arrow_pts$y2[1],
    xend = flux_arrow_pts$x2[2],
    yend = flux_arrow_pts$y2[2]
  )
  uas_az_deg <- 285
  uas_anchor_3d <- data.table(
    x = CONTROL_RADIUS_M * sin(uas_az_deg * pi / 180),
    y = CONTROL_RADIUS_M * cos(uas_az_deg * pi / 180),
    z = 42
  )
  uas_anchor <- iso_project(uas_anchor_3d$x, uas_anchor_3d$y, uas_anchor_3d$z)

  ggplot() +
    geom_polygon(
      data = lake,
      aes(x = x2, y = y2),
      fill = COL_LAKE,
      alpha = 0.58,
      colour = COL_LAKE_DARK,
      linewidth = 0.8
    ) +
    geom_path(
      data = context_rings,
      aes(x = x2, y = y2, group = interaction(z, ring_alpha)),
      colour = COL_RING,
      alpha = 0.25,
      linewidth = 0.8
    ) +
    geom_path(
      data = middle_ring,
      aes(x = x2, y = y2),
      colour = COL_RING,
      alpha = 0.72,
      linewidth = 0.95
    ) +
    geom_segment(
      data = flux_arrow,
      aes(x = x, y = y, xend = xend, yend = yend),
      inherit.aes = FALSE,
      colour = COL_FLUX,
      linewidth = 1.55,
      arrow = arrow(length = unit(0.17, "inches"), type = "closed"),
      lineend = "round"
    ) +
    geom_path(
      data = active_arc,
      aes(x = x2, y = y2),
      colour = COL_ARROW,
      linewidth = 1.35,
      lineend = "round"
    ) +
    geom_segment(
      data = active_tail,
      aes(x = x, y = y, xend = xend, yend = yend),
      inherit.aes = FALSE,
      colour = COL_ARROW,
      linewidth = 1.35,
      arrow = arrow(length = unit(0.18, "inches"), type = "closed")
    ) +
    geom_point(
      data = receptors,
      aes(x = x2, y = y2, alpha = point_alpha, size = point_size, colour = point_colour),
      show.legend = FALSE
    ) +
    uas_symbol_layers(uas_anchor$x2, uas_anchor$y2, size = 8, colour = COL_NORMAL) +
    scale_alpha_identity() +
    scale_size_identity() +
    scale_colour_identity() +
    coord_equal(xlim = c(-130, 125), ylim = c(-65, 178), expand = FALSE) +
    transparent_theme()
}


# -----------------------------------------------------------------------------
# 5. Asset 3: low-resolution control surface cells
# -----------------------------------------------------------------------------
build_lateral_cell_polygons <- function(radius = CONTROL_RADIUS_M,
                                        z_breaks = c(0, 24, 48, 72),
                                        n_azimuth = 6,
                                        n_edge = 10L) {
  az_edges <- seq(0, 360, length.out = n_azimuth + 1L)
  polys <- list()
  idx <- 0L

  for (iz in seq_len(length(z_breaks) - 1L)) {
    for (ia in seq_len(n_azimuth)) {
      az0 <- az_edges[ia]
      az1 <- az_edges[ia + 1L]
      z0 <- z_breaks[iz]
      z1 <- z_breaks[iz + 1L]
      az_forward <- seq(az0, az1, length.out = n_edge)
      az_back <- rev(az_forward)
      x <- c(
        radius * sin(az_forward * pi / 180),
        radius * sin(az_back * pi / 180)
      )
      y <- c(
        radius * cos(az_forward * pi / 180),
        radius * cos(az_back * pi / 180)
      )
      z <- c(rep(z0, length(az_forward)), rep(z1, length(az_back)))
      pr <- iso_project(x, y, z)
      idx <- idx + 1L
      polys[[idx]] <- data.table(
        cell_id = sprintf("lat_%02d_%02d", iz, ia),
        az_mid = (az0 + az1) / 2,
        z_mid = (z0 + z1) / 2,
        depth = mean(x + y),
        x2 = pr$x2,
        y2 = pr$y2
      )
    }
  }

  out <- rbindlist(polys)
  cell_order <- out[, .(depth = unique(depth)), by = cell_id][order(depth)]$cell_id
  out[, cell_id := factor(cell_id, levels = cell_order)]
  setorder(out, cell_id)
  out
}

build_top_cap_polygons <- function(radius = CONTROL_RADIUS_M,
                                   z = 72,
                                   n_azimuth = 6,
                                   n_edge = 12L) {
  az_edges <- seq(0, 360, length.out = n_azimuth + 1L)
  r_inner <- 0
  polys <- vector("list", n_azimuth)

  for (ia in seq_len(n_azimuth)) {
    az0 <- az_edges[ia]
    az1 <- az_edges[ia + 1L]
    az <- seq(az0, az1, length.out = n_edge)
    x <- c(r_inner, radius * sin(az * pi / 180), r_inner)
    y <- c(r_inner, radius * cos(az * pi / 180), r_inner)
    zz <- rep(z, length(x))
    pr <- iso_project(x, y, zz)
    polys[[ia]] <- data.table(
      cell_id = sprintf("top_%02d", ia),
      x2 = pr$x2,
      y2 = pr$y2
    )
  }
  rbindlist(polys)
}

build_normal_arrows <- function(radius = CONTROL_RADIUS_M,
                                z_breaks = c(0, 24, 48, 72),
                                n_azimuth = 6) {
  az_edges <- seq(0, 360, length.out = n_azimuth + 1L)
  arrows <- list()
  idx <- 0L
  for (iz in seq_len(length(z_breaks) - 1L)) {
    z_mid <- mean(z_breaks[iz:(iz + 1L)])
    for (ia in seq_len(n_azimuth)) {
      az_mid <- mean(az_edges[ia:(ia + 1L)])
      nx <- sin(az_mid * pi / 180)
      ny <- cos(az_mid * pi / 180)
      p0 <- iso_project((radius + 2) * nx, (radius + 2) * ny, z_mid)
      p1 <- iso_project((radius + 20) * nx, (radius + 20) * ny, z_mid)
      idx <- idx + 1L
      arrows[[idx]] <- data.table(
        normal_id = sprintf("n_%02d_%02d", iz, ia),
        rear_side = ((radius + 2) * (nx + ny)) > 0,
        x = p0$x2,
        y = p0$y2,
        xend = p1$x2,
        yend = p1$y2
      )
    }
  }
  rbindlist(arrows)
}

build_top_normal_arrow <- function(z = 72) {
  p0 <- iso_project(0, 0, z)
  p1 <- iso_project(0, 0, z + 28)
  data.table(
    x = p0$x2,
    y = p0$y2,
    xend = p1$x2,
    yend = p1$y2
  )
}

build_low_resolution_control_surface_plot <- function() {
  lake <- iso_ring(LAKE_RADIUS_M, z = 0)
  bottom_ring <- iso_ring(CONTROL_RADIUS_M, z = 0)
  top_ring <- iso_ring(CONTROL_RADIUS_M, z = 72)
  lateral_cells <- build_lateral_cell_polygons()
  top_cells <- build_top_cap_polygons()
  normals <- build_normal_arrows()
  rear_normals <- normals[rear_side == TRUE]
  front_normals <- normals[rear_side == FALSE]
  top_normal <- build_top_normal_arrow()
  wind_arrow <- data.table(
    x = -150,
    y = 6,
    xend = -62,
    yend = -48
  )

  ggplot() +
    geom_segment(
      data = rear_normals,
      aes(x = x, y = y, xend = xend, yend = yend),
      inherit.aes = FALSE,
      colour = COL_NORMAL,
      alpha = 1.0,
      linewidth = 0.85,
      arrow = arrow(length = unit(0.11, "inches"), type = "closed")
    ) +
    geom_polygon(
      data = lateral_cells,
      aes(x = x2, y = y2, group = cell_id),
      fill = COL_SURFACE,
      alpha = 0.74,
      colour = COL_SURFACE_EDGE,
      linewidth = 0.55
    ) +
    geom_polygon(
      data = top_cells,
      aes(x = x2, y = y2, group = cell_id),
      fill = COL_SURFACE,
      alpha = 0.68,
      colour = COL_SURFACE_EDGE,
      linewidth = 0.45
    ) +
    geom_path(
      data = bottom_ring,
      aes(x = x2, y = y2),
      colour = COL_SURFACE_EDGE,
      alpha = 0.45,
      linewidth = 0.65
    ) +
    geom_path(
      data = top_ring,
      aes(x = x2, y = y2),
      colour = COL_SURFACE_EDGE,
      alpha = 0.75,
      linewidth = 0.8
    ) +
    geom_polygon(
      data = lake,
      aes(x = x2, y = y2),
      fill = COL_LAKE,
      alpha = 0.58,
      colour = COL_LAKE_DARK,
      linewidth = 0.75
    ) +
    geom_segment(
      data = front_normals,
      aes(x = x, y = y, xend = xend, yend = yend),
      inherit.aes = FALSE,
      colour = COL_NORMAL,
      alpha = 1.0,
      linewidth = 0.85,
      arrow = arrow(length = unit(0.11, "inches"), type = "closed")
    ) +
    geom_segment(
      data = top_normal,
      aes(x = x, y = y, xend = xend, yend = yend),
      inherit.aes = FALSE,
      colour = COL_WIND,
      alpha = 0.92,
      linewidth = 1.25,
      arrow = arrow(length = unit(0.15, "inches"), type = "closed")
    ) +
    geom_segment(
      data = wind_arrow,
      aes(x = x, y = y, xend = xend, yend = yend),
      inherit.aes = FALSE,
      colour = COL_WIND,
      linewidth = 2.2,
      lineend = "round",
      arrow = arrow(length = unit(0.24, "inches"), type = "closed")
    ) +
    coord_equal(xlim = c(-165, 120), ylim = c(-80, 178), expand = FALSE) +
    transparent_theme()
}


# -----------------------------------------------------------------------------
# 6. Build and save
# -----------------------------------------------------------------------------
azimuthal_resolution_plot <- build_azimuthal_resolution_plot()
isometric_retrieval_plot <- build_isometric_retrieval_plot()
low_resolution_control_surface_plot <- build_low_resolution_control_surface_plot()

save_asset(AZIMUTHAL_RESOLUTION_PNG, azimuthal_resolution_plot, width = 5.2, height = 5.2)
save_asset(ISOMETRIC_RETRIEVAL_PNG, isometric_retrieval_plot, width = 6.2, height = 5.1)
save_asset(LOW_RES_CONTROL_SURFACE_PNG, low_resolution_control_surface_plot, width = 6.4, height = 5.1)

message("Diagram asset generation complete.")

