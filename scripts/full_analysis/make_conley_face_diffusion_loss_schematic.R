#!/usr/bin/env Rscript
###############################################################################
# make_conley_face_diffusion_loss_schematic.R
#
# Conceptual top-down schematic showing how lateral diffusion can move plume
# mass out through non-downwind faces of an octagonal Conley control volume.
#
# Figure ingredients:
#   - 100 m lake at the origin.
#   - Very stable analytical plume, disk-integrated over the lake surface.
#   - Identical 8-face control volumes translated downwind.
#   - Colored face overlays show approximate face-normal plume flux:
#       blue   = plume entering the control volume,
#       green  = plume leaving through the opposite/downwind face,
#       orange = plume leaving through non-opposite faces.
###############################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
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
suppressMessages(source(file.path(proj_root, "R", "02b_gaussian_forward.R")))
suppressMessages(source(file.path(proj_root, "R", "06_visualization.R")))


# -----------------------------------------------------------------------------
# 1. User-facing settings
# -----------------------------------------------------------------------------
LAKE_DIAMETER_M <- 100
LAKE_RADIUS_M <- LAKE_DIAMETER_M / 2
CONTROL_RADIUS_M <- LAKE_RADIUS_M + 10
PASQUILL_CLASS <- "F"
U_REF_M_S <- 2
WD_FROM_DEG <- 270
Z_SLICE_M <- 2
FLUX_MG_M2_D <- 50
PLUME_LENGTH_M <- 1000
OCTAGON_CENTERS_M <- c(0, 200, 400, 600, 800)
PLUME_FIELD_X_RANGE <- c(-80, PLUME_LENGTH_M)
PLUME_FIELD_Y_RANGE <- c(-240, 240)
PLUME_FIELD_NX <- 430L
PLUME_FIELD_NY <- 220L
PLUME_VISIBLE_REL_MIN <- 1e-2
FACE_SAMPLES_PER_FACE <- 90L
FACE_SEGMENT_MIN_REL <- 0.008
FACE_ARROW_MIN_REL <- 0.035
FACE_STRIP_MIN_WIDTH_M <- 1.4
FACE_STRIP_MAX_WIDTH_M <- 9.5
FACE_COLORS <- c(
  inflow = "#1C5D99",
  angled_inflow = "#7B4AB2",
  opposite_outflow = "#238B45",
  lateral_outflow = "#C94F1A",
  zero_flux = "#6F6F6F"
)
BAR_CENTERS_M <- OCTAGON_CENTERS_M[2:5]
BAR_BASE_Y <- -170
BAR_MAX_HEIGHT_M <- 48
BAR_WIDTH_M <- 22
BAR_GAP_M <- 8
WIND_ARROW_Y <- 124
WIND_LABEL_Y <- 134

OUTPUT_DIR <- proj_path("output")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
FIG_PNG <- file.path(OUTPUT_DIR, "conley_face_diffusion_loss_schematic.png")


# -----------------------------------------------------------------------------
# 2. Geometry helpers
# -----------------------------------------------------------------------------
circle_path <- function(center_x, center_y, radius_m, n = 240L) {
  theta <- seq(0, 2 * pi, length.out = n)
  data.table(
    x = center_x + radius_m * cos(theta),
    y = center_y + radius_m * sin(theta)
  )
}

octagon_path <- function(center_x, center_y = 0, radius_m = CONTROL_RADIUS_M) {
  # Vertex offset puts flat faces on the upwind and downwind sides.
  theta <- seq(22.5, 22.5 + 360, by = 45) * pi / 180
  data.table(
    x = center_x + radius_m * cos(theta),
    y = center_y + radius_m * sin(theta)
  )
}

octagon_faces <- function(center_x, center_y = 0, radius_m = CONTROL_RADIUS_M) {
  verts <- octagon_path(center_x, center_y, radius_m)
  verts <- verts[1:8]
  nxt <- verts[c(2:8, 1)]
  dx <- nxt$x - verts$x
  dy <- nxt$y - verts$y
  len <- sqrt(dx^2 + dy^2)
  nx <- dy / len
  ny <- -dx / len
  angle_deg <- acos(pmax(pmin(nx, 1), -1)) * 180 / pi

  data.table(
    x = (verts$x + nxt$x) / 2,
    y = (verts$y + nxt$y) / 2,
    x_start = verts$x,
    y_start = verts$y,
    x_end = nxt$x,
    y_end = nxt$y,
    nx = nx,
    ny = ny,
    normal_angle_deg = angle_deg
  )
}

wind_to_unit <- function(wd_from_deg) {
  wd_to_rad <- ((wd_from_deg + 180) %% 360) * pi / 180
  c(dx = sin(wd_to_rad), dy = cos(wd_to_rad))
}

face_role <- function(signed_flux, normal_dot_wind) {
  ifelse(
    signed_flux < 0,
    ifelse(normal_dot_wind < -0.90, "inflow", "angled_inflow"),
    ifelse(
      normal_dot_wind > 0.90,
      "opposite_outflow",
      "lateral_outflow"
    )
  )
}

build_face_samples <- function(centers_m,
                               n_per_face = FACE_SAMPLES_PER_FACE,
                               radius_m = CONTROL_RADIUS_M) {
  rbindlist(lapply(centers_m, function(cx) {
    faces <- octagon_faces(cx, radius_m = radius_m)
    rbindlist(lapply(seq_len(nrow(faces)), function(i) {
      f <- faces[i]
      s <- seq(0, 1, length.out = n_per_face)
      data.table(
        center_x = cx,
        face_index = i,
        point_index = seq_along(s),
        x = f$x_start + s * (f$x_end - f$x_start),
        y = f$y_start + s * (f$y_end - f$y_start),
        nx = f$nx,
        ny = f$ny,
        face_length = sqrt((f$x_end - f$x_start)^2 + (f$y_end - f$y_start)^2)
      )
    }))
  }))
}

build_face_strip_polygons <- function(segments_dt) {
  if (!nrow(segments_dt)) return(data.table())

  dt <- copy(segments_dt)
  dt[, segment_group := .I]
  dt[, strip_width := FACE_STRIP_MIN_WIDTH_M +
       (FACE_STRIP_MAX_WIDTH_M - FACE_STRIP_MIN_WIDTH_M) * sqrt(flux_rel)]
  dt[, strip_side := fifelse(role %in% c("inflow", "angled_inflow"), -1, 1)]
  dt[, `:=`(
    ox = strip_side * nx * strip_width,
    oy = strip_side * ny * strip_width
  )]

  rbindlist(lapply(seq_len(nrow(dt)), function(i) {
    s <- dt[i]
    data.table(
      segment_group = s$segment_group,
      role = s$role,
      x = c(s$x, s$xend, s$xend + s$ox, s$x + s$ox),
      y = c(s$y, s$yend, s$yend + s$oy, s$y + s$oy)
    )
  }))
}

line_rect_segment <- function(b, xmin, xmax, ymin, ymax) {
  pts <- data.table(x = numeric(), y = numeric())

  y_left <- xmin + b
  if (is.finite(y_left) && y_left >= ymin && y_left <= ymax) {
    pts <- rbind(pts, data.table(x = xmin, y = y_left))
  }
  y_right <- xmax + b
  if (is.finite(y_right) && y_right >= ymin && y_right <= ymax) {
    pts <- rbind(pts, data.table(x = xmax, y = y_right))
  }
  x_bottom <- ymin - b
  if (is.finite(x_bottom) && x_bottom >= xmin && x_bottom <= xmax) {
    pts <- rbind(pts, data.table(x = x_bottom, y = ymin))
  }
  x_top <- ymax - b
  if (is.finite(x_top) && x_top >= xmin && x_top <= xmax) {
    pts <- rbind(pts, data.table(x = x_top, y = ymax))
  }

  pts <- unique(pts)
  if (nrow(pts) < 2L) return(NULL)
  data.table(
    x = pts$x[1],
    y = pts$y[1],
    xend = pts$x[2],
    yend = pts$y[2]
  )
}

bar_stripe_polygon <- function(b1, b2, xmin, xmax, ymin, ymax) {
  pts <- data.table(x = numeric(), y = numeric())

  add_boundary_points <- function(b) {
    out <- data.table(x = numeric(), y = numeric())

    y_left <- xmin + b
    if (is.finite(y_left) && y_left >= ymin && y_left <= ymax) {
      out <- rbind(out, data.table(x = xmin, y = y_left))
    }
    y_right <- xmax + b
    if (is.finite(y_right) && y_right >= ymin && y_right <= ymax) {
      out <- rbind(out, data.table(x = xmax, y = y_right))
    }
    x_bottom <- ymin - b
    if (is.finite(x_bottom) && x_bottom >= xmin && x_bottom <= xmax) {
      out <- rbind(out, data.table(x = x_bottom, y = ymin))
    }
    x_top <- ymax - b
    if (is.finite(x_top) && x_top >= xmin && x_top <= xmax) {
      out <- rbind(out, data.table(x = x_top, y = ymax))
    }

    unique(out)
  }

  pts <- rbind(pts, add_boundary_points(b1))
  pts <- rbind(pts, add_boundary_points(b2))

  corners <- data.table(
    x = c(xmin, xmax, xmax, xmin),
    y = c(ymin, ymin, ymax, ymax)
  )
  corner_s <- corners$y - corners$x
  keep_corners <- corner_s >= min(b1, b2) & corner_s <= max(b1, b2)
  pts <- rbind(pts, corners[keep_corners])
  pts <- unique(pts)

  if (nrow(pts) < 3L) return(NULL)

  cx <- mean(pts$x)
  cy <- mean(pts$y)
  pts[, angle := atan2(y - cy, x - cx)]
  setorder(pts, angle)
  pts[, angle := NULL]
  pts
}

build_bar_stripes <- function(bars_dt, spacing_m = 5.5) {
  if (!nrow(bars_dt)) return(data.table())

  rbindlist(lapply(seq_len(nrow(bars_dt)), function(i) {
    bar <- bars_dt[i]
    colours <- if (bar$bar_type == "input") {
      c("inflow", "angled_inflow")
    } else {
      c("opposite_outflow", "lateral_outflow")
    }
    if (!is.finite(bar$ymax) || bar$ymax <= bar$ymin) return(NULL)

    b_vals <- seq(bar$ymin - bar$xmax - spacing_m,
                  bar$ymax - bar$xmin + spacing_m,
                  by = spacing_m)
    stripes <- rbindlist(lapply(seq_len(length(b_vals) - 1L), function(j) {
      poly <- bar_stripe_polygon(
        b1 = b_vals[j],
        b2 = b_vals[j + 1L],
        xmin = bar$xmin,
        xmax = bar$xmax,
        ymin = bar$ymin,
        ymax = bar$ymax
      )
      if (is.null(poly)) return(NULL)
      poly[, `:=`(
        stripe_role = colours[((j - 1L) %% length(colours)) + 1L],
        stripe_fill = FACE_COLORS[colours[((j - 1L) %% length(colours)) + 1L]],
        stripe_group = sprintf("%s_%s_%02d", bar$center_x, bar$bar_type, j)
      )]
      poly
    }), fill = TRUE)
    stripes
  }), fill = TRUE)
}

build_flux_bar_data <- function(face_totals_dt,
                                centers_m = BAR_CENTERS_M) {
  totals <- face_totals_dt[
    center_x %in% centers_m & abs(normal_dot_wind) > 1e-6,
    .(
      input = sum(pmax(-signed_flux_index, 0), na.rm = TRUE),
      output = sum(pmax(signed_flux_index, 0), na.rm = TRUE)
    ),
    by = center_x
  ]

  bars <- melt(
    totals,
    id.vars = "center_x",
    measure.vars = c("input", "output"),
    variable.name = "bar_type",
    value.name = "flux_index"
  )

  max_flux <- max(bars$flux_index, na.rm = TRUE)
  if (is.finite(max_flux) && max_flux > 0) {
    bars[, flux_rel := flux_index / max_flux]
  } else {
    bars[, flux_rel := 0]
  }
  bars[, value_label := sprintf("%d", as.integer(round(100 * flux_rel)))]
  bars[, bar_order := fifelse(bar_type == "input", -0.5, 0.5)]
  bars[, `:=`(
    xmid = center_x + bar_order * (BAR_WIDTH_M + BAR_GAP_M),
    ymin = BAR_BASE_Y,
    ymax = BAR_BASE_Y + BAR_MAX_HEIGHT_M * flux_rel,
    label = fifelse(bar_type == "input", "in", "out")
  )]
  bars[, `:=`(
    xmin = xmid - BAR_WIDTH_M / 2,
    xmax = xmid + BAR_WIDTH_M / 2,
    y_value = pmax(ymax + 8, BAR_BASE_Y + 11),
    y_label = BAR_BASE_Y - 7
  )]

  units <- unique(bars[, .(center_x)])
  units[, `:=`(
    x = center_x,
    y = BAR_BASE_Y - 18,
    label = "relative flux (%)"
  )]
  stripes <- build_bar_stripes(bars)

  list(bars = bars[], units = units[], stripes = stripes[])
}

build_face_flux_layers <- function(samples_dt) {
  wind <- wind_to_unit(WD_FROM_DEG)
  samples <- copy(samples_dt)
  samples[, normal_dot_wind := nx * wind[["dx"]] + ny * wind[["dy"]]]

  conc <- gaussian_plume_disk(
    rec_x = samples$x,
    rec_y = samples$y,
    rec_z = rep(Z_SLICE_M, nrow(samples)),
    source_radius = LAKE_RADIUS_M,
    flux_ug_m2_s = flux_mg_d_to_ug_s(FLUX_MG_M2_D),
    U_ref = U_REF_M_S,
    wd_deg = WD_FROM_DEG,
    pasquill_class = PASQUILL_CLASS,
    n_rings = 18L,
    n_src_az = 40L
  )
  samples[, conc_ppb := ug_m3_to_ppb(conc)]

  setorder(samples, center_x, face_index, point_index)
  samples[, `:=`(
    xend = shift(x, type = "lead"),
    yend = shift(y, type = "lead"),
    conc_next_ppb = shift(conc_ppb, type = "lead")
  ), by = .(center_x, face_index)]

  segments <- samples[!is.na(xend)]
  segments[, `:=`(
    conc_mid_ppb = 0.5 * (conc_ppb + conc_next_ppb),
    ds = face_length / (FACE_SAMPLES_PER_FACE - 1L)
  )]
  segments[, signed_flux_index := conc_mid_ppb * U_REF_M_S * normal_dot_wind * ds]
  segments[, abs_flux_index := abs(signed_flux_index)]

  max_seg_flux <- max(segments$abs_flux_index, na.rm = TRUE)
  if (is.finite(max_seg_flux) && max_seg_flux > 0) {
    segments[, flux_rel := abs_flux_index / max_seg_flux]
  } else {
    segments[, flux_rel := 0]
  }
  segments <- segments[
    flux_rel >= FACE_SEGMENT_MIN_REL &
      abs(normal_dot_wind) > 1e-6 &
      conc_mid_ppb > 0
  ]
  segments[, role := face_role(signed_flux_index, normal_dot_wind)]
  strip_polygons <- build_face_strip_polygons(segments)

  face_totals <- samples[!is.na(xend), .(
    signed_flux_index = sum(
      0.5 * (conc_ppb + conc_next_ppb) *
        U_REF_M_S * normal_dot_wind *
        face_length / (FACE_SAMPLES_PER_FACE - 1L),
      na.rm = TRUE
    ),
    x = weighted.mean(
      0.5 * (x + xend),
      w = pmax(0.5 * (conc_ppb + conc_next_ppb), 1e-12),
      na.rm = TRUE
    ),
    y = weighted.mean(
      0.5 * (y + yend),
      w = pmax(0.5 * (conc_ppb + conc_next_ppb), 1e-12),
      na.rm = TRUE
    ),
    x_geom = mean(0.5 * (x + xend), na.rm = TRUE),
    y_geom = mean(0.5 * (y + yend), na.rm = TRUE),
    nx = nx[1],
    ny = ny[1],
    normal_dot_wind = normal_dot_wind[1]
  ), by = .(center_x, face_index)]
  face_totals[, abs_flux_index := abs(signed_flux_index)]
  max_face_flux <- max(face_totals$abs_flux_index, na.rm = TRUE)
  if (is.finite(max_face_flux) && max_face_flux > 0) {
    face_totals[, flux_rel := abs_flux_index / max_face_flux]
  } else {
    face_totals[, flux_rel := 0]
  }
  face_totals_all <- copy(face_totals)

  zero_arrows <- face_totals[abs(normal_dot_wind) <= 1e-6]
  zero_arrows[, `:=`(x = x_geom, y = y_geom)]
  zero_arrows[, `:=`(
    role = "zero_flux",
    arrow_len = 14
  )]
  zero_arrows[, `:=`(
    x_start = x,
    y_start = y,
    x_end = x + arrow_len * nx,
    y_end = y + arrow_len * ny
  )]
  zero_arrows[, `:=`(
    label = "0",
    x_label = x_end + 4 * nx,
    y_label = y_end + 4 * ny
  )]

  face_totals <- face_totals[
    flux_rel >= FACE_ARROW_MIN_REL &
      abs(normal_dot_wind) > 1e-6
  ]
  face_totals[, role := face_role(signed_flux_index, normal_dot_wind)]
  face_totals[, `:=`(
    arrow_dir_x = sign(signed_flux_index) * nx,
    arrow_dir_y = sign(signed_flux_index) * ny,
    arrow_len = 6 + 26 * sqrt(flux_rel)
  )]
  face_totals[signed_flux_index > 0, `:=`(
    x_start = x,
    y_start = y,
    x_end = x + arrow_len * arrow_dir_x,
    y_end = y + arrow_len * arrow_dir_y
  )]
  face_totals[signed_flux_index < 0, `:=`(
    x_start = x - arrow_len * arrow_dir_x,
    y_start = y - arrow_len * arrow_dir_y,
    x_end = x,
    y_end = y
  )]

  list(
    segments = segments,
    strips = strip_polygons,
    arrows = face_totals,
    zero_arrows = zero_arrows,
    face_totals = face_totals_all
  )
}


# -----------------------------------------------------------------------------
# 3. Build schematic data
# -----------------------------------------------------------------------------
plume_grid <- gaussian_plume_grid(
  source_radius = LAKE_RADIUS_M,
  flux_mg_m2_d = FLUX_MG_M2_D,
  U_ref = U_REF_M_S,
  wd_deg = WD_FROM_DEG,
  pasquill_class = PASQUILL_CLASS,
  plane = "xy",
  z_slice_m = Z_SLICE_M,
  x_range = PLUME_FIELD_X_RANGE,
  y_range = PLUME_FIELD_Y_RANGE,
  nx = PLUME_FIELD_NX,
  ny = PLUME_FIELD_NY,
  n_rings = 18L,
  n_src_az = 40L
)
plume_grid[, conc_rel := conc_ppb / max(conc_ppb, na.rm = TRUE)]
plume_grid[, conc_log_rel := fifelse(
  conc_rel >= PLUME_VISIBLE_REL_MIN,
  log10(pmax(conc_rel, PLUME_VISIBLE_REL_MIN)),
  NA_real_
)]
plume_palette <- grDevices::colorRampPalette(
  c("#E7F3F2", "#B9DAD8", "#62AAA9", "#0D5967")
)(256L)
plume_grid[, plume_fill_scale := pmax(
  0,
  pmin(
    1,
    (conc_log_rel - log10(PLUME_VISIBLE_REL_MIN)) /
      (0 - log10(PLUME_VISIBLE_REL_MIN))
  )
)]
plume_grid[, plume_fill := plume_palette[
  1L + as.integer(floor(plume_fill_scale * 255L))
]]
plume_plot <- plume_grid[conc_rel >= PLUME_VISIBLE_REL_MIN]

lake_dt <- circle_path(0, 0, LAKE_RADIUS_M)
octagon_dt <- rbindlist(lapply(OCTAGON_CENTERS_M, function(cx) {
  dt <- octagon_path(cx)
  dt[, group := sprintf("oct_%g", cx)]
  dt
}))
face_samples <- build_face_samples(OCTAGON_CENTERS_M)
face_flux <- build_face_flux_layers(face_samples)
face_strip_dt <- face_flux$strips
face_arrow_dt <- face_flux$arrows
face_zero_arrow_dt <- face_flux$zero_arrows
flux_bar <- build_flux_bar_data(face_flux$face_totals)
flux_bar_dt <- flux_bar$bars
flux_bar_unit_dt <- flux_bar$units
flux_bar_stripe_dt <- flux_bar$stripes


# -----------------------------------------------------------------------------
# 4. Plot
# -----------------------------------------------------------------------------
p <- ggplot() +
  geom_raster(
    data = plume_plot,
    aes(x = x_down, y = y_cross, fill = plume_fill),
    alpha = 0.52,
    interpolate = TRUE
  ) +
  geom_rect(
    data = flux_bar_dt,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = NA, colour = "#4A4A4A", linewidth = 0.25
  ) +
  geom_polygon(
    data = flux_bar_stripe_dt,
    aes(x = x, y = y, group = stripe_group, fill = stripe_fill),
    colour = NA,
    alpha = 0.95
  ) +
  geom_rect(
    data = flux_bar_dt,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = NA, colour = "#4A4A4A", linewidth = 0.25
  ) +
  geom_text(
    data = flux_bar_dt,
    aes(x = xmid, y = y_value, label = value_label),
    size = 3.9, colour = "#252525"
  ) +
  geom_text(
    data = flux_bar_dt,
    aes(x = xmid, y = y_label, label = label),
    size = 3.5, colour = "#4A4A4A"
  ) +
  geom_text(
    data = flux_bar_unit_dt,
    aes(x = x, y = y, label = label),
    size = 3.2, colour = "#4A4A4A"
  ) +
  geom_polygon(
    data = lake_dt,
    aes(x = x, y = y),
    fill = "#B7D4E8", colour = "#2A5B75", linewidth = 0.7
  ) +
  geom_polygon(
    data = face_strip_dt[role == "inflow"],
    aes(x = x, y = y, group = segment_group),
    fill = FACE_COLORS[["inflow"]], alpha = 0.92, colour = NA
  ) +
  geom_polygon(
    data = face_strip_dt[role == "angled_inflow"],
    aes(x = x, y = y, group = segment_group),
    fill = FACE_COLORS[["angled_inflow"]], alpha = 0.92, colour = NA
  ) +
  geom_polygon(
    data = face_strip_dt[role == "opposite_outflow"],
    aes(x = x, y = y, group = segment_group),
    fill = FACE_COLORS[["opposite_outflow"]], alpha = 0.92, colour = NA
  ) +
  geom_polygon(
    data = face_strip_dt[role == "lateral_outflow"],
    aes(x = x, y = y, group = segment_group),
    fill = FACE_COLORS[["lateral_outflow"]], alpha = 0.92, colour = NA
  ) +
  geom_path(
    data = octagon_dt,
    aes(x = x, y = y, group = group),
    colour = "#252525", linewidth = 0.75
  ) +
  geom_segment(
    data = face_arrow_dt,
    aes(
      x = x_start, y = y_start, xend = x_end, yend = y_end,
      colour = role
    ),
    linewidth = 0.7,
    arrow = grid::arrow(length = grid::unit(0.08, "inches"), type = "closed")
  ) +
  geom_segment(
    data = face_zero_arrow_dt,
    aes(
      x = x_start, y = y_start, xend = x_end, yend = y_end,
      colour = role
    ),
    linewidth = 0.45,
    arrow = grid::arrow(length = grid::unit(0.06, "inches"), type = "closed")
  ) +
  geom_text(
    data = face_zero_arrow_dt,
    aes(x = x_label, y = y_label, label = label),
    colour = FACE_COLORS[["zero_flux"]],
    size = 2.7
  ) +
  geom_segment(
    aes(x = -60, xend = 85, y = WIND_ARROW_Y, yend = WIND_ARROW_Y),
    colour = "#252525", linewidth = 0.7,
    arrow = grid::arrow(length = grid::unit(0.12, "inches"), type = "closed")
  ) +
  geom_text(
    aes(x = 13, y = WIND_LABEL_Y, label = sprintf("wind speed = %.1f m/s", U_REF_M_S)),
    size = 3.0, colour = "#252525"
  ) +
  scale_fill_identity(guide = "none") +
  scale_colour_manual(values = FACE_COLORS, guide = "none") +
  scale_x_continuous(
    breaks = seq(0, 1000, by = 200),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_continuous(
    breaks = seq(-150, 100, by = 50),
    expand = expansion(mult = c(0, 0))
  ) +
  coord_equal(xlim = c(-90, 1050), ylim = c(-205, 145), clip = "off") +
  labs(
    title = "Loss from diffusion of contaminant plume",
    x = "Downwind distance from lake center (m)",
    y = "Crosswind distance (m)"
  ) +
  theme_osse(base_size = 11) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    plot.margin = margin(8, 18, 8, 8)
  )

ggsave_osse(FIG_PNG, p, width = 16, height = 5.4, dpi = 240)
message(sprintf("Saved figure: %s", normalizePath(FIG_PNG, mustWork = FALSE)))

