###############################################################################
# utils.R - Shared helpers, configuration objects, and geometry utilities
#
# This file supports an Observing System Simulation Experiment (OSSE) that
# evaluates UAS-based cylindrical mass-balance retrievals for methane-emitting
# area sources using a Gaussian-plume forward model.
#
# Scientific notes:
#   - Ebullition is treated as part of an effective, lake-averaged surface
#     flux.  We do NOT resolve hotspots or fine-scale spatial structure.
#   - The OSSE is a publication companion for the Gaussian-forward /
#     Conley-inverse manuscript experiments.
###############################################################################

suppressPackageStartupMessages({
  library(data.table)
})

# -----------------------------------------------------------------------------
# Emissions catalog caching (synthetic lake-wide emission fields)
# Bump this when changing how generate_effective_flux_field() (or its inputs)
# affects the generated emissions catalogs.
EMISSIONS_CATALOG_SOFTWARE_VERSION <- "1.0.1"

# -----------------------------------------------------------------------------
# 1.  Lake geometry
# -----------------------------------------------------------------------------
#' Build a square grid covering a circular lake.
#'
#' @param diameter_m  Lake diameter in metres.
#' @param resolution_m  Grid cell size in metres (default 10).
#' @return A list with:
#'   $grid - data.table with columns: x, y, dist_from_center, inside_lake.
#'   $radius - lake radius in metres.
#'   $diameter_m - lake diameter in metres.
define_lake_geometry <- function(diameter_m, resolution_m = 10) {
  radius <- diameter_m / 2
  half_extent <- radius + resolution_m
  seq_vals <- seq(-half_extent, half_extent, by = resolution_m)
  grid <- CJ(x = seq_vals, y = seq_vals)
  grid[, dist_from_center := sqrt(x^2 + y^2)]
  grid[, inside_lake := dist_from_center <= radius]

  list(grid = grid, radius = radius, diameter_m = diameter_m)
}


# -----------------------------------------------------------------------------
# 2.  Effective flux field generator
# -----------------------------------------------------------------------------
#' Generate an effective CH4 surface flux field (mg CH4 m-2 d-1).
#'
#' Returns a constant flux field: same value everywhere on the lake, at all
#' times. No littoral enhancement, spatial noise, or diel variation.
#'
#' @param lake_geom  List from \code{define_lake_geometry} ($grid, $radius).
#' @param base_flux  Lake-mean flux (mg CH4 m-2 d-1), default 200.
#' @param littoral_factor  Unused; retained for API compatibility.
#' @param spatial_corr_length_m  Unused; retained for API compatibility.
#' @param diel_amplitude  Unused; retained for API compatibility.
#' @param diel_phase_hour  Unused; retained for API compatibility.
#' @param times_hour  Numeric vector of hours-of-day at which to evaluate
#'                    (default 0:23).
#' @param seed  Unused; retained for API compatibility.
#' @return A data.table with columns: x, y, hour, flux_mg_m2_d.
generate_effective_flux_field <- function(lake_geom,
                                         base_flux = 200,
                                         littoral_factor = 2.5,
                                         spatial_corr_length_m = 75,
                                         diel_amplitude = 0.30,
                                         diel_phase_hour = 14,
                                         times_hour = 0:23,
                                         seed = NULL) {

  lake <- lake_geom$grid[inside_lake == TRUE]

  out_list <- lapply(times_hour, function(h) {
    data.table(
      x = lake$x,
      y = lake$y,
      hour = h,
      flux_mg_m2_d = base_flux
    )
  })
  rbindlist(out_list)
}


# -----------------------------------------------------------------------------
# 3.  Meteorology time series generator
# -----------------------------------------------------------------------------
#' Generate a synthetic meteorology time series.
#'
#' Uses a small number of archetypal stability regimes rather than a full
#' continuum, because this is a capability study.
#'
#' Regimes:
#'   - Unstable midday:   u* ~ 0.30, L ~ -30
#'   - Near-neutral afternoon: u* ~ 0.25, L ~ -500
#'   - Stable evening:    u* ~ 0.10, L ~ +50
#'
#' Wind speed is derived from u* via log-law for physical consistency:
#'   ws = u* / kappa * ln(z_ref / z0)
#'
#' @param duration_days  Length of simulation (days).
#' @param dt_minutes     Time step in minutes (default 30).
#' @param start_hour     Hour-of-day at which the time series begins (default 0).
#' @param seed           Random seed (default 123).
#' @return A data.table with columns: time_min, hour_of_day, ws_m_s, wd_deg,
#'         ustar_m_s, L_m, z0_m, stability_class.
generate_meteorology_timeseries <- function(duration_days = 1,
                                            dt_minutes = 30,
                                            start_hour = 0,
                                            wd_base_deg = 225,
                                            wd_noise_sd_deg = 15,
                                            wd_diurnal_amp_deg = 0,
                                            ws_reference_height_m = 10,
                                            ws_reference_speed_m_s = NULL,
                                            ws_scale_factor = 1,
                                            ws_noise_sd_m_s = 0,
                                            force_stability = NULL,
                                            seed = 123) {
  n_steps <- as.integer(duration_days * 24 * 60 / dt_minutes)
  if (n_steps < 1) n_steps <- 1L
  time_min <- seq(0, by = dt_minutes, length.out = n_steps)
  hour_of_day <- (start_hour + time_min / 60) %% 24

  regime <- character(n_steps)
  ustar <- numeric(n_steps)
  L_obukhov <- numeric(n_steps)

  set.seed(seed)

  z_ref <- ws_reference_height_m   # reference height for wind speed (m)
  kappa <- 0.4  # von Karman constant

  for (i in seq_len(n_steps)) {
    h <- hour_of_day[i]
    stab <- if (!is.null(force_stability)) {
      force_stability
    } else if (h >= 10 & h < 15) {
      "unstable"
    } else if (h >= 15 & h < 19) {
      "near_neutral"
    } else {
      "stable"
    }
    if (stab == "unstable") {
      regime[i]    <- "unstable"
      ustar[i]     <- 0.30 + rnorm(1, 0, 0.03)
      L_obukhov[i] <- -30 + rnorm(1, 0, 10)
    } else if (stab == "near_neutral") {
      regime[i]    <- "near_neutral"
      ustar[i]     <- 0.25 + rnorm(1, 0, 0.03)
      L_obukhov[i] <- -500 + rnorm(1, 0, 100)
    } else {
      regime[i]    <- "stable"
      ustar[i]     <- 0.10 + rnorm(1, 0, 0.02)
      L_obukhov[i] <- 50 + rnorm(1, 0, 15)
    }
  }

  ustar <- pmax(ustar, 0.03)

  # Charnock relation for open water: z0 = alpha_c * u*^2 / g
  alpha_c <- 0.015
  g <- 9.81
  z0 <- pmax(alpha_c * ustar^2 / g, 1e-5)

  # Derive wind speed from u* via log-law (physically consistent)
  ws <- (ustar / kappa) * log(z_ref / z0)
  # Optional scenario-level target at the reference height. If provided,
  # normalize the profile-derived speeds to this target mean first.
  if (!is.null(ws_reference_speed_m_s) && is.finite(ws_reference_speed_m_s)) {
    ws_mean <- mean(ws, na.rm = TRUE)
    if (is.finite(ws_mean) && ws_mean > 0) {
      ws <- ws * (ws_reference_speed_m_s / ws_mean)
    }
  }
  ws <- ws * ws_scale_factor
  if (ws_noise_sd_m_s > 0) {
    ws <- ws + rnorm(n_steps, 0, ws_noise_sd_m_s)
  }
  ws <- pmax(ws, 0.3)

  # Wind direction truth: by default no diurnal swing; optional amplitude can
  # be supplied by MET_SCENARIOS for experiments.
  base_wd <- wd_base_deg + wd_diurnal_amp_deg * sin(2 * pi * hour_of_day / 24)
  wd <- (base_wd + rnorm(n_steps, 0, wd_noise_sd_deg)) %% 360

  data.table(
    time_min        = time_min,
    hour_of_day     = round(hour_of_day, 2),
    ws_m_s          = round(ws, 3),
    wd_deg          = round(wd, 1),
    ustar_m_s       = round(ustar, 4),
    L_m             = round(L_obukhov, 1),
    z0_m            = round(z0, 6),
    stability_class = regime
  )
}


# -----------------------------------------------------------------------------
# 3b.  Gaussian plume meteorology time series
# -----------------------------------------------------------------------------
#' Generate a Pasquill-class-based met timeseries for the Gaussian plume model.
#'
#' Independent of the MOST-based \code{generate_meteorology_timeseries()}.
#' Uses a fixed Pasquill class and reference wind speed rather than u*, L, z0.
#'
#' @param duration_days  Length of simulation (days).
#' @param dt_minutes  Time step in minutes (default 30).
#' @param start_hour  Hour-of-day at which the time series begins (default 12).
#' @param U_ref_m_s  Reference wind speed at z_ref (m/s).
#' @param wd_base_deg  Base wind-from direction (degrees, met convention).
#' @param wd_noise_sd_deg  SD of wind direction noise (default 0).
#' @param pasquill_class  Pasquill class "A" through "F".
#' @param seed  Random seed (default 123).
#' @return data.table with columns: time_min, hour_of_day, ws_m_s, wd_deg,
#'   ustar_m_s (NA), L_m (NA), z0_m (0.001), stability_class, pasquill_class.
generate_gaussian_met_timeseries <- function(duration_days   = 1,
                                             dt_minutes      = 30,
                                             start_hour      = 12,
                                             U_ref_m_s       = 5,
                                             wd_base_deg     = 225,
                                             wd_noise_sd_deg = 0,
                                             pasquill_class  = "D",
                                             seed            = 123) {
  n_steps <- as.integer(duration_days * 24 * 60 / dt_minutes)
  if (n_steps < 1) n_steps <- 1L
  time_min    <- seq(0, by = dt_minutes, length.out = n_steps)
  hour_of_day <- (start_hour + time_min / 60) %% 24

  set.seed(seed)
  wd <- (wd_base_deg + rnorm(n_steps, 0, wd_noise_sd_deg)) %% 360

  stab_map <- c(A = "unstable", B = "unstable", C = "unstable",
                D = "near_neutral", E = "stable", F = "stable")

  data.table(
    time_min        = time_min,
    hour_of_day     = round(hour_of_day, 2),
    ws_m_s          = round(U_ref_m_s, 3),
    wd_deg          = round(wd, 1),
    ustar_m_s       = NA_real_,
    L_m             = NA_real_,
    z0_m            = 0.001,
    stability_class = unname(stab_map[pasquill_class]),
    pasquill_class  = pasquill_class
  )
}


# -----------------------------------------------------------------------------
# 4.  Receptor (UAV measurement point) configuration
# -----------------------------------------------------------------------------
#' Define receptor coordinates: outer ring plus optional inner rings over the lake.
#'
#' Outer ring is at lake edge + distance_from_edge_m. Optional inner rings
#' are at fractions of lake radius (directly over the lake). Optional
#' background receptors at high altitude (e.g. 150 m over lake center) provide time-varying
#' C_bg for the Conley inverse.
#'
#' @param lake_diameter  Lake diameter in metres.
#' @param n_azimuth  Number of receptor positions around each ring.
#' @param heights_m  Numeric vector of receptor heights above the surface (e.g. c(3, 8, 15)).
#' @param distance_from_edge_m  Radial distance outside the lake edge for outer ring (default 10).
#' @param inner_ring_fracs  Optional numeric vector of radius fractions (0–1) for rings over
#'        the lake, e.g. c(0.5, 0.25) for two inner loops at half and quarter radius.
#' @param background_heights_m  Optional numeric vector of heights (m) for background
#'        receptors. Placed at lake center (0, 0) by default. Use mean of their
#'        values per time step as C_bg for the Conley inverse.
#' @return A data.table with columns: receptor_id, x, y, z, azimuth_deg, ring_type.
define_receptor_config <- function(lake_diameter,
                                   n_azimuth = 8,
                                   heights_m = c(3, 8, 15),
                                   distance_from_edge_m = 10,
                                   inner_ring_fracs = NULL,
                                   background_heights_m = NULL) {
  lake_radius <- lake_diameter / 2
  angles_deg <- seq(0, 360 - 360/n_azimuth, length.out = n_azimuth)

  # Outer ring (outside lake edge)
  radius_outer <- lake_radius + distance_from_edge_m
  outer <- CJ(azimuth_deg = angles_deg, z = heights_m)
  outer[, `:=`(
    x = radius_outer * sin(azimuth_deg * pi / 180),
    y = radius_outer * cos(azimuth_deg * pi / 180),
    ring_type = "outer"
  )]

  rings_list <- list(outer)

  if (!is.null(inner_ring_fracs) && length(inner_ring_fracs) > 0) {
    for (k in seq_along(inner_ring_fracs)) {
      frac <- inner_ring_fracs[k]
      radius_inner <- frac * lake_radius
      inner_k <- CJ(azimuth_deg = angles_deg, z = heights_m)
      inner_k[, `:=`(
        x = radius_inner * sin(azimuth_deg * pi / 180),
        y = radius_inner * cos(azimuth_deg * pi / 180),
        ring_type = sprintf("inner_%g", frac)
      )]
      rings_list[[length(rings_list) + 1L]] <- inner_k
    }
  }

  if (!is.null(background_heights_m) && length(background_heights_m) > 0) {
    bg <- data.table(
      azimuth_deg = 0,
      z = background_heights_m,
      x = 0,
      y = 0,
      ring_type = "background"
    )
    rings_list[[length(rings_list) + 1L]] <- bg
  }

  receptors <- rbindlist(rings_list)
  receptors[, receptor_id := .I]
  setcolorder(receptors, c("receptor_id", "x", "y", "z", "azimuth_deg", "ring_type"))
  receptors[]
}


#' Define a 2-D transect receptor grid along the downwind centreline.
#'
#' Receptors lie on the downwind centreline at the given azimuth (meteorological
#' convention: 0 = north, 90 = east).  Pass \code{(wd_base_deg + 180) \%\% 360}
#' to align with the actual plume direction.
#' Intended purely for plume-structure visualisation; these receptors are NOT
#' used by the Conley inverse.
#'
#' @param x_m                 Numeric vector of downwind distances from lake centre (m).
#' @param heights_m           Numeric vector of receptor heights above surface (m).
#' @param downwind_azimuth_deg Direction the plume travels toward (degrees, met convention).
#' @return data.table with columns receptor_id, x, y, z, azimuth_deg, ring_type.
define_transect_receptors <- function(x_m,
                                      heights_m,
                                      downwind_azimuth_deg = 45) {
  az_rad <- downwind_azimuth_deg * pi / 180
  grid <- CJ(dist = x_m, z = heights_m)
  grid[, `:=`(
    x           = dist * sin(az_rad),
    y           = dist * cos(az_rad),
    azimuth_deg = downwind_azimuth_deg,
    ring_type   = "transect"
  )]
  grid[, receptor_id := paste0("t_d", dist, "_z", z)]
  grid[, dist := NULL]
  setcolorder(grid, c("receptor_id", "x", "y", "z", "azimuth_deg", "ring_type"))
  grid[]
}


#' Default transect configurations for plume cross-section visualisation.
TRANSECT_CONFIGS <- list(
  default = list(
    name      = "default",
    x_m       = seq(2, 126, by = 2),
    heights_m = c(0.1, 0.2, 0.3, 0.5, 0.75, 1, 1.5, 2, 3, 4, 5, 7),
    description = "Downwind cross-section 2-126 m in 2 m steps, 10 heights. Azimuth set at runtime from met_cfg."
  )
)


# -----------------------------------------------------------------------------
# 5.  Configuration objects
# -----------------------------------------------------------------------------
#' Lake-size configurations.
LAKE_CONFIGS <- list(
  small  = list(name = "small",  diameter_m = 100),
  medium = list(name = "medium", diameter_m = 300),
  large  = list(name = "large",  diameter_m = 800),

  # Contamination experiment: small lake (100 m) with a nearby 50 m lake 100 m
  # away in the upwind direction (angle_deg = NULL resolves to wd_base_deg at
  # run time).  Passes resolve_nearby_lake_position() geometry checks by default
  # (ring radius = 60 m; nearest nearby-lake edge = 100 - 25 = 75 m > 60 m).
  small_with_nearby = list(
    name       = "small_with_nearby",
    diameter_m = 100,
    nearby_lake = list(
      diameter_m   = 50,    # nearby lake diameter (m)
      distance_m   = 100,   # centre-to-centre distance from main lake (m)
      angle_deg    = NULL,  # NULL means upwind = wd_base_deg of the met scenario
      flux_mg_m2_d = NULL   # NULL means same base_flux as the main lake
    )
  ),
  small_with_nearby_far = list(
    name       = "small_with_nearby_far",
    diameter_m = 100,
    nearby_lake = list(
      diameter_m   = 50,    # nearby lake diameter (m)
      distance_m   = 1000,   # centre-to-centre distance from main lake (m)
      angle_deg    = NULL,  # NULL means upwind = wd_base_deg of the met scenario
      flux_mg_m2_d = NULL   # NULL means same base_flux as the main lake
    )
  )
)

#' Sensor specifications.
#' Each entry: name, mass_kg, precision_1s_ppb, precision_opt_ppb,
#'             opt_avg_time_s, drift_ppb_hr, power_W.
SENSOR_CONFIGS <- list(
  mira_strato = list(
    name             = "Aeris MIRA Strato",
    short_name       = "mira_strato",
    mass_kg          = 2.0,
    precision_1s_ppb = 1.0,
    precision_opt_ppb = 0.5,
    opt_avg_time_s   = 5,
    drift_ppb_hr     = 0.3,
    power_W          = 18,
    cal_offset_ppb   = 2.0
  ),
  li7810 = list(
    name             = "LI-COR LI-7810",
    short_name       = "li7810",
    mass_kg          = 10.5,
    precision_1s_ppb = 0.60,
    precision_opt_ppb = 0.25,
    opt_avg_time_s   = 5,
    drift_ppb_hr     = 0.05,
    power_W          = 30,
    cal_offset_ppb   = 0.5
  ),
  ugga = list(
    name             = "ABB LGR UGGA",
    short_name       = "ugga",
    mass_kg          = 17.0,
    precision_1s_ppb = 1.4,
    precision_opt_ppb = 0.5,
    opt_avg_time_s   = 10,
    drift_ppb_hr     = 0.15,
    power_W          = 70,
    cal_offset_ppb   = 1.0
  ),
  empa_qcls = list(
    name             = "Empa compact QCLS",
    short_name       = "empa_qcls",
    mass_kg          = 2.1,
    precision_1s_ppb = 0.5,
    precision_opt_ppb = 0.25,
    opt_avg_time_s   = 5,
    drift_ppb_hr     = 0.10,
    power_W          = 25,
    cal_offset_ppb   = 0.8
  )
)

#' Sampling-design configurations.
#' All designs use 5 height layers (1, 3, 5, 7, 10 m) for the Conley
#' cylindrical control volume. Outer ring only (no inner rings).
#'
#' Background subtraction is removed from the Conley inverse, so we do not
#' include a dedicated high-altitude "background" receptor in the forward runs.
#' Plots show one height layer only for clarity.
SAMPLING_CONFIGS <- list(
  design_8az_5h_30min = list(
    name                 = "8az_5h_30min",
    n_azimuth            = 8,
    heights_m            = c(1, 3, 5, 7, 10),
    cycle_min            = 30,
    inner_ring_fracs     = NULL,
    background_heights_m = NULL
  ),
  design_24az_5h_30min = list(
    name                 = "24az_5h_30min",
    n_azimuth            = 24,
    heights_m            = c(1, 3, 5, 7, 10),
    cycle_min            = 30,
    inner_ring_fracs     = NULL,
    background_heights_m = NULL
  ),
  design_36az_5h_30min = list(
    name                 = "36az_5h_30min",
    n_azimuth            = 36,
    heights_m            = c(1, 3, 5, 7, 10),
    cycle_min            = 30,
    inner_ring_fracs     = NULL,
    background_heights_m = NULL
  ),
  design_24az_5h_60min = list(
    name                 = "24az_5h_60min",
    n_azimuth            = 24,
    heights_m            = c(1, 3, 5, 7, 10),
    cycle_min            = 60,
    inner_ring_fracs     = NULL,
    background_heights_m = NULL
  ),
  design_24az_5h_15min = list(
    name                 = "24az_5h_15min",
    n_azimuth            = 24,
    heights_m            = c(1, 3, 5, 7, 10),
    cycle_min            = 15,
    inner_ring_fracs     = NULL,
    background_heights_m = NULL
  ),
  design_36az_8h_15min = list(
    name                 = "36az_8h_15min",
    n_azimuth            = 36,
    heights_m            = c(0.25, 0.5, 1, 1.5, 2, 3, 4, 5),
    cycle_min            = 15,
    inner_ring_fracs     = NULL,
    background_heights_m = NULL
  ),
  design_96az_8h_15min = list(
    name                 = "96az_8h_15min",
    n_azimuth            = 96,
    heights_m            = c(0.25, 0.5, 1, 1.5, 2, 3, 4, 5),
    cycle_min            = 15,
    inner_ring_fracs     = NULL,
    background_heights_m = NULL
  ),
  design_96az_12h_15min = list(
    name                 = "96az_12h_15min",
    n_azimuth            = 96,
    heights_m            = c(0.1, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 2, 2.5, 3, 4, 5),
    cycle_min            = 15,
    inner_ring_fracs     = NULL,
    background_heights_m = NULL
  ),
  design_96az_15h_15min = list(
    name                 = "96az_15h_15min",
    n_azimuth            = 96,
    heights_m            = c(0.01, 0.025, 0.05, 0.1, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 2, 2.5, 3, 4, 5),
    cycle_min            = 15,
    inner_ring_fracs     = NULL,
    background_heights_m = NULL
  ),
  design_96az_15h_17min_high = list(
    name                 = "96az_15h_17min_high",
    n_azimuth            = 96,
    heights_m            = c(0.1, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 2, 2.5, 3, 4, 5, 6, 7, 8, 9, 10),
    cycle_min            = 15,
    inner_ring_fracs     = NULL,
    background_heights_m = NULL
  ),
  design_96az_many_heights_30min = list(
    name                 = "96az_many_heights_30min",
    n_azimuth            = 96,
    heights_m            = c(0.1, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 2, 2.5, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20),
    cycle_min            = 15,
    inner_ring_fracs     = NULL,
    background_heights_m = NULL
  )
)

#' Canonical meteorology scenarios for capability analysis.
#' start_hour controls the beginning of the diurnal window so that
#' each scenario samples the intended stability regime.
MET_SCENARIOS <- list(
  mixed_day = list(
    name          = "mixed_day",
    duration_days = 1,
    dt_minutes    = 30,
    start_hour    = 0,
    wd_base_deg   = 225,
    wd_noise_sd_deg = 15,
    wd_diurnal_amp_deg = 0,
    ws_reference_height_m = 10,
    ws_reference_speed_m_s = NULL,
    ws_scale_factor = 1,
    ws_noise_sd_m_s = 0,
    description   = "Full diurnal cycle with unstable/neutral/stable regimes"
  ),
  unstable_only = list(
    name          = "unstable_only",
    duration_days = 0.25,
    dt_minutes    = 30,
    start_hour    = 10,
    wd_base_deg   = 225,
    wd_noise_sd_deg = 15,
    wd_diurnal_amp_deg = 0,
    ws_reference_height_m = 10,
    ws_reference_speed_m_s = NULL,
    ws_scale_factor = 1,
    ws_noise_sd_m_s = 0,
    description   = "6 hours of convective midday conditions (10:00-16:00)"
  ),
  stable_only = list(
    name          = "stable_only",
    duration_days = 0.25,
    dt_minutes    = 30,
    start_hour    = 20,
    wd_base_deg   = 225,
    wd_noise_sd_deg = 15,
    wd_diurnal_amp_deg = 0,
    ws_reference_height_m = 10,
    ws_reference_speed_m_s = NULL,
    ws_scale_factor = 1,
    ws_noise_sd_m_s = 0,
    description   = "6 hours of stable evening/night conditions (20:00-02:00)"
  ),
  # Full-day forced-stability scenarios.
  stable_full_day = list(
    name             = "stable_full_day",
    duration_days    = 1,
    dt_minutes       = 30,
    start_hour       = 0,
    wd_base_deg      = 225,
    wd_noise_sd_deg  = 15,
    wd_diurnal_amp_deg = 0,
    ws_reference_height_m  = 10,
    ws_reference_speed_m_s = NULL,
    ws_scale_factor  = 1,
    ws_noise_sd_m_s  = 0,
    force_stability  = "stable",
    description      = "Full 24-h day, always stable (nocturnal boundary layer)"
  ),
  unstable_full_day = list(
    name             = "unstable_full_day",
    duration_days    = 1,
    dt_minutes       = 30,
    start_hour       = 0,
    wd_base_deg      = 225,
    wd_noise_sd_deg  = 15,
    wd_diurnal_amp_deg = 0,
    ws_reference_height_m  = 10,
    ws_reference_speed_m_s = NULL,
    ws_scale_factor  = 1,
    ws_noise_sd_m_s  = 0,
    force_stability  = "unstable",
    description      = "Full 24-h day, always unstable (convective boundary layer)"
  ),
  neutral_full_day = list(
    name             = "neutral_full_day",
    duration_days    = 1,
    dt_minutes       = 30,
    start_hour       = 0,
    wd_base_deg      = 225,
    wd_noise_sd_deg  = 15,
    wd_diurnal_amp_deg = 0,
    ws_reference_height_m  = 10,
    ws_reference_speed_m_s = NULL,
    ws_scale_factor  = 1,
    ws_noise_sd_m_s  = 0,
    force_stability  = "near_neutral",
    description      = "Full 24-h day, always near-neutral stability"
  ),
  # Constant-direction scenario.
  constant_from_north = list(
    name             = "constant_from_north",
    duration_days    = 1,
    dt_minutes       = 30,
    start_hour       = 0,
    wd_base_deg      = 0,
    wd_noise_sd_deg  = 0,
    wd_diurnal_amp_deg = 0,
    ws_reference_height_m  = 10,
    ws_reference_speed_m_s = NULL,
    ws_scale_factor  = 1,
    ws_noise_sd_m_s  = 0,
    force_stability  = "near_neutral",
    description      = "Constant wind from north (0 deg), near-neutral, no direction noise"
  )
)

#' Gaussian plume met scenarios (Pasquill class-based, independent of MOST).
#'
#' Each entry specifies a single Pasquill stability class, reference wind
#' speed at 10 m, and wind direction.  Used with
#' \code{generate_gaussian_met_timeseries()} and the Gaussian plume forward
#' model in \code{02b_gaussian_forward.R}.
GAUSSIAN_MET_SCENARIOS <- list(
  unstable_strong = list(
    name = "gp_unstable_A", duration_days = 1, dt_minutes = 30,
    start_hour = 12, U_ref_m_s = 2, wd_base_deg = 225,
    wd_noise_sd_deg = 0, pasquill_class = "A", z_ref_m = 10,
    description = "Very unstable (A), light wind 2 m/s"
  ),
  unstable_moderate = list(
    name = "gp_unstable_B", duration_days = 1, dt_minutes = 30,
    start_hour = 12, U_ref_m_s = 2, wd_base_deg = 225,
    wd_noise_sd_deg = 0, pasquill_class = "B", z_ref_m = 10,
    description = "Moderately unstable (B), 3 m/s"
  ),
  slightly_unstable = list(
    name = "gp_slightly_unstable_C", duration_days = 1, dt_minutes = 30,
    start_hour = 12, U_ref_m_s = 2, wd_base_deg = 225,
    wd_noise_sd_deg = 0, pasquill_class = "C", z_ref_m = 10,
    description = "Slightly unstable (C), 4 m/s"
  ),
  neutral_moderate = list(
    name = "gp_neutral_D", duration_days = 1, dt_minutes = 30,
    start_hour = 12, U_ref_m_s = 2, wd_base_deg = 225,
    wd_noise_sd_deg = 0, pasquill_class = "D", z_ref_m = 10,
    description = "Neutral (D), moderate wind 5 m/s"
  ),
  stable_light = list(
    name = "gp_stable_E", duration_days = 1, dt_minutes = 30,
    start_hour = 0, U_ref_m_s = 2, wd_base_deg = 225,
    wd_noise_sd_deg = 0, pasquill_class = "E", z_ref_m = 10,
    description = "Slightly stable (E), 3 m/s"
  ),
  stable_calm = list(
    name = "gp_stable_F", duration_days = 1, dt_minutes = 30,
    start_hour = 0, U_ref_m_s = 2, wd_base_deg = 225,
    wd_noise_sd_deg = 0, pasquill_class = "F", z_ref_m = 10,
    description = "Moderately stable (F), calm 1.5 m/s"
  )
)

#' Pasquill-Gifford stability classes with Briggs (1973) rural sigma
#' coefficients and power-law wind profile exponents (Hanna et al. 1982).
#'
#' sigma_y = sy_a * x / sqrt(1 + sy_b * x)
#' sigma_z = sz_a * x * (1 + sz_b * x)^sz_p
#'   A-B: sz_p = 0  (unbounded linear)
#'   C-D: sz_p = -0.5
#'   E-F: sz_p = -1
#' 
PASQUILL_CLASSES <- list(
  A = list(name = "A", description = "Very unstable",
           sy_a = 0.22, sy_b = 0.0001,
           sz_a = 0.20, sz_b = 0,      sz_p = 0,
           wind_exp = 0.07),
  B = list(name = "B", description = "Moderately unstable",
           sy_a = 0.16, sy_b = 0.0001,
           sz_a = 0.12, sz_b = 0,      sz_p = 0,
           wind_exp = 0.07),
  C = list(name = "C", description = "Slightly unstable",
           sy_a = 0.11, sy_b = 0.0001,
           sz_a = 0.08, sz_b = 0.0002, sz_p = -0.5,
           wind_exp = 0.10),
  D = list(name = "D", description = "Neutral",
           sy_a = 0.08, sy_b = 0.0001,
           sz_a = 0.06, sz_b = 0.0015, sz_p = -0.5,
           wind_exp = 0.15),
  E = list(name = "E", description = "Slightly stable",
           sy_a = 0.06, sy_b = 0.0001,
           sz_a = 0.03, sz_b = 0.0003, sz_p = -1,
           wind_exp = 0.35),
  F = list(name = "F", description = "Moderately stable",
           sy_a = 0.04, sy_b = 0.0001,
           sz_a = 0.016, sz_b = 0.0003, sz_p = -1,
           wind_exp = 0.55)
)

#' Emission field parameters.
#' Diel amplitude set to 0 for constant (no diel) flux.
EMISSION_DEFAULTS <- list(
  base_flux_range    = c(100, 400),
  base_flux_default  = 200,
  littoral_factor    = 2.5,
  spatial_corr_length_m = 75,
  diel_amplitude     = 0,
  diel_phase_hour    = 14
)


#' Retrieval configurations for the Conley mass-balance inverse.
#'
#' Controls how wind speed is assigned to each receptor in \code{retrieve_flux}.
#'
#' \strong{per_receptor_wind} (default): Each receptor gets an independent
#'   wind-speed observation drawn from \eqn{\mathcal{N}(ws\_at\_z, \sigma_{obs}^2)}.
#'   Represents a drone-mounted anemometer visiting each position separately.
#'
#' \strong{height_avg_wind}: After per-receptor noisy observations are drawn,
#'   wind speed is averaged across all azimuth positions at each height layer
#'   for each time window.  Reduces effective wind noise by ~1/sqrt(n_azimuth),
#'   equivalent to a single fixed-mast measurement that is representative of
#'   the whole ring at that height.
RETRIEVAL_CONFIGS <- list(
  per_receptor_wind = list(
    name                     = "per_receptor_wind",
    use_observed_ws          = TRUE,
    ws_observation_error_m_s = 0.3,   # 1-sigma wind speed error per receptor (m/s)
    ws_uncertainty_m_s       = 0.3,   # used only when use_observed_ws = FALSE
    wd_uncertainty_deg       = 10,    # 1-sigma wind direction error per time window (deg)
    use_height_avg_wind      = FALSE,
    description = "Independent wind observation per receptor (drone-mounted sensor)."
  ),
  height_avg_wind = list(
    name                     = "height_avg_wind",
    use_observed_ws          = TRUE,
    ws_observation_error_m_s = 0.3,
    ws_uncertainty_m_s       = 0.3,
    wd_uncertainty_deg       = 10,
    use_height_avg_wind      = TRUE,
    description = "Wind averaged across all azimuths per height (single-mast equivalent)."
  )
)


# -----------------------------------------------------------------------------
# 5b. Nearby-lake contamination geometry helper
# -----------------------------------------------------------------------------
#' Resolve nearby-lake position and validate geometry for contamination experiments.
#'
#' The nearby lake centre is placed at \code{distance_m} from the main lake
#' centre, along the bearing \code{angle_deg} (meteorological convention: 0 =
#' north, 90 = east, 180 = south, 270 = west; same as receptor azimuth).
#' When \code{angle_deg = NULL} the bearing defaults to \code{wd_base_deg},
#' placing the nearby lake directly upwind.
#'
#' Two geometry checks are enforced (both raise errors if violated):
#' \enumerate{
#'   \item Lakes must not physically overlap:
#'         \code{distance_m >= main_radius_m + nearby_radius_m}.
#'   \item Nearby lake must lie entirely outside the receptor control-volume ring:
#'         \code{distance_m - nearby_radius_m > ring_radius_m}.
#' }
#'
#' @param nearby_lake_cfg  The \code{nearby_lake} sub-list from a LAKE_CONFIG entry.
#'   Required fields: \code{diameter_m}, \code{distance_m}.
#'   Optional: \code{angle_deg} (NULL = upwind), \code{flux_mg_m2_d} (NULL = main flux).
#' @param main_radius_m   Main lake radius (m).
#' @param ring_radius_m   Outer receptor ring radius (m) = main_radius + d_edge (default 10 m).
#' @param wd_base_deg     Mean wind-FROM direction (degrees) for the met scenario.
#'   Used only when \code{angle_deg = NULL}.
#' @param main_base_flux  Main lake base flux (mg CH4 m-2 d-1); used as default
#'   for the nearby lake when \code{flux_mg_m2_d = NULL}.
#' @return Named list: \code{diameter_m}, \code{radius_m}, \code{center_x},
#'   \code{center_y}, \code{flux_mg_m2_d}, \code{angle_deg}, \code{distance_m}.
resolve_nearby_lake_position <- function(nearby_lake_cfg,
                                         main_radius_m,
                                         ring_radius_m,
                                         wd_base_deg,
                                         main_base_flux) {
  nl        <- nearby_lake_cfg
  nearby_r  <- nl$diameter_m / 2
  dist_m    <- nl$distance_m

  # Resolve angle: NULL means use the mean upwind direction (wd_base_deg = FROM direction).
  angle_deg <- if (is.null(nl$angle_deg)) wd_base_deg else nl$angle_deg

  # Resolve flux: NULL means match the main lake base flux.
  flux_mg_m2_d <- if (is.null(nl$flux_mg_m2_d)) main_base_flux else nl$flux_mg_m2_d

  # Cartesian position (x = East, y = North; same convention as receptor positions)
  angle_rad <- angle_deg * pi / 180
  center_x  <- dist_m * sin(angle_rad)
  center_y  <- dist_m * cos(angle_rad)

  # Geometry validation.
  if (dist_m < main_radius_m + nearby_r) {
    stop(sprintf(
      paste0("Nearby lake overlaps the main lake: centre-to-centre distance ",
             "%.1f m < %.1f m (sum of radii = %.1f + %.1f). ",
             "Increase nearby_lake$distance_m or reduce diameter_m."),
      dist_m, main_radius_m + nearby_r, main_radius_m, nearby_r
    ))
  }

  nearest_edge <- dist_m - nearby_r  # nearest nearby-lake edge to the origin
  if (nearest_edge <= ring_radius_m) {
    stop(sprintf(
      paste0("Nearby lake is (partially) inside the receptor control-volume ring: ",
             "nearest edge is %.1f m from the main lake centre, but the ring radius is %.1f m. ",
             "Increase nearby_lake$distance_m or reduce nearby_lake$diameter_m so the ",
             "nearby lake lies entirely outside the control volume."),
      nearest_edge, ring_radius_m
    ))
  }

  list(
    diameter_m   = nl$diameter_m,
    radius_m     = nearby_r,
    center_x     = center_x,
    center_y     = center_y,
    flux_mg_m2_d = flux_mg_m2_d,
    angle_deg    = angle_deg,
    distance_m   = dist_m
  )
}


# -----------------------------------------------------------------------------
# 6.  Background CH4 helpers
# -----------------------------------------------------------------------------
#' Generate a synthetic background CH4 time series.
#'
#' @param n_steps     Number of time steps.
#' @param dt_minutes  Time step duration in minutes.
#' @param mean_ppb    Mean background (default 1920 ppb).
#' @param diel_amp    Diel amplitude (default 0 = no diel).
#' @param synoptic_amp  Synoptic-scale amplitude (default 35 ppb).
#' @param synoptic_tau_h  Autocorrelation timescale in hours (default 6).
#' @param seed        Random seed (if NULL, not set).
#' @return Numeric vector of background CH4 in ppb.
generate_background_ch4 <- function(n_steps, dt_minutes = 30,
                                    mean_ppb = 1920,
                                    diel_amp = 0,
                                    synoptic_amp = 35,
                                    synoptic_tau_h = 6,
                                    start_hour = 0,
                                    seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  hours <- seq(0, by = dt_minutes / 60, length.out = n_steps)
  hour_of_day <- (start_hour + hours) %% 24
  diel <- diel_amp * cos(2 * pi * (hour_of_day - 6) / 24)

  # AR(1) synoptic component
  alpha <- exp(-dt_minutes / 60 / synoptic_tau_h)
  synoptic <- numeric(n_steps)
  synoptic[1] <- rnorm(1, 0, synoptic_amp)
  if (n_steps > 1) {
    for (i in 2:n_steps) {
      synoptic[i] <- alpha * synoptic[i - 1] +
        sqrt(1 - alpha^2) * rnorm(1, 0, synoptic_amp)
    }
  }

  mean_ppb + diel + synoptic
}


# -----------------------------------------------------------------------------
# 6b. Monin-Obukhov wind profile (inverse / diagnostics)
# -----------------------------------------------------------------------------
#' Monin-Obukhov wind speed profile (Businger-Dyer stability corrections)
#'
#' @param z_m     Numeric vector of heights above surface (m). Must be > z0_m.
#' @param ustar   Friction velocity (m/s). Scalar.
#' @param L_m     Obukhov length (m). Scalar. Use large magnitude (e.g. -1e6) for neutral.
#' @param z0_m    Aerodynamic roughness length (m). Scalar.
#' @param kappa   von Karman constant (default 0.4).
#' @return Numeric vector of wind speeds (m/s), same length as z_m. Clamped to >= 0.1 m/s.
wind_at_height <- function(z_m, ustar, L_m, z0_m, kappa = 0.4) {
  stopifnot(length(ustar) == 1, length(L_m) == 1, length(z0_m) == 1)
  # Neutral when L is missing or too small (avoid division by zero)
  if (!is.finite(L_m) || abs(L_m) < 0.1) L_m <- -1e6
  # Guard: avoid log(0) or negative z
  z_m  <- pmax(z_m, z0_m * 1.01)
  zeta <- z_m / L_m

  psi <- numeric(length(zeta))

  # Unstable (zeta < 0): Businger-Dyer
  u_idx <- zeta < 0
  if (any(u_idx)) {
    x           <- (1 - 16 * zeta[u_idx])^0.25
    psi[u_idx]  <- 2 * log((1 + x) / 2) +
                       log((1 + x^2) / 2) -
                       2 * atan(x) + pi / 2
  }

  # Stable (zeta >= 0): Businger-Dyer linear form
  s_idx       <- zeta >= 0
  psi[s_idx]  <- -5 * zeta[s_idx]

  ws <- (ustar / kappa) * (log(z_m / z0_m) - psi)
  pmax(ws, 0.1)   # physical floor
}


# 8.  Flux unit conversions
# -----------------------------------------------------------------------------
#' Molar volume of air (L mol-1) at given temperature and pressure.
#' @param temp_C  Temperature in Celsius (default 10, typical boreal).
#' @param p_kPa   Pressure in kPa (default 101.325).
molar_volume_L <- function(temp_C = 10, p_kPa = 101.325) {
  R <- 8.31446  # J mol-1 K-1
  T_K <- temp_C + 273.15
  R * T_K / (p_kPa * 1000) * 1000  # m3/mol -> L/mol
}

#' Convert mg CH4 m-2 d-1 to ug m-2 s-1.
flux_mg_d_to_ug_s <- function(flux_mg_m2_d) {
  flux_mg_m2_d * 1000 / 86400
}

#' Convert ug m-2 s-1 to mg CH4 m-2 d-1.
flux_ug_s_to_mg_d <- function(flux_ug_m2_s) {
  flux_ug_m2_s * 86400 / 1000
}

#' Convert ppb CH4 to ug m-3.
#' Molar mass CH4 = 16.04 g/mol.
#' @param ppb  CH4 mixing ratio in ppb.
#' @param temp_C  Temperature (default 10C for boreal conditions).
ppb_to_ug_m3 <- function(ppb, temp_C = 10) {
  mv <- molar_volume_L(temp_C)
  ppb * 16.04 / mv
}

#' Convert ug m-3 CH4 to ppb.
#' @param ug_m3  CH4 mass concentration in ug m-3.
#' @param temp_C  Temperature (default 10C for boreal conditions).
ug_m3_to_ppb <- function(ug_m3, temp_C = 10) {
  mv <- molar_volume_L(temp_C)
  ug_m3 * mv / 16.04
}


# -----------------------------------------------------------------------------
# 9.  Project path helpers
# -----------------------------------------------------------------------------
#' Return project root, robust to Rscript and interactive use.
project_root <- function() {
  root <- getOption("uas_gauss_osse_root", default = NULL)
  if (!is.null(root)) return(root)
  if (file.exists("R/utils.R")) return(normalizePath("."))
  if (file.exists("../R/utils.R")) return(normalizePath(".."))
  getwd()
}

proj_path <- function(...) {
  file.path(project_root(), ...)
}


#' Save ggplot2 figures to PNG with Windows-safe fallbacks.
#'
#' Default bitmap \code{png()} can fail on Windows (path length, or device).
#' Tries Cairo PNG, then \code{ragg::agg_png}, then default; optionally writes
#' via a short path and copies into place.
#'
#' @param filename Output path (typically \code{.png}).
#' @param plot ggplot object.
#' @param ... Passed to \code{ggplot2::ggsave} (e.g. \code{width}, \code{height}, \code{dpi}).
#' @return Invisibly, the path that was written (may differ if a short-path copy was used).
ggsave_osse <- function(filename, plot, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required for ggsave_osse")
  }
  fn <- as.character(filename)[1L]
  d <- dirname(fn)
  if (!dir.exists(d)) {
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
  }

  fn_abs <- tryCatch(
    normalizePath(fn, winslash = "/", mustWork = FALSE),
    error = function(e) fn
  )

  use_short <- (.Platform$OS.type == "windows" && nchar(fn_abs) > 220L)
  short_base <- "C:/osse_fig"
  write_path <- if (use_short) {
    if (!dir.exists(short_base)) {
      dir.create(short_base, recursive = TRUE, showWarnings = FALSE)
    }
    file.path(short_base, basename(fn_abs))
  } else {
    fn_abs
  }

  do_save <- function(path) {
    ggplot2::ggsave(path, plot, ..., type = "cairo")
  }
  do_ragg <- function(path) {
    ggplot2::ggsave(path, plot, ..., device = ragg::agg_png)
  }
  do_default <- function(path) {
    ggplot2::ggsave(path, plot, ...)
  }

  ok <- tryCatch({ do_save(write_path); TRUE }, error = function(e) FALSE)
  if (!ok && requireNamespace("ragg", quietly = TRUE)) {
    ok <- tryCatch({ do_ragg(write_path); TRUE }, error = function(e) FALSE)
  }
  if (!ok) ok <- tryCatch({ do_default(write_path); TRUE }, error = function(e) FALSE)

  if (!ok) {
    stop(
      "Unable to save figure (png device failed): ", fn_abs,
      ". On Windows, try install.packages(\"ragg\") or shorten the output path."
    )
  }

  if (use_short) {
    copied <- tryCatch(
      file.copy(write_path, fn_abs, overwrite = TRUE),
      error = function(e) FALSE
    )
    if (!isTRUE(copied)) {
      warning(
        "Figure saved to short path (copy to project figures dir failed): ",
        write_path,
        call. = FALSE
      )
      return(invisible(write_path))
    }
  }

  invisible(fn_abs)
}


# -----------------------------------------------------------------------------
# 10. Scenario key helper
# -----------------------------------------------------------------------------
#' Build a compact scenario identifier string.
scenario_key <- function(lake_name, sensor_name, design_name,
                         met_name = "mixed_day") {
  paste(lake_name, sensor_name, design_name, met_name, sep = "__")
}

#' Parse a scenario key back to components.
parse_scenario_key <- function(key) {
  parts <- strsplit(key, "__", fixed = TRUE)[[1]]
  list(lake = parts[1], sensor = parts[2], design = parts[3],
       met = if (length(parts) >= 4) parts[4] else "mixed_day")
}


# -----------------------------------------------------------------------------
# 11. Progress-bar helper (works in Rscript and interactive)
# -----------------------------------------------------------------------------
#' Create a progress bar with elapsed time and ETA.
#'
#' Uses a fixed-width single-line bar that updates in place. On completion,
#' prints elapsed time. Works in Rscript and interactive sessions.
#'
#' @param total  Total number of items.
#' @param label  Character label for the bar.
#' @return A list with \code{$tick(msg = "")} and \code{$done()} methods.
make_progress <- function(total, label = "Progress") {
  count <- 0L
  start_time <- Sys.time()
  next_pct <- 5
  total <- max(1L, as.integer(total))
  bar_width <- 54L
  last_render_nchar <- 0L

  fmt_duration_short <- function(seconds) {
    if (!is.finite(seconds) || is.na(seconds) || seconds < 0) return("--")
    seconds <- as.integer(round(seconds))
    hrs <- seconds %/% 3600L
    mins <- (seconds %% 3600L) %/% 60L
    secs <- seconds %% 60L
    if (hrs > 0L) {
      sprintf("%dh %02dm %02ds", hrs, mins, secs)
    } else if (mins > 0L) {
      sprintf("%dm %02ds", mins, secs)
    } else {
      sprintf("%ds", secs)
    }
  }

  render <- function(final = FALSE) {
    pct <- 100 * count / total
    filled <- min(bar_width, max(0L, as.integer(round(bar_width * count / total))))
    empty <- max(0L, bar_width - filled)
    bar <- paste0(strrep("=", filled), strrep(" ", empty))
    elapsed_sec <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    eta_sec <- if (count > 0L && count < total) elapsed_sec / count * (total - count) else 0
    eta_label <- if (count == 0L) {
      "--"
    } else if (count < total) {
      fmt_duration_short(eta_sec)
    } else {
      "0s"
    }
    line <- sprintf(
      "\r%s |%s| %3d%% ETA %s",
      label,
      bar,
      as.integer(round(pct)),
      eta_label
    )
    pad <- max(0L, last_render_nchar - nchar(line, type = "width"))
    cat(line, strrep(" ", pad), sep = "")
    flush.console()
    last_render_nchar <<- nchar(line, type = "width")
    if (isTRUE(final)) cat("\n")
  }

  render()

  list(
    tick = function(msg = "") {
      count <<- count + 1L
      pct <- 100 * count / total
      if (pct >= next_pct || count == total) {
        render(final = (count == total))
        if (count < total) next_pct <<- min(100, next_pct + 5)
      }
      if (count == total) {
        elapsed <- difftime(Sys.time(), start_time, units = "secs")
        message(sprintf("  %s: %s", label, format(elapsed, digits = 2)))
      }
    },
    done = function() {
      if (count < total) cat("\n")
    }
  )
}


# -----------------------------------------------------------------------------
# 12. Deterministic seed helper
# -----------------------------------------------------------------------------
#' Deterministic seed from a character string.
#'
#' Uses positional weighting so that anagrams and strings with the same
#' character-code sum produce different seeds.
#' @param label Character string.
#' @return Integer seed.
string_to_seed <- function(label) {
  bytes <- utf8ToInt(as.character(label))
  hash <- sum(bytes * seq_along(bytes)) + sum(bytes^2)
  (hash %% 100000L) + 1L
}

#' Evaluate an expression with an optional seed, restoring RNG state afterward.
#'
#' Useful for deterministic helper draws that should not perturb the caller's
#' RNG stream.
#'
#' @param expr Expression to evaluate.
#' @param seed Optional integer seed. When NULL, the expression is evaluated
#'   without reseeding.
#' @return Result of \code{expr}.
with_preserved_seed <- function(expr, seed = NULL) {
  rng_exists <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  rng_state <- if (rng_exists) {
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  } else {
    NULL
  }

  on.exit({
    if (rng_exists) {
      assign(".Random.seed", rng_state, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)

  if (!is.null(seed)) set.seed(as.integer(seed))
  force(expr)
}

message("utils.R loaded successfully.")

