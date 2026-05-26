###############################################################################
# 04_conley_inverse_flux.R - Gaussian-surface (mass balance) inverse for
#                             whole-lake CH4 flux retrieval
#
# Method: Conley et al. (2017) "Application of Gauss's theorem to quantify
#   localized surface emissions from airborne measurements of wind and trace
#   gases", Atmospheric Measurement Techniques, 10(9), 3345-3358.
#
# Control volume:
#   - Vertical cylinder centred on the lake.
#   - Lateral surface sampled by the outer receptor ring.
#   - Top flux = 0 (well-mixed boundary layer assumption).
#   - Bottom flux = lake surface emission (what we retrieve).
#
# Divergence theorem:
#   integral(F_surface dA) = closed_integral(C * v_n * dA) over the lateral wall
#
# Discrete form:
#   Q_total [ug s-1] = sum_i(C_i * v_r,i * dA_i)
#
#   where:
#     C_i    = CH4 concentration at receptor i [ug m-3]
#     v_r,i  = outward radial wind at receptor i [m s-1]
#               = ws * cos(azimuth_i - wd_TO)
#               wd_TO = (wd_FROM + 180) %% 360
#     dA_i   = R_outer * dtheta * dz [m2], midpoint height bands
#     R_outer = lake_radius + distance_from_edge_m  (from design geometry)
#
# Lake-averaged flux [ug m-2 s-1] = Q_total / A_lake, then converted to
# mg m-2 d-1.
#
# Inputs:
#   - Simulated measurements from 03_drone_and_sensor_sim.R
#     (must include receptor metadata such as ring_type and azimuth_deg)
#
# Outputs:
#   - output/retrieved_flux_<scenario>.csv when save_output = TRUE
###############################################################################

if (!exists("project_root", mode = "function")) {
  proj_root <- if (file.exists("R/utils.R")) "." else ".."
  source(file.path(proj_root, "R", "utils.R"))
}


# -----------------------------------------------------------------------------
# 1.  Vertical height-band areas (mid-point rule)
# -----------------------------------------------------------------------------
#' Compute dz for each measurement height using the midpoint rule.
#'
#' Lower bound of the bottom height is 0 (ground). Upper bound of the top
#' height mirrors the spacing below it (symmetric extrapolation), so the
#' top band is not artificially truncated.  For a single-height design the
#' band runs from 0 to 2 × height.
#'
#' Example for heights_m = c(3, 6, 10, 15, 20):
#'   z= 3: lo=0,    hi=4.5,  dz=4.5
#'   z= 6: lo=4.5,  hi=8,    dz=3.5
#'   z=10: lo=8,    hi=12.5, dz=4.5
#'   z=15: lo=12.5, hi=17.5, dz=5.0
#'   z=20: lo=17.5, hi=22.5, dz=5.0   (mirrors gap below: 20-15=5)
#'   Total cylinder height: 22.5 m
#'
#' @param heights_m  Numeric vector of unique measurement heights (m).
#' @return data.table with columns: z (height), dz (band width in m).
compute_height_bands <- function(heights_m) {
  h <- sort(unique(heights_m))
  n <- length(h)
  dz <- numeric(n)

  for (k in seq_len(n)) {
    lo <- if (k == 1L) 0 else (h[k - 1L] + h[k]) / 2
    hi <- if (k == n && n > 1L) h[k] + (h[k] - h[k - 1L]) / 2 else
          if (k == n && n == 1L) h[k] * 2              else
          (h[k] + h[k + 1L]) / 2
    dz[k] <- hi - lo
  }

  data.table(z = h, dz = dz)
}


# 2.  Conley mass-balance inverse
# -----------------------------------------------------------------------------
#' Retrieve whole-lake CH4 flux using the Conley Gaussian-surface method.
#'
#' Applies Gauss's theorem over a vertical cylinder:
#'   Q = sum(C_i * v_r,i * dA_i)          [ug s-1]
#'   F = Q / A_lake                       [ug m-2 s-1] -> mg m-2 d-1
#'
#' Only outer-ring receptors (ring_type == "outer") are used.  Per-window
#' calculations are fully vectorised via data.table grouping.
#'
#' @param meas_dt  data.table from simulate_drone_measurements(), with
#'   ring_type, azimuth_deg, ws_m_s, wd_deg, z, ustar_m_s, L_m,
#'   stability_class, ch4_measured_ppb, flux_true_mg_m2_d.
#' @param lake_diameter  Lake diameter in metres.
#' @param bg_percentile  Ignored (background subtraction removed).
#' @param bg_dt  Ignored (background subtraction removed). Kept for API
#'   compatibility with earlier saved analysis calls.
#' @param distance_from_edge_m  Distance of the outer receptor ring beyond
#'   the lake edge in metres (default 10, matching define_receptor_config).
#'   Used to set R_outer = lake_radius + distance_from_edge_m from design
#'   geometry rather than from averaged receptor positions.
#' @param n_azimuth  Number of azimuth positions in the outer ring.  When
#'   NULL (default) it is inferred from the first time window in boundary_dt.
#'   Pass design_cfg$n_azimuth explicitly for robustness.
#' @param wd_uncertainty_deg  SD of per-step wind direction noise applied to
#'   wd_deg before computing v_r (default 0; set > 0 for realistic runs).
#' @param ws_uncertainty_m_s  SD of per-step wind speed error (m s-1) applied
#'   to ws(z) after MOST profile (ignored when use_observed_ws is TRUE).
#' @param use_observed_ws  If TRUE, use ws_at_z_m_s from meas_dt with error
#'   ws_observation_error_m_s; no MOST profile or z0 used for ws.
#' @param ws_observation_error_m_s  SD of error (m s-1) on observed ws when
#'   use_observed_ws is TRUE (default 0.5).
#' @param z0_multiplier  Multiplier on z0_m for wind profile (ignored when
#'   use_observed_ws is TRUE).
#' @param min_ws  Minimum wind speed (m s-1) below which a time window is
#'   flagged "low_ws" and excluded from performance metrics (default 0.5).
#' @param scenario_label  Label for output files.
#' @param save_output  Logical.
#' @return data.table with retrieved flux per time window plus QC diagnostics.
retrieve_flux <- function(meas_dt,
                          lake_diameter,
                          bg_percentile        = 0.05,
                          bg_dt                = NULL,
                          distance_from_edge_m = 10,
                          n_azimuth            = NULL,
                          wd_uncertainty_deg   = 0,
                          ws_uncertainty_m_s   = 0.3,
                          use_observed_ws      = TRUE,
                          ws_observation_error_m_s = 0.5,
                          z0_multiplier        = 1.0,
                          use_height_avg_wind  = FALSE,
                          min_ws               = 0.5,
                          scenario_label       = "default",
                          save_output          = TRUE) {

  lake_radius <- lake_diameter / 2
  lake_area   <- pi * lake_radius^2                    # m2
  r_outer     <- lake_radius + distance_from_edge_m   # from design geometry

  dt <- copy(meas_dt)

  # Identify outer ring receptors.
  if ("ring_type" %in% names(dt)) {
    boundary_dt <- dt[ring_type == "outer"]
  } else {
    dt[, radius_m := sqrt(x^2 + y^2)]
    r_max <- max(dt$radius_m, na.rm = TRUE)
    boundary_dt <- dt[abs(radius_m - r_max) < 0.5]
    dt[, radius_m := NULL]
  }

  if (nrow(boundary_dt) == 0L) {
    stop("No outer-ring receptors found. Check that ring_type == 'outer' exists in meas_dt.")
  }

  # Cylinder geometry.
  # Azimuth spacing: from explicit design parameter or first window's unique
  # azimuths (design is fixed across all time windows).
  if (!is.null(n_azimuth)) {
    n_az <- as.integer(n_azimuth)
  } else {
    first_win <- sort(unique(boundary_dt$time_min))[1L]
    n_az      <- length(unique(boundary_dt[time_min == first_win]$azimuth_deg))
  }
  d_theta <- 2 * pi / n_az  # radians between adjacent azimuths

  # Vertical height bands (midpoint rule)
  height_bands <- compute_height_bands(unique(boundary_dt$z))
  boundary_dt  <- merge(boundary_dt, height_bands, by = "z", all.x = TRUE)

  # Area element per receptor: R * dtheta * dz [m2]
  boundary_dt[, dA := r_outer * d_theta * dz]

  # Background concentration is not subtracted in this publication workflow.
  # We intentionally do NOT estimate or subtract a background term.
  # In the discrete inversion integral, treat C_bg = 0.
  bg_source <- "C_bg=0 (background subtraction removed)"

  # Optional per-step wind direction noise.
  if (wd_uncertainty_deg > 0) {
    set.seed(string_to_seed(scenario_label) + 7777L)
    times    <- sort(unique(boundary_dt$time_min))
    wd_noise <- data.table(time_min = times,
                           wd_noise = rnorm(length(times), 0, wd_uncertainty_deg))
    boundary_dt <- merge(boundary_dt, wd_noise, by = "time_min", all.x = TRUE)
    boundary_dt[, wd_deg := (wd_deg + wd_noise) %% 360]
    boundary_dt[, wd_noise := NULL]
  }

  # Wind speed at each receptor: observed (with error) or MOST profile.
  if (use_observed_ws && "ws_at_z_m_s" %in% names(boundary_dt)) {
    set.seed(string_to_seed(scenario_label) + 7778L)
    boundary_dt[, ws_z := pmax(0.1, ws_at_z_m_s + rnorm(.N, 0, ws_observation_error_m_s))]
  } else {
    if (!"z0_m" %in% names(boundary_dt)) boundary_dt[, z0_m := 0.001]
    boundary_dt[, z0_eff := z0_m * z0_multiplier]
    boundary_dt[, ws_z := wind_at_height(z_m = z, ustar = ustar_m_s[1L],
                                          L_m = L_m[1L], z0_m = z0_eff[1L]),
                by = time_min]
    if (ws_uncertainty_m_s > 0) {
      set.seed(string_to_seed(scenario_label) + 7778L)
      times    <- sort(unique(boundary_dt$time_min))
      ws_noise <- data.table(time_min = times,
                             ws_noise = rnorm(length(times), 0, ws_uncertainty_m_s))
      boundary_dt <- merge(boundary_dt, ws_noise, by = "time_min", all.x = TRUE)
      boundary_dt[, ws_z := pmax(0.1, ws_z + ws_noise)]
      boundary_dt[, ws_noise := NULL]
    }
  }

  # Optional height-layer wind averaging (single-mast equivalent).
  # Replaces per-receptor noisy ws_z with the mean across all azimuth positions
  # at the same height and time window. Reduces effective wind noise by
  # ~1/sqrt(n_azimuth), equivalent to a fixed reference anemometer whose
  # reading applies uniformly around the ring at each height.
  if (use_height_avg_wind) {
    boundary_dt[, ws_z := mean(ws_z), by = .(time_min, z)]
  }

  # Outward radial wind component.
  # wd_deg is wind-FROM direction; wd_to is the direction the wind blows TOWARD.
  # v_r = ws(z) * cos(azimuth - wd_to); positive is outward through cylinder.
  boundary_dt[, wd_to := (wd_deg + 180) %% 360]
  boundary_dt[, v_r   := ws_z * cos((azimuth_deg - wd_to) * pi / 180)]

  # Per-receptor flux contribution.
  boundary_dt[, enh_ug_m3    := ppb_to_ug_m3(ch4_measured_ppb)]
  boundary_dt[, flux_contrib := enh_ug_m3 * v_r * dA]

  # Per-window mass balance.
  # All time windows are grouped simultaneously; data.table handles the
  # aggregation.
  flux_dt <- boundary_dt[, .(
    flux_true            = flux_true_mg_m2_d[1L],
    n_receptors          = .N,
    total_Q_ug_s         = sum(flux_contrib, na.rm = TRUE),
    net_vr_dA_m3_s       = sum(v_r * dA, na.rm = TRUE),
    mean_enhancement_ppb = round(mean(ch4_measured_ppb, na.rm = TRUE), 4),
    bg_ppb               = 0,
    mean_abs_vr_m_s      = round(mean(abs(v_r), na.rm = TRUE), 4),
    ustar                = ustar_m_s[1L],
    L                    = L_m[1L],
    ws                   = ws_m_s[1L],
    ws_at_z_mean_m_s     = round(mean(ws_z, na.rm = TRUE), 4),
    stability            = stability_class[1L]
  ), by = time_min][order(time_min)]

  # Derived flux and QC flags.
  # Compute flux_mg_d for every window; negative retrieved fluxes are kept and
  # included in metrics (no separate negative_flux flag).
  flux_dt[, .flux_raw := flux_ug_s_to_mg_d(total_Q_ug_s / lake_area)]
  flux_dt[, qc_flag := fcase(
    !is.finite(ws) | ws < min_ws, "low_ws",
    !is.finite(.flux_raw),        "invalid",
    default =                     "ok"
  )]
  flux_dt[, flux_retrieved := fcase(
    qc_flag %in% c("low_ws", "invalid"), NA_real_,
    default = round(.flux_raw, 4)
  )]
  flux_dt[, .flux_raw := NULL]

  # Reorder columns for readability
  setcolorder(flux_dt, c("time_min", "flux_retrieved", "flux_true",
                          "n_receptors", "total_Q_ug_s",
                          "net_vr_dA_m3_s",
                          "mean_enhancement_ppb", "bg_ppb",
                          "mean_abs_vr_m_s",
                          "ustar", "L", "ws", "ws_at_z_mean_m_s", "stability", "qc_flag"))

  # Performance metrics.
  valid_dt <- flux_dt[qc_flag == "ok" & is.finite(flux_retrieved)]
  if (nrow(valid_dt) > 0L) {
    rel_err  <- (valid_dt$flux_retrieved - valid_dt$flux_true) / valid_dt$flux_true
    bias_pct <- 100 * mean(rel_err)
    rmse_pct <- 100 * sqrt(mean(rel_err^2))
    corr     <- if (nrow(valid_dt) > 2L)
      cor(valid_dt$flux_retrieved, valid_dt$flux_true) else NA_real_

    flux_dt[, `:=`(
      bias_pct       = round(bias_pct, 2),
      rmse_pct       = round(rmse_pct, 2),
      correlation    = round(corr, 4),
      valid_fraction = round(nrow(valid_dt) / nrow(flux_dt), 4)
    )]

    message(sprintf(
      "  Conley inverse [bg=%s]: bias=%.1f%%, RMSE=%.1f%%, r=%.3f, valid=%.0f%%",
      bg_source, bias_pct, rmse_pct,
      ifelse(is.na(corr), 0, corr),
      100 * nrow(valid_dt) / nrow(flux_dt)
    ))
  }

  if (save_output) {
    outfile <- proj_path("data", "raw",
                         sprintf("retrieved_flux_%s.csv", scenario_label))
    fwrite(flux_dt, outfile)
  }

  flux_dt
}


# -----------------------------------------------------------------------------
# 3.  Capability metrics
# -----------------------------------------------------------------------------
#' Summarise retrieval performance for a single scenario.
#'
#' @param flux_dt  data.table from retrieve_flux().
#' @param scenario_label  Character label.
#' @param uncertainty_target_pct  RMSE threshold for "meets target" (default 25).
#' @return One-row data.table of capability metrics.
compute_capability_metrics <- function(flux_dt,
                                       scenario_label = "default",
                                       uncertainty_target_pct = 25) {

  valid_dt <- flux_dt[qc_flag == "ok" & is.finite(flux_retrieved)]
  n_total  <- nrow(flux_dt)
  n_valid  <- nrow(valid_dt)

  if (n_valid < 2L) {
    return(data.table(
      scenario       = scenario_label,
      bias_pct       = NA_real_,
      rmse_pct       = NA_real_,
      correlation    = NA_real_,
      valid_fraction = n_valid / max(n_total, 1L),
      meets_target   = FALSE,
      n_windows      = n_total,
      n_valid        = n_valid
    ))
  }

  rel_err  <- (valid_dt$flux_retrieved - valid_dt$flux_true) / valid_dt$flux_true
  bias_pct <- 100 * mean(rel_err)
  rmse_pct <- 100 * sqrt(mean(rel_err^2))
  corr     <- cor(valid_dt$flux_retrieved, valid_dt$flux_true)

  data.table(
    scenario       = scenario_label,
    bias_pct       = round(bias_pct, 2),
    rmse_pct       = round(rmse_pct, 2),
    correlation    = round(corr, 4),
    valid_fraction = round(n_valid / n_total, 4),
    meets_target   = rmse_pct < uncertainty_target_pct,
    n_windows      = n_total,
    n_valid        = n_valid
  )
}



