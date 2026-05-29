#!/usr/bin/env Rscript
###############################################################################
# analyze_temporal_measurement_changing_atmosphere.R
#
# Full-analysis sweep: how does Conley retrieval change when one drone samples
# a cylindrical measurement surface over time while the atmosphere transitions?
#
# Experiment Configuration:
#   - 100 m lake, 1 m lake grid, 50 mg CH4 m^-2 d^-1
#   - LI-COR LI-7810 sensor
#   - 12 hour stability transitions: B-D, C-D, C-F, F-D, D-B
#   - mission durations from 30 min to 12 h
#   - repeated full-profile cadences of 30, 60, 90, and 120 min
#   - 20 Monte Carlo repeats per transition-duration configuration by default
#   - Conley inverse treats each mission as one snapshot (time_min = 0),
#     while concentrations and winds evolve during the sequential mission
#   - repeated-profile tests reuse the same evolving atmosphere, retrieve one
#     flux per full profile, and store the mean retrieval across profiles
#
# Outputs (under output/):
#   - temporal_measurement_changing_atmosphere_runs.csv
#   - temporal_measurement_changing_atmosphere_summary.csv
#   - temporal_measurement_changing_atmosphere_composite.png
#   - temporal_measurement_changing_atmosphere_additional_transitions_2x2.png
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
suppressMessages(source(file.path(proj_root, "R", "01_synthetic_emissions.R")))
suppressMessages(source(file.path(proj_root, "R", "02b_gaussian_forward.R")))
suppressMessages(source(file.path(proj_root, "R", "04_conley_inverse_flux.R")))
suppressMessages(source(file.path(proj_root, "R", "06_visualization.R")))
suppressMessages(source(file.path(proj_root, "scripts", "full_analysis", "results_reporting_helpers.R")))

if (!requireNamespace("cowplot", quietly = TRUE)) {
  stop("Package `cowplot` is required for the composite figure.")
}


# -----------------------------------------------------------------------------
# 1. CLI arguments
# -----------------------------------------------------------------------------
parse_named_args <- function(args) {
  out <- list()
  for (arg in args) {
    if (!startsWith(arg, "--")) next
    parts <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1]]
    key <- parts[1]
    val <- if (length(parts) >= 2L) paste(parts[-1], collapse = "=") else "TRUE"
    out[[key]] <- val
  }
  out
}

parse_num_vec <- function(x) {
  if (is.null(x) || !nzchar(x)) return(numeric(0))
  as.numeric(strsplit(x, ",", fixed = TRUE)[[1]])
}

parse_bool <- function(x) {
  if (is.null(x) || !nzchar(x)) return(FALSE)
  tolower(trimws(x)) %in% c("true", "t", "1", "yes", "y")
}

parse_transition_vec <- function(x) {
  if (is.null(x) || !nzchar(x)) return(list())
  raw_vals <- trimws(strsplit(x, ",", fixed = TRUE)[[1]])
  lapply(raw_vals, function(val) {
    parts <- trimws(strsplit(val, "-", fixed = TRUE)[[1]])
    if (length(parts) != 2L) stop("Transitions must look like B-D,C-D,C-F,F-D,D-B.")
    parts
  })
}

transition_label <- function(tr) paste(tr, collapse = "-")

cli <- parse_named_args(commandArgs(trailingOnly = TRUE))

LAKE_NAME <- "small"
DESIGN_NAME <- if (!is.null(cli$design)) cli$design else "design_96az_many_heights_30min"
SENSOR_NAME <- if (!is.null(cli$sensor)) cli$sensor else "li7810"
BASE_FLUX_MG_D <- if (!is.null(cli$flux_mg_m2_d)) as.numeric(cli$flux_mg_m2_d) else 50
LAKE_GRID_RESOLUTION_M <- if (!is.null(cli$lake_grid_m)) as.numeric(cli$lake_grid_m) else 1
TRANSITION_DURATION_MIN <- if (!is.null(cli$transition_duration_min)) as.numeric(cli$transition_duration_min) else 720
MISSION_DURATIONS_MIN <- {
  vals <- parse_num_vec(cli$mission_durations_min)
  if (length(vals)) vals else c(30, 90, 150, 210, 270, 330, 390, 450, 510, 570, 630, 720)
}
REPEAT_PROFILE_DURATIONS_MIN <- {
  vals <- parse_num_vec(cli$repeat_profile_min)
  if (length(vals)) vals else c(30, 60, 90, 120)
}
PRIMARY_TRANSITIONS <- list(c("B", "D"), c("C", "D"), c("C", "F"))
ADDITIONAL_TRANSITIONS <- list(c("F", "D"), c("D", "B"))
DEFAULT_TRANSITIONS <- c(PRIMARY_TRANSITIONS, ADDITIONAL_TRANSITIONS)
PRIMARY_TRANSITION_LABELS <- vapply(PRIMARY_TRANSITIONS, transition_label, character(1))
ADDITIONAL_TRANSITION_LABELS <- vapply(ADDITIONAL_TRANSITIONS, transition_label, character(1))
TRANSITIONS <- {
  vals <- parse_transition_vec(cli$transitions)
  if (length(vals)) vals else DEFAULT_TRANSITIONS
}
TRANSITION_LABELS <- vapply(TRANSITIONS, transition_label, character(1))
N_MC <- if (!is.null(cli$n_mc)) as.integer(cli$n_mc) else 20L
N_WORKERS_REQUESTED <- if (!is.null(cli$n_workers)) as.integer(cli$n_workers) else NA_integer_
U_MEAN_M_S <- 2
WD_MEAN_DEG <- if (!is.null(cli$wd_mean_deg)) as.numeric(cli$wd_mean_deg) else 225
WS_FLOOR_M_S <- if (!is.null(cli$ws_floor_m_s)) as.numeric(cli$ws_floor_m_s) else 0.5
MET_DT_MIN <- if (!is.null(cli$met_dt_min)) as.numeric(cli$met_dt_min) else 30
SEED_BASE <- if (!is.null(cli$seed)) as.integer(cli$seed) else 20260414L
PAPER_STATS_ONLY <- parse_bool(cli$paper_stats_only)
OUTPUT_STEM <- if (!is.null(cli$output_stem) && nzchar(cli$output_stem)) {
  cli$output_stem
} else {
  "temporal_measurement_changing_atmosphere"
}

stopifnot(LAKE_NAME %in% names(LAKE_CONFIGS))
stopifnot(DESIGN_NAME %in% names(SAMPLING_CONFIGS))
stopifnot(SENSOR_NAME %in% names(SENSOR_CONFIGS))
stopifnot(all(vapply(TRANSITIONS, function(tr) {
  length(tr) == 2L && all(tr %in% names(PASQUILL_CLASSES))
}, logical(1))))
stopifnot(all(is.finite(MISSION_DURATIONS_MIN)), all(MISSION_DURATIONS_MIN > 0))
stopifnot(all(is.finite(REPEAT_PROFILE_DURATIONS_MIN)), all(REPEAT_PROFILE_DURATIONS_MIN > 0))
stopifnot(all(REPEAT_PROFILE_DURATIONS_MIN <= TRANSITION_DURATION_MIN + 1e-9))

BESTFIT_PATH <- file.path(proj_root, "output", "ameriflux_wind_variability_bestfit.csv")
if (!file.exists(BESTFIT_PATH)) {
  stop("AmeriFlux best-fit CSV not found: ", BESTFIT_PATH,
       "\nRun process_ameriflux_wind_variability.R first.")
}

OUTPUT_DIR <- proj_path("output")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

RUNS_CSV <- file.path(OUTPUT_DIR, sprintf("%s_runs.csv", OUTPUT_STEM))
SUMMARY_CSV <- file.path(OUTPUT_DIR, sprintf("%s_summary.csv", OUTPUT_STEM))
FIG_PNG <- file.path(OUTPUT_DIR, sprintf("%s_composite.png", OUTPUT_STEM))
ADDITIONAL_FIG_PNG <- file.path(OUTPUT_DIR, sprintf("%s_additional_transitions_2x2.png", OUTPUT_STEM))
KEY_CSV <- file.path(OUTPUT_DIR, sprintf("%s_key_table.csv", OUTPUT_STEM))
PAPER_STATS_CSV <- file.path(OUTPUT_DIR, sprintf("%s_paper_stats.csv", OUTPUT_STEM))

FIG_BASE_SIZE <- 12
FIG_POINT_SIZE <- 2.4
FIG_LINEWIDTH <- 1.0


# -----------------------------------------------------------------------------
# 2. Helpers
# -----------------------------------------------------------------------------
PASQUILL_ZETA_ANCHORS <- c(A = -2.0, B = -1.0, C = -0.5,
                           D =  0.0, E =  0.5, F =  2.0)

quietly <- function(expr) {
  tf_out <- tempfile()
  tf_msg <- tempfile()
  out_con <- file(tf_out, open = "wt")
  msg_con <- file(tf_msg, open = "wt")
  on.exit({
    try(sink(type = "message"), silent = TRUE)
    try(sink(), silent = TRUE)
    try(close(msg_con), silent = TRUE)
    try(close(out_con), silent = TRUE)
    unlink(c(tf_out, tf_msg), force = TRUE)
  }, add = TRUE)

  sink(out_con)
  sink(msg_con, type = "message")
  force(expr)
}

configure_main_lake_emissions <- function(lake_name, base_flux_mg_m2_d, lake_grid_resolution_m) {
  emission_result <- quietly(generate_single_lake_emissions(
    lake_name = lake_name,
    save_output = FALSE,
    base_flux_override = base_flux_mg_m2_d
  ))

  emission_times_hour <- if ("hour" %in% names(emission_result$flux_field)) {
    sort(unique(emission_result$flux_field$hour))
  } else {
    0:23
  }

  emission_result$lake_geom <- define_lake_geometry(
    diameter_m = emission_result$diameter_m,
    resolution_m = lake_grid_resolution_m
  )
  emission_result$flux_field <- generate_effective_flux_field(
    lake_geom = emission_result$lake_geom,
    base_flux = emission_result$base_flux,
    littoral_factor = EMISSION_DEFAULTS$littoral_factor,
    spatial_corr_length_m = EMISSION_DEFAULTS$spatial_corr_length_m,
    diel_amplitude = EMISSION_DEFAULTS$diel_amplitude,
    diel_phase_hour = EMISSION_DEFAULTS$diel_phase_hour,
    times_hour = emission_times_hour
  )
  emission_result$lake_mean <- emission_result$flux_field[, .(
    lake_mean_flux = mean(flux_mg_m2_d),
    flux_sd = stats::sd(flux_mg_m2_d),
    flux_cv = stats::sd(flux_mg_m2_d) / mean(flux_mg_m2_d),
    n_cells = .N
  ), by = hour]
  emission_result$lake_mean[, `:=`(
    lake_name = lake_name,
    diameter_m = emission_result$diameter_m,
    base_flux = emission_result$base_flux
  )]

  emission_result
}

interpolate_briggs <- function(t_min, transition_duration_min, class_start, class_end) {
  w <- pmin(pmax(t_min / transition_duration_min, 0), 1)
  cs <- PASQUILL_CLASSES[[class_start]]
  ce <- PASQUILL_CLASSES[[class_end]]
  list(
    sy_a     = cs$sy_a     + w * (ce$sy_a     - cs$sy_a),
    sy_b     = cs$sy_b     + w * (ce$sy_b     - cs$sy_b),
    sz_a     = cs$sz_a     + w * (ce$sz_a     - cs$sz_a),
    sz_b     = cs$sz_b     + w * (ce$sz_b     - cs$sz_b),
    sz_p     = cs$sz_p     + w * (ce$sz_p     - cs$sz_p),
    wind_exp = cs$wind_exp + w * (ce$wind_exp - cs$wind_exp)
  )
}

zeta_at_time <- function(t_min, transition_duration_min, class_start, class_end) {
  w <- pmin(pmax(t_min / transition_duration_min, 0), 1)
  z0 <- PASQUILL_ZETA_ANCHORS[[class_start]]
  z1 <- PASQUILL_ZETA_ANCHORS[[class_end]]
  z0 + w * (z1 - z0)
}

load_bestfit_interpolators <- function(path) {
  bf <- fread(path)
  bf <- bf[order(ZETA_BIN_CENTER)]
  list(
    sigma_u = approxfun(bf$ZETA_BIN_CENTER, bf$U_SIGMA_POLY_FIT, rule = 2),
    sigma_theta = approxfun(bf$ZETA_BIN_CENTER, bf$SIGMA_THETA_POLY_FIT, rule = 2)
  )
}

generate_wind_timeline <- function(transition_duration_min, met_dt_min,
                                   u_mean_m_s, wd_mean_deg, ws_floor_m_s,
                                   bestfit_interp, class_start, class_end,
                                   seed) {
  set.seed(seed)
  boundaries <- seq(0, transition_duration_min, by = met_dt_min)
  if (tail(boundaries, 1) < transition_duration_min) {
    boundaries <- c(boundaries, transition_duration_min)
  }
  n <- length(boundaries)

  ws_draws <- numeric(n)
  wd_draws <- numeric(n)
  sigma_u_vec <- numeric(n)
  sigma_theta_vec <- numeric(n)
  zeta_vec <- numeric(n)

  for (i in seq_len(n)) {
    zeta_now <- zeta_at_time(boundaries[i], transition_duration_min, class_start, class_end)
    sigma_u_now <- bestfit_interp$sigma_u(zeta_now)
    sigma_theta_now <- bestfit_interp$sigma_theta(zeta_now)

    zeta_vec[i] <- zeta_now
    sigma_u_vec[i] <- sigma_u_now
    sigma_theta_vec[i] <- sigma_theta_now
    ws_draws[i] <- u_mean_m_s
    wd_draws[i] <- (wd_mean_deg + rnorm(1, 0, sigma_theta_now)) %% 360
  }

  data.table(
    time_min = boundaries,
    zeta = zeta_vec,
    sigma_u = sigma_u_vec,
    sigma_theta = sigma_theta_vec,
    ws = ws_draws,
    wd = wd_draws
  )
}

ws_at_time <- function(t_min, wind_tl) {
  approx(wind_tl$time_min, wind_tl$ws, xout = t_min, rule = 2)$y
}

wd_at_time <- function(t_min, wind_tl) {
  approx(wind_tl$time_min, wind_tl$wd, xout = t_min, rule = 2)$y
}

pg_sigma_y_interp <- function(x_m, coefs) {
  coefs$sy_a * x_m / sqrt(1 + coefs$sy_b * x_m)
}

pg_sigma_z_interp <- function(x_m, coefs) {
  coefs$sz_a * x_m * (1 + coefs$sz_b * x_m)^coefs$sz_p
}

gaussian_plume_point_interp <- function(x_down, y_cross, z_rec, q, u_ref, coefs) {
  valid <- x_down > 1e-3
  conc <- numeric(length(x_down))
  if (!any(valid)) return(conc)

  xv <- x_down[valid]
  yv <- y_cross[valid]
  zv <- z_rec[valid]

  sy <- pmax(pg_sigma_y_interp(xv, coefs), 0.01)
  sz <- pmax(pg_sigma_z_interp(xv, coefs), 0.01)

  conc[valid] <- (q / (pi * sy * sz * u_ref)) *
    exp(-yv^2 / (2 * sy^2)) *
    exp(-zv^2 / (2 * sz^2))
  conc
}

gaussian_plume_disk_interp <- function(rec_x, rec_y, rec_z,
                                       source_radius, flux_ug_m2_s,
                                       u_ref, wd_deg, coefs,
                                       n_rings = 20L, n_src_az = 36L) {
  conc_total <- numeric(length(rec_x))

  wd_to_rad <- ((wd_deg + 180) %% 360) * pi / 180
  wind_dx <- sin(wd_to_rad)
  wind_dy <- cos(wd_to_rad)

  r_edges <- seq(0, source_radius, length.out = n_rings + 1L)
  r_mids <- (r_edges[-1L] + r_edges[-(n_rings + 1L)]) / 2
  theta_mids <- seq(0, 2 * pi * (1 - 1 / n_src_az), length.out = n_src_az)
  d_theta <- 2 * pi / n_src_az

  for (kr in seq_len(n_rings)) {
    r_lo <- r_edges[kr]
    r_hi <- r_edges[kr + 1L]
    r_mid <- r_mids[kr]
    dA_cell <- (r_hi^2 - r_lo^2) / 2 * d_theta
    q_cell <- flux_ug_m2_s * dA_cell

    for (ka in seq_len(n_src_az)) {
      theta_now <- theta_mids[ka]
      sx <- r_mid * cos(theta_now)
      sy <- r_mid * sin(theta_now)

      dx <- rec_x - sx
      dy <- rec_y - sy
      x_down <- dx * wind_dx + dy * wind_dy
      y_cross <- dx * wind_dy - dy * wind_dx

      conc_total <- conc_total + gaussian_plume_point_interp(
        x_down = x_down,
        y_cross = y_cross,
        z_rec = rec_z,
        q = q_cell,
        u_ref = u_ref,
        coefs = coefs
      )
    }
  }

  conc_total
}

build_profile_truth_and_measurements <- function(total_flight_min,
                                                 design_cfg,
                                                 lake_diameter,
                                                 flux_mg_m2_d,
                                                 wind_tl,
                                                 transition_duration_min,
                                                 class_start,
                                                 class_end,
                                                 sensor_cfg,
                                                 seed,
                                                 window_time_min = 0,
                                                 start_time_min = 0) {
  lake_radius <- lake_diameter / 2
  receptors <- define_receptor_config(
    lake_diameter = lake_diameter,
    n_azimuth = design_cfg$n_azimuth,
    heights_m = design_cfg$heights_m,
    distance_from_edge_m = 10,
    inner_ring_fracs = NULL,
    background_heights_m = NULL
  )
  setorder(receptors, z, azimuth_deg)

  n_samples <- nrow(receptors)
  dwell_min <- total_flight_min / n_samples
  dwell_s <- dwell_min * 60
  flux_ug_s <- flux_mg_d_to_ug_s(flux_mg_m2_d)

  set.seed(seed)
  precision_at_dwell <- sensor_cfg$precision_1s_ppb / sqrt(max(dwell_s, 1))
  precision_at_dwell <- max(precision_at_dwell, sensor_cfg$precision_opt_ppb)

  rows <- vector("list", n_samples)

  for (k in seq_len(n_samples)) {
    rec <- receptors[k]
    t_k <- start_time_min + (k - 1L) * dwell_min
    ws_k <- ws_at_time(t_k, wind_tl)
    wd_k <- wd_at_time(t_k, wind_tl)
    zeta_k <- zeta_at_time(t_k, transition_duration_min, class_start, class_end)
    coefs_k <- interpolate_briggs(t_k, transition_duration_min, class_start, class_end)

    enh_ug_m3 <- gaussian_plume_disk_interp(
      rec_x = rec$x,
      rec_y = rec$y,
      rec_z = rec$z,
      source_radius = lake_radius,
      flux_ug_m2_s = flux_ug_s,
      u_ref = ws_k,
      wd_deg = wd_k,
      coefs = coefs_k
    )
    enh_ppb <- ug_m3_to_ppb(enh_ug_m3)

    noise_ppb <- rnorm(1, 0, precision_at_dwell)
    drift_ppb <- sensor_cfg$drift_ppb_hr * (t_k / 60)
    offset_ppb <- sensor_cfg$cal_offset_ppb
    measured_ppb <- enh_ppb + noise_ppb + drift_ppb + offset_ppb

    rows[[k]] <- data.table(
      receptor_id = rec$receptor_id,
      x = rec$x,
      y = rec$y,
      z = rec$z,
      azimuth_deg = rec$azimuth_deg,
      ring_type = rec$ring_type,
      time_min = window_time_min,
      sample_time_min = t_k,
      ch4_bg_ppb = 0,
      ch4_enhancement_ppb = round(enh_ppb, 4),
      ch4_total_ppb = round(enh_ppb, 4),
      ch4_measured_ppb = round(measured_ppb, 4),
      noise_ppb = round(noise_ppb, 4),
      drift_ppb = round(drift_ppb, 4),
      offset_ppb = offset_ppb,
      precision_ppb = round(precision_at_dwell, 4),
      dwell_s = dwell_s,
      ws_m_s = ws_k,
      ws_at_z_m_s = ws_k,
      wd_deg = wd_k,
      flux_true_mg_m2_d = flux_mg_m2_d,
      ustar_m_s = NA_real_,
      L_m = NA_real_,
      z0_m = 0.001,
      zeta = zeta_k,
      class_start = class_start,
      class_end = class_end,
      stability_class = "transitional"
    )
  }

  dt <- rbindlist(rows)
  list(
    truth_dt = copy(dt)[, ch4_measured_ppb := NULL],
    measurements_dt = dt,
    receptors = receptors
  )
}

build_repeated_profile_truth_and_measurements <- function(profile_cycle_min,
                                                          total_window_min,
                                                          design_cfg,
                                                          lake_diameter,
                                                          flux_mg_m2_d,
                                                          wind_tl,
                                                          transition_duration_min,
                                                          class_start,
                                                          class_end,
                                                          sensor_cfg,
                                                          seed) {
  n_full_suites <- max(1L, floor(total_window_min / profile_cycle_min + 1e-9))
  suite_starts_min <- (seq_len(n_full_suites) - 1L) * profile_cycle_min

  suite_list <- lapply(seq_along(suite_starts_min), function(idx) {
    suite_start <- suite_starts_min[idx]
    build_profile_truth_and_measurements(
      total_flight_min = profile_cycle_min,
      design_cfg = design_cfg,
      lake_diameter = lake_diameter,
      flux_mg_m2_d = flux_mg_m2_d,
      wind_tl = wind_tl,
      transition_duration_min = transition_duration_min,
      class_start = class_start,
      class_end = class_end,
      sensor_cfg = sensor_cfg,
      seed = seed + idx - 1L,
      window_time_min = suite_start,
      start_time_min = suite_start
    )
  })

  list(
    truth_dt = rbindlist(lapply(suite_list, `[[`, "truth_dt"), fill = TRUE),
    measurements_dt = rbindlist(lapply(suite_list, `[[`, "measurements_dt"), fill = TRUE),
    receptors = copy(suite_list[[1]]$receptors),
    n_full_suites = n_full_suites
  )
}

summarise_retrieval_windows <- function(flux_dt) {
  valid_idx <- flux_dt$qc_flag == "ok" & is.finite(flux_dt$flux_retrieved)
  valid_flux <- flux_dt$flux_retrieved[valid_idx]
  valid_enh <- flux_dt$mean_enhancement_ppb[valid_idx]
  valid_ws <- flux_dt$ws_at_z_mean_m_s[valid_idx]

  qc_summary <- if (any(valid_idx)) {
    "ok"
  } else if (all(flux_dt$qc_flag == "low_ws")) {
    "low_ws"
  } else {
    "invalid"
  }

  list(
    qc_flag = qc_summary,
    n_retrieval_windows = nrow(flux_dt),
    n_valid_retrieval_windows = sum(valid_idx),
    retrieval_window_valid_fraction = mean(valid_idx),
    flux_retrieved_mg_m2_d = if (length(valid_flux)) mean(valid_flux, na.rm = TRUE) else NA_real_,
    within_run_flux_sd_mg_m2_d = safe_sd(valid_flux),
    mean_enhancement_ppb = if (length(valid_enh)) mean(valid_enh, na.rm = TRUE) else NA_real_,
    ws_at_z_mean_m_s = if (length(valid_ws)) mean(valid_ws, na.rm = TRUE) else NA_real_,
    flux_true_mg_m2_d = flux_dt$flux_true[1]
  )
}

run_transition_mc <- function(mc_idx, tr_row, tr_idx) {
  wind_seed <- SEED_BASE + tr_idx * 100000L + mc_idx
  wind_tl <- generate_wind_timeline(
    transition_duration_min = TRANSITION_DURATION_MIN,
    met_dt_min = MET_DT_MIN,
    u_mean_m_s = U_MEAN_M_S,
    wd_mean_deg = WD_MEAN_DEG,
    ws_floor_m_s = WS_FLOOR_M_S,
    bestfit_interp = bestfit_interp,
    class_start = tr_row$class_start,
    class_end = tr_row$class_end,
    seed = wind_seed
  )

  n_cases_local <- length(MISSION_DURATIONS_MIN) + length(REPEAT_PROFILE_DURATIONS_MIN)
  run_rows_local <- vector("list", n_cases_local)
  row_idx_local <- 0L
  representative_local <- NULL

  for (dur_idx in seq_along(MISSION_DURATIONS_MIN)) {
    mission_duration_min <- MISSION_DURATIONS_MIN[dur_idx]
    scenario_label <- sprintf(
      "temporal_%s_%03dmin_mc%03d",
      tr_row$transition,
      as.integer(round(mission_duration_min)),
      mc_idx
    )

    flight_seed <- SEED_BASE + tr_idx * 100000L + mc_idx * 1000L + dur_idx
    flight <- build_profile_truth_and_measurements(
      total_flight_min = mission_duration_min,
      design_cfg = design_cfg,
      lake_diameter = lk$diameter_m,
      flux_mg_m2_d = BASE_FLUX_MG_D,
      wind_tl = wind_tl,
      transition_duration_min = TRANSITION_DURATION_MIN,
      class_start = tr_row$class_start,
      class_end = tr_row$class_end,
      sensor_cfg = sensor_cfg,
      seed = flight_seed
    )

    flux_dt <- quietly(retrieve_flux(
      meas_dt = flight$measurements_dt,
      lake_diameter = lk$diameter_m,
      distance_from_edge_m = 10,
      n_azimuth = design_cfg$n_azimuth,
      use_observed_ws = TRUE,
      ws_observation_error_m_s = 0,
      wd_uncertainty_deg = 0,
      use_height_avg_wind = TRUE,
      scenario_label = scenario_label,
      save_output = FALSE
    ))
    retrieval_summary <- summarise_retrieval_windows(flux_dt)

    row_idx_local <- row_idx_local + 1L
    run_rows_local[[row_idx_local]] <- data.table(
      transition = tr_row$transition,
      class_start = tr_row$class_start,
      class_end = tr_row$class_end,
      mc_id = mc_idx,
      experiment_family = "single_mission",
      experiment_label = "Single sequential mission",
      duration_axis_label = "Mission duration",
      mission_duration_min = mission_duration_min,
      observation_window_min = mission_duration_min,
      profile_cycle_min = NA_real_,
      analysis_duration_min = mission_duration_min,
      analysis_duration_h = mission_duration_min / 60,
      n_azimuth = design_cfg$n_azimuth,
      n_height_levels = length(design_cfg$heights_m),
      n_samples = nrow(flight$receptors),
      n_full_suites = 1L,
      n_retrieval_windows = retrieval_summary$n_retrieval_windows,
      n_valid_retrieval_windows = retrieval_summary$n_valid_retrieval_windows,
      retrieval_window_valid_fraction = retrieval_summary$retrieval_window_valid_fraction,
      dwell_s = unique(flight$measurements_dt$dwell_s)[1],
      precision_ppb = unique(flight$measurements_dt$precision_ppb)[1],
      flux_retrieved_mg_m2_d = retrieval_summary$flux_retrieved_mg_m2_d,
      within_run_flux_sd_mg_m2_d = retrieval_summary$within_run_flux_sd_mg_m2_d,
      flux_true_mg_m2_d = retrieval_summary$flux_true_mg_m2_d,
      qc_flag = retrieval_summary$qc_flag,
      mean_enhancement_ppb = retrieval_summary$mean_enhancement_ppb,
      ws_at_z_mean_m_s = retrieval_summary$ws_at_z_mean_m_s
    )

    if (is.null(representative_local) &&
        mc_idx == 1L &&
        mission_duration_min == max(MISSION_DURATIONS_MIN)) {
      representative_local <- list(
        measurements = copy(flight$measurements_dt),
        truth = copy(flight$truth_dt),
        receptors = copy(flight$receptors),
        scenario_label = scenario_label,
        mission_duration_min = mission_duration_min
      )
    }
  }

  for (cad_idx in seq_along(REPEAT_PROFILE_DURATIONS_MIN)) {
    profile_cycle_min <- REPEAT_PROFILE_DURATIONS_MIN[cad_idx]
    scenario_label <- sprintf(
      "temporal_repeat_%s_%03dmin_mc%03d",
      tr_row$transition,
      as.integer(round(profile_cycle_min)),
      mc_idx
    )

    repeat_seed <- SEED_BASE + tr_idx * 100000L + mc_idx * 1000L + 500L + cad_idx
    repeated_flight <- build_repeated_profile_truth_and_measurements(
      profile_cycle_min = profile_cycle_min,
      total_window_min = TRANSITION_DURATION_MIN,
      design_cfg = design_cfg,
      lake_diameter = lk$diameter_m,
      flux_mg_m2_d = BASE_FLUX_MG_D,
      wind_tl = wind_tl,
      transition_duration_min = TRANSITION_DURATION_MIN,
      class_start = tr_row$class_start,
      class_end = tr_row$class_end,
      sensor_cfg = sensor_cfg,
      seed = repeat_seed
    )

    flux_dt <- quietly(retrieve_flux(
      meas_dt = repeated_flight$measurements_dt,
      lake_diameter = lk$diameter_m,
      distance_from_edge_m = 10,
      n_azimuth = design_cfg$n_azimuth,
      use_observed_ws = TRUE,
      ws_observation_error_m_s = 0,
      wd_uncertainty_deg = 0,
      use_height_avg_wind = TRUE,
      scenario_label = scenario_label,
      save_output = FALSE
    ))
    retrieval_summary <- summarise_retrieval_windows(flux_dt)

    row_idx_local <- row_idx_local + 1L
    run_rows_local[[row_idx_local]] <- data.table(
      transition = tr_row$transition,
      class_start = tr_row$class_start,
      class_end = tr_row$class_end,
      mc_id = mc_idx,
      experiment_family = "repeat_profile",
      experiment_label = "Repeated full-profile mean",
      duration_axis_label = "Full-profile duration",
      mission_duration_min = NA_real_,
      observation_window_min = TRANSITION_DURATION_MIN,
      profile_cycle_min = profile_cycle_min,
      analysis_duration_min = profile_cycle_min,
      analysis_duration_h = profile_cycle_min / 60,
      n_azimuth = design_cfg$n_azimuth,
      n_height_levels = length(design_cfg$heights_m),
      n_samples = nrow(repeated_flight$receptors),
      n_full_suites = repeated_flight$n_full_suites,
      n_retrieval_windows = retrieval_summary$n_retrieval_windows,
      n_valid_retrieval_windows = retrieval_summary$n_valid_retrieval_windows,
      retrieval_window_valid_fraction = retrieval_summary$retrieval_window_valid_fraction,
      dwell_s = unique(repeated_flight$measurements_dt$dwell_s)[1],
      precision_ppb = unique(repeated_flight$measurements_dt$precision_ppb)[1],
      flux_retrieved_mg_m2_d = retrieval_summary$flux_retrieved_mg_m2_d,
      within_run_flux_sd_mg_m2_d = retrieval_summary$within_run_flux_sd_mg_m2_d,
      flux_true_mg_m2_d = retrieval_summary$flux_true_mg_m2_d,
      qc_flag = retrieval_summary$qc_flag,
      mean_enhancement_ppb = retrieval_summary$mean_enhancement_ppb,
      ws_at_z_mean_m_s = retrieval_summary$ws_at_z_mean_m_s
    )
  }

  list(
    run_rows = rbindlist(run_rows_local, fill = TRUE),
    representative = representative_local
  )
}

compute_summary <- function(run_dt, true_flux_mg_m2_d) {
  run_dt[, valid := qc_flag == "ok" & is.finite(flux_retrieved_mg_m2_d)]
  run_dt[, retrieval_bias_pct := 100 * (flux_retrieved_mg_m2_d - true_flux_mg_m2_d) / true_flux_mg_m2_d]

  run_dt[, {
    valid_flux <- flux_retrieved_mg_m2_d[valid]
    valid_bias <- retrieval_bias_pct[valid]
    scenario_seed <- string_to_seed(sprintf(
      "temporal|%s|%s|%.4f",
      transition[1], experiment_family[1], analysis_duration_h[1]
    ))
    flux_stats <- distribution_detail_cols(
      valid_flux,
      base_name = "retrieved_flux",
      suffix = "_mg_m2_d",
      seed = scenario_seed + 13L,
      include_abs_mean = FALSE
    )
    bias_stats <- distribution_detail_cols(
      valid_bias,
      base_name = "retrieval_bias",
      suffix = "_pct",
      seed = scenario_seed + 31L,
      include_abs_mean = TRUE
    )
    rmse_now <- if (length(valid_bias)) sqrt(mean(valid_bias^2, na.rm = TRUE)) else NA_real_
    bias_mean_now <- if (length(valid_bias)) mean(valid_bias, na.rm = TRUE) else NA_real_

    c(list(
      flux_true_mg_m2_d = true_flux_mg_m2_d,
      mc_runs = .N,
      n_valid = sum(valid),
      valid_fraction = mean(valid),
      mean_retrieval_window_valid_fraction = mean(retrieval_window_valid_fraction, na.rm = TRUE),
      retrieved_flux_mean_mg_m2_d = if (length(valid_flux)) mean(valid_flux, na.rm = TRUE) else NA_real_,
      retrieved_flux_sd_mg_m2_d = safe_sd(valid_flux),
      retrieval_bias_mean_pct = bias_mean_now,
      retrieval_bias_sd_pct = safe_sd(valid_bias),
      retrieval_rmse_pct = rmse_now,
      retrieval_bias_direction = bias_direction_label(bias_mean_now),
      rmse_target_met = is.finite(rmse_now) && rmse_now < RESULTS_RMSE_TARGET_PCT,
      abs_bias_target_met = is.finite(bias_mean_now) && abs(bias_mean_now) < RESULTS_ABS_BIAS_TARGET_PCT,
      mean_n_full_suites = mean(n_full_suites, na.rm = TRUE),
      mean_n_valid_retrieval_windows = mean(n_valid_retrieval_windows, na.rm = TRUE),
      mean_dwell_s = mean(dwell_s, na.rm = TRUE),
      mean_precision_ppb = mean(precision_ppb, na.rm = TRUE),
      mean_within_run_flux_sd_mg_m2_d = mean(within_run_flux_sd_mg_m2_d, na.rm = TRUE)
    ), flux_stats, bias_stats)
  }, by = .(
    transition,
    class_start,
    class_end,
    experiment_family,
    experiment_label,
    duration_axis_label,
    mission_duration_min,
    observation_window_min,
    profile_cycle_min,
    analysis_duration_min,
    analysis_duration_h,
    n_azimuth,
    n_height_levels,
    n_samples
  )][order(transition, factor(experiment_family, levels = c("single_mission", "repeat_profile")), analysis_duration_min)]
}

build_temporal_key_table <- function(summary_dt) {
  summary_cases <- unique(summary_dt[, .(transition, experiment_family)])
  setorder(summary_cases, transition)
  summary_cases[, family_order := match(experiment_family, c("single_mission", "repeat_profile"))]
  setorder(summary_cases, transition, family_order)
  summary_cases[, family_order := NULL]

  rbindlist(lapply(seq_len(nrow(summary_cases)), function(idx) {
    case_row <- summary_cases[idx]
    dt_case <- copy(summary_dt[
      transition == case_row$transition &
        experiment_family == case_row$experiment_family
    ])[order(analysis_duration_h)]
    shortest <- dt_case[1]
    longest <- dt_case[.N]
    rmse_target_row <- pick_threshold_row(
      dt_case,
      x_col = "analysis_duration_h",
      metric_col = "retrieval_rmse_pct",
      threshold = RESULTS_RMSE_TARGET_PCT,
      direction = "largest",
      abs_metric = FALSE,
      bias_col = "retrieval_bias_mean_pct"
    )
    bias_target_row <- pick_threshold_row(
      dt_case,
      x_col = "analysis_duration_h",
      metric_col = "retrieval_bias_mean_pct",
      threshold = RESULTS_ABS_BIAS_TARGET_PCT,
      direction = "largest",
      abs_metric = TRUE,
      bias_col = "retrieval_bias_mean_pct"
    )
    best_rmse_idx <- safe_which_min(dt_case$retrieval_rmse_pct)
    best_rmse_row <- if (is.finite(best_rmse_idx)) dt_case[best_rmse_idx] else dt_case[1]
    bias_change <- longest$retrieval_bias_mean_pct[1] - shortest$retrieval_bias_mean_pct[1]
    mean_sd <- mean(c(shortest$retrieval_bias_sd_pct[1], longest$retrieval_bias_sd_pct[1]), na.rm = TRUE)

    data.table(
      transition = case_row$transition,
      class_start = shortest$class_start[1],
      class_end = shortest$class_end[1],
      experiment_family = case_row$experiment_family,
      experiment_label = shortest$experiment_label[1],
      duration_axis_label = shortest$duration_axis_label[1],
      observation_window_h = longest$observation_window_min[1] / 60,
      shortest_duration_h = shortest$analysis_duration_h[1],
      longest_duration_h = longest$analysis_duration_h[1],
      bias_at_shortest_pct = shortest$retrieval_bias_mean_pct[1],
      bias_at_longest_pct = longest$retrieval_bias_mean_pct[1],
      rmse_at_shortest_pct = shortest$retrieval_rmse_pct[1],
      rmse_at_longest_pct = longest$retrieval_rmse_pct[1],
      valid_fraction_at_shortest = shortest$valid_fraction[1],
      valid_fraction_at_longest = longest$valid_fraction[1],
      mean_precision_shortest_ppb = shortest$mean_precision_ppb[1],
      mean_precision_longest_ppb = longest$mean_precision_ppb[1],
      most_accurate_duration_h = best_rmse_row$analysis_duration_h[1],
      min_rmse_pct = best_rmse_row$retrieval_rmse_pct[1],
      max_duration_rmse25_h = rmse_target_row$analysis_duration_h[1],
      max_duration_abs_bias10_h = bias_target_row$analysis_duration_h[1],
      bias_growth_pct_per_h = endpoint_slope(dt_case$analysis_duration_h, dt_case$retrieval_bias_mean_pct, log_x = FALSE),
      rmse_change_pct_points = longest$retrieval_rmse_pct[1] - shortest$retrieval_rmse_pct[1],
      valid_fraction_change = longest$valid_fraction[1] - shortest$valid_fraction[1],
      bias_change_to_sd_ratio = if (is.finite(mean_sd) && mean_sd > 0) abs(bias_change) / mean_sd else NA_real_,
      mean_n_full_suites_at_shortest = shortest$mean_n_full_suites[1],
      mean_n_full_suites_at_longest = longest$mean_n_full_suites[1]
    )
  }), fill = TRUE)
}

build_temporal_paired_bias <- function(summary_dt, transition_name) {
  single_dt <- summary_dt[
    transition == transition_name & experiment_family == "single_mission",
    .(
      transition,
      class_start,
      class_end,
      individual_duration_min = mission_duration_min,
      flux_true_mg_m2_d,
      single_mc_runs = mc_runs,
      single_flux_mean_mg_m2_d = retrieved_flux_mean_mg_m2_d,
      single_flux_sd_mg_m2_d = retrieved_flux_sd_mg_m2_d,
      single_flux_mean_pct_true = 100 * retrieved_flux_mean_mg_m2_d / flux_true_mg_m2_d,
      single_flux_sd_pct_true = 100 * retrieved_flux_sd_mg_m2_d / flux_true_mg_m2_d,
      single_bias_pct = retrieval_bias_mean_pct,
      single_bias_sd_pct = retrieval_bias_sd_pct
    )
  ]
  repeated_dt <- summary_dt[
    transition == transition_name & experiment_family == "repeat_profile",
    .(
      transition,
      individual_duration_min = profile_cycle_min,
      repeated_observation_window_min = observation_window_min,
      repeated_mean_n_full_suites = mean_n_full_suites,
      repeated_mc_runs = mc_runs,
      repeated_flux_mean_mg_m2_d = retrieved_flux_mean_mg_m2_d,
      repeated_flux_sd_mg_m2_d = retrieved_flux_sd_mg_m2_d,
      repeated_flux_mean_pct_true = 100 * retrieved_flux_mean_mg_m2_d / flux_true_mg_m2_d,
      repeated_flux_sd_pct_true = 100 * retrieved_flux_sd_mg_m2_d / flux_true_mg_m2_d,
      repeated_bias_pct = retrieval_bias_mean_pct,
      repeated_bias_sd_pct = retrieval_bias_sd_pct
    )
  ]

  merge(single_dt, repeated_dt, by = c("transition", "individual_duration_min"))
}

build_temporal_paper_stats <- function(summary_dt) {
  dt <- copy(summary_dt)
  primary_transitions <- c("B-D", "C-D", "C-F")
  rows <- list()

  add_rows <- function(x) {
    if (nrow(x)) rows[[length(rows) + 1L]] <<- x
  }

  single_30_dt <- dt[
    transition %in% primary_transitions &
      experiment_family == "single_mission" &
      mission_duration_min == 30
  ][order(transition)]
  add_rows(single_30_dt[, .(
    stat_id = "temporal_01_retrieved_flux_30min_single_flights",
    stat_description = "Retrieved flux for 30-minute single flights for each primary stability transition.",
    transition,
    class_start,
    class_end,
    experiment_family,
    mission_duration_min,
    analysis_duration_h,
    n_azimuth,
    n_height_levels,
    n_samples,
    flux_true_mg_m2_d,
    mc_runs,
    n_valid,
    retrieved_flux_mean_mg_m2_d,
    retrieved_flux_sd_mg_m2_d,
    retrieved_flux_mean_pct_true = 100 * retrieved_flux_mean_mg_m2_d / flux_true_mg_m2_d,
    retrieved_flux_sd_pct_true = 100 * retrieved_flux_sd_mg_m2_d / flux_true_mg_m2_d,
    retrieval_bias_mean_pct,
    retrieval_bias_sd_pct,
    value_text = sprintf(
      "%.2f +/- %.2f mg CH4 m^-2 d^-1; %.2f +/- %.2f%% of true flux",
      retrieved_flux_mean_mg_m2_d,
      retrieved_flux_sd_mg_m2_d,
      100 * retrieved_flux_mean_mg_m2_d / flux_true_mg_m2_d,
      100 * retrieved_flux_sd_mg_m2_d / flux_true_mg_m2_d
    ),
    note = "mean +/- std dev across Monte Carlo replicates"
  )])

  longest_single_dt <- dt[
    transition %in% primary_transitions & experiment_family == "single_mission"
  ][order(transition, -mission_duration_min), .SD[1], by = transition]
  add_rows(longest_single_dt[, .(
    stat_id = "temporal_02_bias_at_longest_single_flight",
    stat_description = "Accuracy loss at the longest single-flight duration tested.",
    transition,
    class_start,
    class_end,
    experiment_family,
    mission_duration_min,
    analysis_duration_h,
    flux_true_mg_m2_d,
    mc_runs,
    retrieval_bias_mean_pct,
    retrieval_bias_sd_pct,
    retrieved_flux_mean_mg_m2_d,
    retrieved_flux_sd_mg_m2_d,
    value_text = sprintf("%.2f%% bias at %.1f h", retrieval_bias_mean_pct, analysis_duration_h),
    note = "bias is 100 * (retrieved flux - true flux) / true flux"
  )])

  single_30_cd <- single_30_dt[transition == "C-D"]
  repeated_30_cd <- dt[
    transition == "C-D" &
      experiment_family == "repeat_profile" &
      profile_cycle_min == 30
  ]
  if (nrow(single_30_cd)) {
    srow <- single_30_cd[1]
    add_rows(data.table(
      stat_id = "temporal_03_precision_c_d_repeated_profile_improvement",
      stat_description = "Precision improvement with repeated 30-minute profiles for the C-D transition.",
      transition = "C-D",
      experiment_family = "single_mission",
      comparison_label = "single_30min_flight",
      mission_duration_min = 30,
      profile_cycle_min = NA_real_,
      total_profile_time_min = 30,
      n_profiles = 1,
      flux_true_mg_m2_d = srow$flux_true_mg_m2_d[1],
      mc_runs = srow$mc_runs[1],
      std_dev_mg_m2_d = srow$retrieved_flux_sd_mg_m2_d[1],
      std_dev_pct_true = 100 * srow$retrieved_flux_sd_mg_m2_d[1] / srow$flux_true_mg_m2_d[1],
      value_text = sprintf(
        "%.2f mg CH4 m^-2 d^-1 (%.2f%% of true flux)",
        srow$retrieved_flux_sd_mg_m2_d[1],
        100 * srow$retrieved_flux_sd_mg_m2_d[1] / srow$flux_true_mg_m2_d[1]
      ),
      note = "direct value from the 30-minute single-flight summary row"
    ))
    add_rows(data.table(
      stat_id = "temporal_03_precision_c_d_repeated_profile_improvement",
      stat_description = "Precision improvement with repeated 30-minute profiles for the C-D transition.",
      transition = "C-D",
      experiment_family = "repeat_profile_estimated",
      comparison_label = "five_30min_profiles_2.5h",
      mission_duration_min = NA_real_,
      profile_cycle_min = 30,
      total_profile_time_min = 150,
      n_profiles = 5,
      flux_true_mg_m2_d = srow$flux_true_mg_m2_d[1],
      mc_runs = srow$mc_runs[1],
      std_dev_mg_m2_d = srow$retrieved_flux_sd_mg_m2_d[1] / sqrt(5),
      std_dev_pct_true = 100 * (srow$retrieved_flux_sd_mg_m2_d[1] / sqrt(5)) / srow$flux_true_mg_m2_d[1],
      value_text = sprintf(
        "%.2f mg CH4 m^-2 d^-1 (%.2f%% of true flux)",
        srow$retrieved_flux_sd_mg_m2_d[1] / sqrt(5),
        100 * (srow$retrieved_flux_sd_mg_m2_d[1] / sqrt(5)) / srow$flux_true_mg_m2_d[1]
      ),
      note = "estimated as single 30-minute std dev divided by sqrt(5); the cached summary does not contain a direct 2.5-hour repeated-profile scenario"
    ))
  }
  if (nrow(repeated_30_cd)) {
    rrow <- repeated_30_cd[1]
    add_rows(data.table(
      stat_id = "temporal_03_precision_c_d_repeated_profile_improvement",
      stat_description = "Precision improvement with repeated 30-minute profiles for the C-D transition.",
      transition = "C-D",
      experiment_family = "repeat_profile",
      comparison_label = "twenty_four_30min_profiles_12h",
      mission_duration_min = NA_real_,
      profile_cycle_min = 30,
      observation_window_min = rrow$observation_window_min[1],
      total_profile_time_min = rrow$observation_window_min[1],
      n_profiles = rrow$mean_n_full_suites[1],
      flux_true_mg_m2_d = rrow$flux_true_mg_m2_d[1],
      mc_runs = rrow$mc_runs[1],
      std_dev_mg_m2_d = rrow$retrieved_flux_sd_mg_m2_d[1],
      std_dev_pct_true = 100 * rrow$retrieved_flux_sd_mg_m2_d[1] / rrow$flux_true_mg_m2_d[1],
      value_text = sprintf(
        "%.2f mg CH4 m^-2 d^-1 (%.2f%% of true flux)",
        rrow$retrieved_flux_sd_mg_m2_d[1],
        100 * rrow$retrieved_flux_sd_mg_m2_d[1] / rrow$flux_true_mg_m2_d[1]
      ),
      note = "direct value from the 30-minute repeated-profile summary row over the full 12-hour transition"
    ))
  }

  for (tr in primary_transitions) {
    pair_dt <- build_temporal_paired_bias(dt, tr)
    if (!nrow(pair_dt)) next
    pair_dt[, abs_single_bias_pct := abs(single_bias_pct)]
    best_pair <- pair_dt[which.min(abs_single_bias_pct)]
    add_rows(best_pair[, .(
      stat_id = "temporal_04_bias_reduction_single_to_repeated_at_min_bias_duration",
      stat_description = "Mean bias reduction when switching from single to repeated profiles at the common individual duration that minimizes single-flight bias.",
      transition,
      class_start,
      class_end,
      individual_duration_min,
      flux_true_mg_m2_d,
      single_bias_pct,
      repeated_bias_pct,
      signed_bias_change_pct = repeated_bias_pct - single_bias_pct,
      abs_bias_reduction_pct = abs(single_bias_pct) - abs(repeated_bias_pct),
      single_flux_mean_mg_m2_d,
      repeated_flux_mean_mg_m2_d,
      value_text = sprintf(
        "single %.2f%%; repeated %.2f%%; abs-bias reduction %.2f percentage points",
        single_bias_pct,
        repeated_bias_pct,
        abs(single_bias_pct) - abs(repeated_bias_pct)
      ),
      note = "duration selected from durations present in both single and repeated summaries"
    )])

    ordered_pair <- pair_dt[order(individual_duration_min)]
    selected_idx <- unique(c(1L, nrow(ordered_pair)))
    selected_pair <- ordered_pair[selected_idx]
    selected_pair[, comparison_label := ifelse(
      individual_duration_min == min(ordered_pair$individual_duration_min),
      "shortest_individual_flight",
      "longest_individual_flight"
    )]
    add_rows(selected_pair[, .(
      stat_id = "temporal_05_pairwise_bias_single_vs_repeated_shortest_longest",
      stat_description = "Pairwise bias comparison for single versus repeated profiles at the shortest and longest common individual flight durations.",
      transition,
      class_start,
      class_end,
      comparison_label,
      individual_duration_min,
      flux_true_mg_m2_d,
      single_bias_pct,
      repeated_bias_pct,
      signed_bias_change_pct = repeated_bias_pct - single_bias_pct,
      abs_bias_reduction_pct = abs(single_bias_pct) - abs(repeated_bias_pct),
      value_text = sprintf(
        "%s %.0f min: single %.2f%%; repeated %.2f%%; change %.2f percentage points",
        comparison_label,
        individual_duration_min,
        single_bias_pct,
        repeated_bias_pct,
        repeated_bias_pct - single_bias_pct
      ),
      note = "common durations are the individual durations available in both single-mission and repeated-profile summaries"
    )])

    pair_30 <- pair_dt[individual_duration_min == 30]
    add_rows(pair_30[, .(
      stat_id = "temporal_06_marginal_accuracy_gain_multiple_30min_profiles",
      stat_description = "Marginal accuracy gain from flying multiple 30-minute profiles versus a single 30-minute flight.",
      transition,
      class_start,
      class_end,
      individual_duration_min,
      repeated_observation_window_min,
      repeated_mean_n_full_suites,
      flux_true_mg_m2_d,
      single_bias_pct,
      repeated_bias_pct,
      signed_bias_change_pct = repeated_bias_pct - single_bias_pct,
      abs_bias_reduction_pct = abs(single_bias_pct) - abs(repeated_bias_pct),
      value_text = sprintf(
        "single %.2f%%; repeated %.2f%%; abs-bias reduction %.2f percentage points",
        single_bias_pct,
        repeated_bias_pct,
        abs(single_bias_pct) - abs(repeated_bias_pct)
      ),
      note = "repeated value is the 30-minute full-profile cadence averaged over the full transition window"
    )])
  }

  repeated_pool_dt <- dt[
    transition %in% primary_transitions &
      experiment_family == "repeat_profile" &
      profile_cycle_min >= 30 &
      profile_cycle_min <= 120
  ]
  add_rows(repeated_pool_dt[, .(
    stat_id = "temporal_07_average_sd_repeated_profile_strategies",
    stat_description = "Average std dev of retrieved flux across repeated-profile strategies from 30 min to 2 h, pooled across primary transitions.",
    n_strategy_transition_rows = .N,
    profile_cycle_min_min = min(profile_cycle_min, na.rm = TRUE),
    profile_cycle_min_max = max(profile_cycle_min, na.rm = TRUE),
    avg_std_dev_mg_m2_d = mean(retrieved_flux_sd_mg_m2_d, na.rm = TRUE),
    sd_of_std_dev_mg_m2_d = stats::sd(retrieved_flux_sd_mg_m2_d, na.rm = TRUE),
    avg_std_dev_pct_true = mean(100 * retrieved_flux_sd_mg_m2_d / flux_true_mg_m2_d, na.rm = TRUE),
    value_text = sprintf(
      "%.2f mg CH4 m^-2 d^-1 average std dev (%.2f%% of true flux)",
      mean(retrieved_flux_sd_mg_m2_d, na.rm = TRUE),
      mean(100 * retrieved_flux_sd_mg_m2_d / flux_true_mg_m2_d, na.rm = TRUE)
    ),
    note = "pooled across B-D, C-D, and C-F repeated-profile rows with 30-, 60-, 90-, and 120-minute full-profile durations"
  )])

  radius_m <- 60
  n_azimuth <- 96
  total_receptors <- 2592
  circumference_m <- 2 * pi * radius_m
  dwell_assumed_s <- 17
  dwell_implied_s <- 30 * 60 / total_receptors
  dwell_values <- c(dwell_assumed_s, dwell_implied_s)
  assumption_labels <- c("stated_17s_dwell", "implied_30min_over_2592_receptors")
  orbit_time_s <- n_azimuth * dwell_values
  ground_speed_kmh <- circumference_m / orbit_time_s * 3.6
  add_rows(data.table(
    stat_id = "temporal_08_uas_ground_speed_one_orbit",
    stat_description = "UAS ground speed required to complete one full orbit at one altitude level.",
    assumption_label = assumption_labels,
    radius_m = radius_m,
    n_azimuth = n_azimuth,
    total_receptors = total_receptors,
    dwell_s = dwell_values,
    orbit_time_s = orbit_time_s,
    ground_speed_kmh = ground_speed_kmh,
    value_text = sprintf("%.2f km/h", ground_speed_kmh),
    note = c(
      "computed from 17 s per receptor: speed = 2*pi*radius / (96*dwell)",
      "computed from the stated 30-min mission divided by 2592 total receptors; this implies 0.694 s per receptor, not 17 s"
    )
  ))

  rbindlist(rows, fill = TRUE)
}

write_temporal_paper_stats <- function(summary_dt, output_csv = PAPER_STATS_CSV) {
  paper_stats_dt <- build_temporal_paper_stats(summary_dt)
  fwrite(paper_stats_dt, output_csv)
  message(sprintf("Saved paper stats CSV: %s", normalizePath(output_csv, mustWork = FALSE)))
  invisible(paper_stats_dt)
}

if (isTRUE(PAPER_STATS_ONLY)) {
  if (!file.exists(SUMMARY_CSV)) {
    stop("Paper-stats-only mode requested, but summary CSV was not found: ", SUMMARY_CSV)
  }
  summary_dt <- fread(SUMMARY_CSV)
  write_temporal_paper_stats(summary_dt, PAPER_STATS_CSV)
  quit(save = "no", status = 0)
}


# -----------------------------------------------------------------------------
# 3. Experiment setup
# -----------------------------------------------------------------------------
lk <- LAKE_CONFIGS[[LAKE_NAME]]
design_cfg <- SAMPLING_CONFIGS[[DESIGN_NAME]]
sensor_cfg <- SENSOR_CONFIGS[[SENSOR_NAME]]
sensor_cfg$drift_ppb_hr <- 0
sensor_cfg$cal_offset_ppb <- 0
transition_dt <- rbindlist(lapply(TRANSITIONS, function(tr) {
  data.table(
    class_start = tr[1],
    class_end = tr[2],
    transition = paste(tr, collapse = "-")
  )
}))
total_mc_jobs <- nrow(transition_dt) * N_MC
n_samples_total <- design_cfg$n_azimuth * length(design_cfg$heights_m)
all_profile_durations_min <- sort(unique(c(MISSION_DURATIONS_MIN, REPEAT_PROFILE_DURATIONS_MIN)))
dwell_range_s <- sort(all_profile_durations_min * 60 / n_samples_total)
available_cores <- parallel::detectCores(logical = TRUE)
if (!is.finite(available_cores) || available_cores < 1) available_cores <- 1L
PARALLEL_WORKER_MULTIPLIER <- 2L
N_WORKERS_BASE <- if (is.finite(N_WORKERS_REQUESTED) && !is.na(N_WORKERS_REQUESTED)) {
  max(1L, as.integer(N_WORKERS_REQUESTED))
} else {
  max(1L, available_cores - 1L)
}
N_WORKERS <- max(1L, PARALLEL_WORKER_MULTIPLIER * N_WORKERS_BASE)
N_WORKERS <- min(N_WORKERS, total_mc_jobs)
PROGRESS_WAVES_TARGET <- min(total_mc_jobs, 4L)
MC_BATCH_SIZE <- min(N_WORKERS, max(1L, ceiling(total_mc_jobs / PROGRESS_WAVES_TARGET)))

message("========================================")
message(" Temporal Measurement Analysis")
message("========================================")
message(sprintf(" Lake:        %s (%d m diameter)", LAKE_NAME, lk$diameter_m))
message(sprintf(" Lake grid:   %.2f m", LAKE_GRID_RESOLUTION_M))
message(sprintf(" Flux:        %.1f mg CH4 m^-2 d^-1", BASE_FLUX_MG_D))
message(sprintf(" Sensor:      %s", sensor_cfg$name))
message(sprintf(" Wind speed:  %.1f m/s (fixed true wind)", U_MEAN_M_S))
message(" Sensor err:  drift = 0 ppb/hr, calibration offset = 0 ppb")
message(sprintf(" Sampling:    %d azimuths x %d heights = %d receptors", design_cfg$n_azimuth, length(design_cfg$heights_m), n_samples_total))
message(sprintf(" Heights:     %.2f to %.2f m", min(design_cfg$heights_m), max(design_cfg$heights_m)))
message(sprintf(" Dwell:       %.2f to %.2f s per receptor across requested tests", min(dwell_range_s), max(dwell_range_s)))
message(sprintf(" Transition:  %.1f h", TRANSITION_DURATION_MIN / 60))
message(sprintf(" Missions:    %s min", paste(MISSION_DURATIONS_MIN, collapse = ", ")))
message(sprintf(" Repeats:     %s min per full profile", paste(REPEAT_PROFILE_DURATIONS_MIN, collapse = ", ")))
message(sprintf(" Monte Carlo: %d repeats", N_MC))
message(sprintf(" Parallel:    %d worker%s (%dx base request of %d) from %d core%s | %d MC run%s per batch",
                N_WORKERS, ifelse(N_WORKERS == 1L, "", "s"),
                PARALLEL_WORKER_MULTIPLIER, N_WORKERS_BASE,
                available_cores, ifelse(available_cores == 1L, "", "s"),
                MC_BATCH_SIZE, ifelse(MC_BATCH_SIZE == 1L, "", "s")))
message(sprintf(" Cases:       %s", paste(TRANSITION_LABELS, collapse = ", ")))
message("========================================")

emission_result <- configure_main_lake_emissions(
  lake_name = LAKE_NAME,
  base_flux_mg_m2_d = BASE_FLUX_MG_D,
  lake_grid_resolution_m = LAKE_GRID_RESOLUTION_M
)

bestfit_interp <- load_bestfit_interpolators(BESTFIT_PATH)

cl <- NULL
if (N_WORKERS > 1L) {
  worker_outfile <- tempfile(pattern = "temporal_parallel_workers_", fileext = ".log")
  cl <- parallel::makePSOCKcluster(
    N_WORKERS,
    outfile = worker_outfile,
    rscript_args = c("--vanilla")
  )
  on.exit(parallel::stopCluster(cl), add = TRUE)
  on.exit(unlink(worker_outfile, force = TRUE), add = TRUE)
  parallel::clusterExport(cl, varlist = c("proj_root"), envir = environment())
  parallel::clusterEvalQ(cl, {
    options(uas_gauss_osse_root = proj_root)
    suppressMessages(source(file.path(proj_root, "R", "utils.R")))
    suppressMessages(source(file.path(proj_root, "R", "01_synthetic_emissions.R")))
    suppressMessages(source(file.path(proj_root, "R", "02b_gaussian_forward.R")))
    suppressMessages(source(file.path(proj_root, "R", "04_conley_inverse_flux.R")))
    suppressMessages(source(file.path(proj_root, "scripts", "full_analysis", "results_reporting_helpers.R")))
    NULL
  })
  parallel::clusterExport(
    cl,
    varlist = c(
      "PASQUILL_ZETA_ANCHORS",
      "quietly",
      "interpolate_briggs",
      "zeta_at_time",
      "generate_wind_timeline",
      "ws_at_time",
      "wd_at_time",
      "pg_sigma_y_interp",
      "pg_sigma_z_interp",
      "gaussian_plume_point_interp",
      "gaussian_plume_disk_interp",
      "build_profile_truth_and_measurements",
      "build_repeated_profile_truth_and_measurements",
      "summarise_retrieval_windows",
      "safe_sd",
      "run_transition_mc",
      "MISSION_DURATIONS_MIN",
      "REPEAT_PROFILE_DURATIONS_MIN",
      "design_cfg",
      "lk",
      "sensor_cfg",
      "BASE_FLUX_MG_D",
      "TRANSITION_DURATION_MIN",
      "MET_DT_MIN",
      "U_MEAN_M_S",
      "WD_MEAN_DEG",
      "WS_FLOOR_M_S",
      "SEED_BASE",
      "bestfit_interp"
    ),
    envir = environment()
  )
}

job_dt <- rbindlist(lapply(seq_len(nrow(transition_dt)), function(tr_idx) {
  cbind(
    data.table(tr_idx = tr_idx, mc_idx = seq_len(N_MC)),
    transition_dt[tr_idx]
  )
}), fill = TRUE)
job_chunks <- split(seq_len(nrow(job_dt)), ceiling(seq_len(nrow(job_dt)) / MC_BATCH_SIZE))
run_pb <- make_progress(
  nrow(transition_dt) * N_MC * (length(MISSION_DURATIONS_MIN) + length(REPEAT_PROFILE_DURATIONS_MIN)),
  label = sprintf("Temporal sweep (%d worker%s, %d-run batches)", N_WORKERS, ifelse(N_WORKERS == 1L, "", "s"), MC_BATCH_SIZE)
)
run_batches <- vector("list", length(job_chunks))
batch_idx <- 0L

representative <- NULL


# -----------------------------------------------------------------------------
# 4. Run sweep
# -----------------------------------------------------------------------------
for (job_chunk in job_chunks) {
  chunk_specs <- lapply(job_chunk, function(job_idx) as.list(job_dt[job_idx]))
  chunk_results <- if (!is.null(cl)) {
    parallel::parLapplyLB(
      cl,
      chunk_specs,
      function(job_spec) {
        tr_row <- list(
          class_start = job_spec$class_start,
          class_end = job_spec$class_end,
          transition = job_spec$transition
        )
        run_transition_mc(job_spec$mc_idx, tr_row, job_spec$tr_idx)
      }
    )
  } else {
    lapply(chunk_specs, function(job_spec) {
      tr_row <- list(
        class_start = job_spec$class_start,
        class_end = job_spec$class_end,
        transition = job_spec$transition
      )
      run_transition_mc(job_spec$mc_idx, tr_row, job_spec$tr_idx)
    })
  }

  batch_idx <- batch_idx + 1L
  batch_rows <- rbindlist(lapply(chunk_results, `[[`, "run_rows"), fill = TRUE)
  run_batches[[batch_idx]] <- batch_rows

  if (is.null(representative)) {
    rep_candidates <- Filter(Negate(is.null), lapply(chunk_results, `[[`, "representative"))
    if (length(rep_candidates)) {
      representative <- rep_candidates[[1]]
    }
  }

  for (tick_idx in seq_len(nrow(batch_rows))) {
    run_pb$tick()
  }
}

run_dt <- rbindlist(run_batches, fill = TRUE)
summary_dt <- compute_summary(run_dt, true_flux_mg_m2_d = BASE_FLUX_MG_D)
key_dt <- build_temporal_key_table(summary_dt)
paper_stats_dt <- write_temporal_paper_stats(summary_dt, PAPER_STATS_CSV)


# -----------------------------------------------------------------------------
# 5. Save tables
# -----------------------------------------------------------------------------
fwrite(run_dt, RUNS_CSV)
fwrite(summary_dt, SUMMARY_CSV)
fwrite(key_dt, KEY_CSV)
message(sprintf("Saved run-level CSV: %s", normalizePath(RUNS_CSV, mustWork = FALSE)))
message(sprintf("Saved summary CSV:   %s", normalizePath(SUMMARY_CSV, mustWork = FALSE)))
message(sprintf("Saved key table:     %s", normalizePath(KEY_CSV, mustWork = FALSE)))
if (isTRUE(maybe_write_full_analysis_synthesis(OUTPUT_DIR))) {
  message("Updated cross-analysis synthesis table.")
}


# -----------------------------------------------------------------------------
# 6. Build composite figure
# -----------------------------------------------------------------------------
if (is.null(representative)) {
  stop("Representative run was not captured; cannot build composite figure.")
}

iso_limits <- compute_isometric_scale_limits(
  truth_dt = representative$truth,
  measurements_dt = representative$measurements,
  receptors = representative$receptors,
  flux_field = emission_result$flux_field,
  scenario_label = representative$scenario_label,
  ws_observation_error_m_s = 0
)

representative_plot_measurements <- representative$measurements[, .(
  receptor_id,
  time_min,
  ch4_bg_ppb,
  ch4_measured_ppb,
  x_measured = x,
  y_measured = y,
  z_measured = z
)]

left_plot <- plot_measured_isometric_over_emission_map(
  measurements_dt = representative_plot_measurements,
  receptors = representative$receptors,
  flux_field = emission_result$flux_field,
  lake_diameter = lk$diameter_m,
  time_min = 0,
  scenario_label = representative$scenario_label,
  value_mode = "enhancement",
  fill_limits = iso_limits$fill,
  colour_limits = iso_limits$measured_enhancement,
  emission_underlay_mode = "polar",
  polar_n_rings = 20L,
  polar_n_azimuth = 36L,
  scale_point_alpha_by_value = TRUE,
  point_alpha_floor = 0.28,
  point_alpha_ceiling = 0.98,
  point_size = 3.0,
  base_size = FIG_BASE_SIZE,
  save = FALSE
) +
  labs(
    title = "a) Representative measured methane enhancement",
    subtitle = sprintf(
      "%s transition | %.1f h mission",
      sub("^temporal_", "", sub("_\\d{3}min_mc\\d{3}$", "", representative$scenario_label)),
      representative$mission_duration_min / 60
    )
  ) +
  guides(fill = "none", shape = "none") +
  theme(
    legend.position = "bottom",
    legend.box.margin = margin(0, 0, 0, 0),
    legend.margin = margin(0, 0, 0, 0),
    legend.key.height = unit(0.28, "cm"),
    legend.key.width = unit(0.62, "cm"),
    legend.spacing.x = unit(0.04, "cm"),
    legend.spacing.y = unit(0.01, "cm")
  )

plot_dt <- copy(summary_dt)
plot_dt[, transition_facet := factor(
  transition,
  levels = TRANSITION_LABELS
)]
composite_transition_labels <- intersect(PRIMARY_TRANSITION_LABELS, TRANSITION_LABELS)
if (!length(composite_transition_labels)) {
  composite_transition_labels <- TRANSITION_LABELS
}
composite_plot_dt <- copy(plot_dt[transition %in% composite_transition_labels])
composite_plot_dt[, transition_facet := factor(transition, levels = composite_transition_labels)]

temporal_series_colour <- "#4E5D6C"
temporal_ribbon_fill <- "#C7D1DB"

build_temporal_flux_plot <- function(plot_data,
                                     title,
                                     subtitle,
                                     x_label,
                                     x_breaks) {
  ggplot(
    plot_data,
    aes(x = analysis_duration_h, y = retrieved_flux_mean_mg_m2_d)
  ) +
    geom_ribbon(
      aes(ymin = retrieved_flux_mean_mg_m2_d - retrieved_flux_sd_mg_m2_d,
          ymax = retrieved_flux_mean_mg_m2_d + retrieved_flux_sd_mg_m2_d),
      fill = temporal_ribbon_fill,
      alpha = 0.32, linewidth = 0, show.legend = FALSE
    ) +
    geom_hline(yintercept = BASE_FLUX_MG_D, linetype = "dashed", colour = "black", linewidth = 0.8) +
    geom_line(linewidth = FIG_LINEWIDTH, colour = temporal_series_colour, show.legend = FALSE) +
    geom_point(size = FIG_POINT_SIZE, colour = temporal_series_colour, show.legend = FALSE) +
    facet_wrap(~transition_facet, ncol = 1) +
    labs(
      title = title,
      subtitle = subtitle,
      x = x_label,
      y = expression("Retrieved flux (mg CH"[4]~m^-2~d^-1*")")
    ) +
    scale_x_continuous(breaks = x_breaks) +
    theme_osse(base_size = FIG_BASE_SIZE) +
    theme(
      legend.position = "none",
      strip.text = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )
}

build_repeated_temporal_flux_plot <- function(plot_data,
                                              title,
                                              subtitle,
                                              x_label,
                                              x_breaks) {
  plot_data <- copy(plot_data)
  plot_data[, suite_count_label := factor(
    sprintf("%d profiles", as.integer(round(mean_n_full_suites))),
    levels = sprintf(
      "%d profiles",
      sort(unique(as.integer(round(mean_n_full_suites))))
    )
  )]

  ggplot(
    plot_data,
    aes(x = analysis_duration_h, y = retrieved_flux_mean_mg_m2_d)
  ) +
    geom_ribbon(
      aes(ymin = retrieved_flux_mean_mg_m2_d - retrieved_flux_sd_mg_m2_d,
          ymax = retrieved_flux_mean_mg_m2_d + retrieved_flux_sd_mg_m2_d),
      fill = temporal_ribbon_fill,
      alpha = 0.32, linewidth = 0, show.legend = FALSE
    ) +
    geom_hline(yintercept = BASE_FLUX_MG_D, linetype = "dashed", colour = "black", linewidth = 0.8) +
    geom_line(linewidth = FIG_LINEWIDTH, colour = temporal_series_colour, show.legend = FALSE) +
    geom_point(
      aes(shape = suite_count_label),
      size = FIG_POINT_SIZE + 0.2,
      stroke = 0.55,
      colour = temporal_series_colour
    ) +
    facet_wrap(~transition_facet, ncol = 1) +
    labs(
      title = title,
      subtitle = subtitle,
      x = x_label,
      y = expression("Retrieved flux (mg CH"[4]~m^-2~d^-1*")"),
      shape = "Profiles in mean"
    ) +
    scale_x_continuous(breaks = x_breaks) +
    theme_osse(base_size = FIG_BASE_SIZE) +
    theme(
      legend.position = "right",
      strip.text = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    ) +
    guides(
      shape = guide_legend(override.aes = list(colour = "#333333", fill = "white"))
    )
}

single_plot_dt <- composite_plot_dt[experiment_family == "single_mission"]
repeat_plot_dt <- composite_plot_dt[experiment_family == "repeat_profile"]
column_subtitle <- sprintf(
  "%d azimuths x %d heights | %.1f h transitions | mean +/- 1 stdev across %d repetitions",
  design_cfg$n_azimuth,
  length(design_cfg$heights_m),
  TRANSITION_DURATION_MIN / 60,
  N_MC
)

middle_plot <- build_temporal_flux_plot(
  plot_data = single_plot_dt,
  title = "b) Retrieved flux vs mission duration",
  subtitle = column_subtitle,
  x_label = "Mission duration (hours)",
  x_breaks = seq(min(MISSION_DURATIONS_MIN) / 60, max(MISSION_DURATIONS_MIN) / 60, by = 1)
)

right_plot <- build_repeated_temporal_flux_plot(
  plot_data = repeat_plot_dt,
  title = "c) Retrieved flux vs repeated full-profile duration",
  subtitle = column_subtitle,
  x_label = "Time per full profile (hours)",
  x_breaks = sort(unique(REPEAT_PROFILE_DURATIONS_MIN / 60))
)

left_legend <- cowplot::get_legend(left_plot)
left_main <- left_plot + theme(legend.position = "none")

left_column <- if (!is.null(left_legend)) {
  cowplot::plot_grid(
    left_main,
    left_legend,
    ncol = 1,
    rel_heights = c(1, 0.09)
  )
} else {
  left_main
}

composite_plot <- cowplot::plot_grid(
  left_column,
  middle_plot,
  right_plot,
  ncol = 3,
  rel_widths = c(1.02, 1.15, 1.15)
)

ggsave_osse(FIG_PNG, composite_plot, width = 20, height = 11, dpi = 220)
message(sprintf("Saved composite figure: %s", normalizePath(FIG_PNG, mustWork = FALSE)))

additional_transition_labels <- intersect(ADDITIONAL_TRANSITION_LABELS, TRANSITION_LABELS)
additional_plot_dt <- copy(plot_dt[transition %in% additional_transition_labels])
if (nrow(additional_plot_dt)) {
  additional_plot_dt[, transition_facet := factor(transition, levels = additional_transition_labels)]
  single_x_breaks <- seq(
    min(MISSION_DURATIONS_MIN) / 60,
    max(MISSION_DURATIONS_MIN) / 60,
    by = 1
  )
  repeat_x_breaks <- sort(unique(REPEAT_PROFILE_DURATIONS_MIN / 60))

  additional_single_plot <- build_temporal_flux_plot(
    plot_data = additional_plot_dt[experiment_family == "single_mission"],
    title = "a) Retrieved flux vs mission duration",
    subtitle = column_subtitle,
    x_label = "Mission duration (hours)",
    x_breaks = single_x_breaks
  )

  additional_repeat_plot <- build_repeated_temporal_flux_plot(
    plot_data = additional_plot_dt[experiment_family == "repeat_profile"],
    title = "b) Retrieved flux vs repeated full-profile duration",
    subtitle = column_subtitle,
    x_label = "Time per full profile (hours)",
    x_breaks = repeat_x_breaks
  )

  additional_transition_plot <- cowplot::plot_grid(
    additional_single_plot,
    additional_repeat_plot,
    ncol = 2,
    align = "hv"
  )

  ggsave_osse(ADDITIONAL_FIG_PNG, additional_transition_plot, width = 14, height = 8.5, dpi = 220)
  message(sprintf("Saved additional transition figure: %s", normalizePath(ADDITIONAL_FIG_PNG, mustWork = FALSE)))
} else {
  message("No F-D or D-B results found; skipped additional transition figure.")
}

message("Analysis complete.")
