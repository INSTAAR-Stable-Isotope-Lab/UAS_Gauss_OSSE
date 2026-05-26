###############################################################################
# 03_drone_and_sensor_sim.R - Simulate drone flight patterns + sensor noise,
#                              drift, GPS error
#
# Scientific notes:
#   - Ebullition is treated as part of an effective, lake-averaged surface
#     flux.  We do NOT resolve hotspots or fine-scale spatial structure.
#   - The OSSE is a publication companion for Gaussian-forward /
#     Conley-inverse manuscript experiments.
#   - Instrument models are idealised, focusing on precision and drift.
#     More complex behaviours (thermal transients, pressure effects) are
#     neglected.
#
# Inputs:
#   - "Truth" concentrations from 02b_gaussian_forward.R
#   - Sensor specifications from utils.R (SENSOR_CONFIGS)
#
# Optional outputs when save_output = TRUE:
#   - data/raw/simulated_measurements_<scenario>.csv
###############################################################################

if (!exists("project_root", mode = "function")) {
  proj_root <- if (file.exists("R/utils.R")) "." else ".."
  source(file.path(proj_root, "R", "utils.R"))
}

# -----------------------------------------------------------------------------
# Core: apply sensor model to truth concentrations
# -----------------------------------------------------------------------------
#' Simulate drone-measured CH4 concentrations from truth concentrations.
#'
#' Applies:
#'   1. Measurement noise (Gaussian, scaled by precision and averaging time).
#'   2. Linear drift over the elapsed time.
#'   3. Fixed calibration offset.
#'   4. GPS position errors (horizontal and vertical), optional.
#'
#' @param truth_dt  data.table of truth concentrations (from run_forward_model).
#' @param sensor_cfg  Sensor configuration list (from SENSOR_CONFIGS).
#' @param design_cfg  Sampling design list (from SAMPLING_CONFIGS).
#' @param add_position_uncertainty  If FALSE (default), use true positions
#'        (x_measured = x, etc.) so inverse uses same positions as forward.
#'        If TRUE, add GPS errors to simulate realistic drone position uncertainty.
#' @param gps_h_error_m  Horizontal GPS error std dev (default 1.5 m), used when
#'        add_position_uncertainty is TRUE.
#' @param gps_v_error_m  Vertical GPS error std dev (default 2.5 m), used when
#'        add_position_uncertainty is TRUE.
#' @param seed  Random seed for reproducibility.
#' @return data.table with simulated measurements and error columns.
apply_sensor_model <- function(truth_dt,
                               sensor_cfg,
                               design_cfg,
                               add_position_uncertainty = FALSE,
                               gps_h_error_m = 1.5,
                               gps_v_error_m = 2.5,
                               seed = 42) {

  set.seed(seed)
  dt <- copy(truth_dt)
  n <- nrow(dt)

  # --- 1. Measurement noise ------------------------------------------------
  # Effective precision depends on dwell time at each receptor. Use actual
  # receptor count (includes outer, inner, and background rings) so dwell_s
  # reflects the full flight cycle.
  n_receptors <- n
  dwell_s <- design_cfg$cycle_min * 60 / n_receptors
  dwell_s <- max(dwell_s, 1)

  # Allan-variance scaling: precision improves as sqrt(averaging_time)
  precision_at_dwell <- sensor_cfg$precision_1s_ppb / sqrt(dwell_s)
  precision_at_dwell <- max(precision_at_dwell,
                            sensor_cfg$precision_opt_ppb)

  noise_ppb <- rnorm(n, 0, precision_at_dwell)

  # --- 2. Linear drift over elapsed time -----------------------------------
  elapsed_h <- dt$time_min / 60
  drift_ppb <- sensor_cfg$drift_ppb_hr * elapsed_h

  # --- 3. Calibration offset ------------------------------------------------
  offset_ppb <- sensor_cfg$cal_offset_ppb

  # --- 4. Total measured concentration --------------------------------------
  dt[, `:=`(
    ch4_measured_ppb = ch4_total_ppb + noise_ppb + drift_ppb + offset_ppb,
    noise_ppb        = round(noise_ppb, 4),
    drift_ppb        = round(drift_ppb, 4),
    offset_ppb       = offset_ppb,
    precision_ppb    = round(precision_at_dwell, 4),
    dwell_s          = dwell_s
  )]

  # --- 5. GPS position errors (optional) -------------------------------------
  if (add_position_uncertainty) {
    dt[, `:=`(
      x_measured = x + rnorm(.N, 0, gps_h_error_m),
      y_measured = y + rnorm(.N, 0, gps_h_error_m),
      z_measured = pmax(0.5, z + rnorm(.N, 0, gps_v_error_m)),
      gps_h_error_m = gps_h_error_m,
      gps_v_error_m = gps_v_error_m
    )]
  } else {
    dt[, `:=`(
      x_measured = x,
      y_measured = y,
      z_measured = z,
      gps_h_error_m = 0,
      gps_v_error_m = 0
    )]
  }

  dt[, ch4_measured_ppb := round(ch4_measured_ppb, 4)]
  dt[]
}


# -----------------------------------------------------------------------------
# Convenience: simulate measurements for a full scenario
# -----------------------------------------------------------------------------
#' Simulate drone measurements for a given truth dataset and sensor.
#'
#' @param truth_dt  data.table of truth concentrations.
#' @param sensor_name  Key into SENSOR_CONFIGS.
#' @param design_cfg  Sampling design list.
#' @param scenario_label  Label for output file.
#' @param save_output  Logical.
#' @param add_position_uncertainty  If FALSE (default), use true positions for
#'        inverse; if TRUE, add GPS errors. Passed to apply_sensor_model.
#' @return data.table of simulated measurements.
simulate_drone_measurements <- function(truth_dt,
                                        sensor_name,
                                        design_cfg,
                                        scenario_label = "default",
                                        save_output = TRUE,
                                        add_position_uncertainty = FALSE) {

  sensor_cfg <- SENSOR_CONFIGS[[sensor_name]]
  if (is.null(sensor_cfg)) {
    stop("Unknown sensor: ", sensor_name,
         ". Choose from: ", paste(names(SENSOR_CONFIGS), collapse = ", "))
  }

  message(sprintf("  Simulating sensor: %s (%s)",
                  sensor_cfg$name, sensor_cfg$short_name))

  meas_dt <- apply_sensor_model(
    truth_dt                 = truth_dt,
    sensor_cfg               = sensor_cfg,
    design_cfg               = design_cfg,
    add_position_uncertainty = add_position_uncertainty,
    seed                     = string_to_seed(scenario_label)
  )

  meas_dt[, sensor := sensor_cfg$short_name]

  if (save_output) {
    outfile <- proj_path("data", "raw",
                         sprintf("simulated_measurements_%s.csv", scenario_label))
    dir.create(dirname(outfile), recursive = TRUE, showWarnings = FALSE)
    fwrite(meas_dt, outfile)
    message(sprintf("  Measurements saved: %s (%d rows)",
                    basename(outfile), nrow(meas_dt)))
  }

  meas_dt
}


# -----------------------------------------------------------------------------
# Batch: loop over multiple sensor types for the same truth
# -----------------------------------------------------------------------------
#' Simulate measurements for multiple sensors against the same truth.
#'
#' @param truth_dt  data.table of truth concentrations.
#' @param sensor_names  Character vector of SENSOR_CONFIGS keys.
#' @param design_cfg  Sampling design list.
#' @param lake_name  Lake name for labelling.
#' @param met_name  Met scenario name for labelling.
#' @param save_output  Logical.
#' @return Named list of data.tables.
simulate_all_sensors <- function(truth_dt,
                                 sensor_names = names(SENSOR_CONFIGS),
                                 design_cfg,
                                 lake_name = "small",
                                 met_name = "mixed_day",
                                 save_output = TRUE) {

  pb <- make_progress(length(sensor_names), label = "Sensors")
  results <- list()

  for (sn in sensor_names) {
    label <- scenario_key(lake_name, sn, design_cfg$name, met_name)
    results[[sn]] <- simulate_drone_measurements(
      truth_dt       = truth_dt,
      sensor_name    = sn,
      design_cfg     = design_cfg,
      scenario_label = label,
      save_output    = save_output
    )
    pb$tick(msg = sn)
  }

  results
}


# -----------------------------------------------------------------------------
# Save sensor specs table
# -----------------------------------------------------------------------------
#' Write sensor specification summary to CSV.
#'
#' @param outdir  Output directory (default data/processed).
save_sensor_specs <- function(outdir = NULL) {
  if (is.null(outdir)) outdir <- proj_path("data", "processed")
  specs <- rbindlist(lapply(SENSOR_CONFIGS, as.data.table))
  fwrite(specs, file.path(outdir, "sensor_specs.csv"))
  message("  Sensor specs saved to sensor_specs.csv")
  invisible(specs)
}



