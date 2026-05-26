###############################################################################
# 02b_gaussian_forward.R - Gaussian plume forward model (analytical)
#
# Uses Pasquill-Gifford dispersion with Briggs (1973) rural sigma
# parameterizations and disk-integrated area source for circular lakes.
#
# This is the forward model used by the publication companion analyses.
#
# Inputs:
#   - Emission fields from 01_synthetic_emissions.R
#   - Gaussian met timeseries (Pasquill class A-F, wind speed/direction)
#   - Receptor configurations from define_receptor_config()
#
# Outputs:
#   - output/truth_concentrations_<scenario>.csv when save_output = TRUE
###############################################################################

if (!exists("project_root", mode = "function")) {
  proj_root <- if (file.exists("R/utils.R")) "." else ".."
  source(file.path(proj_root, "R", "utils.R"))
}


# -----------------------------------------------------------------------------
# 1. Pasquill-Gifford dispersion coefficients (Briggs 1973, rural)
# -----------------------------------------------------------------------------
#' Lateral dispersion coefficient sigma_y (m).
#' Briggs (1973) rural: sigma_y = a * x / sqrt(1 + b * x).
#' @param x_m  Downwind distance (m). Must be > 0.
#' @param pasquill_class  Character "A" through "F".
#' @return sigma_y in metres.
pg_sigma_y <- function(x_m, pasquill_class) {
  coefs <- PASQUILL_CLASSES[[pasquill_class]]
  if (is.null(coefs)) stop("Unknown Pasquill class: ", pasquill_class)
  coefs$sy_a * x_m / sqrt(1 + coefs$sy_b * x_m)
}

#' Vertical dispersion coefficient sigma_z (m).
#' Briggs (1973) rural. Functional form varies by class:
#'   A-B: a*x (unbounded), C-D: a*x*(1+b*x)^{-1/2}, E-F: a*x*(1+b*x)^{-1}.
#' @param x_m  Downwind distance (m). Must be > 0.
#' @param pasquill_class  Character "A" through "F".
#' @return sigma_z in metres.
pg_sigma_z <- function(x_m, pasquill_class) {
  coefs <- PASQUILL_CLASSES[[pasquill_class]]
  if (is.null(coefs)) stop("Unknown Pasquill class: ", pasquill_class)
  coefs$sz_a * x_m * (1 + coefs$sz_b * x_m)^coefs$sz_p
}


# -----------------------------------------------------------------------------
# 2. Gaussian plume point source (ground-level, with ground reflection)
# -----------------------------------------------------------------------------
#' Concentration from a ground-level point source (H = 0).
#'
#' C = (Q / (pi * sy * sz * U)) * exp(-y^2/(2*sy^2)) * exp(-z^2/(2*sz^2))
#'
#' Ground reflection doubles the ground-level contribution (factor of 2
#' absorbed into denominator: pi instead of 2*pi).
#'
#' @param x_down  Downwind distance (m). Vector.
#' @param y_cross  Crosswind distance (m). Same length as x_down.
#' @param z_rec  Receptor height (m). Same length as x_down.
#' @param q  Source emission rate (ug/s).
#' @param U  Wind speed (m/s).
#' @param pasquill_class  Pasquill stability class "A"-"F".
#' @return Concentration in ug/m3.
gaussian_plume_point <- function(x_down, y_cross, z_rec, q, U, pasquill_class) {
  valid <- x_down > 1e-3
  conc <- numeric(length(x_down))
  if (!any(valid)) return(conc)

  xv <- x_down[valid]
  yv <- y_cross[valid]
  zv <- z_rec[valid]

  sy <- pmax(pg_sigma_y(xv, pasquill_class), 0.01)
  sz <- pmax(pg_sigma_z(xv, pasquill_class), 0.01)

  conc[valid] <- (q / (pi * sy * sz * U)) *
    exp(-yv^2 / (2 * sy^2)) *
    exp(-zv^2 / (2 * sz^2))

  conc
}


# -----------------------------------------------------------------------------
# 3. Disk-integrated area source
# -----------------------------------------------------------------------------
#' Concentration at receptors from a circular area source via polar quadrature.
#'
#' Each area element on the lake disk is treated as a point source with
#' emission rate F * dA, then summed at each receptor.
#'
#' @param rec_x,rec_y,rec_z  Receptor coordinates (m). Vectors of equal length.
#' @param source_radius  Lake radius (m).
#' @param flux_ug_m2_s  Flux in ug m-2 s-1 (uniform over disk).
#' @param U_ref  Reference wind speed (m/s).
#' @param wd_deg  Wind-from direction (meteorological convention, degrees).
#' @param pasquill_class  Pasquill stability class "A"-"F".
#' @param n_rings  Radial quadrature rings (default 20).
#' @param n_src_az  Azimuthal quadrature points (default 36).
#' @param source_center  Source centre c(x, y), default c(0, 0).
#' @return Concentration in ug/m3 at each receptor.
gaussian_plume_disk <- function(rec_x, rec_y, rec_z,
                                source_radius, flux_ug_m2_s,
                                U_ref, wd_deg, pasquill_class,
                                n_rings = 20L, n_src_az = 36L,
                                source_center = c(0, 0)) {

  n_rec <- length(rec_x)
  conc_total <- numeric(n_rec)

  wd_to_rad <- ((wd_deg + 180) %% 360) * pi / 180
  wind_dx <- sin(wd_to_rad)
  wind_dy <- cos(wd_to_rad)

  r_edges <- seq(0, source_radius, length.out = n_rings + 1L)
  r_mids  <- (r_edges[-1L] + r_edges[-(n_rings + 1L)]) / 2
  theta_mids <- seq(0, 2 * pi * (1 - 1 / n_src_az), length.out = n_src_az)
  d_theta <- 2 * pi / n_src_az

  for (kr in seq_len(n_rings)) {
    r_lo <- r_edges[kr]
    r_hi <- r_edges[kr + 1L]
    r_m  <- r_mids[kr]
    dA_cell <- (r_hi^2 - r_lo^2) / 2 * d_theta
    q_cell  <- flux_ug_m2_s * dA_cell

    for (ka in seq_len(n_src_az)) {
      th <- theta_mids[ka]
      sx <- source_center[1] + r_m * cos(th)
      sy <- source_center[2] + r_m * sin(th)

      dx <- rec_x - sx
      dy <- rec_y - sy

      x_down  <- dx * wind_dx + dy * wind_dy
      y_cross <- dx * wind_dy - dy * wind_dx

      conc_total <- conc_total + gaussian_plume_point(
        x_down, y_cross, rec_z, q_cell, U_ref, pasquill_class
      )
    }
  }

  conc_total
}


# -----------------------------------------------------------------------------
# 4. Vectorised field computation (all receptors, one time step)
# -----------------------------------------------------------------------------
#' Gaussian plume concentration field at all receptor positions.
#'
#' @param receptors  data.table with x, y, z columns.
#' @param source_radius  Lake radius (m).
#' @param flux_ug_m2_s  Flux in ug m-2 s-1.
#' @param U_ref  Reference wind speed at z_ref (m/s).
#' @param wd_deg  Wind-from direction (degrees, met convention).
#' @param pasquill_class  Pasquill stability class "A"-"F".
#' @param n_rings  Quadrature rings (default 20).
#' @param n_src_az  Quadrature azimuth points (default 36).
#' @param source_center  Source centre c(x, y).
#' @return Numeric vector of concentrations (ug/m3), one per receptor.
gaussian_plume_field <- function(receptors, source_radius, flux_ug_m2_s,
                                 U_ref, wd_deg, pasquill_class,
                                 n_rings = 20L, n_src_az = 36L,
                                 source_center = c(0, 0)) {
  gaussian_plume_disk(
    rec_x          = as.numeric(receptors$x),
    rec_y          = as.numeric(receptors$y),
    rec_z          = as.numeric(receptors$z),
    source_radius  = source_radius,
    flux_ug_m2_s   = flux_ug_m2_s,
    U_ref          = U_ref,
    wd_deg         = wd_deg,
    pasquill_class = pasquill_class,
    n_rings        = n_rings,
    n_src_az       = n_src_az,
    source_center  = source_center
  )
}


#' Dense regular grid of Gaussian enhancement for isometric plots
#'
#' The OSSE forward model only stores concentrations at receptor locations.
#' This evaluates the same analytical plume (disk-integrated source) on a
#' horizontal \code{grid_n × grid_n} lattice at each height in \code{heights_m},
#' so visualizations can show a smooth field.  Quadrature defaults are slightly
#' coarser than the default forward model (\code{n_rings}/\code{n_src_az}) for
#' speed; the field is for visualization only.
#'
#' @param lake_diameter_m  Lake diameter (m).
#' @param flux_mg_m2_d  Area-mean flux (mg m-2 d-1).
#' @param ws_m_s  Reference wind speed (m/s), same as plume dilution.
#' @param wd_deg  Wind-from direction (deg).
#' @param pasquill_class  Pasquill class "A"–"F".
#' @param heights_m  Numeric vector of z levels (m).
#' @param grid_n  Grid points per axis on the domain (each axis).
#' @param pad_m  Extra extent beyond lake + outer ring (m).
#' @param n_rings  Disk quadrature rings (visualization; default 16).
#' @param n_src_az  Disk quadrature azimuth points (default 32).
#' @param nearby_lake_config  Optional list with \code{center_x}, \code{center_y},
#'   \code{radius_m}, \code{flux_mg_m2_d} (same as \code{resolve_nearby_lake_position()}
#'   output / isometric \code{nearby_lake_config}). When non-NULL, adds the second
#'   disk-integrated source and expands the horizontal domain to include both plumes.
#' @return \code{data.table} with columns \code{x}, \code{y}, \code{z},
#'   \code{ch4_enhancement_ppb}.
dense_gaussian_enhancement_for_isometric <- function(lake_diameter_m,
                                                     flux_mg_m2_d,
                                                     ws_m_s,
                                                     wd_deg,
                                                     pasquill_class,
                                                     heights_m,
                                                     grid_n    = 40L,
                                                     pad_m     = 30,
                                                     n_rings   = 16L,
                                                     n_src_az  = 32L,
                                                     nearby_lake_config = NULL) {
  lake_r    <- lake_diameter_m / 2
  flux_ug_s <- flux_mg_d_to_ug_s(flux_mg_m2_d)
  half_main <- lake_r + 10 + pad_m

  if (is.null(nearby_lake_config)) {
    r_extent <- half_main
    xg <- seq(-r_extent, r_extent, length.out = grid_n)
    yg <- seq(-r_extent, r_extent, length.out = grid_n)
  } else {
    cx <- nearby_lake_config$center_x
    cy <- nearby_lake_config$center_y
    nr <- nearby_lake_config$radius_m
    xmin <- min(-half_main, cx - nr) - pad_m
    xmax <- max(half_main, cx + nr) + pad_m
    ymin <- min(-half_main, cy - nr) - pad_m
    ymax <- max(half_main, cy + nr) + pad_m
    xg <- seq(xmin, xmax, length.out = grid_n)
    yg <- seq(ymin, ymax, length.out = grid_n)
  }

  rows <- vector("list", length(heights_m))
  for (ih in seq_along(heights_m)) {
    zz <- heights_m[ih]
    gr <- CJ(x = xg, y = yg)
    gr[, z := zz]
    gr[, receptor_id := seq_len(.N)]
    conc <- gaussian_plume_field(
      gr, lake_r, flux_ug_s, ws_m_s, wd_deg, pasquill_class,
      n_rings = n_rings, n_src_az = n_src_az
    )
    if (!is.null(nearby_lake_config)) {
      nearby_flux_ug_s <- flux_mg_d_to_ug_s(nearby_lake_config$flux_mg_m2_d)
      conc <- conc + gaussian_plume_field(
        gr, nearby_lake_config$radius_m, nearby_flux_ug_s,
        ws_m_s, wd_deg, pasquill_class,
        n_rings = n_rings, n_src_az = n_src_az,
        source_center = c(nearby_lake_config$center_x, nearby_lake_config$center_y)
      )
    }
    gr[, ch4_enhancement_ppb := ug_m3_to_ppb(conc)]
    rows[[ih]] <- gr[, .(x, y, z, ch4_enhancement_ppb)]
  }
  rbindlist(rows)
}


#' Wind speed at height z using power-law profile.
#'
#' U(z) = U_ref * (z / z_ref)^p, where p depends on Pasquill class and
#' terrain (Hanna et al. 1982, rural).
#'
#' @param z_m  Heights (m). Vector.
#' @param U_ref  Reference wind speed (m/s).
#' @param z_ref  Reference height (m), default 10.
#' @param pasquill_class  Pasquill stability class "A"-"F".
#' @return Wind speed at each height (m/s). Clamped >= 0.1.
wind_at_height_powerlaw <- function(z_m, U_ref, z_ref = 10,
                                    pasquill_class = "D") {
  p <- PASQUILL_CLASSES[[pasquill_class]]$wind_exp
  if (is.null(p)) p <- 0.15
  z_m <- pmax(z_m, 0.01)
  pmax(U_ref * (z_m / z_ref)^p, 0.1)
}


# -----------------------------------------------------------------------------
# 5. Pipeline-compatible forward model wrapper
# -----------------------------------------------------------------------------
#' Run Gaussian plume forward model producing the manuscript truth table.
#'
#' Drop-in replacement for \code{run_forward_model()}.  Uses analytical
#' Gaussian plume with Pasquill-Gifford dispersion and disk-integrated area
#' source.  MOST columns (\code{ustar_m_s}, \code{L_m}) are set to NA.
#' \code{ws_at_z_m_s} is the same constant reference wind as \code{ws_m_s} at
#' every receptor height, matching the transport wind used in the plume
#' dilution term so the Conley inverse samples forward-model wind truth.
#'
#' @param emission_result  Output from generate_single_lake_emissions().
#' @param met_ts  data.table from \code{generate_gaussian_met_timeseries()}.
#'   Must contain \code{pasquill_class}.
#' @param receptors  data.table from define_receptor_config().
#' @param scenario_label  Character label for output files.
#' @param save_output  Logical; write CSV outputs?
#' @param bg_seed  Random seed for background CH4 generation.
#' @param bg_diel_amp  Background diel amplitude in ppb (default 0).
#' @param bg_synoptic_amp  Background synoptic amplitude in ppb (default 35).
#' @param flux_spatial_cv  CV for per-receptor flux noise (0 = homogeneous).
#' @param z_ref_m  Reference height for \code{U_ref} in met (default 10 m);
#'   not used to vary \code{ws_at_z_m_s}, which is constant for Gaussian plume.
#' @param n_rings  Disk quadrature rings (default 20).
#' @param n_src_az  Disk quadrature azimuth points (default 36).
#' @param nearby_lake_config  Optional nearby-lake contamination config.
#' @return data.table with receptor-level methane enhancement and meteorology.
run_forward_model_gaussian <- function(emission_result,
                                       met_ts,
                                       receptors,
                                       scenario_label  = "default",
                                       save_output     = TRUE,
                                       bg_seed         = NULL,
                                       bg_diel_amp     = 0,
                                       bg_synoptic_amp = 35,
                                       flux_spatial_cv = 0,
                                       z_ref_m         = 10,
                                       n_rings         = 20L,
                                       n_src_az        = 36L,
                                       nearby_lake_config = NULL) {

  n_steps     <- nrow(met_ts)
  lake_radius <- emission_result$diameter_m / 2

  if (is.null(bg_seed)) bg_seed <- string_to_seed(scenario_label)

  lake_mean_dt <- emission_result$lake_mean

  bg_ch4 <- generate_background_ch4(
    n_steps,
    dt_minutes   = if (n_steps > 1) met_ts$time_min[2] - met_ts$time_min[1] else 30,
    diel_amp     = bg_diel_amp,
    synoptic_amp = bg_synoptic_amp,
    start_hour   = met_ts$hour_of_day[1],
    seed         = bg_seed
  )

  if (flux_spatial_cv > 0) {
    rng_state <- if (exists(".Random.seed", envir = .GlobalEnv))
      get(".Random.seed", envir = .GlobalEnv) else NULL
    set.seed(bg_seed + 7777L)
    receptor_flux_multiplier <- pmax(1 + rnorm(nrow(receptors), 0, flux_spatial_cv), 0.01)
    if (!is.null(rng_state)) assign(".Random.seed", rng_state, envir = .GlobalEnv)
  } else {
    receptor_flux_multiplier <- rep(1, nrow(receptors))
  }

  results <- vector("list", n_steps)
  pb <- make_progress(n_steps, label = "Gaussian forward")

  for (t_idx in seq_len(n_steps)) {
    hour_now <- met_ts$hour_of_day[t_idx]
    flux_idx <- which.min(abs(lake_mean_dt$hour - hour_now))
    flux_true <- lake_mean_dt$lake_mean_flux[flux_idx]

    ws_now <- met_ts$ws_m_s[t_idx]
    wd_now <- met_ts$wd_deg[t_idx]
    pc_now <- met_ts$pasquill_class[t_idx]

    flux_ug_s <- flux_mg_d_to_ug_s(flux_true)

    enhancement_ug_m3 <- gaussian_plume_field(
      receptors, lake_radius, flux_ug_s, ws_now, wd_now, pc_now,
      n_rings = n_rings, n_src_az = n_src_az
    ) * receptor_flux_multiplier

    ce_vals <- if (flux_ug_s > 0) {
      enhancement_ug_m3 / (flux_ug_s * receptor_flux_multiplier)
    } else {
      rep(0, nrow(receptors))
    }

    contamination_ug_m3   <- rep(0, nrow(receptors))
    ce_nearby             <- rep(0, nrow(receptors))
    contamination_ppb_vec <- rep(0, nrow(receptors))
    if (!is.null(nearby_lake_config)) {
      nearby_flux_ug_s <- flux_mg_d_to_ug_s(nearby_lake_config$flux_mg_m2_d)
      contamination_ug_m3 <- gaussian_plume_field(
        receptors, nearby_lake_config$radius_m, nearby_flux_ug_s,
        ws_now, wd_now, pc_now,
        n_rings = n_rings, n_src_az = n_src_az,
        source_center = c(nearby_lake_config$center_x, nearby_lake_config$center_y)
      )
      ce_nearby <- if (nearby_flux_ug_s > 0) {
        contamination_ug_m3 / nearby_flux_ug_s
      } else {
        rep(0, nrow(receptors))
      }
      contamination_ppb_vec <- round(ug_m3_to_ppb(contamination_ug_m3), 4)
    }

    total_enhancement_ug_m3 <- enhancement_ug_m3 + contamination_ug_m3
    enhancement_ppb <- ug_m3_to_ppb(total_enhancement_ug_m3)

    # Constant wind at all heights: same U as in gaussian_plume_* dilution (ws_now).
    # Conley inverse uses ws_at_z_m_s as observed wind; it must match forward truth.
    ws_at_z <- rep(ws_now, nrow(receptors))

    ustar_now <- if ("ustar_m_s" %in% names(met_ts)) met_ts$ustar_m_s[t_idx] else NA_real_
    L_now     <- if ("L_m" %in% names(met_ts)) met_ts$L_m[t_idx] else NA_real_
    z0_now    <- if ("z0_m" %in% names(met_ts)) met_ts$z0_m[t_idx] else 0.001

    results[[t_idx]] <- data.table(
      time_min              = met_ts$time_min[t_idx],
      hour_of_day           = hour_now,
      receptor_id           = receptors$receptor_id,
      x                     = receptors$x,
      y                     = receptors$y,
      z                     = receptors$z,
      ch4_bg_ppb            = bg_ch4[t_idx],
      ch4_enhancement_ppb   = round(enhancement_ppb, 4),
      ch4_contamination_ppb = contamination_ppb_vec,
      ch4_total_ppb         = round(bg_ch4[t_idx] + enhancement_ppb, 4),
      flux_true_mg_m2_d     = round(flux_true, 2),
      CE_ratio              = round(ce_vals, 6),
      CE_ratio_nearby       = round(ce_nearby, 6),
      ustar_m_s             = ustar_now,
      L_m                   = L_now,
      ws_m_s                = ws_now,
      ws_at_z_m_s           = round(ws_at_z, 4),
      wd_deg                = wd_now,
      z0_m                  = z0_now,
      stability_class       = met_ts$stability_class[t_idx],
      pasquill_class        = as.character(met_ts$pasquill_class[t_idx])
    )

    pb$tick()
  }

  truth_dt <- rbindlist(results)

  if (save_output) {
    outfile <- proj_path("data", "raw",
                         sprintf("truth_concentrations_%s.csv", scenario_label))
    fwrite(truth_dt, outfile)
    message(sprintf("  Truth concentrations saved: %s (%d rows)",
                    basename(outfile), nrow(truth_dt)))
  }

  truth_dt
}


#' Convenience wrapper: Gaussian forward model for a standard scenario.
#'
#' Analogous to \code{run_forward_scenario()} but uses the Gaussian plume.
#'
#' @param lake_name  "small", "medium", or "large".
#' @param design_cfg  Sampling design list (from SAMPLING_CONFIGS).
#' @param met_cfg  Gaussian met config (from GAUSSIAN_MET_SCENARIOS).
#' @param sensor_name  Sensor name for labelling (default "truth").
#' @param emission_result  Pre-computed emission result (optional).
#' @param save_output  Logical.
#' @param met_seed  Seed for met generation (default 123).
#' @return List with $truth_dt and $receptors.
run_forward_scenario_gaussian <- function(lake_name,
                                          design_cfg,
                                          met_cfg,
                                          sensor_name     = "truth",
                                          emission_result = NULL,
                                          save_output     = TRUE,
                                          met_seed        = 123) {
  lk <- LAKE_CONFIGS[[lake_name]]

  if (is.null(emission_result)) {
    emission_result <- generate_single_lake_emissions(lake_name,
                                                      save_output = save_output)
  }

  nearby_lake_config <- if (!is.null(lk$nearby_lake)) {
    wd_base <- if (!is.null(met_cfg$wd_base_deg)) met_cfg$wd_base_deg else 225
    resolve_nearby_lake_position(
      nearby_lake_cfg = lk$nearby_lake,
      main_radius_m   = lk$diameter_m / 2,
      ring_radius_m   = lk$diameter_m / 2 + 10,
      wd_base_deg     = wd_base,
      main_base_flux  = emission_result$base_flux
    )
  } else {
    NULL
  }

  met_ts <- generate_gaussian_met_timeseries(
    duration_days   = met_cfg$duration_days,
    dt_minutes      = met_cfg$dt_minutes,
    start_hour      = if (!is.null(met_cfg$start_hour)) met_cfg$start_hour else 12,
    U_ref_m_s       = met_cfg$U_ref_m_s,
    wd_base_deg     = if (!is.null(met_cfg$wd_base_deg)) met_cfg$wd_base_deg else 225,
    wd_noise_sd_deg = if (!is.null(met_cfg$wd_noise_sd_deg)) met_cfg$wd_noise_sd_deg else 0,
    pasquill_class  = met_cfg$pasquill_class,
    seed            = met_seed
  )

  receptors <- define_receptor_config(
    lake_diameter        = lk$diameter_m,
    n_azimuth            = design_cfg$n_azimuth,
    heights_m            = design_cfg$heights_m,
    distance_from_edge_m = 10,
    inner_ring_fracs     = design_cfg$inner_ring_fracs,
    background_heights_m = design_cfg$background_heights_m
  )

  label   <- scenario_key(lake_name, sensor_name, design_cfg$name, met_cfg$name)
  bg_seed <- string_to_seed(
    paste0(lake_name, "__", design_cfg$name, "__", met_cfg$name, "__metSeed", met_seed)
  )
  z_ref_m <- if (!is.null(met_cfg$z_ref_m)) met_cfg$z_ref_m else 10

  truth_dt <- run_forward_model_gaussian(
    emission_result    = emission_result,
    met_ts             = met_ts,
    receptors          = receptors,
    scenario_label     = label,
    save_output        = save_output,
    bg_seed            = bg_seed,
    z_ref_m            = z_ref_m,
    nearby_lake_config = nearby_lake_config
  )

  list(truth_dt = truth_dt, receptors = receptors)
}


# -----------------------------------------------------------------------------
# 6. Grid evaluation for visualisation
# -----------------------------------------------------------------------------
#' Evaluate Gaussian plume on a fine grid for visualisation.
#'
#' Returns a data.table suitable for ggplot2 \code{geom_tile} / contour plots.
#'
#' @param source_radius  Lake radius (m).
#' @param flux_mg_m2_d  Flux in mg CH4 m-2 d-1.
#' @param U_ref  Reference wind speed at z_ref (m/s).
#' @param wd_deg  Wind-from direction (degrees).
#' @param pasquill_class  Pasquill class "A"-"F".
#' @param plane  "xz" (downwind vs height, y=0) or "xy" (downwind vs crosswind, fixed z).
#' @param z_slice_m  Height for "xy" plane (m, default 1).
#' @param x_range  Downwind range (m), default c(1, 500).
#' @param y_range  Crosswind range for "xy" (m), default c(-200, 200).
#' @param z_range  Height range for "xz" (m), default c(0.1, 30).
#' @param nx  Grid points in x (default 200).
#' @param ny  Grid points in y or z (default 100).
#' @param n_rings  Quadrature rings.
#' @param n_src_az  Quadrature azimuth points.
#' @return data.table with x_down, z or y_cross, conc_ppb, conc_ug_m3.
gaussian_plume_grid <- function(source_radius, flux_mg_m2_d,
                                U_ref, wd_deg, pasquill_class,
                                plane     = "xz",
                                z_slice_m = 1,
                                x_range   = c(1, 500),
                                y_range   = c(-200, 200),
                                z_range   = c(0.1, 30),
                                nx = 200L, ny = 100L,
                                n_rings = 20L, n_src_az = 36L) {

  flux_ug_s <- flux_mg_d_to_ug_s(flux_mg_m2_d)

  wd_to_rad <- ((wd_deg + 180) %% 360) * pi / 180
  wind_dx <- sin(wd_to_rad)
  wind_dy <- cos(wd_to_rad)

  x_seq <- seq(x_range[1], x_range[2], length.out = nx)

  if (plane == "xz") {
    z_seq <- seq(z_range[1], z_range[2], length.out = ny)
    grid  <- CJ(x_down = x_seq, z = z_seq)
    grid[, `:=`(abs_x = x_down * wind_dx, abs_y = x_down * wind_dy)]

    conc <- gaussian_plume_disk(
      grid$abs_x, grid$abs_y, grid$z,
      source_radius, flux_ug_s, U_ref, wd_deg, pasquill_class,
      n_rings, n_src_az
    )
    grid[, `:=`(conc_ug_m3 = conc, conc_ppb = ug_m3_to_ppb(conc))]
    grid[, c("abs_x", "abs_y") := NULL]

  } else if (plane == "xy") {
    y_seq <- seq(y_range[1], y_range[2], length.out = ny)
    grid  <- CJ(x_down = x_seq, y_cross = y_seq)

    cross_dx <- wind_dy
    cross_dy <- -wind_dx
    grid[, `:=`(
      abs_x = x_down * wind_dx + y_cross * cross_dx,
      abs_y = x_down * wind_dy + y_cross * cross_dy
    )]

    conc <- gaussian_plume_disk(
      grid$abs_x, grid$abs_y, rep(z_slice_m, nrow(grid)),
      source_radius, flux_ug_s, U_ref, wd_deg, pasquill_class,
      n_rings, n_src_az
    )
    grid[, `:=`(conc_ug_m3 = conc, conc_ppb = ug_m3_to_ppb(conc))]
    grid[, c("abs_x", "abs_y") := NULL]

  } else {
    stop("plane must be 'xz' or 'xy'")
  }

  grid[]
}


# -----------------------------------------------------------------------------
# Standalone execution
# -----------------------------------------------------------------------------
if (sys.nframe() == 0) {
  proj_root <- project_root()
  source(file.path(proj_root, "R", "01_synthetic_emissions.R"))
  message("=== 02b: Running Gaussian plume forward model ===")
  fwd <- run_forward_scenario_gaussian(
    "small",
    design_cfg = SAMPLING_CONFIGS$design_24az_5h_30min,
    met_cfg    = GAUSSIAN_MET_SCENARIOS$neutral_moderate
  )
  message("Done.  Rows: ", nrow(fwd$truth_dt))
}

