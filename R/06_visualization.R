###############################################################################
# 06_visualization.R - Plotting code for paper-style figures
#
# Scientific notes:
#   - Ebullition is treated as part of an effective, lake-averaged surface
#     flux.  We do NOT resolve hotspots or fine-scale spatial structure.
#   - The OSSE is a publication companion for Gaussian-forward /
#     Conley-inverse manuscript experiments.
#
# Generates:
#   Capability-oriented figures:
#     - Flux RMSE vs sensor type and sampling configuration
#     - Receptor design analysis (n_azimuth, vertical levels, cycle time)
#     - Flux error vs meteorological stability class
#     - Uncertainty budget bar chart
#     - Sensitivity tornado (LHS first-order Si) and Sobol Si/STi bar chart
#   Diagnostic figures:
#     - Effective flux field maps + diel time series
#     - True vs retrieved flux time series with uncertainty bands
#   CSV summary tables:
#     - capability_metrics.csv, scenarios_summary.csv
###############################################################################

if (!exists("project_root", mode = "function")) {
  proj_root <- if (file.exists("R/utils.R")) "." else ".."
  source(file.path(proj_root, "R", "utils.R"))
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(viridis)
})

cowplot_available <- requireNamespace("cowplot", quietly = TRUE)
if (cowplot_available) library(cowplot)

# Figure directory: use option if set so callers can write figures beside their
# analysis outputs.
get_figure_dir <- function() {
  getOption("uas_gauss_osse_figures", proj_path("figures"))
}
DATA_PROC  <- proj_path("data", "processed")

# Consistent theme for publication figures
theme_osse <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey92", colour = NA),
      legend.position  = "bottom",
      plot.title       = element_text(face = "bold", size = base_size + 1)
    )
}


# -----------------------------------------------------------------------------
# 1. Flux RMSE vs sensor type and configuration
# -----------------------------------------------------------------------------
#' Plot flux RMSE by sensor type for each lake size.
#'
#' @param metrics_dt  data.table with columns: scenario, lake, sensor,
#'                    design, rmse_pct, meets_target.
#' @param target_pct  Uncertainty target line (default 25).
#' @param save  Logical.
plot_flux_rmse_vs_sensor <- function(metrics_dt, target_pct = 25, save = TRUE) {
  if (nrow(metrics_dt) == 0) { message("  No data for sensor plot"); return(invisible()) }

  p <- ggplot(metrics_dt, aes(x = sensor, y = rmse_pct, fill = sensor)) +
    geom_col(width = 0.7, alpha = 0.85) +
    geom_hline(yintercept = target_pct, linetype = "dashed", colour = "red",
               linewidth = 0.7) +
    facet_grid(design ~ lake, scales = "free_x") +
    scale_fill_viridis_d(option = "D", end = 0.85) +
    labs(
      title    = "Flux RMSE by sensor and sampling configuration",
      subtitle = sprintf("Dashed line = %d%% uncertainty target", target_pct),
      x = "Sensor", y = "Flux RMSE (%)",
      fill = "Sensor"
    ) +
    theme_osse(base_size = base_size) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))

  if (save) {
    ggsave_osse(file.path(get_figure_dir(), "flux_rmse_vs_sensor.png"), p,
           width = 10, height = 7, dpi = 200)
    message("  Saved: flux_rmse_vs_sensor.png")
  }
  invisible(p)
}


# -----------------------------------------------------------------------------
# 2. Receptor design: RMSE vs receptors, levels, and cycle time
# -----------------------------------------------------------------------------
#' Plot flux RMSE vs receptor count, vertical levels, and cycle time.
#'
#' @param metrics_dt  data.table with design details parsed.
#' @param save  Logical.
plot_receptor_design <- function(metrics_dt, save = TRUE) {
  if (nrow(metrics_dt) == 0) { message("  No data for receptor plot"); return(invisible()) }

  dt <- copy(metrics_dt)

  # Parse design components from name: e.g. "8az_3h_30min"
  dt[, n_azimuth  := as.integer(sub("^(\\d+)az.*", "\\1", design))]
  dt[, n_levels   := as.integer(sub(".*?(\\d+)h_.*", "\\1", design))]
  dt[, cycle_min  := as.integer(sub(".*_(\\d+)min$", "\\1", design))]
  dt[is.na(n_levels), n_levels := 1L]
  dt[is.na(cycle_min), cycle_min := 30L]

  # Create a combined label for colour mapping
  dt[, design_label := sprintf("%d levels, %d min", n_levels, cycle_min)]

  p <- ggplot(dt, aes(x = factor(n_azimuth), y = rmse_pct,
                      colour = design_label, shape = design_label,
                      group = design_label)) +
    geom_point(size = 3, position = position_dodge(width = 0.3)) +
    geom_line(linewidth = 0.8, position = position_dodge(width = 0.3)) +
    facet_wrap(~ lake, nrow = 1) +
    geom_hline(yintercept = 25, linetype = "dashed", colour = "red") +
    scale_colour_viridis_d(option = "C", end = 0.85) +
    labs(
      title  = "Receptor design: RMSE vs azimuthal count, levels, and cycle time",
      x      = "Number of azimuthal receptors",
      y      = "Flux RMSE (%)",
      colour = "Design",
      shape  = "Design"
    ) +
    theme_osse()

  if (save) {
    ggsave_osse(file.path(get_figure_dir(), "receptor_design.png"), p,
           width = 10, height = 5, dpi = 200)
    message("  Saved: receptor_design.png")
  }
  invisible(p)
}


# -----------------------------------------------------------------------------
# 3. Flux error vs meteorological stability
# -----------------------------------------------------------------------------
#' Plot flux error by stability class.
#'
#' @param flux_results  List of flux data.tables (from retrieve_flux),
#'                      or a single combined data.table.
#' @param save  Logical.
plot_flux_vs_stability <- function(flux_results, save = TRUE) {
  if (is.data.table(flux_results)) {
    dt <- flux_results
  } else {
    dt <- rbindlist(flux_results, fill = TRUE)
  }

  dt_valid <- dt[qc_flag == "ok" & is.finite(flux_retrieved)]
  if (nrow(dt_valid) == 0) { message("  No valid data for stability plot"); return(invisible()) }

  dt_valid[, rel_error_pct := 100 * (flux_retrieved - flux_true) / flux_true]

  p <- ggplot(dt_valid, aes(x = stability, y = rel_error_pct, fill = stability)) +
    geom_boxplot(alpha = 0.8, outlier.size = 1) +
    geom_hline(yintercept = 0, linetype = "solid", colour = "grey40") +
    geom_hline(yintercept = c(-25, 25), linetype = "dashed", colour = "red") +
    scale_fill_viridis_d(option = "B", begin = 0.2, end = 0.8) +
    labs(
      title = "Flux retrieval error by meteorological stability class",
      x     = "Stability class",
      y     = "Relative flux error (%)",
      fill  = "Stability"
    ) +
    theme_osse()

  if (save) {
    ggsave_osse(file.path(get_figure_dir(), "flux_vs_stability.png"), p,
           width = 7, height = 5, dpi = 200)
    message("  Saved: flux_vs_stability.png")
  }
  invisible(p)
}


# -----------------------------------------------------------------------------
# 4. Uncertainty budget bar chart
# -----------------------------------------------------------------------------
#' Plot variance-decomposition uncertainty budget.
#'
#' @param budget_dt  data.table from compute_variance_budget().
#' @param lake_name  Lake name for reading from file if budget_dt is NULL.
#' @param processed_dir  Optional: directory (relative to project root) for CSV
#'   lookup; e.g. \code{"data/processed/run_20240301"} to use a previous run's data.
#' @param save  Logical.
plot_uncertainty_budget <- function(budget_dt = NULL, lake_name = "small",
                                    processed_dir = NULL, save = TRUE) {
  if (is.null(budget_dt)) {
    data_dir <- if (is.null(processed_dir)) DATA_PROC else file.path(project_root(), processed_dir)
    fpath <- file.path(data_dir, sprintf("variance_budget_%s.csv", lake_name))
    if (!file.exists(fpath)) {
      fpath <- file.path(data_dir, "variance_budget.csv")
    }
    if (!file.exists(fpath)) { message("  No variance budget data"); return(invisible()) }
    budget_dt <- fread(fpath)
  }

  if (nrow(budget_dt) == 0) return(invisible())

  budget_dt[, category := factor(category, levels = unique(category))]

  p <- ggplot(budget_dt, aes(x = category, y = variance_share_norm, fill = category)) +
    geom_col(width = 0.7, alpha = 0.85) +
    scale_fill_viridis_d(option = "E", end = 0.85) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(
      title = "Uncertainty budget: variance contributions by source",
      x     = NULL,
      y     = "Fraction of total RMSE variance",
      fill  = NULL
    ) +
    theme_osse() +
    theme(axis.text.x = element_text(angle = 25, hjust = 1),
          legend.position = "none")

  if (save) {
    ggsave_osse(file.path(get_figure_dir(), "uncertainty_budget.png"), p,
           width = 7, height = 5, dpi = 200)
    message("  Saved: uncertainty_budget.png")
  }
  invisible(p)
}


# -----------------------------------------------------------------------------
# 5. Effective flux field maps + diel time series
# -----------------------------------------------------------------------------
#' Plot spatial map and diel cycle for one lake's emission field.
#'
#' @param flux_field  data.table from generate_effective_flux_field.
#' @param lake_mean   data.table with hour and lake_mean_flux.
#' @param lake_name   Character label.
#' @param save  Logical.
plot_emission_field <- function(flux_field, lake_mean, lake_name = "small",
                                save = TRUE) {

  # Spatial map at hour = 12 (or any hour; flux is constant with no diel)
  snap <- flux_field[hour == 12]

  p_map <- ggplot(snap, aes(x = x, y = y, fill = flux_mg_m2_d)) +
    geom_tile() +
    scale_fill_viridis_c(option = "A",
                         name = expression(paste("Flux (mg CH"[4], " m"^{-2}, " d"^{-1}, ")"))) +
    coord_equal() +
    labs(title = sprintf("%s lake - effective flux", lake_name)) +
    theme_osse()

  # Lake-mean flux (constant when no diel)
  p_diel <- ggplot(lake_mean, aes(x = hour, y = lake_mean_flux)) +
    geom_line(linewidth = 1, colour = viridis::viridis(3)[2]) +
    geom_point(size = 2, colour = viridis::viridis(3)[2]) +
    labs(
      title = sprintf("%s lake - lake-mean flux", lake_name),
      x     = "Hour of day",
      y     = expression(paste("Lake-mean flux (mg CH"[4], " m"^{-2}, " d"^{-1}, ")"))
    ) +
    theme_osse()

  if (save) {
    ggsave_osse(file.path(get_figure_dir(), sprintf("emission_map_%s.png", lake_name)),
           p_map, width = 7, height = 6, dpi = 200)
    ggsave_osse(file.path(get_figure_dir(), sprintf("emission_diel_%s.png", lake_name)),
           p_diel, width = 6, height = 4, dpi = 200)
    message(sprintf("  Saved: emission_map_%s.png, emission_diel_%s.png",
                    lake_name, lake_name))
  }

  if (cowplot_available) {
    combined <- cowplot::plot_grid(p_map, p_diel, ncol = 2, rel_widths = c(1.3, 1))
    if (save) {
      ggsave_osse(file.path(get_figure_dir(), sprintf("emission_combined_%s.png", lake_name)),
             combined, width = 13, height = 5, dpi = 200)
    }
    return(invisible(combined))
  }

  invisible(list(map = p_map, diel = p_diel))
}


# -----------------------------------------------------------------------------
# 5b. Measurement locations (one height layer)
# -----------------------------------------------------------------------------
#' Plot measurement locations for one height layer (x-y map).
#'
#' Shows the lake circle and receptor positions at a single height. The design
#' may have multiple heights; this plot shows only one layer for clarity.
#'
#' @param receptors  data.table from define_receptor_config (x, y, z, ring_type).
#' @param lake_diameter  Lake diameter in metres.
#' @param lake_name  Character label.
#' @param show_height_m  Height (m) to display. If NULL, use the first height in receptors.
#' @param save  Logical.
plot_measurement_locations <- function(receptors, lake_diameter, lake_name = "small",
                                       show_height_m = NULL, save = TRUE) {
  if (is.null(show_height_m)) show_height_m <- if (8 %in% receptors$z) 8 else receptors[, z[1]]
  one_layer <- receptors[z == show_height_m]
  if (nrow(one_layer) == 0) one_layer <- receptors[z == receptors[, min(z)]]

  lake_radius <- lake_diameter / 2
  theta <- seq(0, 2 * pi, length.out = 100)
  lake_df <- data.frame(
    x = lake_radius * cos(theta),
    y = lake_radius * sin(theta)
  )

  p <- ggplot() +
    geom_polygon(data = lake_df, aes(x = x, y = y), fill = "grey88", colour = "grey40") +
    geom_point(data = one_layer, aes(x = x, y = y, colour = ring_type), size = 2.5, alpha = 0.9) +
    scale_colour_viridis_d(option = "C", end = 0.9, name = "Ring") +
    coord_equal() +
    labs(
      title = sprintf("%s lake - measurement locations (%g m height)", lake_name, show_height_m),
      x     = "x (m)", y = "y (m)"
    ) +
    theme_osse()

  if (save) {
    ggsave_osse(file.path(get_figure_dir(), sprintf("measurement_locations_%s.png", lake_name)),
           p, width = 7, height = 6, dpi = 200)
    message(sprintf("  Saved: measurement_locations_%s.png", lake_name))
  }
  invisible(p)
}


# -----------------------------------------------------------------------------
# 5c. Forward model plume visualization (x-z cross section and z profile)
# -----------------------------------------------------------------------------
#' Plot CH4 enhancement in x-z cross section (slice |y| < y_slice_m).
#'
#' Shows receptor positions coloured by enhancement. Useful to see plume
#' vertical extent and whether the background receptor (150 m) is above the plume.
#'
#' @param truth_dt  data.table from run_forward_model (receptor_id, x, y, z,
#'        ch4_enhancement_ppb, time_min).
#' @param receptors  data.table from define_receptor_config (receptor_id, ring_type).
#' @param lake_diameter  Lake diameter in metres.
#' @param time_min  Time step to plot. If NULL, use first.
#' @param y_slice_m  Include receptors with |y| < this (default 15).
#' @param scenario_label  Label for title and filename.
#' @param save  Logical.
plot_plume_xz <- function(truth_dt, receptors, lake_diameter,
                          time_min = NULL, y_slice_m = 15,
                          scenario_label = "default", save = TRUE) {
  dt <- merge(truth_dt, receptors[, .(receptor_id, ring_type)], by = "receptor_id")
  t_use <- if (is.null(time_min)) dt[, min(time_min)] else time_min
  dt <- dt[time_min == t_use]
  dt <- dt[abs(y) < y_slice_m]
  if (nrow(dt) == 0) { message("  No receptors in x-z slice"); return(invisible()) }

  lake_radius <- lake_diameter / 2

  p <- ggplot(dt, aes(x = x, y = z, colour = ch4_enhancement_ppb, shape = ring_type)) +
    geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.5) +
    geom_segment(aes(x = -lake_radius, xend = lake_radius, y = 0, yend = 0),
                 data = data.frame(lake_radius = lake_radius),
                 colour = "grey40", linewidth = 3, inherit.aes = FALSE) +
    geom_point(size = 3, alpha = 0.9) +
    scale_colour_viridis_c(option = "plasma", name = "Enhancement (ppb)") +
    scale_shape_manual(values = c(outer = 16, inner_0.5 = 17, inner_0.25 = 15, background = 8),
                       name = "Ring", drop = TRUE) +
    coord_cartesian(ylim = c(0, NA)) +
    labs(
      title = paste("Plume x-z cross section:", scenario_label),
      subtitle = sprintf("Time %g min, slice |y| < %g m", t_use, y_slice_m),
      x = "x (m)", y = "Height z (m)"
    ) +
    theme_osse()

  if (save) {
    fname <- sprintf("plume_xz_%s.png", gsub("[^a-zA-Z0-9_]", "_", scenario_label))
    ggsave_osse(file.path(get_figure_dir(), fname), p, width = 9, height = 6, dpi = 200)
    message(sprintf("  Saved: %s", fname))
  }
  invisible(p)
}


#' Plot vertical profile of CH4 enhancement (mean by height).
#'
#' @param truth_dt  data.table from run_forward_model.
#' @param receptors  data.table from define_receptor_config.
#' @param time_min  Time step. If NULL, use first.
#' @param scenario_label  Label for title and filename.
#' @param save  Logical.
plot_plume_z_profile <- function(truth_dt, receptors,
                                 time_min = NULL, scenario_label = "default", save = TRUE) {
  dt <- merge(truth_dt, receptors[, .(receptor_id, ring_type)], by = "receptor_id")
  t_use <- if (is.null(time_min)) dt[, min(time_min)] else time_min
  dt <- dt[time_min == t_use]

  prof <- dt[, .(
    mean_enh_ppb = mean(ch4_enhancement_ppb, na.rm = TRUE),
    sd_enh_ppb   = sd(ch4_enhancement_ppb, na.rm = TRUE)
  ), by = .(z, ring_type)]
  if (nrow(prof) == 0) { message("  No data for z profile"); return(invisible()) }

  p <- ggplot(prof, aes(x = mean_enh_ppb, y = z, colour = ring_type)) +
    geom_path(linewidth = 0.8) +
    geom_point(size = 2.5) +
    geom_errorbar(aes(xmin = pmax(0, mean_enh_ppb - sd_enh_ppb),
                      xmax = mean_enh_ppb + sd_enh_ppb),
                  orientation = "y", width = 0, linewidth = 0.5, alpha = 0.6) +
    scale_colour_viridis_d(option = "C", end = 0.9, name = "Ring") +
    labs(
      title = paste("Plume vertical profile:", scenario_label),
      subtitle = sprintf("Time %g min - mean enhancement by height", t_use),
      x = expression(paste("CH"[4], " enhancement (ppb)")),
      y = "Height z (m)"
    ) +
    theme_osse()

  if (save) {
    fname <- sprintf("plume_z_profile_%s.png", gsub("[^a-zA-Z0-9_]", "_", scenario_label))
    ggsave_osse(file.path(get_figure_dir(), fname), p, width = 8, height = 6, dpi = 200)
    message(sprintf("  Saved: %s", fname))
  }
  invisible(p)
}


# -----------------------------------------------------------------------------
# 6b2. Isometric plot colour limits (global min/max over time for static + anim)
# -----------------------------------------------------------------------------
#' @keywords internal
.safe_colour_limits_pair <- function(lo, hi) {
  if (!is.finite(lo) || !is.finite(hi)) return(NULL)
  if (lo > hi) {
    tmp <- lo
    lo <- hi
    hi <- tmp
  }
  if (lo == hi) {
    eps <- max(abs(lo) * 1e-6, 1e-9)
    lo <- lo - eps
    hi <- hi + eps
  }
  c(lo, hi)
}

#' Compute fill (emission flux) and per–plot-type colour limits over all time steps.
#'
#' Used so static isometric PNGs and frame-by-frame animations share the same
#' colour scale (min/max across the full scenario time range).
#'
#' @param truth_dt,measurements_dt,receptors,flux_field As in isometric plots.
#' @param scenario_label Passed to wind-perturbation RNG (must match plots).
#' @param ws_observation_error_m_s Observation noise on measured wind (m/s).
#' @return A list with \code{fill}, \code{truth_enhancement}, \code{measured_enhancement},
#'   \code{measured_wind}, \code{true_wind}, \code{methane_diff}, \code{wind_diff},
#'   \code{fluxvector_measured}, \code{fluxvector_true}; each NULL or numeric c(lo, hi).
compute_isometric_scale_limits <- function(truth_dt,
                                             measurements_dt,
                                             receptors,
                                             flux_field,
                                             scenario_label = "default",
                                             ws_observation_error_m_s = 0) {
  out <- list(
    fill                  = NULL,
    truth_enhancement     = NULL,
    measured_enhancement  = NULL,
    measured_wind         = NULL,
    true_wind             = NULL,
    methane_diff          = NULL,
    wind_diff             = NULL,
    fluxvector_measured   = NULL,
    fluxvector_true       = NULL,
    contam_fv             = NULL
  )

  if ("flux_mg_m2_d" %in% names(flux_field)) {
    out$fill <- .safe_colour_limits_pair(
      min(flux_field$flux_mg_m2_d, na.rm = TRUE),
      max(flux_field$flux_mg_m2_d, na.rm = TRUE)
    )
  }

  if ("ch4_enhancement_ppb" %in% names(truth_dt)) {
    out$truth_enhancement <- .safe_colour_limits_pair(
      min(truth_dt$ch4_enhancement_ppb, na.rm = TRUE),
      max(truth_dt$ch4_enhancement_ppb, na.rm = TRUE)
    )
  }

  if ("ch4_measured_ppb" %in% names(measurements_dt)) {
    if ("ch4_bg_ppb" %in% names(measurements_dt)) {
      v <- measurements_dt$ch4_measured_ppb - measurements_dt$ch4_bg_ppb
    } else {
      v <- measurements_dt$ch4_measured_ppb
    }
    out$measured_enhancement <- .safe_colour_limits_pair(min(v, na.rm = TRUE), max(v, na.rm = TRUE))
  }

  meas_ws_col <- if ("ws_at_z_m_s" %in% names(measurements_dt)) {
    "ws_at_z_m_s"
  } else if ("ws_m_s" %in% names(measurements_dt)) {
    "ws_m_s"
  } else {
    NA_character_
  }
  if (!is.na(meas_ws_col)) {
    mp <- copy(measurements_dt)
    col <- meas_ws_col
    if (ws_observation_error_m_s > 0) {
      rng_state <- if (exists(".Random.seed", envir = .GlobalEnv)) {
        get(".Random.seed", envir = .GlobalEnv)
      } else {
        NULL
      }
      set.seed(string_to_seed(scenario_label) + 7778L)
      mp[, ws_obs_plot := pmax(0.1, get(meas_ws_col) + rnorm(.N, 0, ws_observation_error_m_s))]
      if (!is.null(rng_state)) assign(".Random.seed", rng_state, envir = .GlobalEnv)
      col <- "ws_obs_plot"
    }
    out$measured_wind <- .safe_colour_limits_pair(
      min(mp[[col]], na.rm = TRUE),
      max(mp[[col]], na.rm = TRUE)
    )
  }

  truth_ws_col <- if ("ws_at_z_m_s" %in% names(truth_dt)) {
    "ws_at_z_m_s"
  } else if ("ws_m_s" %in% names(truth_dt)) {
    "ws_m_s"
  } else {
    NA_character_
  }
  if (!is.na(truth_ws_col)) {
    out$true_wind <- .safe_colour_limits_pair(
      min(truth_dt[[truth_ws_col]], na.rm = TRUE),
      max(truth_dt[[truth_ws_col]], na.rm = TRUE)
    )
  }

  if ("ch4_total_ppb" %in% names(truth_dt) && "ch4_measured_ppb" %in% names(measurements_dt)) {
    meas_cols <- c("receptor_id", "time_min", "ch4_measured_ppb")
    truth_cols <- c("receptor_id", "time_min", "ch4_total_ppb")
    md <- merge(
      measurements_dt[, ..meas_cols],
      truth_dt[, ..truth_cols],
      by = c("receptor_id", "time_min"),
      all = FALSE
    )
    if (nrow(md) > 0L) {
      md[, delta_ch4_ppb := ch4_measured_ppb - ch4_total_ppb]
      out$methane_diff <- .safe_colour_limits_pair(
        min(md$delta_ch4_ppb, na.rm = TRUE),
        max(md$delta_ch4_ppb, na.rm = TRUE)
      )
    }
  }

  if (!is.na(meas_ws_col) && !is.na(truth_ws_col)) {
    meas_cols <- c("receptor_id", "time_min", meas_ws_col)
    truth_cols <- c("receptor_id", "time_min", truth_ws_col)
    dt_all <- merge(
      measurements_dt[, ..meas_cols],
      truth_dt[, ..truth_cols],
      by = c("receptor_id", "time_min"),
      all = FALSE,
      suffixes = c("_meas", "_true")
    )
    ws_meas_name <- paste0(meas_ws_col, "_meas")
    ws_true_name <- paste0(truth_ws_col, "_true")
    wind_diff_vals <- numeric(0)
    for (t in sort(unique(dt_all$time_min))) {
      d <- dt_all[time_min == t]
      if (nrow(d) == 0L) next
      if (ws_observation_error_m_s > 0) {
        rng_state <- if (exists(".Random.seed", envir = .GlobalEnv)) {
          get(".Random.seed", envir = .GlobalEnv)
        } else {
          NULL
        }
        set.seed(string_to_seed(scenario_label) + 7778L)
        d[, ws_meas_for_diff := pmax(0.1, get(ws_meas_name) + rnorm(.N, 0, ws_observation_error_m_s))]
        if (!is.null(rng_state)) assign(".Random.seed", rng_state, envir = .GlobalEnv)
        ws_m <- "ws_meas_for_diff"
      } else {
        ws_m <- ws_meas_name
      }
      d[, delta_ws_m_s := get(ws_m) - get(ws_true_name)]
      wind_diff_vals <- c(wind_diff_vals, d$delta_ws_m_s)
    }
    if (length(wind_diff_vals) > 0L) {
      out$wind_diff <- .safe_colour_limits_pair(
        min(wind_diff_vals, na.rm = TRUE),
        max(wind_diff_vals, na.rm = TRUE)
      )
    }
  }

  if (all(c("ch4_measured_ppb", "ws_at_z_m_s", "wd_deg") %in% names(measurements_dt))) {
    dt <- copy(measurements_dt)
    if (!("azimuth_deg" %in% names(dt))) {
      dt <- merge(dt, receptors[, .(receptor_id, x, y, azimuth_deg)], by = "receptor_id", all.x = TRUE)
    }
    if (!("azimuth_deg" %in% names(dt)) || any(is.na(dt$azimuth_deg))) {
      if (!all(c("x", "y") %in% names(dt))) {
        dt <- merge(dt, receptors[, .(receptor_id, x, y)], by = "receptor_id", all.x = TRUE)
      }
      dt[, azimuth_deg := (atan2(x, y) * 180 / pi) %% 360]
    }
    dt[, wd_to := (wd_deg + 180) %% 360]
    dt[, v_r := ws_at_z_m_s * cos((azimuth_deg - wd_to) * pi / 180)]
    dt[, flux_vector_sum := ch4_measured_ppb * v_r]
    out$fluxvector_measured <- .safe_colour_limits_pair(
      min(dt$flux_vector_sum, na.rm = TRUE),
      max(dt$flux_vector_sum, na.rm = TRUE)
    )
  }

  if (all(c("ch4_total_ppb", "ws_at_z_m_s", "wd_deg") %in% names(truth_dt))) {
    dt <- copy(truth_dt)
    if (!("azimuth_deg" %in% names(dt))) {
      dt <- merge(dt, receptors[, .(receptor_id, x, y, azimuth_deg)], by = "receptor_id", all.x = TRUE)
    }
    if (!("azimuth_deg" %in% names(dt)) || any(is.na(dt$azimuth_deg))) {
      if (!all(c("x", "y") %in% names(dt))) {
        dt <- merge(dt, receptors[, .(receptor_id, x, y)], by = "receptor_id", all.x = TRUE)
      }
      dt[, azimuth_deg := (atan2(x, y) * 180 / pi) %% 360]
    }
    dt[, wd_to := (wd_deg + 180) %% 360]
    dt[, v_r := ws_at_z_m_s * cos((azimuth_deg - wd_to) * pi / 180)]
    dt[, flux_vector_sum := ch4_total_ppb * v_r]
    out$fluxvector_true <- .safe_colour_limits_pair(
      min(dt$flux_vector_sum, na.rm = TRUE),
      max(dt$flux_vector_sum, na.rm = TRUE)
    )
  }

  if (all(c("ch4_contamination_ppb", "ws_at_z_m_s", "wd_deg") %in% names(truth_dt)) &&
      any(truth_dt$ch4_contamination_ppb != 0, na.rm = TRUE)) {
    dt <- copy(truth_dt)
    if (!("azimuth_deg" %in% names(dt))) {
      dt <- merge(dt, receptors[, .(receptor_id, x, y, azimuth_deg)],
                  by = "receptor_id", all.x = TRUE)
    }
    dt[, wd_to := (wd_deg + 180) %% 360]
    dt[, v_r := ws_at_z_m_s * cos((azimuth_deg - wd_to) * pi / 180)]
    dt[, contam_fv := ch4_contamination_ppb * v_r]
    out$contam_fv <- .safe_colour_limits_pair(
      min(dt$contam_fv, na.rm = TRUE),
      max(dt$contam_fv, na.rm = TRUE)
    )
  }

  out
}


# -----------------------------------------------------------------------------
# 6c. Isometric plume view over emission map
# -----------------------------------------------------------------------------
#
#' Plot an isometric-looking 2D view of receptors colored by plume enhancement,
#' overlaid on the emission map as a ground plane.
#'
#' This uses a simple isometric-like projection onto screen coordinates:
#'   x2d = (x - y) * cos(iso_angle_deg)
#'   y2d = z * z_scale + (x + y) * sin(iso_angle_deg)
#' where iso_angle_deg is currently set to 15 degrees (30 degrees lower
#' viewing angle than the previous 45-degree setup).
#'
#' @param truth_dt data.table from run_forward_model.
#' @param receptors data.table from define_receptor_config.
#' @param flux_field data.table from generate_effective_flux_field().
#' @param lake_diameter Lake diameter in metres.
#' @param time_min Time step to plot (defaults to first).
#' @param scenario_label Character label.
#' @param z_scale Vertical exaggeration for z (m -> projected units).
#' @param gaussian_surface_grid Optional integer: horizontal resolution (per
#'   axis) for a semi-transparent analytical Gaussian plume cloud under the
#'   receptor points. Only when \code{truth_dt} contains \code{pasquill_class}
#'   and this is \code{>= 10} (e.g. \code{40}). Default \code{NULL} = \code{0}:
#'   plot only receptors. When \code{nearby_lake_config} is
#'   non-NULL, the dense field includes both disk sources (same sum as forward).
#' @param save Logical.
plot_plume_isometric_over_emission_map <- function(truth_dt,
                                                     receptors,
                                                     flux_field,
                                                     lake_diameter,
                                                     time_min = NULL,
                                                     scenario_label = "default",
                                                     z_scale = 0.35,
                                                     fill_limits = NULL,
                                                     colour_limits = NULL,
                                                     nearby_lake_config = NULL,
                                                     gaussian_surface_grid = NULL,
                                                     save = TRUE) {
  # Choose a time slice.
  t_use <- if (is.null(time_min)) min(truth_dt$time_min) else time_min
  tr_slice <- truth_dt[time_min == t_use]
  if (nrow(tr_slice) == 0) {
    message("  No data for isometric plume plot")
    return(invisible())
  }

  cols_truth <- c("receptor_id", "ch4_enhancement_ppb")
  for (cn in c("pasquill_class", "ws_m_s", "wd_deg", "flux_true_mg_m2_d")) {
    if (cn %in% names(tr_slice)) cols_truth <- c(cols_truth, cn)
  }
  cols_truth <- unique(cols_truth)

  # Use ring_type + receptor geometry from receptors table.
  dt_t <- merge(
    unique(tr_slice[, ..cols_truth]),
    receptors[, .(receptor_id, x, y, z, ring_type)],
    by = "receptor_id",
    all.x = TRUE
  )

  iso_angle_deg <- 176
  cos_iso <- cos(iso_angle_deg * pi / 180)
  sin_iso <- sin(iso_angle_deg * pi / 180)

  # Project receptor points.
  dt_t[, `:=`(
    x2d = (x - y) * cos_iso,
    y2d = (x + y) * sin_iso + z * z_scale
  )]

  # Optional dense analytical Gaussian field. By default, only receptors are plotted.
  gsg <- gaussian_surface_grid
  if (is.null(gsg)) {
    gsg <- 0L
  } else {
    gsg <- suppressWarnings(as.integer(gsg))
    if (length(gsg) != 1L || is.na(gsg)) gsg <- 0L
  }

  dense_dt <- NULL
  if (gsg >= 10L &&
      all(c("pasquill_class", "ws_m_s", "wd_deg", "flux_true_mg_m2_d") %in% names(dt_t))) {
    if (!exists("dense_gaussian_enhancement_for_isometric", mode = "function")) {
      source(file.path(project_root(), "R", "02b_gaussian_forward.R"))
    }
    dense_dt <- tryCatch(
      dense_gaussian_enhancement_for_isometric(
        lake_diameter_m      = lake_diameter,
        flux_mg_m2_d         = dt_t$flux_true_mg_m2_d[1],
        ws_m_s               = dt_t$ws_m_s[1],
        wd_deg               = dt_t$wd_deg[1],
        pasquill_class       = as.character(dt_t$pasquill_class[1]),
        heights_m            = sort(unique(dt_t$z)),
        grid_n               = gsg,
        nearby_lake_config   = nearby_lake_config
      ),
      error = function(e) {
        message("  Dense Gaussian isometric cloud skipped: ", conditionMessage(e))
        NULL
      }
    )
    if (!is.null(dense_dt) && nrow(dense_dt) > 0L) {
      dense_dt[, `:=`(
        x2d = (x - y) * cos_iso,
        y2d = (x + y) * sin_iso + z * z_scale
      )]
    }
  }

  # Choose an emission-map snapshot hour (flux is constant without diel, but
  # we keep this for robustness).
  snap_hour <- NULL
  if ("hour" %in% names(flux_field)) {
    if (any(flux_field$hour == 12)) snap_hour <- 12 else snap_hour <- unique(flux_field$hour)[1]
  }
  snap <- if (is.null(snap_hour)) flux_field else flux_field[hour == snap_hour]

  if (nrow(snap) == 0) {
    message("  No data for isometric emission-map underlay")
    return(invisible())
  }

  # Project emission map ground plane (z=0).
  if (!("flux_mg_m2_d" %in% names(snap))) {
    message("  Emission map underlay missing `flux_mg_m2_d`")
    return(invisible())
  }

  snap_xy <- copy(snap)[, .(
    x = x,
    y = y,
    flux = flux_mg_m2_d
  )]
  snap_xy[, `:=`(
    x2d = (x - y) * cos_iso,
    y2d = (x + y) * sin_iso
  )]

  # Estimate underlying emission grid cell size in the original (x,y)
  # coordinate system. We then project each cell's 4 corners so the
  # underlay is continuous (no checkerboard gaps from axis-aligned tiles).
  ux0 <- sort(unique(snap_xy$x))
  uy0 <- sort(unique(snap_xy$y))
  dx0 <- if (length(ux0) >= 2L) median(diff(ux0)) else NA_real_
  dy0 <- if (length(uy0) >= 2L) median(diff(uy0)) else NA_real_
  half_dx <- if (is.finite(dx0)) dx0 / 2 else 1
  half_dy <- if (is.finite(dy0)) dy0 / 2 else 1

  # Append nearby lake emission cells to the ground-plane underlay.
  if (!is.null(nearby_lake_config) && is.finite(dx0) && is.finite(dy0)) {
    nl_xs <- seq(nearby_lake_config$center_x - nearby_lake_config$radius_m,
                 nearby_lake_config$center_x + nearby_lake_config$radius_m,
                 by = dx0)
    nl_ys <- seq(nearby_lake_config$center_y - nearby_lake_config$radius_m,
                 nearby_lake_config$center_y + nearby_lake_config$radius_m,
                 by = dy0)
    nl_grid <- CJ(x = nl_xs, y = nl_ys)
    nl_grid <- nl_grid[
      (x - nearby_lake_config$center_x)^2 +
        (y - nearby_lake_config$center_y)^2 <= nearby_lake_config$radius_m^2
    ]
    if (nrow(nl_grid) > 0L) {
      nl_flux_val <- if (!is.null(nearby_lake_config$flux_mg_m2_d)) nearby_lake_config$flux_mg_m2_d else median(snap_xy$flux, na.rm = TRUE)
      nl_grid[, flux := nl_flux_val]
      nl_grid[, `:=`(x2d = (x - y) * cos_iso, y2d = (x + y) * sin_iso)]
      snap_xy <- rbind(snap_xy, nl_grid)
    }
  }

  # Project the lake boundary circle (z=0).
  lake_radius <- lake_diameter / 2
  theta <- seq(0, 2 * pi, length.out = 220)
  lake_x <- lake_radius * cos(theta)
  lake_y <- lake_radius * sin(theta)
  lake_proj <- data.table(
    x2d = (lake_x - lake_y) * cos_iso,
    y2d = (lake_x + lake_y) * sin_iso
  )

  # Project nearby lake boundary circle if present.
  nearby_lake_proj <- NULL
  if (!is.null(nearby_lake_config)) {
    nl_x <- nearby_lake_config$center_x + nearby_lake_config$radius_m * cos(theta)
    nl_y <- nearby_lake_config$center_y + nearby_lake_config$radius_m * sin(theta)
    nearby_lake_proj <- data.table(
      x2d = (nl_x - nl_y) * cos_iso,
      y2d = (nl_x + nl_y) * sin_iso
    )
  }

  poly_data <- {
    n_cells <- nrow(snap_xy)
    if (n_cells == 0L) NULL else {
      cell_id <- seq_len(n_cells)
      x_left   <- snap_xy$x - half_dx
      x_right  <- snap_xy$x + half_dx
      y_bottom <- snap_xy$y - half_dy
      y_top    <- snap_xy$y + half_dy
      x_corner <- as.vector(t(cbind(x_left, x_right, x_right, x_left)))
      y_corner <- as.vector(t(cbind(y_bottom, y_bottom, y_top, y_top)))
      flux_vec <- rep(snap_xy$flux, each = 4L)
      id_vec   <- rep(cell_id, each = 4L)
      data.table(
        cell_id = id_vec,
        flux = flux_vec,
        x2d = (x_corner - y_corner) * cos_iso,
        y2d = (x_corner + y_corner) * sin_iso
      )
    }
  }

  p <- ggplot() +
    geom_polygon(
      data = poly_data,
      aes(x = x2d, y = y2d, group = cell_id, fill = flux),
      color = NA
    ) +
    geom_path(
      data = lake_proj,
      aes(x = x2d, y = y2d),
      colour = "grey20",
      linewidth = 0.6
    )

  if (!is.null(dense_dt) && nrow(dense_dt) > 0L) {
    p <- p + geom_point(
      data = dense_dt,
      aes(x = x2d, y = y2d, colour = ch4_enhancement_ppb),
      size = 0.35,
      alpha = 0.2,
      shape = 16
    )
  }

  p <- p +
    geom_point(
      data = dt_t,
      aes(x = x2d, y = y2d, colour = ch4_enhancement_ppb, shape = ring_type),
      size = 2.6,
      alpha = 0.92
    ) +
    scale_fill_viridis_c(
      option = "A",
      name = "Emission Flux (mg CH4 m^-2 d^-1)",
      limits = fill_limits,
      oob = if (!is.null(fill_limits)) scales::squish else ggplot2::waiver()
    ) +
    scale_colour_viridis_c(
      option = "plasma",
      name = "Methane Enhancement (ppb)",
      limits = colour_limits,
      oob = if (!is.null(colour_limits)) scales::squish else ggplot2::waiver()
    ) +
    labs(
      title = sprintf("Isometric plume over emission map: %s", scenario_label),
      subtitle = sprintf("Time %g min (projection z_scale=%g)", t_use, z_scale),
      x = NULL,
      y = NULL
    ) +
    theme_osse() +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      legend.position = "right"
    ) +
    guides(
      fill = guide_colorbar(order = 1),
      colour = guide_colorbar(order = 2),
      shape = "none"
    )

  if (!is.null(nearby_lake_proj)) {
    p <- p + geom_path(
      data = nearby_lake_proj,
      aes(x = x2d, y = y2d),
      colour = "grey20",
      linewidth = 0.6,
      linetype = "dashed"
    )
  }

  if (save) {
    fname <- sprintf(
      "plume_isometric_emission_map_%s.png",
      gsub("[^a-zA-Z0-9_]", "_", scenario_label)
    )
    ggsave_osse(file.path(get_figure_dir(), fname), p, width = 10, height = 7, dpi = 200)
    message(sprintf("  Saved: %s", fname))
  }

  invisible(p)
}


# -----------------------------------------------------------------------------
# 6d. Isometric measured values view over emission map
# -----------------------------------------------------------------------------
#
#' Plot an isometric-looking 2D view of *measured* receptor values over the
#' emission map ground plane (same projection as the truth isometric plot).
#'
#' @param measurements_dt data.table from simulate_drone_measurements().
#' @param receptors data.table from define_receptor_config().
#' @param flux_field data.table from generate_effective_flux_field().
#' @param lake_diameter Lake diameter in metres.
#' @param time_min Time step to plot (defaults to first).
#' @param scenario_label Character label.
#' @param z_scale Vertical exaggeration for z (m -> projected units).
#' @param value_mode Which measured value to plot:
#'        - "enhancement" (default): ch4_measured_ppb - ch4_bg_ppb (if present)
#'        - "total": ch4_measured_ppb
#'        - "custom": uses value_col
#' @param value_col Column name when value_mode="custom".
#' @param save Logical.
plot_measured_isometric_over_emission_map <- function(measurements_dt,
                                                      receptors,
                                                      flux_field,
                                                      lake_diameter,
                                                      time_min = NULL,
                                                      scenario_label = "default",
                                                      z_scale = 0.35,
                                                      value_mode = "enhancement",
                                                      value_col = NULL,
                                                      value_label_override = NULL,
                                                      fill_limits = NULL,
                                                      colour_limits = NULL,
                                                      nearby_lake_config = NULL,
                                                      emission_underlay_mode = c("cartesian", "polar"),
                                                      polar_n_rings = NULL,
                                                      polar_n_azimuth = NULL,
                                                      polar_radial_step_m = NULL,
                                                      scale_point_alpha_by_value = FALSE,
                                                      point_alpha_floor = 0.25,
                                                      point_alpha_ceiling = 0.95,
                                                      point_size = 2.6,
                                                      base_size = 11,
                                                      save = TRUE) {
  emission_underlay_mode <- match.arg(emission_underlay_mode)
  # Choose a time slice.
  t_use <- if (is.null(time_min)) min(measurements_dt$time_min) else time_min
  meas_t <- measurements_dt[time_min == t_use]
  if (nrow(meas_t) == 0) {
    message("  No data for measured isometric plume plot")
    return(invisible())
  }

  iso_angle_deg <- 176
  cos_iso <- cos(iso_angle_deg * pi / 180)
  sin_iso <- sin(iso_angle_deg * pi / 180)

  # Join ring_type and receptor geometry.
  receptor_join <- receptors
  if (!("ring_type" %in% names(receptor_join))) {
    receptor_join[, ring_type := "unknown"]
  }
  meas_t <- merge(
    meas_t,
    receptor_join[, .(receptor_id, x, y, z, ring_type)],
    by = "receptor_id",
    all.x = TRUE
  )

  # Final safety net: ensure ring_type exists for geom_point aesthetics.
  if (!("ring_type" %in% names(meas_t))) {
    meas_t[, ring_type := "unknown"]
  }

  # Prefer measured coordinates if present (GPS error option).
  if ("x_measured" %in% names(meas_t)) meas_t[, x_plot := x_measured]
  else meas_t[, x_plot := x]
  if ("y_measured" %in% names(meas_t)) meas_t[, y_plot := y_measured]
  else meas_t[, y_plot := y]
  if ("z_measured" %in% names(meas_t)) meas_t[, z_plot := z_measured]
  else meas_t[, z_plot := z]

  # Pick measured value for coloring.
  if (value_mode == "enhancement") {
    if ("ch4_bg_ppb" %in% names(meas_t) && "ch4_measured_ppb" %in% names(meas_t)) {
      meas_t[, value_plot := ch4_measured_ppb - ch4_bg_ppb]
      value_label <- expression(CH[4]~"enhancement (ppb)")
    } else if ("ch4_measured_ppb" %in% names(meas_t)) {
      meas_t[, value_plot := ch4_measured_ppb]
      value_label <- "Measured total (ppb)"
    } else {
      stop("Measured value columns not found: need ch4_measured_ppb (and ch4_bg_ppb for enhancement).")
    }
  } else if (value_mode == "total") {
    if (!("ch4_measured_ppb" %in% names(meas_t))) {
      stop("Column `ch4_measured_ppb` not found in measurements_dt.")
    }
    meas_t[, value_plot := ch4_measured_ppb]
    value_label <- "Measured total (ppb)"
  } else if (value_mode == "custom") {
    if (is.null(value_col) || !(value_col %in% names(meas_t))) {
      stop("value_mode='custom' requires value_col to be set and present in measurements_dt.")
    }
    meas_t[, value_plot := get(value_col)]
    value_label <- value_col
  } else {
    stop("Unknown value_mode: ", value_mode)
  }
  if (!is.null(value_label_override) && nzchar(value_label_override)) {
    value_label <- value_label_override
  }

  if (isTRUE(scale_point_alpha_by_value)) {
    mag <- abs(meas_t$value_plot)
    finite_mag <- mag[is.finite(mag)]
    if (length(finite_mag) == 0L) {
      meas_t[, alpha_plot := point_alpha_floor]
    } else {
      mag_min <- min(finite_mag)
      mag_max <- max(finite_mag)
      if (!is.finite(mag_min) || !is.finite(mag_max) || mag_max <= mag_min) {
        meas_t[, alpha_plot := point_alpha_ceiling]
      } else {
        meas_t[, alpha_plot := point_alpha_floor +
                 (pmin(pmax(abs(value_plot), mag_min), mag_max) - mag_min) /
                 (mag_max - mag_min) *
                 (point_alpha_ceiling - point_alpha_floor)]
      }
    }
  } else {
    meas_t[, alpha_plot := point_alpha_ceiling]
  }

  # Project receptor points.
  meas_t[, `:=`(
    x2d = (x_plot - y_plot) * cos_iso,
    y2d = (x_plot + y_plot) * sin_iso + z_plot * z_scale
  )]

  # Choose an emission-map snapshot hour.
  snap_hour <- NULL
  if ("hour" %in% names(flux_field)) {
    if (any(flux_field$hour == 12)) snap_hour <- 12 else snap_hour <- unique(flux_field$hour)[1]
  }
  snap <- if (is.null(snap_hour)) flux_field else flux_field[hour == snap_hour]

  if (nrow(snap) == 0 || !("flux_mg_m2_d" %in% names(snap))) {
    message("  No data for isometric emission-map underlay")
    return(invisible())
  }

  snap_xy <- copy(snap)[, .(
    x = x,
    y = y,
    flux = flux_mg_m2_d
  )]

  snap_xy[, `:=`(
    x2d = (x - y) * cos_iso,
    y2d = (x + y) * sin_iso
  )]

  # Estimate underlying emission grid cell half-sizes in original (x,y).
  ux0 <- sort(unique(snap_xy$x))
  uy0 <- sort(unique(snap_xy$y))
  dx0 <- if (length(ux0) >= 2L) median(diff(ux0)) else NA_real_
  dy0 <- if (length(uy0) >= 2L) median(diff(uy0)) else NA_real_
  half_dx <- if (is.finite(dx0)) dx0 / 2 else 1
  half_dy <- if (is.finite(dy0)) dy0 / 2 else 1

  # Append nearby lake emission cells to the ground-plane underlay.
  if (!is.null(nearby_lake_config) && is.finite(dx0) && is.finite(dy0)) {
    nl_xs <- seq(nearby_lake_config$center_x - nearby_lake_config$radius_m,
                 nearby_lake_config$center_x + nearby_lake_config$radius_m,
                 by = dx0)
    nl_ys <- seq(nearby_lake_config$center_y - nearby_lake_config$radius_m,
                 nearby_lake_config$center_y + nearby_lake_config$radius_m,
                 by = dy0)
    nl_grid <- CJ(x = nl_xs, y = nl_ys)
    nl_grid <- nl_grid[
      (x - nearby_lake_config$center_x)^2 +
        (y - nearby_lake_config$center_y)^2 <= nearby_lake_config$radius_m^2
    ]
    if (nrow(nl_grid) > 0L) {
      nl_flux_val <- if (!is.null(nearby_lake_config$flux_mg_m2_d)) nearby_lake_config$flux_mg_m2_d else median(snap_xy$flux, na.rm = TRUE)
      nl_grid[, flux := nl_flux_val]
      nl_grid[, `:=`(x2d = (x - y) * cos_iso, y2d = (x + y) * sin_iso)]
      snap_xy <- rbind(snap_xy, nl_grid)
    }
  }

  # Project the lake boundary circle (z=0).
  lake_radius <- lake_diameter / 2
  theta <- seq(0, 2 * pi, length.out = 220)
  lake_x <- lake_radius * cos(theta)
  lake_y <- lake_radius * sin(theta)
  lake_proj <- data.table(
    x2d = (lake_x - lake_y) * cos_iso,
    y2d = (lake_x + lake_y) * sin_iso
  )

  # Project nearby lake boundary circle if present.
  nearby_lake_proj <- NULL
  if (!is.null(nearby_lake_config)) {
    nl_x <- nearby_lake_config$center_x + nearby_lake_config$radius_m * cos(theta)
    nl_y <- nearby_lake_config$center_y + nearby_lake_config$radius_m * sin(theta)
    nearby_lake_proj <- data.table(
      x2d = (nl_x - nl_y) * cos_iso,
      y2d = (nl_x + nl_y) * sin_iso
    )
  }

  emission_polygons <- if (emission_underlay_mode == "polar") {
    n_rings <- if (is.null(polar_n_rings) || !is.finite(polar_n_rings) || polar_n_rings < 1) {
      NA_integer_
    } else {
      as.integer(round(polar_n_rings))
    }
    radial_step <- if (is.na(n_rings)) {
      if (is.null(polar_radial_step_m) || !is.finite(polar_radial_step_m) || polar_radial_step_m <= 0) {
        if (is.finite(dx0)) dx0 else 5
      } else {
        polar_radial_step_m
      }
    } else {
      lake_radius / n_rings
    }
    n_az <- if (is.null(polar_n_azimuth) || !is.finite(polar_n_azimuth) || polar_n_azimuth < 3) {
      24L
    } else {
      as.integer(round(polar_n_azimuth))
    }
    if (!is.na(n_rings)) {
      r_edges <- seq(0, lake_radius, length.out = n_rings + 1L)
    } else {
      r_edges <- seq(0, lake_radius, by = radial_step)
      if (tail(r_edges, 1) < lake_radius) r_edges <- c(r_edges, lake_radius)
      if (length(r_edges) < 2L) r_edges <- c(0, lake_radius)
    }
    theta_edges <- seq(0, 2 * pi, length.out = n_az + 1L)
    flux_val <- median(snap_xy$flux, na.rm = TRUE)
    n_arc <- 14L
    poly_rows <- vector("list", (length(r_edges) - 1L) * n_az)
    idx <- 0L
    for (ir in seq_len(length(r_edges) - 1L)) {
      r_in <- r_edges[ir]
      r_out <- r_edges[ir + 1L]
      for (ia in seq_len(n_az)) {
        th0 <- theta_edges[ia]
        th1 <- theta_edges[ia + 1L]
        th_seq <- seq(th0, th1, length.out = n_arc)
        x_outer <- r_out * cos(th_seq)
        y_outer <- r_out * sin(th_seq)
        x_inner <- rev(r_in * cos(th_seq))
        y_inner <- rev(r_in * sin(th_seq))
        x_corner <- c(x_outer, x_inner)
        y_corner <- c(y_outer, y_inner)
        idx <- idx + 1L
        poly_rows[[idx]] <- data.table(
          cell_id = idx,
          flux = flux_val,
          x2d = (x_corner - y_corner) * cos_iso,
          y2d = (x_corner + y_corner) * sin_iso
        )
      }
    }
    rbindlist(poly_rows, fill = TRUE)
  } else {
    n_cells <- nrow(snap_xy)
    if (n_cells == 0L) {
      NULL
    } else {
      cell_id <- seq_len(n_cells)

      x_left   <- snap_xy$x - half_dx
      x_right  <- snap_xy$x + half_dx
      y_bottom <- snap_xy$y - half_dy
      y_top    <- snap_xy$y + half_dy

      x_corner <- as.vector(t(cbind(x_left, x_right, x_right, x_left)))
      y_corner <- as.vector(t(cbind(y_bottom, y_bottom, y_top, y_top)))

      flux_vec <- rep(snap_xy$flux, each = 4L)
      id_vec   <- rep(cell_id, each = 4L)

      data.table(
        cell_id = id_vec,
        flux = flux_vec,
        x2d = (x_corner - y_corner) * cos_iso,
        y2d = (x_corner + y_corner) * sin_iso
      )
    }
  }

  p <- ggplot() +
    geom_polygon(
      data = emission_polygons,
      aes(x = x2d, y = y2d, group = cell_id, fill = flux),
      colour = if (emission_underlay_mode == "polar") alpha("grey100", 0.65) else NA,
      linewidth = if (emission_underlay_mode == "polar") 0.16 else 0
    ) +
    geom_path(
      data = lake_proj,
      aes(x = x2d, y = y2d),
      colour = "grey20",
      linewidth = 0.6
    ) +
    geom_point(
      data = meas_t,
      aes(x = x2d, y = y2d, colour = value_plot, shape = ring_type, alpha = alpha_plot),
      size = point_size
    ) +
    scale_fill_viridis_c(
      option = "A",
      name = "Emission Flux (mg CH4 m^-2 d^-1)",
      limits = fill_limits,
      oob = if (!is.null(fill_limits)) scales::squish else ggplot2::waiver()
    ) +
    scale_colour_viridis_c(
      option = "plasma",
      name = value_label,
      limits = colour_limits,
      oob = if (!is.null(colour_limits)) scales::squish else ggplot2::waiver()
    ) +
    scale_alpha_identity() +
    labs(
      title = sprintf("Isometric measured overlay: %s", scenario_label),
      subtitle = sprintf("Time %g min", t_use),
      x = NULL,
      y = NULL
    ) +
    theme_osse() +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      legend.position = "right"
    ) +
    guides(
      fill = guide_colorbar(order = 1),
      colour = guide_colorbar(order = 2),
      shape = "none"
    )

  if (!is.null(nearby_lake_proj)) {
    p <- p + geom_path(
      data = nearby_lake_proj,
      aes(x = x2d, y = y2d),
      colour = "grey20",
      linewidth = 0.6,
      linetype = "dashed"
    )
  }

  if (save) {
    fname <- sprintf(
      "plume_isometric_measured_emission_map_%s.png",
      gsub("[^a-zA-Z0-9_]", "_", scenario_label)
    )
    ggsave_osse(file.path(get_figure_dir(), fname), p, width = 10, height = 7, dpi = 200)
    message(sprintf("  Saved: %s", fname))
  }

  invisible(p)
}


# -----------------------------------------------------------------------------
# 6e. Isometric measured wind-speed view over emission map
# -----------------------------------------------------------------------------
plot_windspeed_isometric_measured_over_emission_map <- function(measurements_dt,
                                                                receptors,
                                                                flux_field,
                                                                lake_diameter,
                                                                time_min = NULL,
                                                                scenario_label = "default",
                                                                z_scale = 0.35,
                                                                ws_observation_error_m_s = 0,
                                                                fill_limits = NULL,
                                                                colour_limits = NULL,
                                                                nearby_lake_config = NULL,
                                                                save = TRUE) {
  if ("ws_at_z_m_s" %in% names(measurements_dt)) {
    value_col <- "ws_at_z_m_s"
  } else if ("ws_m_s" %in% names(measurements_dt)) {
    value_col <- "ws_m_s"
  } else {
    stop("No wind-speed column found in measurements_dt (`ws_at_z_m_s` or `ws_m_s`).")
  }

  measurements_plot <- copy(measurements_dt)
  if (ws_observation_error_m_s > 0) {
    rng_state <- if (exists(".Random.seed", envir = .GlobalEnv))
      get(".Random.seed", envir = .GlobalEnv) else NULL
    set.seed(string_to_seed(scenario_label) + 7778L)
    measurements_plot[, ws_obs_plot := pmax(0.1, get(value_col) + rnorm(.N, 0, ws_observation_error_m_s))]
    if (!is.null(rng_state)) assign(".Random.seed", rng_state, envir = .GlobalEnv)
    value_col <- "ws_obs_plot"
  }

  p <- plot_measured_isometric_over_emission_map(
    measurements_dt = measurements_plot,
    receptors = receptors,
    flux_field = flux_field,
    lake_diameter = lake_diameter,
    time_min = time_min,
    scenario_label = scenario_label,
    z_scale = z_scale,
    value_mode = "custom",
    value_col = value_col,
    value_label_override = "Measured Wind Speed (m s^-1)",
    fill_limits = fill_limits,
    colour_limits = colour_limits,
    nearby_lake_config = nearby_lake_config,
    save = FALSE
  )

  p <- p + labs(
    title = sprintf("Isometric measured wind-speed overlay: %s", scenario_label),
    subtitle = sprintf("Time %g min", if (is.null(time_min)) min(measurements_dt$time_min) else time_min)
  )

  if (save) {
    fname <- sprintf(
      "plume_isometric_windspeed_measured_emission_map_%s.png",
      gsub("[^a-zA-Z0-9_]", "_", scenario_label)
    )
    ggsave_osse(file.path(get_figure_dir(), fname), p, width = 10, height = 7, dpi = 200)
    message(sprintf("  Saved: %s", fname))
  }

  invisible(p)
}


# -----------------------------------------------------------------------------
# 6f. Isometric true wind-speed view over emission map
# -----------------------------------------------------------------------------
plot_windspeed_isometric_true_over_emission_map <- function(truth_dt,
                                                            receptors,
                                                            flux_field,
                                                            lake_diameter,
                                                            time_min = NULL,
                                                            scenario_label = "default",
                                                            z_scale = 0.35,
                                                            fill_limits = NULL,
                                                            colour_limits = NULL,
                                                            nearby_lake_config = NULL,
                                                            save = TRUE) {
  t_use <- if (is.null(time_min)) min(truth_dt$time_min) else time_min
  truth_t <- truth_dt[time_min == t_use]
  if (nrow(truth_t) == 0) {
    message("  No data for true wind-speed isometric plot")
    return(invisible())
  }

  if ("ws_at_z_m_s" %in% names(truth_t)) {
    truth_t[, ws_true_plot := ws_at_z_m_s]
  } else if ("ws_m_s" %in% names(truth_t)) {
    truth_t[, ws_true_plot := ws_m_s]
  } else {
    stop("No wind-speed column found in truth_dt (`ws_at_z_m_s` or `ws_m_s`).")
  }

  # Reuse measured-isometric infrastructure with a minimal table so geometry
  # always comes from `receptors` during merge.
  truth_like_meas <- truth_t[, .(receptor_id, time_min, ws_true_plot)]
  truth_like_meas[, ch4_measured_ppb := ws_true_plot]
  truth_like_meas[, ch4_bg_ppb := 0]

  p <- plot_measured_isometric_over_emission_map(
    measurements_dt = truth_like_meas,
    receptors = receptors,
    flux_field = flux_field,
    lake_diameter = lake_diameter,
    time_min = t_use,
    scenario_label = scenario_label,
    z_scale = z_scale,
    value_mode = "custom",
    value_col = "ws_true_plot",
    value_label_override = "True Wind Speed (m s^-1)",
    fill_limits = fill_limits,
    colour_limits = colour_limits,
    nearby_lake_config = nearby_lake_config,
    save = FALSE
  )

  p <- p + labs(
    title = sprintf("Isometric true wind-speed overlay: %s", scenario_label),
    subtitle = sprintf("Time %g min", t_use)
  )

  if (save) {
    fname <- sprintf(
      "plume_isometric_windspeed_true_emission_map_%s.png",
      gsub("[^a-zA-Z0-9_]", "_", scenario_label)
    )
    ggsave_osse(file.path(get_figure_dir(), fname), p, width = 10, height = 7, dpi = 200)
    message(sprintf("  Saved: %s", fname))
  }

  invisible(p)
}


# -----------------------------------------------------------------------------
# 6g. Isometric methane difference (measured - true) over emission map
# -----------------------------------------------------------------------------
plot_methane_difference_isometric_over_emission_map <- function(measurements_dt,
                                                                truth_dt,
                                                                receptors,
                                                                flux_field,
                                                                lake_diameter,
                                                                time_min = NULL,
                                                                scenario_label = "default",
                                                                z_scale = 0.35,
                                                                fill_limits = NULL,
                                                                colour_limits = NULL,
                                                                nearby_lake_config = NULL,
                                                                save = TRUE) {
  t_use <- if (is.null(time_min)) min(measurements_dt$time_min) else time_min

  meas_cols <- c("receptor_id", "time_min", "ch4_measured_ppb")
  if ("x_measured" %in% names(measurements_dt)) meas_cols <- c(meas_cols, "x_measured")
  if ("y_measured" %in% names(measurements_dt)) meas_cols <- c(meas_cols, "y_measured")
  if ("z_measured" %in% names(measurements_dt)) meas_cols <- c(meas_cols, "z_measured")

  truth_cols <- c("receptor_id", "time_min")
  if ("ch4_total_ppb" %in% names(truth_dt)) truth_cols <- c(truth_cols, "ch4_total_ppb")
  else stop("truth_dt is missing `ch4_total_ppb` for methane-difference plotting.")

  dt <- merge(
    measurements_dt[, ..meas_cols],
    truth_dt[, ..truth_cols],
    by = c("receptor_id", "time_min"),
    all = FALSE
  )
  dt <- dt[time_min == t_use]
  if (nrow(dt) == 0L) {
    message("  No rows available for methane-difference isometric plot")
    return(invisible())
  }

  dt[, delta_ch4_ppb := ch4_measured_ppb - ch4_total_ppb]

  p <- plot_measured_isometric_over_emission_map(
    measurements_dt = dt,
    receptors = receptors,
    flux_field = flux_field,
    lake_diameter = lake_diameter,
    time_min = t_use,
    scenario_label = scenario_label,
    z_scale = z_scale,
    value_mode = "custom",
    value_col = "delta_ch4_ppb",
    value_label_override = "Methane Difference (Measured - True, ppb)",
    fill_limits = fill_limits,
    colour_limits = colour_limits,
    nearby_lake_config = nearby_lake_config,
    save = FALSE
  )

  p <- p + labs(
    title = sprintf("Isometric methane difference (measured - true): %s", scenario_label),
    subtitle = sprintf("Time %g min", t_use)
  )

  if (save) {
    fname <- sprintf(
      "plume_isometric_ch4_difference_emission_map_%s.png",
      gsub("[^a-zA-Z0-9_]", "_", scenario_label)
    )
    ggsave_osse(file.path(get_figure_dir(), fname), p, width = 10, height = 7, dpi = 200)
    message(sprintf("  Saved: %s", fname))
  }

  invisible(p)
}


# -----------------------------------------------------------------------------
# 6h. Isometric wind-speed difference (measured - true) over emission map
# -----------------------------------------------------------------------------
plot_windspeed_difference_isometric_over_emission_map <- function(measurements_dt,
                                                                  truth_dt,
                                                                  receptors,
                                                                  flux_field,
                                                                  lake_diameter,
                                                                  time_min = NULL,
                                                                  scenario_label = "default",
                                                                  z_scale = 0.35,
                                                                  ws_observation_error_m_s = 0,
                                                                  fill_limits = NULL,
                                                                  colour_limits = NULL,
                                                                  nearby_lake_config = NULL,
                                                                  save = TRUE) {
  t_use <- if (is.null(time_min)) min(measurements_dt$time_min) else time_min

  meas_ws_col <- if ("ws_at_z_m_s" %in% names(measurements_dt)) "ws_at_z_m_s" else
    if ("ws_m_s" %in% names(measurements_dt)) "ws_m_s" else NA_character_
  truth_ws_col <- if ("ws_at_z_m_s" %in% names(truth_dt)) "ws_at_z_m_s" else
    if ("ws_m_s" %in% names(truth_dt)) "ws_m_s" else NA_character_

  if (!is.finite(nchar(meas_ws_col)) || !is.finite(nchar(truth_ws_col))) {
    stop("No compatible wind-speed columns found for wind-difference plotting.")
  }

  meas_cols <- c("receptor_id", "time_min", meas_ws_col)
  if ("x_measured" %in% names(measurements_dt)) meas_cols <- c(meas_cols, "x_measured")
  if ("y_measured" %in% names(measurements_dt)) meas_cols <- c(meas_cols, "y_measured")
  if ("z_measured" %in% names(measurements_dt)) meas_cols <- c(meas_cols, "z_measured")

  truth_cols <- c("receptor_id", "time_min", truth_ws_col)

  dt <- merge(
    measurements_dt[, ..meas_cols],
    truth_dt[, ..truth_cols],
    by = c("receptor_id", "time_min"),
    all = FALSE,
    suffixes = c("_meas", "_true")
  )
  dt <- dt[time_min == t_use]
  if (nrow(dt) == 0L) {
    message("  No rows available for wind-difference isometric plot")
    return(invisible())
  }

  ws_meas_name <- paste0(meas_ws_col, "_meas")
  ws_true_name <- paste0(truth_ws_col, "_true")
  if (ws_observation_error_m_s > 0) {
    rng_state <- if (exists(".Random.seed", envir = .GlobalEnv))
      get(".Random.seed", envir = .GlobalEnv) else NULL
    set.seed(string_to_seed(scenario_label) + 7778L)
    dt[, ws_meas_for_diff := pmax(0.1, get(ws_meas_name) + rnorm(.N, 0, ws_observation_error_m_s))]
    if (!is.null(rng_state)) assign(".Random.seed", rng_state, envir = .GlobalEnv)
    ws_meas_name <- "ws_meas_for_diff"
  }
  dt[, delta_ws_m_s := get(ws_meas_name) - get(ws_true_name)]

  p <- plot_measured_isometric_over_emission_map(
    measurements_dt = dt,
    receptors = receptors,
    flux_field = flux_field,
    lake_diameter = lake_diameter,
    time_min = t_use,
    scenario_label = scenario_label,
    z_scale = z_scale,
    value_mode = "custom",
    value_col = "delta_ws_m_s",
    value_label_override = "Wind-Speed Difference (Measured - True, m s^-1)",
    fill_limits = fill_limits,
    colour_limits = colour_limits,
    nearby_lake_config = nearby_lake_config,
    save = FALSE
  )

  p <- p + labs(
    title = sprintf("Isometric wind-speed difference (measured - true): %s", scenario_label),
    subtitle = sprintf("Time %g min", t_use)
  )

  if (save) {
    fname <- sprintf(
      "plume_isometric_windspeed_difference_emission_map_%s.png",
      gsub("[^a-zA-Z0-9_]", "_", scenario_label)
    )
    ggsave_osse(file.path(get_figure_dir(), fname), p, width = 10, height = 7, dpi = 200)
    message(sprintf("  Saved: %s", fname))
  }

  invisible(p)
}


# -----------------------------------------------------------------------------
# 6i. Isometric flux-vector summary (measured)
# -----------------------------------------------------------------------------
plot_fluxvector_isometric_measured_over_emission_map <- function(measurements_dt,
                                                                 receptors,
                                                                 flux_field,
                                                                 lake_diameter,
                                                                 time_min = NULL,
                                                                 scenario_label = "default",
                                                                 z_scale = 0.35,
                                                                 fill_limits = NULL,
                                                                 colour_limits = NULL,
                                                                 nearby_lake_config = NULL,
                                                                 save = TRUE) {
  t_use <- if (is.null(time_min)) min(measurements_dt$time_min) else time_min
  dt <- measurements_dt[time_min == t_use]
  if (nrow(dt) == 0L) {
    message("  No rows available for measured flux-vector isometric plot")
    return(invisible())
  }

  required_cols <- c("ch4_measured_ppb", "ws_at_z_m_s", "wd_deg")
  miss <- setdiff(required_cols, names(dt))
  if (length(miss) > 0L) {
    stop("Missing columns for measured flux-vector plot: ",
         paste(miss, collapse = ", "))
  }

  # Match inversion sign convention:
  # v_r = ws * cos(azimuth - wd_to), where v_r > 0 is outward transport.
  # Build azimuth from receptor geometry so "inside/outside" direction is explicit.
  if (!("azimuth_deg" %in% names(dt))) {
    receptor_meta <- receptors[, .(receptor_id, x, y, azimuth_deg)]
    dt <- merge(dt, receptor_meta, by = "receptor_id", all.x = TRUE)
  }
  if (!("azimuth_deg" %in% names(dt)) || any(is.na(dt$azimuth_deg))) {
    if (!all(c("x", "y") %in% names(dt))) {
      receptor_xy <- receptors[, .(receptor_id, x, y)]
      dt <- merge(dt, receptor_xy, by = "receptor_id", all.x = TRUE)
    }
    # Keep project azimuth convention used by define_receptor_config():
    # x = r*sin(az), y = r*cos(az)  => az = atan2(x, y)
    dt[, azimuth_deg := (atan2(x, y) * 180 / pi) %% 360]
  }

  # Wind direction is FROM; convert to TO.
  dt[, wd_to := (wd_deg + 180) %% 360]
  dt[, v_r := ws_at_z_m_s * cos((azimuth_deg - wd_to) * pi / 180)]
  dt[, flux_vector_sum := ch4_measured_ppb * v_r]
  dt[, direction := ifelse(flux_vector_sum >= 0, "outside (+)", "inside (-)")]

  p <- plot_measured_isometric_over_emission_map(
    measurements_dt = dt,
    receptors = receptors,
    flux_field = flux_field,
    lake_diameter = lake_diameter,
    time_min = t_use,
    scenario_label = scenario_label,
    z_scale = z_scale,
    value_mode = "custom",
    value_col = "flux_vector_sum",
    value_label_override = "Measured Flux-Vector Sum (ppb·m s^-1)",
    fill_limits = fill_limits,
    colour_limits = colour_limits,
    nearby_lake_config = nearby_lake_config,
    save = FALSE
  )

  n_out <- dt[flux_vector_sum >= 0, .N]
  n_in <- dt[flux_vector_sum < 0, .N]
  p <- p + labs(
    title = sprintf("Isometric measured flux-vector sum: %s", scenario_label),
    subtitle = sprintf("Time %g min | outside (+): %d, inside (-): %d", t_use, n_out, n_in)
  )

  if (save) {
    fname <- sprintf(
      "plume_isometric_fluxvector_measured_emission_map_%s.png",
      gsub("[^a-zA-Z0-9_]", "_", scenario_label)
    )
    ggsave_osse(file.path(get_figure_dir(), fname), p, width = 10, height = 7, dpi = 200)
    message(sprintf("  Saved: %s", fname))
  }

  invisible(p)
}


# -----------------------------------------------------------------------------
# 6j. Isometric flux-vector summary (true)
# -----------------------------------------------------------------------------
plot_fluxvector_isometric_true_over_emission_map <- function(truth_dt,
                                                             receptors,
                                                             flux_field,
                                                             lake_diameter,
                                                             time_min = NULL,
                                                             scenario_label = "default",
                                                             z_scale = 0.35,
                                                             fill_limits = NULL,
                                                             colour_limits = NULL,
                                                             nearby_lake_config = NULL,
                                                             save = TRUE) {
  t_use <- if (is.null(time_min)) min(truth_dt$time_min) else time_min
  dt <- truth_dt[time_min == t_use]
  if (nrow(dt) == 0L) {
    message("  No rows available for true flux-vector isometric plot")
    return(invisible())
  }

  required_cols <- c("ch4_total_ppb", "ws_at_z_m_s", "wd_deg")
  miss <- setdiff(required_cols, names(dt))
  if (length(miss) > 0L) {
    stop("Missing columns for true flux-vector plot: ",
         paste(miss, collapse = ", "))
  }

  if (!("azimuth_deg" %in% names(dt))) {
    receptor_meta <- receptors[, .(receptor_id, x, y, azimuth_deg)]
    dt <- merge(dt, receptor_meta, by = "receptor_id", all.x = TRUE)
  }
  if (!("azimuth_deg" %in% names(dt)) || any(is.na(dt$azimuth_deg))) {
    if (!all(c("x", "y") %in% names(dt))) {
      receptor_xy <- receptors[, .(receptor_id, x, y)]
      dt <- merge(dt, receptor_xy, by = "receptor_id", all.x = TRUE)
    }
    # Keep project azimuth convention used by define_receptor_config():
    # x = r*sin(az), y = r*cos(az)  => az = atan2(x, y)
    dt[, azimuth_deg := (atan2(x, y) * 180 / pi) %% 360]
  }

  dt[, wd_to := (wd_deg + 180) %% 360]
  dt[, v_r := ws_at_z_m_s * cos((azimuth_deg - wd_to) * pi / 180)]
  dt[, flux_vector_sum := ch4_total_ppb * v_r]
  dt[, direction := ifelse(flux_vector_sum >= 0, "outside (+)", "inside (-)")]

  # Reuse measured-isometric infrastructure with minimal columns so geometry
  # comes from `receptors`.
  truth_like_meas <- dt[, .(receptor_id, time_min, flux_vector_sum)]
  truth_like_meas[, `:=`(ch4_measured_ppb = flux_vector_sum, ch4_bg_ppb = 0)]

  p <- plot_measured_isometric_over_emission_map(
    measurements_dt = truth_like_meas,
    receptors = receptors,
    flux_field = flux_field,
    lake_diameter = lake_diameter,
    time_min = t_use,
    scenario_label = scenario_label,
    z_scale = z_scale,
    value_mode = "custom",
    value_col = "flux_vector_sum",
    value_label_override = "True Flux-Vector Sum (ppb·m s^-1)",
    fill_limits = fill_limits,
    colour_limits = colour_limits,
    nearby_lake_config = nearby_lake_config,
    save = FALSE
  )

  n_out <- dt[flux_vector_sum >= 0, .N]
  n_in <- dt[flux_vector_sum < 0, .N]
  p <- p + labs(
    title = sprintf("Isometric true flux-vector sum: %s", scenario_label),
    subtitle = sprintf("Time %g min | outside (+): %d, inside (-): %d", t_use, n_out, n_in)
  )

  if (save) {
    fname <- sprintf(
      "plume_isometric_fluxvector_true_emission_map_%s.png",
      gsub("[^a-zA-Z0-9_]", "_", scenario_label)
    )
    ggsave_osse(file.path(get_figure_dir(), fname), p, width = 10, height = 7, dpi = 200)
    message(sprintf("  Saved: %s", fname))
  }

  invisible(p)
}


# -----------------------------------------------------------------------------
# 6k. Isometric contamination flux-vector (multi-lake only)
# -----------------------------------------------------------------------------
plot_contamination_fv_isometric_over_emission_map <- function(truth_dt,
                                                               receptors,
                                                               flux_field,
                                                               lake_diameter,
                                                               time_min = NULL,
                                                               scenario_label = "default",
                                                               z_scale = 0.35,
                                                               fill_limits = NULL,
                                                               colour_limits = NULL,
                                                               nearby_lake_config = NULL,
                                                               save = TRUE) {
  if (!("ch4_contamination_ppb" %in% names(truth_dt)) ||
      all(truth_dt$ch4_contamination_ppb == 0, na.rm = TRUE)) {
    return(invisible(NULL))
  }

  t_use <- if (is.null(time_min)) min(truth_dt$time_min) else time_min
  dt <- truth_dt[time_min == t_use]
  if (nrow(dt) == 0L) {
    message("  No rows for contamination flux-vector isometric plot")
    return(invisible())
  }

  if (!("azimuth_deg" %in% names(dt))) {
    dt <- merge(dt, receptors[, .(receptor_id, x, y, azimuth_deg)],
                by = "receptor_id", all.x = TRUE)
  }
  if (!("azimuth_deg" %in% names(dt)) || any(is.na(dt$azimuth_deg))) {
    if (!all(c("x", "y") %in% names(dt))) {
      dt <- merge(dt, receptors[, .(receptor_id, x, y)],
                  by = "receptor_id", all.x = TRUE)
    }
    dt[, azimuth_deg := (atan2(x, y) * 180 / pi) %% 360]
  }

  dt[, wd_to := (wd_deg + 180) %% 360]
  dt[, v_r := ws_at_z_m_s * cos((azimuth_deg - wd_to) * pi / 180)]
  dt[, contam_fv := ch4_contamination_ppb * v_r]

  truth_like_meas <- dt[, .(receptor_id, time_min, contam_fv)]
  truth_like_meas[, `:=`(ch4_measured_ppb = contam_fv, ch4_bg_ppb = 0)]

  p <- plot_measured_isometric_over_emission_map(
    measurements_dt = truth_like_meas,
    receptors = receptors,
    flux_field = flux_field,
    lake_diameter = lake_diameter,
    time_min = t_use,
    scenario_label = scenario_label,
    z_scale = z_scale,
    value_mode = "custom",
    value_col = "contam_fv",
    value_label_override = expression("Contamination Flux-Vector (ppb" %.% "m s"^-1 * ")"),
    fill_limits = fill_limits,
    colour_limits = colour_limits,
    nearby_lake_config = nearby_lake_config,
    save = FALSE
  )

  # Override with diverging red-white-blue scale (white = 0)
  suppressWarnings(p <- p + scale_colour_gradient2(
    low      = "#2166ac",
    mid      = "white",
    high     = "#d6604d",
    midpoint = 0,
    name     = expression("Contam. FV (ppb" %.% "m s"^-1 * ")"),
    limits   = colour_limits,
    oob      = scales::squish
  ))

  n_pos <- dt[contam_fv > 0, .N]
  n_neg <- dt[contam_fv < 0, .N]
  net   <- round(sum(dt$contam_fv), 3)
  p <- p + labs(
    title = sprintf("Isometric contamination flux-vector: %s", scenario_label),
    subtitle = sprintf("Time %g min | outward (+): %d, inward (-): %d | net: %s",
                       t_use, n_pos, n_neg, net)
  )

  if (save) {
    fname <- sprintf(
      "iso_contam_fv_%s.png",
      gsub("[^a-zA-Z0-9_]", "_", scenario_label)
    )
    ggsave_osse(file.path(get_figure_dir(), fname), p, width = 10, height = 7, dpi = 200)
    message(sprintf("  Saved: %s", fname))
  }

  invisible(p)
}


# -----------------------------------------------------------------------------
# 6l. Contamination flux-vector time series (multi-lake only)
# -----------------------------------------------------------------------------
.topdown_circle_path <- function(center_x, center_y, radius_m, n = 240L) {
  theta <- seq(0, 2 * pi, length.out = n)
  data.table(
    x = center_x + radius_m * cos(theta),
    y = center_y + radius_m * sin(theta)
  )
}

.rotate_topdown_xy <- function(x, y, angle_deg = 0) {
  if (is.null(angle_deg) || !is.finite(angle_deg) || abs(angle_deg) < 1e-12) {
    return(list(x = x, y = y))
  }
  angle_rad <- angle_deg * pi / 180
  list(
    x = x * cos(angle_rad) - y * sin(angle_rad),
    y = x * sin(angle_rad) + y * cos(angle_rad)
  )
}

.rotate_topdown_dt <- function(dt, angle_deg = 0) {
  if (is.null(dt) || !nrow(dt) || is.null(angle_deg) ||
      !is.finite(angle_deg) || abs(angle_deg) < 1e-12) {
    return(dt)
  }
  out <- copy(dt)
  rot <- .rotate_topdown_xy(out$x, out$y, angle_deg)
  out[, `:=`(x = rot$x, y = rot$y)]
  out
}

.smooth_topdown_field <- function(grid_dt,
                                  value_col = "contamination_ppb",
                                  kernel = matrix(c(1, 2, 1,
                                                    2, 4, 2,
                                                    1, 2, 1), nrow = 3, byrow = TRUE),
                                  passes = 2L) {
  x_vals <- sort(unique(grid_dt$x))
  y_vals <- sort(unique(grid_dt$y))
  if (!length(x_vals) || !length(y_vals)) return(grid_dt)

  mat <- matrix(grid_dt[[value_col]], nrow = length(x_vals), ncol = length(y_vals), byrow = TRUE)
  kernel <- kernel / sum(kernel)
  kx <- nrow(kernel)
  ky <- ncol(kernel)
  hx <- (kx - 1L) %/% 2L
  hy <- (ky - 1L) %/% 2L

  smooth_once <- function(m) {
    out <- matrix(NA_real_, nrow = nrow(m), ncol = ncol(m))
    for (ix in seq_len(nrow(m))) {
      x_idx <- seq.int(max(1L, ix - hx), min(nrow(m), ix + hx))
      kx_idx <- x_idx - ix + hx + 1L
      for (iy in seq_len(ncol(m))) {
        y_idx <- seq.int(max(1L, iy - hy), min(ncol(m), iy + hy))
        ky_idx <- y_idx - iy + hy + 1L
        vals <- m[x_idx, y_idx, drop = FALSE]
        wts <- kernel[kx_idx, ky_idx, drop = FALSE]
        valid <- is.finite(vals)
        if (any(valid)) {
          out[ix, iy] <- sum(vals[valid] * wts[valid]) / sum(wts[valid])
        }
      }
    }
    out
  }

  for (i in seq_len(as.integer(passes))) {
    mat <- smooth_once(mat)
  }

  grid_dt[, (value_col) := as.vector(t(mat))]
  grid_dt
}

.contamination_topdown_grid <- function(lake_diameter,
                                        nearby_lake_config,
                                        met_cfg,
                                        z_m,
                                        ring_radius_m = lake_diameter / 2 + 10,
                                        pad_m = 80,
                                        grid_n = 220L,
                                        n_rings = 48L,
                                        n_src_az = 144L,
                                        smooth_passes = 2L,
                                        display_rotate_deg = 0) {
  if (is.null(nearby_lake_config)) {
    stop("nearby_lake_config is required for top-down contamination plan views.")
  }

  nearby_center_plot <- .rotate_topdown_xy(
    nearby_lake_config$center_x,
    nearby_lake_config$center_y,
    display_rotate_deg
  )

  x_min <- min(-ring_radius_m,
               -lake_diameter / 2,
               nearby_center_plot$x - nearby_lake_config$radius_m) - pad_m
  x_max <- max(ring_radius_m,
               lake_diameter / 2,
               nearby_center_plot$x + nearby_lake_config$radius_m) + pad_m
  y_min <- min(-ring_radius_m,
               -lake_diameter / 2,
               nearby_center_plot$y - nearby_lake_config$radius_m) - pad_m
  y_max <- max(ring_radius_m,
               lake_diameter / 2,
               nearby_center_plot$y + nearby_lake_config$radius_m) + pad_m

  grid_dt <- CJ(
    x = seq(x_min, x_max, length.out = as.integer(grid_n)),
    y = seq(y_min, y_max, length.out = as.integer(grid_n))
  )
  grid_dt[, `:=`(z = z_m, receptor_id = seq_len(.N))]
  inv_rot <- .rotate_topdown_xy(grid_dt$x, grid_dt$y, -display_rotate_deg)
  grid_eval <- grid_dt[, .(receptor_id, z)]
  grid_eval[, `:=`(x = inv_rot$x, y = inv_rot$y)]

  contamination_ug_m3 <- gaussian_plume_field(
    receptors = grid_eval,
    source_radius = nearby_lake_config$radius_m,
    flux_ug_m2_s = flux_mg_d_to_ug_s(nearby_lake_config$flux_mg_m2_d),
    U_ref = met_cfg$U_ref_m_s,
    wd_deg = met_cfg$wd_base_deg,
    pasquill_class = met_cfg$pasquill_class,
    n_rings = n_rings,
    n_src_az = n_src_az,
    source_center = c(nearby_lake_config$center_x, nearby_lake_config$center_y)
  )
  grid_dt[, contamination_ppb := ug_m3_to_ppb(contamination_ug_m3)]
  grid_dt <- .smooth_topdown_field(
    grid_dt = grid_dt,
    value_col = "contamination_ppb",
    passes = smooth_passes
  )
  grid_dt[]
}

.build_contamination_topdown_base_plot <- function(grid_dt,
                                                   lake_diameter,
                                                   nearby_lake_config,
                                                   ring_radius_m = lake_diameter / 2 + 10,
                                                   title = NULL,
                                                   subtitle = NULL,
                                                   show_fill_legend = TRUE,
                                                   display_rotate_deg = 0,
                                                   base_size = 11) {
  is_rotated_display <- is.finite(display_rotate_deg) && abs(display_rotate_deg) > 1e-12
  main_lake_path <- .rotate_topdown_dt(
    .topdown_circle_path(0, 0, lake_diameter / 2),
    display_rotate_deg
  )
  nearby_lake_path <- .rotate_topdown_dt(.topdown_circle_path(
    nearby_lake_config$center_x,
    nearby_lake_config$center_y,
    nearby_lake_config$radius_m
  ), display_rotate_deg)
  ring_path <- .rotate_topdown_dt(
    .topdown_circle_path(0, 0, ring_radius_m),
    display_rotate_deg
  )
  fill_max <- max(grid_dt$contamination_ppb, na.rm = TRUE)
  if (!is.finite(fill_max) || fill_max <= 0) fill_max <- 1

  ggplot(grid_dt, aes(x = x, y = y)) +
    geom_raster(aes(fill = contamination_ppb), interpolate = TRUE) +
    geom_path(
      data = ring_path,
      aes(x = x, y = y),
      inherit.aes = FALSE,
      colour = "grey35",
      linewidth = 0.55,
      linetype = "dotted"
    ) +
    geom_path(
      data = main_lake_path,
      aes(x = x, y = y),
      inherit.aes = FALSE,
      colour = "grey15",
      linewidth = 0.75,
      linetype = "dashed"
    ) +
    geom_path(
      data = nearby_lake_path,
      aes(x = x, y = y),
      inherit.aes = FALSE,
      colour = "grey15",
      linewidth = 0.75,
      linetype = "dashed"
    ) +
    coord_equal() +
    scale_fill_gradientn(
      colours = c("white", "#fee8c8", "#fdbb84", "#e34a33", "#7f0000"),
      limits = c(0, fill_max),
      oob = scales::squish,
      na.value = "white",
      name = "Contam. (ppb)"
    ) +
    labs(
      title = title,
      subtitle = subtitle,
      x = if (is_rotated_display) "Along-lake separation (m)" else "East-West position (m)",
      y = if (is_rotated_display) "Cross-separation position (m)" else "North-South position (m)"
    ) +
    theme_osse(base_size = base_size) +
    theme(
      panel.grid = element_blank(),
      legend.position = if (show_fill_legend) "right" else "none",
      legend.margin = margin(0, 0, 0, 0),
      legend.box.margin = margin(0, 0, 0, 0)
    )
}

.contamination_fv_receptor_dt <- function(truth_dt,
                                          receptors,
                                          z_m,
                                          time_min = NULL) {
  t_use <- if (is.null(time_min)) min(truth_dt$time_min) else time_min
  dt <- copy(truth_dt[time_min == t_use])
  dt <- merge(
    dt,
    receptors[, .(receptor_id, azimuth_deg, ring_type)],
    by = "receptor_id",
    all.x = TRUE
  )
  dt <- dt[ring_type == "outer" & abs(z - z_m) < 1e-8]
  if (!nrow(dt)) return(dt)

  if (!("azimuth_deg" %in% names(dt)) || any(is.na(dt$azimuth_deg))) {
    dt[, azimuth_deg := (atan2(x, y) * 180 / pi) %% 360]
  }

  dt[, wd_to := (wd_deg + 180) %% 360]
  dt[, v_r := ws_at_z_m_s * cos((azimuth_deg - wd_to) * pi / 180)]
  dt[, contam_fv := ch4_contamination_ppb * v_r]
  dt
}

#' Top-down plan view of nearby-lake contamination at a fixed height.
plot_contamination_topdown_planview <- function(lake_diameter,
                                                nearby_lake_config,
                                                met_cfg,
                                                z_m,
                                                ring_radius_m = lake_diameter / 2 + 10,
                                                pad_m = 80,
                                                grid_n = 220L,
                                                n_rings = 48L,
                                                n_src_az = 144L,
                                                smooth_passes = 2L,
                                                title = NULL,
                                                subtitle = NULL,
                                                show_fill_legend = TRUE,
                                                display_rotate_deg = 0,
                                                base_size = 11) {
  grid_dt <- .contamination_topdown_grid(
    lake_diameter = lake_diameter,
    nearby_lake_config = nearby_lake_config,
    met_cfg = met_cfg,
    z_m = z_m,
    ring_radius_m = ring_radius_m,
    pad_m = pad_m,
    grid_n = grid_n,
    n_rings = n_rings,
    n_src_az = n_src_az,
    smooth_passes = smooth_passes,
    display_rotate_deg = display_rotate_deg
  )

  .build_contamination_topdown_base_plot(
    grid_dt = grid_dt,
    lake_diameter = lake_diameter,
    nearby_lake_config = nearby_lake_config,
    ring_radius_m = ring_radius_m,
    title = title,
    subtitle = subtitle,
    show_fill_legend = show_fill_legend,
    display_rotate_deg = display_rotate_deg,
    base_size = base_size
  )
}

#' Top-down plan view of contamination flux-vector contribution at receptors.
plot_contamination_fv_topdown_planview <- function(truth_dt,
                                                   receptors,
                                                   lake_diameter,
                                                   nearby_lake_config,
                                                   met_cfg,
                                                   z_m,
                                                   time_min = NULL,
                                                   ring_radius_m = lake_diameter / 2 + 10,
                                                   pad_m = 80,
                                                   grid_n = 220L,
                                                   n_rings = 48L,
                                                   n_src_az = 144L,
                                                   smooth_passes = 2L,
                                                   title = NULL,
                                                   subtitle = NULL,
                                                   show_fill_legend = FALSE,
                                                   display_rotate_deg = 0,
                                                   receptor_point_size = 2.2,
                                                   base_size = 11) {
  grid_dt <- .contamination_topdown_grid(
    lake_diameter = lake_diameter,
    nearby_lake_config = nearby_lake_config,
    met_cfg = met_cfg,
    z_m = z_m,
    ring_radius_m = ring_radius_m,
    pad_m = pad_m,
    grid_n = grid_n,
    n_rings = n_rings,
    n_src_az = n_src_az,
    smooth_passes = smooth_passes,
    display_rotate_deg = display_rotate_deg
  )
  receptor_dt <- .contamination_fv_receptor_dt(
    truth_dt = truth_dt,
    receptors = receptors,
    z_m = z_m,
    time_min = time_min
  )
  max_abs_fv <- max(abs(receptor_dt$contam_fv), na.rm = TRUE)
  if (!is.finite(max_abs_fv) || max_abs_fv <= 0) max_abs_fv <- 1
  receptor_dt <- .rotate_topdown_dt(receptor_dt, display_rotate_deg)

  .build_contamination_topdown_base_plot(
    grid_dt = grid_dt,
    lake_diameter = lake_diameter,
    nearby_lake_config = nearby_lake_config,
    ring_radius_m = ring_radius_m,
    title = title,
    subtitle = subtitle,
    show_fill_legend = show_fill_legend,
    display_rotate_deg = display_rotate_deg,
    base_size = base_size
  ) +
    geom_point(
      data = receptor_dt,
      aes(x = x, y = y, colour = contam_fv),
      inherit.aes = FALSE,
      size = receptor_point_size
    ) +
    scale_colour_gradient2(
      low = "#2166ac",
      mid = "white",
      high = "#d6604d",
      midpoint = 0,
      limits = c(-max_abs_fv, max_abs_fv),
      oob = scales::squish,
      name = expression("Contam. FV (ppb" %.% "m s"^-1 * ")")
    ) +
    guides(
      fill = if (show_fill_legend) guide_colorbar(order = 1) else "none",
      colour = guide_colorbar(order = 1 + as.integer(show_fill_legend))
    ) +
    theme(legend.position = "right")
}

plot_contamination_fv_timeseries <- function(truth_dt,
                                             receptors,
                                             scenario_label = "default",
                                             save = TRUE) {
  if (!("ch4_contamination_ppb" %in% names(truth_dt)) ||
      all(truth_dt$ch4_contamination_ppb == 0, na.rm = TRUE)) {
    return(invisible(NULL))
  }

  dt <- copy(truth_dt)
  if (!("azimuth_deg" %in% names(dt))) {
    dt <- merge(dt, receptors[, .(receptor_id, x, y, azimuth_deg)],
                by = "receptor_id", all.x = TRUE)
  }
  if (!("azimuth_deg" %in% names(dt)) || any(is.na(dt$azimuth_deg))) {
    if (!all(c("x", "y") %in% names(dt))) {
      dt <- merge(dt, receptors[, .(receptor_id, x, y)],
                  by = "receptor_id", all.x = TRUE)
    }
    dt[, azimuth_deg := (atan2(x, y) * 180 / pi) %% 360]
  }

  dt[, wd_to := (wd_deg + 180) %% 360]
  dt[, v_r := ws_at_z_m_s * cos((azimuth_deg - wd_to) * pi / 180)]
  dt[, contam_fv := ch4_contamination_ppb * v_r]
  dt[, total_fv  := ch4_total_ppb * v_r]

  agg <- dt[, .(
    net_contam_fv  = sum(contam_fv),
    net_total_fv   = sum(total_fv),
    contam_inflow  = sum(contam_fv[contam_fv < 0]),
    contam_outflow = sum(contam_fv[contam_fv > 0])
  ), by = time_min]
  agg[, contam_pct := ifelse(abs(net_total_fv) > 1e-9,
                             100 * net_contam_fv / net_total_fv,
                             NA_real_)]
  agg[, exit_pct := ifelse(abs(contam_inflow) > 1e-12,
                           100 * contam_outflow / abs(contam_inflow),
                           NA_real_)]
  agg[, hour_of_day := (time_min / 60) %% 24]

  p1 <- ggplot(agg, aes(x = hour_of_day, y = net_contam_fv)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_line(colour = "firebrick", linewidth = 0.8) +
    geom_point(colour = "firebrick", size = 2) +
    labs(
      x = "Hour of day",
      y = expression("Net contamination flux-vector (ppb" %.% "m s"^-1 * ")")
    ) +
    theme_osse()

  p2 <- ggplot(agg, aes(x = hour_of_day, y = contam_pct)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_line(colour = "steelblue", linewidth = 0.8) +
    geom_point(colour = "steelblue", size = 2) +
    labs(
      x = "Hour of day",
      y = "Contamination bias (%)"
    ) +
    theme_osse()

  p3 <- ggplot(agg, aes(x = hour_of_day, y = exit_pct)) +
    geom_hline(yintercept = 100, linetype = "dashed", colour = "grey50") +
    geom_line(colour = "darkgreen", linewidth = 0.8) +
    geom_point(colour = "darkgreen", size = 2) +
    scale_y_continuous(limits = c(0, NA)) +
    labs(
      x = "Hour of day",
      y = "Lateral exit fraction (%)"
    ) +
    theme_osse()

  if (cowplot_available) {
    p <- cowplot::plot_grid(p1, p2, p3, ncol = 1, align = "v",
                           labels = c("A", "B", "C"))
    title_gg <- cowplot::ggdraw() +
      cowplot::draw_label(
        sprintf("Contamination flux-vector diagnostics: %s", scenario_label),
        fontface = "bold", size = 13
      )
    p <- cowplot::plot_grid(title_gg, p, ncol = 1, rel_heights = c(0.05, 1))
  } else {
    p <- p1 + labs(
      title = sprintf("Contamination flux-vector: %s", scenario_label)
    )
  }

  if (save) {
    fname <- sprintf(
      "contam_fv_timeseries_%s.png",
      gsub("[^a-zA-Z0-9_]", "_", scenario_label)
    )
    ggsave_osse(file.path(get_figure_dir(), fname), p, width = 9, height = 7, dpi = 200)
    message(sprintf("  Saved: %s", fname))
  }

  invisible(p)
}


# -----------------------------------------------------------------------------
# 6. True vs retrieved flux time series
# -----------------------------------------------------------------------------
#' Plot time series of true vs retrieved whole-lake flux.
#'
#' @param flux_dt  data.table from retrieve_flux().
#' @param scenario_label  Character label.
#' @param save  Logical.
plot_flux_timeseries <- function(flux_dt, scenario_label = "default", save = TRUE) {
  dt <- flux_dt[qc_flag == "ok" & is.finite(flux_retrieved)]
  if (nrow(dt) == 0) { message("  No valid data for time series plot"); return(invisible()) }

  dt[, hour_of_day := (time_min / 60) %% 24]

  p <- ggplot(dt, aes(x = hour_of_day)) +
    geom_ribbon(aes(ymin = flux_true * 0.75, ymax = flux_true * 1.25),
                fill = "grey85", alpha = 0.6) +
    geom_line(aes(y = flux_true, colour = "True"), linewidth = 1) +
    geom_point(aes(y = flux_retrieved, colour = "Retrieved"),
               size = 2, alpha = 0.7) +
    geom_line(aes(y = flux_retrieved, colour = "Retrieved"),
              linewidth = 0.6, linetype = "dashed") +
    scale_colour_manual(
      values = c("True" = "black", "Retrieved" = viridis::viridis(5)[3]),
      name   = NULL
    ) +
    labs(
      title = paste("True vs retrieved flux:", scenario_label),
      subtitle = "Shaded band = +/- 25% of true flux",
      x     = "Hour of day",
      y     = expression(paste("Flux (mg CH"[4], " m"^{-2}, " d"^{-1}, ")"))
    ) +
    theme_osse()

  if (save) {
    fname <- sprintf("flux_timeseries_%s.png",
                     gsub("[^a-zA-Z0-9_]", "_", scenario_label))
    ggsave_osse(file.path(get_figure_dir(), fname), p, width = 8, height = 5, dpi = 200)
    message(sprintf("  Saved: %s", fname))
  }
  invisible(p)
}


# -----------------------------------------------------------------------------
# 7. Sensitivity tornado chart
# -----------------------------------------------------------------------------
#' Plot sensitivity indices as a horizontal bar (tornado) chart.
#'
#' @param si_dt  data.table from compute_sensitivity_indices().
#' @param lake_name  Lake name for reading from file if si_dt is NULL.
#' @param processed_dir  Optional: directory (relative to project root) for CSV
#'   lookup; e.g. \code{"data/processed/run_20240301"} to use a previous run's data.
#' @param save  Logical.
plot_sensitivity_tornado <- function(si_dt = NULL, lake_name = "small",
                                     processed_dir = NULL, save = TRUE) {
  if (is.null(si_dt)) {
    data_dir <- if (is.null(processed_dir)) DATA_PROC else file.path(project_root(), processed_dir)
    fpath <- file.path(data_dir, sprintf("sensitivity_indices_%s.csv", lake_name))
    if (!file.exists(fpath)) {
      fpath <- file.path(data_dir, "sensitivity_indices.csv")
    }
    if (!file.exists(fpath)) { message("  No sensitivity data"); return(invisible()) }
    si_dt <- fread(fpath)
  }

  if (nrow(si_dt) == 0) return(invisible())

  si_dt[, parameter := factor(parameter, levels = rev(parameter))]

  p <- ggplot(si_dt, aes(x = parameter, y = Si, fill = Si)) +
    geom_col(width = 0.7) +
    coord_flip() +
    scale_fill_viridis_c(option = "D", direction = -1) +
    labs(
      title = "First-order sensitivity indices (variance-based screening)",
      x     = NULL,
      y     = "First-order index (Si)",
      fill  = "Si"
    ) +
    theme_osse() +
    theme(legend.position = "right")

  if (save) {
    ggsave_osse(file.path(get_figure_dir(), "sensitivity_tornado.png"), p,
           width = 8, height = 5, dpi = 200)
    message("  Saved: sensitivity_tornado.png")
  }
  invisible(p)
}


# -----------------------------------------------------------------------------
# 7b. Sobol sensitivity (first-order and total-effect)
# -----------------------------------------------------------------------------
#' Plot Sobol first-order (Si) and total-effect (STi) indices.
#'
#' Reads \code{sensitivity_sobol_fast_<lake>.csv} from \code{run_fast_ensemble_sobol()}
#' and plots Si and STi as grouped horizontal bars. STi - Si indicates interaction share.
#'
#' @param sobol_dt  data.table from run_fast_ensemble_sobol() (columns: parameter, Si, STi).
#' @param lake_name  Lake name for reading from file if sobol_dt is NULL.
#' @param processed_dir  Optional: directory (relative to project root) for CSV
#'   lookup; e.g. \code{"data/processed/run_20240301"} to use a previous run's data.
#' @param save  Logical. If TRUE, saves to figures/sobol_sensitivity_<lake>.png.
#' @return The ggplot object (invisibly), or invisible() if no data.
plot_sobol_indices <- function(sobol_dt = NULL, lake_name = "small",
                               processed_dir = NULL, save = TRUE) {
  if (is.null(sobol_dt)) {
    data_dir <- if (is.null(processed_dir)) DATA_PROC else file.path(project_root(), processed_dir)
    fpath <- file.path(data_dir, sprintf("sensitivity_sobol_fast_%s.csv", lake_name))
    if (!file.exists(fpath)) {
      message("  No Sobol data at ", fpath)
      return(invisible())
    }
    sobol_dt <- fread(fpath)
  }

  if (nrow(sobol_dt) == 0L || !all(c("parameter", "Si", "STi") %in% names(sobol_dt))) {
    message("  Sobol data missing required columns (parameter, Si, STi)")
    return(invisible())
  }

  sobol_dt <- sobol_dt[is.finite(Si) & is.finite(STi)]
  if (nrow(sobol_dt) == 0L) return(invisible())

  # Order parameters by total-effect (STi) descending
  sobol_dt[, parameter := factor(parameter, levels = parameter[order(-STi)])]

  # Long format for grouped bars
  long <- melt(sobol_dt, id.vars = "parameter",
               measure.vars = c("Si", "STi"),
               variable.name = "index", value.name = "value")
  long[, index := factor(index, levels = c("STi", "Si"))]

  p <- ggplot(long, aes(x = parameter, y = value, fill = index)) +
    geom_col(position = position_dodge(width = 0.75), width = 0.7, alpha = 0.9) +
    coord_flip() +
    scale_fill_manual(
      values = c("Si" = "#21918c", "STi" = "#440154"),
      labels = c("Si" = "First-order (Si)", "STi" = "Total-effect (STi)")
    ) +
    labs(
      title = "Sobol sensitivity indices (fast stage)",
      subtitle = sprintf("Lake: %s; response: rmse_pct", lake_name),
      x = NULL,
      y = "Sobol index",
      fill = NULL
    ) +
    theme_osse() +
    theme(legend.position = "bottom")

  if (save) {
    dir.create(get_figure_dir(), showWarnings = FALSE, recursive = TRUE)
    outfile <- file.path(get_figure_dir(), sprintf("sobol_sensitivity_%s.png", lake_name))
    ggsave_osse(outfile, p, width = 8, height = 5, dpi = 200)
    message("  Saved: ", basename(outfile))
  }
  invisible(p)
}


# -----------------------------------------------------------------------------
# 8. Save summary tables
# -----------------------------------------------------------------------------
#' Save capability metrics CSV.
save_capability_metrics <- function(all_metrics, outdir = NULL) {
  if (is.null(outdir)) outdir <- DATA_PROC
  if (is.list(all_metrics) && !is.data.table(all_metrics)) {
    all_metrics <- rbindlist(all_metrics, fill = TRUE)
  }
  fwrite(all_metrics, file.path(outdir, "capability_metrics.csv"))
  message("  Saved: capability_metrics.csv")
}

#' Save scenarios summary CSV.
save_scenarios_summary <- function(all_metrics, outdir = NULL) {
  if (is.null(outdir)) outdir <- DATA_PROC
  if (is.list(all_metrics) && !is.data.table(all_metrics)) {
    all_metrics <- rbindlist(all_metrics, fill = TRUE)
  }

  dt <- copy(all_metrics)
  dt[, c("lake", "sensor", "design", "met") :=
       as.list(parse_scenario_key(scenario)), by = scenario]

  summary_dt <- dt[, .(lake, sensor, design, met, rmse_pct, bias_pct,
                        correlation, valid_fraction, meets_target)]
  fwrite(summary_dt, file.path(outdir, "scenarios_summary.csv"))
  message("  Saved: scenarios_summary.csv")
}


# -----------------------------------------------------------------------------
# 9. Master plotting function
# -----------------------------------------------------------------------------
#' Generate all figures and summary tables from processed data.
#'
#' @param all_metrics  List or data.table of metrics from all OSSE scenarios.
#' @param emission_results  Named list from generate_all_emissions().
#' @param flux_results  Named list of flux data.tables from retrieve_flux().
#' @param processed_dir  Optional: directory (relative to project root) for
#'   sensitivity/budget CSVs; e.g. \code{"data/processed/run_20240301"} to plot
#'   from a previous run instead of the default \code{data/processed}.
#' @param sensitivity_lake  Lake name for sensitivity and budget plots (default \code{"small"}).
generate_all_figures <- function(all_metrics = NULL,
                                 emission_results = NULL,
                                 flux_results = NULL,
                                 processed_dir = NULL,
                                 sensitivity_lake = "small") {

  message("\n=== Generating figures ===")
  dir.create(get_figure_dir(), showWarnings = FALSE, recursive = TRUE)

  if (!is.null(all_metrics)) {
    if (is.list(all_metrics) && !is.data.table(all_metrics)) {
      metrics_dt <- rbindlist(all_metrics, fill = TRUE)
    } else {
      metrics_dt <- copy(all_metrics)
    }

    if (!"lake" %in% names(metrics_dt)) {
      metrics_dt[, c("lake", "sensor", "design", "met") :=
                   as.list(parse_scenario_key(scenario)), by = scenario]
    }

    plot_flux_rmse_vs_sensor(metrics_dt)
    plot_receptor_design(metrics_dt)
    save_capability_metrics(metrics_dt)
    save_scenarios_summary(metrics_dt)
  }

  if (!is.null(emission_results)) {
    for (lk_name in names(emission_results)) {
      em <- emission_results[[lk_name]]
      plot_emission_field(em$flux_field, em$lake_mean, lk_name)
    }
    # Measurement locations (one height layer) for first lake with first design
    first_lk <- names(emission_results)[1]
    first_design <- SAMPLING_CONFIGS[[1]]
    diam <- emission_results[[first_lk]]$diameter_m
    receptors <- define_receptor_config(
      lake_diameter         = diam,
      n_azimuth             = first_design$n_azimuth,
      heights_m             = first_design$heights_m,
      distance_from_edge_m  = 10,
      inner_ring_fracs      = first_design$inner_ring_fracs,
      background_heights_m  = first_design$background_heights_m
    )
    plot_measurement_locations(receptors, diam, first_lk)
  }

  # Individual flux time series are plotted selectively by the caller to avoid
  # generating one PNG per scenario.

  if (!is.null(flux_results)) {
    combined_flux <- rbindlist(flux_results, fill = TRUE)
    if (nrow(combined_flux) > 0) plot_flux_vs_stability(combined_flux)
  }

  plot_uncertainty_budget(lake_name = sensitivity_lake, processed_dir = processed_dir)
  plot_sensitivity_tornado(lake_name = sensitivity_lake, processed_dir = processed_dir)
  plot_sobol_indices(lake_name = sensitivity_lake, processed_dir = processed_dir)

  message("=== Figure generation complete ===")
}


# -----------------------------------------------------------------------------
# Standalone execution
# -----------------------------------------------------------------------------
if (sys.nframe() == 0) {
  proj_root <- project_root()
  message("=== 06: Generating visualisation (from saved data) ===")

  met_file <- file.path(DATA_PROC, "capability_metrics.csv")
  if (file.exists(met_file)) {
    metrics_dt <- fread(met_file)
    if (!"lake" %in% names(metrics_dt)) {
      metrics_dt[, c("lake", "sensor", "design", "met") :=
                   as.list(parse_scenario_key(scenario)), by = scenario]
    }
    plot_flux_rmse_vs_sensor(metrics_dt)
    plot_receptor_design(metrics_dt)
  }
  plot_uncertainty_budget()
  plot_sensitivity_tornado()
  plot_sobol_indices(lake_name = "small")
  message("Done.")
}


# -----------------------------------------------------------------------------
# 8. Plume transect heatmap (x-z cross-section visualization)
# -----------------------------------------------------------------------------
#' Heatmap of CH4 enhancement in the downwind x-z cross-section.
#'
#' Plots concentration at transect receptors (y = 0, varying x and z) as a
#' filled tile plot for a single time step or an animated GIF over all steps.
#'
#' @param truth_dt          data.table from run_forward_model / build_truth_from_ce,
#'                          filtered to transect receptors.
#' @param transect_receptors data.table from define_transect_receptors().
#' @param scenario_label    Character label for titles and filenames.
#' @param time_min          Single time step to plot; NULL = animate all steps.
#' @param fig_dir           Directory for output files (default: get_figure_dir()).
#' @param fps               Frames per second for animation (default 4).
#' @param save              Logical; save to disk?
#' @return Invisibly: the last ggplot (single frame) or NULL (animation).
plot_plume_transect_heatmap <- function(truth_dt,
                                        transect_receptors,
                                        scenario_label = "default",
                                        time_min = NULL,
                                        fig_dir = NULL,
                                        fps = 4,
                                        save = TRUE) {

  if (is.null(fig_dir)) fig_dir <- get_figure_dir()

  # Join azimuth_deg onto truth (x, y, z are already in truth_dt from build_truth_from_ce)
  rec_pos <- transect_receptors[, .(receptor_id, azimuth_deg)]
  dt <- merge(truth_dt, rec_pos, by = "receptor_id")

  if (!"ch4_enhancement_ppb" %in% names(dt)) {
    if ("ch4_total_ppb" %in% names(dt) && "ch4_bg_ppb" %in% names(dt)) {
      dt[, ch4_enhancement_ppb := ch4_total_ppb - ch4_bg_ppb]
    } else {
      stop("truth_dt must contain ch4_enhancement_ppb or ch4_total_ppb + ch4_bg_ppb.")
    }
  }

  all_times <- sort(unique(dt$time_min))
  if (!is.null(time_min)) all_times <- intersect(all_times, time_min)
  if (length(all_times) == 0L) {
    message("plot_plume_transect_heatmap: no matching time steps.")
    return(invisible(NULL))
  }

  enh_max <- quantile(dt$ch4_enhancement_ppb, 0.99, na.rm = TRUE)
  enh_max <- max(enh_max, 1e-3)

  safe_label <- gsub("[^a-zA-Z0-9_]", "_", scenario_label)

  make_frame <- function(t_use) {
    snap <- dt[time_min == t_use]
    hour_str <- sprintf("%04.1f h", (t_use / 60) %% 24)
    az <- if ("azimuth_deg" %in% names(snap) && nrow(snap) > 0)
      unique(snap$azimuth_deg)[1] else NA_real_

    # Compute downwind distance along the transect centreline
    snap[, dist_m := sqrt(x^2 + y^2)]

    ggplot(snap, aes(x = dist_m, y = z, fill = ch4_enhancement_ppb)) +
      geom_tile() +
      scale_fill_gradientn(
        colours = c("#f7fbff", "#c6dbef", "#6baed6", "#2171b5",
                    "#084594", "#ffeda0", "#feb24c", "#f03b20"),
        name   = "Enhancement\n(ppb)",
        limits = c(0, enh_max),
        oob    = scales::squish
      ) +
      scale_x_continuous(name = "Downwind distance (m)") +
      scale_y_continuous(name = "Height z (m)") +
      labs(
        title    = sprintf("Plume transect - %s", scenario_label),
        subtitle = sprintf(
          "Time %g min (%s) | downwind azimuth %.0f\u00b0",
          t_use, hour_str, if (is.finite(az)) az else 0
        )
      ) +
      theme_bw(base_size = 11) +
      theme(
        panel.grid.minor = element_blank(),
        legend.position  = "right"
      )
  }

  if (length(all_times) == 1L) {
    p <- make_frame(all_times[1L])
    if (save) {
      fname <- file.path(fig_dir,
                         sprintf("transect_heatmap_%s_t%g.png", safe_label, all_times[1L]))
      ggsave_osse(fname, p, width = 9, height = 5, dpi = 200)
      message("  Saved: ", basename(fname))
    }
    return(invisible(p))
  }

  # Multi-frame: animate
  if (!requireNamespace("magick", quietly = TRUE)) {
    message("Skipping transect animation: package `magick` not installed.")
    return(invisible(NULL))
  }

  tmp_dir <- file.path(tempdir(), paste0("transect_frames_", safe_label))
  dir.create(tmp_dir, showWarnings = FALSE)
  frame_files <- character(0)

  for (t_i in all_times) {
    p <- make_frame(t_i)
    fp <- file.path(tmp_dir, sprintf("frame_%06g.png", t_i))
    ok <- tryCatch({
      ggsave_osse(fp, p, width = 9, height = 5, dpi = 150)
      TRUE
    }, error = function(e) FALSE)
    if (ok && file.exists(fp)) frame_files <- c(frame_files, fp)
  }

  if (length(frame_files) == 0L) {
    message("Skipping transect animation: no frames rendered.")
    return(invisible(NULL))
  }

  img   <- magick::image_read(frame_files)
  anim  <- magick::image_animate(img, fps = fps)
  out_gif <- file.path(fig_dir,
                       sprintf("transect_heatmap_anim_%s.gif", safe_label))
  if (save) {
    magick::image_write(anim, path = out_gif, format = "gif")
    message("  Saved transect animation: ", basename(out_gif))
  }
  invisible(NULL)
}


# -----------------------------------------------------------------------------
# Gaussian plume visualisation
# -----------------------------------------------------------------------------
.ensure_gaussian_fwd <- function() {
  if (!exists("gaussian_plume_grid", mode = "function")) {
    source(file.path(project_root(), "R", "02b_gaussian_forward.R"))
  }
}

#' Downwind-distance vs height cross-section (x-z, y = 0).
#'
#' Evaluates the Gaussian plume on a fine grid through the plume centreline
#' and plots concentration as a filled heatmap.
#'
#' @param source_radius  Lake radius (m).
#' @param flux_mg_m2_d  Flux (mg CH4 m-2 d-1), default 200.
#' @param U_ref  Reference wind speed (m/s), default 5.
#' @param pasquill_class  Pasquill class "A"-"F", default "D".
#' @param wd_deg  Wind-from direction (degrees), default 225.
#' @param x_range  Downwind range (m), default c(1, 500).
#' @param z_range  Height range (m), default c(0.1, 30).
#' @param nx,ny  Grid resolution (default 250 x 120).
#' @param receptor_heights_m  Optional vector of heights to overlay as dashed lines.
#' @param save  Logical (default TRUE).
#' @param fig_dir  Output directory (default from get_figure_dir()).
#' @return Invisibly, the ggplot object.
plot_gaussian_plume_xz <- function(source_radius, flux_mg_m2_d = 200,
                                   U_ref = 5, pasquill_class = "D",
                                   wd_deg = 225,
                                   x_range = c(1, 500), z_range = c(0.1, 30),
                                   nx = 250L, ny = 120L,
                                   receptor_heights_m = NULL,
                                   save = TRUE, fig_dir = get_figure_dir()) {
  .ensure_gaussian_fwd()

  grid <- gaussian_plume_grid(
    source_radius, flux_mg_m2_d, U_ref, wd_deg, pasquill_class,
    plane = "xz", x_range = x_range, z_range = z_range, nx = nx, ny = ny
  )
  grid[conc_ppb < 1e-4, conc_ppb := NA_real_]

  p <- ggplot(grid, aes(x = x_down, y = z, fill = conc_ppb)) +
    geom_tile() +
    scale_fill_viridis_c(option = "C", na.value = "white", trans = "log10",
                         name = expression(Delta*CH[4]~"(ppb)")) +
    labs(
      title = sprintf("Gaussian plume (class %s): downwind x-z cross-section",
                      pasquill_class),
      subtitle = sprintf("F = %g mg m-2 d-1, U = %g m/s, R = %g m",
                         flux_mg_m2_d, U_ref, source_radius),
      x = "Downwind distance (m)", y = "Height (m)"
    ) +
    theme_osse() +
    theme(legend.position = "right")

  if (!is.null(receptor_heights_m)) {
    for (hz in receptor_heights_m) {
      p <- p + geom_hline(yintercept = hz, linetype = "dashed",
                          colour = "grey40", linewidth = 0.3)
    }
  }

  if (save) {
    fp <- file.path(fig_dir, sprintf("gaussian_plume_xz_%s.png", pasquill_class))
    ggsave_osse(fp, p, width = 10, height = 5, dpi = 200)
    message("  Saved: ", basename(fp))
  }
  invisible(p)
}


#' Plan view (x-y) at a fixed height.
#'
#' @param source_radius  Lake radius (m).
#' @param flux_mg_m2_d  Flux (mg CH4 m-2 d-1).
#' @param U_ref  Reference wind speed (m/s).
#' @param pasquill_class  Pasquill class "A"-"F".
#' @param wd_deg  Wind-from direction (degrees).
#' @param z_slice_m  Height of the horizontal slice (m), default 1.
#' @param x_range  Downwind range (m).
#' @param y_range  Crosswind range (m).
#' @param nx,ny  Grid resolution.
#' @param lake_diameter  Lake diameter for outline (m). NULL = 2*source_radius.
#' @param save  Logical.
#' @param fig_dir  Output directory.
#' @return Invisibly, the ggplot object.
plot_gaussian_plume_xy <- function(source_radius, flux_mg_m2_d = 200,
                                   U_ref = 5, pasquill_class = "D",
                                   wd_deg = 225, z_slice_m = 1,
                                   x_range = c(-50, 500),
                                   y_range = c(-250, 250),
                                   nx = 250L, ny = 200L,
                                   lake_diameter = NULL,
                                   save = TRUE, fig_dir = get_figure_dir()) {
  .ensure_gaussian_fwd()

  grid <- gaussian_plume_grid(
    source_radius, flux_mg_m2_d, U_ref, wd_deg, pasquill_class,
    plane = "xy", z_slice_m = z_slice_m,
    x_range = x_range, y_range = y_range, nx = nx, ny = ny
  )
  grid[conc_ppb < 1e-4, conc_ppb := NA_real_]

  if (is.null(lake_diameter)) lake_diameter <- 2 * source_radius

  lake_circle <- data.table(
    theta = seq(0, 2 * pi, length.out = 100)
  )
  wd_to_rad <- ((wd_deg + 180) %% 360) * pi / 180
  wind_dx   <- sin(wd_to_rad)
  wind_dy   <- cos(wd_to_rad)
  cross_dx  <- wind_dy
  cross_dy  <- -wind_dx
  lake_circle[, `:=`(
    lx = (lake_diameter / 2) * cos(theta),
    ly = (lake_diameter / 2) * sin(theta)
  )]
  lake_circle[, `:=`(
    x_down  = lx * wind_dx + ly * wind_dy,
    y_cross = lx * cross_dx + ly * cross_dy
  )]

  p <- ggplot(grid, aes(x = x_down, y = y_cross, fill = conc_ppb)) +
    geom_tile() +
    geom_path(data = lake_circle, aes(x = x_down, y = y_cross),
              inherit.aes = FALSE, colour = "white", linewidth = 0.8) +
    scale_fill_viridis_c(option = "C", na.value = "white", trans = "log10",
                         name = expression(Delta*CH[4]~"(ppb)")) +
    coord_equal() +
    labs(
      title = sprintf("Gaussian plume (class %s): plan view at z = %g m",
                      pasquill_class, z_slice_m),
      subtitle = sprintf("F = %g mg m-2 d-1, U = %g m/s, wind from %g deg",
                         flux_mg_m2_d, U_ref, wd_deg),
      x = "Downwind distance (m)", y = "Crosswind distance (m)"
    ) +
    theme_osse() +
    theme(legend.position = "right")

  if (save) {
    fp <- file.path(fig_dir, sprintf("gaussian_plume_xy_%s_z%gm.png",
                                     pasquill_class, z_slice_m))
    ggsave_osse(fp, p, width = 10, height = 7, dpi = 200)
    message("  Saved: ", basename(fp))
  }
  invisible(p)
}


#' Faceted stability comparison panel (A through F).
#'
#' Creates a 2x3 faceted plot showing the x-z plume cross-section for each
#' Pasquill class at the same flux and wind speed.
#'
#' @param source_radius  Lake radius (m).
#' @param flux_mg_m2_d  Flux (mg CH4 m-2 d-1).
#' @param U_ref  Reference wind speed (m/s).
#' @param wd_deg  Wind-from direction (degrees).
#' @param classes  Character vector of Pasquill classes (default A-F).
#' @param x_range  Downwind range (m).
#' @param z_range  Height range (m).
#' @param nx,ny  Grid resolution per panel.
#' @param save  Logical.
#' @param fig_dir  Output directory.
#' @return Invisibly, the ggplot object.
plot_gaussian_plume_stability_comparison <- function(
    source_radius, flux_mg_m2_d = 200, U_ref = 5, wd_deg = 225,
    classes = c("A", "B", "C", "D", "E", "F"),
    x_range = c(1, 400), z_range = c(0.1, 25),
    nx = 200L, ny = 80L,
    save = TRUE, fig_dir = get_figure_dir()) {

  .ensure_gaussian_fwd()

  grids <- lapply(classes, function(pc) {
    g <- gaussian_plume_grid(
      source_radius, flux_mg_m2_d, U_ref, wd_deg, pc,
      plane = "xz", x_range = x_range, z_range = z_range, nx = nx, ny = ny
    )
    g[, pasquill_class := pc]
    g
  })
  dt <- rbindlist(grids)
  dt[conc_ppb < 1e-4, conc_ppb := NA_real_]
  dt[, class_label := sprintf("Class %s (%s)", pasquill_class,
                              PASQUILL_CLASSES[[pasquill_class]]$description),
     by = pasquill_class]
  dt[, class_label := factor(class_label,
                             levels = sapply(classes, function(pc)
                               sprintf("Class %s (%s)", pc,
                                       PASQUILL_CLASSES[[pc]]$description)))]

  p <- ggplot(dt, aes(x = x_down, y = z, fill = conc_ppb)) +
    geom_tile() +
    facet_wrap(~class_label, ncol = 3L, scales = "fixed") +
    scale_fill_viridis_c(option = "C", na.value = "white", trans = "log10",
                         name = expression(Delta*CH[4]~"(ppb)")) +
    labs(
      title = "Gaussian plume: Pasquill stability class comparison (A-F)",
      subtitle = sprintf("F = %g mg m-2 d-1, U = %g m/s, R = %g m, centreline x-z",
                         flux_mg_m2_d, U_ref, source_radius),
      x = "Downwind distance (m)", y = "Height (m)"
    ) +
    theme_osse(base_size = 10) +
    theme(legend.position = "right",
          strip.text = element_text(size = 9))

  if (save) {
    fp <- file.path(fig_dir, "gaussian_plume_stability_comparison.png")
    ggsave_osse(fp, p, width = 14, height = 8, dpi = 200)
    message("  Saved: ", basename(fp))
  }
  invisible(p)
}


#' Dispersion coefficient curves: sigma_y and sigma_z vs downwind distance.
#'
#' Reference plot showing how the Briggs (1973) rural parameterisation
#' varies across Pasquill classes.
#'
#' @param classes  Character vector (default A-F).
#' @param x_range  Downwind range (m), default c(1, 2000).
#' @param nx  Number of points (default 500).
#' @param save  Logical.
#' @param fig_dir  Output directory.
#' @return Invisibly, the ggplot object.
plot_gaussian_plume_sigma_curves <- function(classes = c("A", "B", "C", "D", "E", "F"),
                                             x_range = c(1, 2000), nx = 500L,
                                             save = TRUE,
                                             fig_dir = get_figure_dir()) {
  .ensure_gaussian_fwd()

  x_seq <- seq(x_range[1], x_range[2], length.out = nx)
  rows <- lapply(classes, function(pc) {
    data.table(
      x_m            = x_seq,
      sigma_y        = pg_sigma_y(x_seq, pc),
      sigma_z        = pg_sigma_z(x_seq, pc),
      pasquill_class = pc
    )
  })
  dt <- rbindlist(rows)
  dt_long <- melt(dt, id.vars = c("x_m", "pasquill_class"),
                  measure.vars = c("sigma_y", "sigma_z"),
                  variable.name = "sigma_type", value.name = "sigma_m")

  dt_long[, sigma_label := factor(sigma_type,
                                  levels = c("sigma_y", "sigma_z"),
                                  labels = c("sigma[y]", "sigma[z]"))]

  p <- ggplot(dt_long, aes(x = x_m, y = sigma_m,
                            colour = pasquill_class)) +
    geom_line(linewidth = 0.7) +
    facet_wrap(~sigma_label, scales = "free_y", labeller = label_parsed) +
    scale_colour_viridis_d(option = "D", end = 0.9, name = "Pasquill\nclass") +
    scale_x_log10() + scale_y_log10() +
    labs(
      title = "Pasquill-Gifford dispersion coefficients (Briggs 1973, rural)",
      x = "Downwind distance (m)", y = "Dispersion coefficient (m)"
    ) +
    theme_osse() +
    theme(legend.position = "right")

  if (save) {
    fp <- file.path(fig_dir, "gaussian_plume_sigma_curves.png")
    ggsave_osse(fp, p, width = 11, height = 5, dpi = 200)
    message("  Saved: ", basename(fp))
  }
  invisible(p)
}

