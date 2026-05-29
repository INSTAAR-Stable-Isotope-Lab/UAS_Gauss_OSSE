#!/usr/bin/env Rscript
###############################################################################
# analyze_nearby_lake_contamination_distance.R
#
# Full-analysis sweep: how does contamination from a nearby lake affect the
# Conley inverse for the main lake?
#
# Requested experiment:
#   - main lake: 100 m diameter, 1 m emission grid, 50 mg CH4 m^-2 d^-1
#   - nearby lake: 50 m diameter, same areal flux by default
#   - nearby-lake center moved from 200 m to 10,000 m from the main lake
#   - one Gaussian-plume time step
#   - Pasquill classes: C, D, F
#   - LI-7810 sensor model
#   - Conley inverse with wind speed + direction error and height-averaged wind
#   - 100 Monte Carlo repeats per class-distance pair
#
# Outputs (under output/):
#   - nearby_lake_contamination_distance_runs.csv
#   - nearby_lake_contamination_distance_summary.csv
#   - nearby_lake_contamination_distance.png
#   - nearby_lake_contamination_distance_10km.png
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
suppressMessages(source(file.path(proj_root, "R", "03_drone_and_sensor_sim.R")))
suppressMessages(source(file.path(proj_root, "R", "04_conley_inverse_flux.R")))
suppressMessages(source(file.path(proj_root, "R", "06_visualization.R")))
suppressMessages(source(file.path(proj_root, "scripts", "full_analysis", "results_reporting_helpers.R")))

if (!requireNamespace("cowplot", quietly = TRUE)) {
  stop("Package `cowplot` is required for the stacked figure.")
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

parse_chr_vec <- function(x) {
  if (is.null(x) || !nzchar(x)) return(character(0))
  trimws(strsplit(x, ",", fixed = TRUE)[[1]])
}

parse_bool <- function(x) {
  if (is.null(x) || !nzchar(x)) return(FALSE)
  tolower(trimws(x)) %in% c("true", "t", "1", "yes", "y")
}

cli <- parse_named_args(commandArgs(trailingOnly = TRUE))

LAKE_NAME <- "small"
SENSOR_NAME <- "li7810"
RETRIEVAL_NAME <- if (!is.null(cli$retrieval)) cli$retrieval else "height_avg_wind"
DESIGN_NAME <- if (!is.null(cli$design)) cli$design else "design_96az_many_heights_30min"
TARGET_N_AZIMUTH <- if (!is.null(cli$n_azimuth)) as.integer(cli$n_azimuth) else 96L
BASE_FLUX_MG_D <- if (!is.null(cli$flux_mg_m2_d)) as.numeric(cli$flux_mg_m2_d) else 50
NEARBY_DIAMETER_M <- if (!is.null(cli$nearby_diameter_m)) as.numeric(cli$nearby_diameter_m) else 50
NEARBY_FLUX_MG_D <- if (!is.null(cli$nearby_flux_mg_m2_d)) as.numeric(cli$nearby_flux_mg_m2_d) else BASE_FLUX_MG_D
DISTANCES_M <- {
  vals <- parse_num_vec(cli$distances_m)
  if (length(vals)) vals else c(seq(200, 1000, by = 100), seq(2000, 10000, by = 1000))
}
PASQUILL_CLASSES_TO_RUN <- {
  vals <- parse_chr_vec(cli$classes)
  if (length(vals)) vals else c("C", "D", "F")
}
N_MC <- if (!is.null(cli$n_mc)) as.integer(cli$n_mc) else 100L
U_REF_M_S <- 2
WD_BASE_DEG <- if (!is.null(cli$wd_base_deg)) as.numeric(cli$wd_base_deg) else 225
SEED_BASE <- if (!is.null(cli$seed)) as.integer(cli$seed) else 20260414L
LAKE_GRID_RESOLUTION_M <- if (!is.null(cli$lake_grid_m)) as.numeric(cli$lake_grid_m) else 1
PLOT_ONLY <- parse_bool(cli$plot_only) || parse_bool(cli$plots_only)
PAPER_STATS_ONLY <- parse_bool(cli$paper_stats_only)
OUTPUT_STEM <- if (!is.null(cli$output_stem) && nzchar(cli$output_stem)) {
  cli$output_stem
} else {
  "nearby_lake_contamination_distance"
}

stopifnot(LAKE_NAME %in% names(LAKE_CONFIGS))
stopifnot(SENSOR_NAME %in% names(SENSOR_CONFIGS))
stopifnot(RETRIEVAL_NAME %in% names(RETRIEVAL_CONFIGS))
stopifnot(DESIGN_NAME %in% names(SAMPLING_CONFIGS))
stopifnot(all(PASQUILL_CLASSES_TO_RUN %in% names(PASQUILL_CLASSES)))

OUTPUT_DIR <- proj_path("output")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

RUNS_CSV <- file.path(OUTPUT_DIR, sprintf("%s_runs.csv", OUTPUT_STEM))
SUMMARY_CSV <- file.path(OUTPUT_DIR, sprintf("%s_summary.csv", OUTPUT_STEM))
FIG_PNG <- file.path(OUTPUT_DIR, sprintf("%s.png", OUTPUT_STEM))
FIG_10KM_PNG <- file.path(OUTPUT_DIR, sprintf("%s_10km.png", OUTPUT_STEM))
KEY_CSV <- file.path(OUTPUT_DIR, sprintf("%s_key_table.csv", OUTPUT_STEM))
BASELINE_CSV <- file.path(OUTPUT_DIR, sprintf("%s_clean_baseline_table.csv", OUTPUT_STEM))
PAPER_STATS_CSV <- file.path(OUTPUT_DIR, sprintf("%s_paper_stats.csv", OUTPUT_STEM))

FIG_BASE_SIZE <- 12
FIG_POINT_SIZE <- 2.5
FIG_LINEWIDTH <- 1.0
FIG_LABEL_SIZE <- 13


# -----------------------------------------------------------------------------
# 2. Helpers
# -----------------------------------------------------------------------------
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

build_one_step_met_cfg <- function(pasquill_class, u_ref_m_s, wd_base_deg) {
  start_hour <- if (pasquill_class %in% c("E", "F")) 0 else 12
  list(
    name = sprintf("contam_class_%s_one_step", pasquill_class),
    duration_days = 1 / 48,
    dt_minutes = 30,
    start_hour = start_hour,
    U_ref_m_s = u_ref_m_s,
    wd_base_deg = wd_base_deg,
    wd_noise_sd_deg = 0,
    pasquill_class = pasquill_class,
    z_ref_m = 10
  )
}

prepare_measurements_for_inverse <- function(meas_dt, truth_dt, receptors) {
  dt <- copy(meas_dt)
  if (!"azimuth_deg" %in% names(dt) || !"ring_type" %in% names(dt)) {
    dt <- merge(
      dt,
      receptors[, .(receptor_id, azimuth_deg, ring_type)],
      by = "receptor_id",
      all.x = TRUE
    )
  }
  if (!"z0_m" %in% names(dt) && "z0_m" %in% names(truth_dt)) {
    dt <- merge(
      dt,
      unique(truth_dt[, .(time_min, z0_m)]),
      by = "time_min",
      all.x = TRUE
    )
  } else if (!"z0_m" %in% names(dt)) {
    dt[, z0_m := 0.001]
  }
  dt
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

compute_contamination_metrics <- function(truth_dt,
                                          receptors,
                                          scenario_label,
                                          ret_cfg) {
  if (!("ch4_contamination_ppb" %in% names(truth_dt))) {
    return(data.table(
      flux_vector_contamination_pct = NA_real_,
      lateral_exit_fraction_pct = NA_real_
    ))
  }

  dt <- copy(truth_dt)
  dt <- merge(
    dt,
    receptors[, .(receptor_id, azimuth_deg, ring_type)],
    by = "receptor_id",
    all.x = TRUE
  )
  dt <- dt[ring_type == "outer"]

  if (nrow(dt) == 0L) {
    return(data.table(
      flux_vector_contamination_pct = NA_real_,
      lateral_exit_fraction_pct = NA_real_
    ))
  }

  n_az <- length(unique(dt$azimuth_deg))
  r_outer <- median(sqrt(dt$x^2 + dt$y^2), na.rm = TRUE)
  if (!is.finite(n_az) || n_az < 1L || !is.finite(r_outer) || r_outer <= 0) {
    return(data.table(
      flux_vector_contamination_pct = NA_real_,
      lateral_exit_fraction_pct = NA_real_
    ))
  }

  d_theta <- 2 * pi / as.integer(n_az)
  height_bands <- compute_height_bands(unique(dt$z))
  dt <- merge(dt, height_bands, by = "z", all.x = TRUE)
  dt[, dA := r_outer * d_theta * dz]

  if (ret_cfg$wd_uncertainty_deg > 0) {
    set.seed(string_to_seed(scenario_label) + 7777L)
    wd_noise <- data.table(
      time_min = sort(unique(dt$time_min)),
      wd_noise = rnorm(length(unique(dt$time_min)), 0, ret_cfg$wd_uncertainty_deg)
    )
    dt <- merge(dt, wd_noise, by = "time_min", all.x = TRUE)
    dt[, wd_deg := (wd_deg + wd_noise) %% 360]
    dt[, wd_noise := NULL]
  }

  if (ret_cfg$use_observed_ws && "ws_at_z_m_s" %in% names(dt)) {
    set.seed(string_to_seed(scenario_label) + 7778L)
    dt[, ws_z := pmax(0.1, ws_at_z_m_s + rnorm(.N, 0, ret_cfg$ws_observation_error_m_s))]
  } else {
    if (!"z0_m" %in% names(dt)) dt[, z0_m := 0.001]
    dt[, ws_z := wind_at_height(z_m = z, ustar = ustar_m_s[1L],
                                L_m = L_m[1L], z0_m = z0_m[1L]),
       by = time_min]
    if (ret_cfg$ws_uncertainty_m_s > 0) {
      set.seed(string_to_seed(scenario_label) + 7778L)
      ws_noise <- data.table(
        time_min = sort(unique(dt$time_min)),
        ws_noise = rnorm(length(unique(dt$time_min)), 0, ret_cfg$ws_uncertainty_m_s)
      )
      dt <- merge(dt, ws_noise, by = "time_min", all.x = TRUE)
      dt[, ws_z := pmax(0.1, ws_z + ws_noise)]
      dt[, ws_noise := NULL]
    }
  }

  if (ret_cfg$use_height_avg_wind) {
    dt[, ws_z := mean(ws_z), by = .(time_min, z)]
  }

  dt[, wd_to := (wd_deg + 180) %% 360]
  dt[, v_r := ws_z * cos((azimuth_deg - wd_to) * pi / 180)]
  dt[, contam_fv := ch4_contamination_ppb * v_r * dA]
  dt[, total_fv := ch4_total_ppb * v_r * dA]

  agg <- dt[, .(
    net_contam_fv = sum(contam_fv, na.rm = TRUE),
    net_total_fv = sum(total_fv, na.rm = TRUE),
    contam_inflow = sum(contam_fv[contam_fv < 0], na.rm = TRUE),
    contam_outflow = sum(contam_fv[contam_fv > 0], na.rm = TRUE)
  ), by = time_min]

  agg[, flux_vector_contamination_pct := ifelse(abs(net_total_fv) > 1e-9,
                                                100 * net_contam_fv / net_total_fv,
                                                NA_real_)]
  agg[, lateral_exit_fraction_pct := ifelse(abs(contam_inflow) > 1e-12,
                                            100 * contam_outflow / abs(contam_inflow),
                                            NA_real_)]

  agg[, .(
    flux_vector_contamination_pct = mean(flux_vector_contamination_pct, na.rm = TRUE),
    lateral_exit_fraction_pct = mean(lateral_exit_fraction_pct, na.rm = TRUE)
  )]
}

compute_summary <- function(run_dt, true_flux_mg_m2_d) {
  run_dt[, valid := qc_flag == "ok" & is.finite(flux_retrieved_mg_m2_d)]
  run_dt[, retrieval_bias_pct := 100 * (flux_retrieved_mg_m2_d - true_flux_mg_m2_d) / true_flux_mg_m2_d]
  run_dt[, baseline_bias_pct := 100 * (baseline_flux_retrieved_mg_m2_d - true_flux_mg_m2_d) / true_flux_mg_m2_d]

  run_dt[, {
    valid_flux <- flux_retrieved_mg_m2_d[valid]
    valid_retrieval_bias <- retrieval_bias_pct[valid]
    valid_contam_bias <- contamination_bias_pct[is.finite(contamination_bias_pct)]
    scenario_seed <- string_to_seed(sprintf("nearby|%s|%.2f", pasquill_class[1], distance_m[1]))

    flux_stats <- distribution_detail_cols(
      valid_flux,
      base_name = "retrieved_flux",
      suffix = "_mg_m2_d",
      seed = scenario_seed + 7L,
      include_abs_mean = FALSE
    )
    retrieval_bias_stats <- distribution_detail_cols(
      valid_retrieval_bias,
      base_name = "retrieval_bias",
      suffix = "_pct",
      seed = scenario_seed + 19L,
      include_abs_mean = TRUE
    )
    contamination_bias_stats <- distribution_detail_cols(
      valid_contam_bias,
      base_name = "contamination_bias",
      suffix = "_pct",
      seed = scenario_seed + 37L,
      include_abs_mean = TRUE
    )

    retrieval_rmse_now <- if (length(valid_retrieval_bias)) sqrt(mean(valid_retrieval_bias^2, na.rm = TRUE)) else NA_real_
    contamination_rmse_now <- if (length(valid_contam_bias)) sqrt(mean(valid_contam_bias^2, na.rm = TRUE)) else NA_real_
    retrieval_bias_mean_now <- if (length(valid_retrieval_bias)) mean(valid_retrieval_bias, na.rm = TRUE) else NA_real_

    c(list(
      flux_true_mg_m2_d = true_flux_mg_m2_d,
      mc_runs = .N,
      n_valid = sum(valid),
      valid_fraction = mean(valid),
      retrieved_flux_mean_mg_m2_d = if (length(valid_flux)) mean(valid_flux, na.rm = TRUE) else NA_real_,
      retrieved_flux_sd_mg_m2_d = safe_sd(valid_flux),
      retrieval_bias_mean_pct = retrieval_bias_mean_now,
      retrieval_bias_sd_pct = safe_sd(valid_retrieval_bias),
      retrieval_rmse_pct = retrieval_rmse_now,
      retrieval_bias_direction = bias_direction_label(retrieval_bias_mean_now),
      contamination_bias_mean_pct = if (length(valid_contam_bias)) mean(valid_contam_bias, na.rm = TRUE) else NA_real_,
      contamination_bias_sd_pct = safe_sd(valid_contam_bias),
      contamination_rmse_pct = contamination_rmse_now,
      flux_vector_contamination_mean_pct = mean(flux_vector_contamination_pct, na.rm = TRUE),
      flux_vector_contamination_sd_pct = safe_sd(flux_vector_contamination_pct),
      lateral_exit_fraction_mean_pct = mean(lateral_exit_fraction_pct, na.rm = TRUE),
      lateral_exit_fraction_sd_pct = safe_sd(lateral_exit_fraction_pct),
      rmse_target_met = is.finite(retrieval_rmse_now) && retrieval_rmse_now < RESULTS_RMSE_TARGET_PCT,
      mean_dwell_s = mean(dwell_s, na.rm = TRUE)
    ), flux_stats, retrieval_bias_stats, contamination_bias_stats)
  }, by = .(
    pasquill_class,
    distance_m,
    nearby_diameter_m,
    nearby_flux_mg_m2_d
  )][order(pasquill_class, distance_m)]
}

build_metric_panel <- function(plot_dt, y_col, y_sd_col, y_label, show_x = TRUE, show_true_flux = FALSE) {
  class_colors <- c(C = "#B26A00", D = "#2C7FB8", F = "#238B45")
  class_levels <- unique(as.character(plot_dt$pasquill_class))

  p <- ggplot(plot_dt, aes(x = distance_m, y = .data[[y_col]], colour = pasquill_class, fill = pasquill_class)) +
    geom_ribbon(aes(ymin = .data[[y_col]] - .data[[y_sd_col]], ymax = .data[[y_col]] + .data[[y_sd_col]]),
                alpha = 0.2, linewidth = 0, show.legend = FALSE, na.rm = TRUE) +
    geom_line(linewidth = FIG_LINEWIDTH) +
    geom_point(size = FIG_POINT_SIZE) +
    scale_colour_manual(values = class_colors[class_levels],
      labels = paste("Pasquill", class_levels),
      name = "Stability") +
    scale_fill_manual(values = class_colors[class_levels], guide = "none") +
    labs(
      x = if (show_x) "Distance between lake centers (m)" else NULL,
      y = y_label,
      colour = "Stability"
    ) +
    theme_osse(base_size = FIG_BASE_SIZE) +
    theme(
      legend.position = "bottom",
      axis.title.x = if (show_x) element_text() else element_blank(),
      axis.text.x = if (show_x) element_text() else element_blank(),
      axis.ticks.x = if (show_x) element_line() else element_blank()
    )

  if (show_true_flux) {
    p <- p + geom_hline(yintercept = BASE_FLUX_MG_D, linetype = "dashed", colour = "black", linewidth = 0.8)
  }

  p
}

build_composite_figure <- function(plot_dt,
                                   rep_truth_dt,
                                   rep_nearby_lake_config,
                                   rep_met_cfg,
                                   rep_subtitle,
                                   distance_label) {
  topdown_zoom_pad_m <- 15
  topdown_display_rotate_deg <- -atan2(
    rep_nearby_lake_config$center_y,
    rep_nearby_lake_config$center_x
  ) * 180 / pi
  rep_display_distance_m <- sqrt(
    rep_nearby_lake_config$center_x^2 +
      rep_nearby_lake_config$center_y^2
  )
  topdown_x_limits <- range(
    -ring_radius_m,
    ring_radius_m,
    rep_display_distance_m - rep_nearby_lake_config$radius_m,
    rep_display_distance_m + rep_nearby_lake_config$radius_m
  ) + c(-topdown_zoom_pad_m, topdown_zoom_pad_m)
  topdown_y_limits <- range(
    -ring_radius_m,
    ring_radius_m,
    -rep_nearby_lake_config$radius_m,
    rep_nearby_lake_config$radius_m
  ) + c(-topdown_zoom_pad_m, topdown_zoom_pad_m)

  p_topdown_fv <- plot_contamination_fv_topdown_planview(
    truth_dt = rep_truth_dt,
    receptors = receptors,
    lake_diameter = lk$diameter_m,
    nearby_lake_config = rep_nearby_lake_config,
    met_cfg = rep_met_cfg,
    z_m = rep_lowest_height_m,
    ring_radius_m = ring_radius_m,
    title = "Representative contamination flux-vector at receptors",
    subtitle = rep_subtitle,
    show_fill_legend = TRUE,
    display_rotate_deg = topdown_display_rotate_deg,
    receptor_point_size = 2.7,
    base_size = FIG_BASE_SIZE
  ) +
    coord_equal(
      xlim = topdown_x_limits,
      ylim = topdown_y_limits,
      expand = FALSE
    ) +
    theme(
      plot.margin = margin(5.5, 5.5, 5.5, 5.5),
      axis.line.x.bottom = element_line(colour = "grey30", linewidth = 0.45),
      axis.line.y.left = element_line(colour = "grey30", linewidth = 0.45),
      axis.ticks = element_line(colour = "grey30", linewidth = 0.35),
      legend.margin = margin(0, 0, 0, 0),
      legend.box.margin = margin(0, 0, 0, 0)
    )

  p_contam <- build_metric_panel(
    plot_dt = plot_dt,
    y_col = "contamination_bias_mean_pct",
    y_sd_col = "contamination_bias_sd_pct",
    y_label = "Contamination bias (%)",
    show_x = TRUE,
    show_true_flux = FALSE
  ) +
    labs(
      title = "Nearby-lake contamination impact on cylindrical-flux method",
      subtitle = sprintf(
        "%d azimuths x %d heights | 100 m main lake, 50 m nearby lake | %s",
        design_cfg$n_azimuth,
        length(design_cfg$heights_m),
        distance_label
      )
    ) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "black", linewidth = 0.7) +
    theme(plot.margin = margin(5.5, 5.5, 5.5, 0))

  p_exit <- build_metric_panel(
    plot_dt = plot_dt,
    y_col = "lateral_exit_fraction_mean_pct",
    y_sd_col = "lateral_exit_fraction_sd_pct",
    y_label = "Lateral exit fraction (%)",
    show_x = TRUE,
    show_true_flux = FALSE
  ) +
    geom_hline(yintercept = 100, linetype = "dashed", colour = "black", linewidth = 0.7) +
    theme(plot.margin = margin(5.5, 5.5, 5.5, 0))

  legend_plot <- p_contam +
    theme(
      legend.position = "bottom",
      legend.justification = "center",
      legend.margin = margin(0, 0, 0, 0),
      legend.box.margin = margin(0, 0, 0, 0),
      plot.margin = margin(0, 0, 0, 0)
    ) +
    guides(colour = guide_legend(nrow = 1))
  legend <- cowplot::get_legend(legend_plot)

  top_row <- cowplot::plot_grid(
    p_topdown_fv,
    ncol = 1,
    labels = "a",
    label_size = FIG_LABEL_SIZE
  )

  bottom_panels <- cowplot::plot_grid(
    p_contam + theme(legend.position = "none"),
    p_exit + theme(legend.position = "none"),
    ncol = 2,
    align = "hv",
    axis = "tb",
    labels = c("b", "c"),
    label_size = FIG_LABEL_SIZE,
    rel_widths = c(1, 1)
  )

  bottom_row <- if (!is.null(legend)) {
    cowplot::plot_grid(
      bottom_panels,
      legend,
      ncol = 1,
      rel_heights = c(1, 0.12)
    )
  } else {
    bottom_panels
  }

  cowplot::plot_grid(
    top_row,
    bottom_row,
    ncol = 1,
    rel_heights = c(1.35, 1),
    align = "v",
    axis = "lr"
  )
}

build_distance_response_figure <- function(plot_dt,
                                           distance_label) {
  p_contam <- build_metric_panel(
    plot_dt = plot_dt,
    y_col = "contamination_bias_mean_pct",
    y_sd_col = "contamination_bias_sd_pct",
    y_label = "Contamination bias (%)",
    show_x = TRUE,
    show_true_flux = FALSE
  ) +
    labs(
      title = "Nearby-lake contamination impact on cylindrical-flux retrieval",
      subtitle = sprintf(
        "%d azimuths x %d heights | 100 m main lake, 50 m nearby lake | %s",
        design_cfg$n_azimuth,
        length(design_cfg$heights_m),
        distance_label
      )
    ) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "black", linewidth = 0.7)

  p_exit <- build_metric_panel(
    plot_dt = plot_dt,
    y_col = "lateral_exit_fraction_mean_pct",
    y_sd_col = "lateral_exit_fraction_sd_pct",
    y_label = "Lateral exit fraction (%)",
    show_x = TRUE,
    show_true_flux = FALSE
  ) +
    geom_hline(yintercept = 100, linetype = "dashed", colour = "black", linewidth = 0.7)

  legend_plot <- p_contam +
    theme(
      legend.position = "bottom",
      legend.justification = "center",
      legend.margin = margin(0, 0, 0, 0),
      legend.box.margin = margin(0, 0, 0, 0),
      plot.margin = margin(0, 0, 0, 0)
    ) +
    guides(colour = guide_legend(nrow = 1))
  legend <- cowplot::get_legend(legend_plot)

  panels <- cowplot::plot_grid(
    p_contam + theme(legend.position = "none"),
    p_exit + theme(legend.position = "none"),
    ncol = 2,
    align = "hv",
    axis = "tb",
    labels = c("a", "b"),
    label_size = FIG_LABEL_SIZE
  )

  if (!is.null(legend)) {
    cowplot::plot_grid(
      panels,
      legend,
      ncol = 1,
      rel_heights = c(1, 0.11)
    )
  } else {
    panels
  }
}

build_nearby_paper_stats <- function(summary_dt) {
  dt <- copy(summary_dt)
  target_dt <- dt[
    pasquill_class %in% c("C", "D", "F") &
      distance_m %in% c(200, 1000)
  ][order(pasquill_class, distance_m)]

  target_dt[, .(
    stat_id = "nearby_01_contamination_bias_at_200m_and_1000m",
    stat_description = "Contamination bias at closest separation and at 1000 m for Pasquill C, D, and F.",
    pasquill_class,
    distance_m,
    nearby_diameter_m,
    nearby_flux_mg_m2_d,
    flux_true_mg_m2_d,
    mc_runs,
    contamination_bias_mean_pct,
    contamination_bias_sd_pct,
    contamination_rmse_pct,
    value_text = sprintf("%.2f +/- %.2f%%", contamination_bias_mean_pct, contamination_bias_sd_pct),
    note = "contamination_bias_pct = 100 * (retrieved flux with nearby lake - clean baseline retrieved flux) / true flux"
  )]
}

write_nearby_paper_stats <- function(summary_dt, output_csv = PAPER_STATS_CSV) {
  paper_stats_dt <- build_nearby_paper_stats(summary_dt)
  fwrite(paper_stats_dt, output_csv)
  message(sprintf("Saved paper stats CSV: %s", normalizePath(output_csv, mustWork = FALSE)))
  invisible(paper_stats_dt)
}

if (isTRUE(PAPER_STATS_ONLY)) {
  if (!file.exists(SUMMARY_CSV)) {
    stop("Paper-stats-only mode requested, but summary CSV was not found: ", SUMMARY_CSV)
  }
  summary_dt <- fread(SUMMARY_CSV)
  write_nearby_paper_stats(summary_dt, PAPER_STATS_CSV)
  quit(save = "no", status = 0)
}


# -----------------------------------------------------------------------------
# 3. Experiment setup
# -----------------------------------------------------------------------------
lk <- LAKE_CONFIGS[[LAKE_NAME]]
design_cfg <- SAMPLING_CONFIGS[[DESIGN_NAME]]
ret_cfg <- RETRIEVAL_CONFIGS[[RETRIEVAL_NAME]]
sensor_cfg <- SENSOR_CONFIGS[[SENSOR_NAME]]
sensor_cfg$drift_ppb_hr <- 0
sensor_cfg$cal_offset_ppb <- 0
if (!isTRUE(all.equal(design_cfg$n_azimuth, TARGET_N_AZIMUTH))) {
  design_cfg$n_azimuth <- TARGET_N_AZIMUTH
  design_cfg$name <- sprintf("%s_%daz", design_cfg$name, TARGET_N_AZIMUTH)
}

message("========================================")
message(" Nearby-Lake Contamination Analysis")
message("========================================")
message(sprintf(" Main lake:   %s (%d m diameter)", LAKE_NAME, lk$diameter_m))
message(sprintf(" Main flux:   %.1f mg CH4 m^-2 d^-1", BASE_FLUX_MG_D))
message(sprintf(" Lake grid:   %.2f m", LAKE_GRID_RESOLUTION_M))
message(sprintf(" Nearby lake: %.1f m diameter, %.1f mg CH4 m^-2 d^-1", NEARBY_DIAMETER_M, NEARBY_FLUX_MG_D))
message(sprintf(" Sensor:      %s", sensor_cfg$name))
message(sprintf(" Wind speed:  %.1f m/s (fixed true wind)", U_REF_M_S))
message(" Sensor err:  drift = 0 ppb/hr, calibration offset = 0 ppb")
message(sprintf(" Design:      %s", design_cfg$name))
message(sprintf(" Azimuths:    %d", design_cfg$n_azimuth))
message(sprintf(" Retrieval:   %s", ret_cfg$name))
message(sprintf(" Classes:     %s", paste(PASQUILL_CLASSES_TO_RUN, collapse = ", ")))
message(sprintf(" Distances:   %s", paste(DISTANCES_M, collapse = ", ")))
message(sprintf(" Monte Carlo: %d repeats", N_MC))
message("========================================")

emission_result <- configure_main_lake_emissions(
  lake_name = LAKE_NAME,
  base_flux_mg_m2_d = BASE_FLUX_MG_D,
  lake_grid_resolution_m = LAKE_GRID_RESOLUTION_M
)

receptors <- define_receptor_config(
  lake_diameter = lk$diameter_m,
  n_azimuth = design_cfg$n_azimuth,
  heights_m = design_cfg$heights_m,
  distance_from_edge_m = 10,
  inner_ring_fracs = design_cfg$inner_ring_fracs,
  background_heights_m = design_cfg$background_heights_m
)

scenario_grid <- CJ(
  pasquill_class = PASQUILL_CLASSES_TO_RUN,
  distance_m = sort(unique(DISTANCES_M))
)

if (isTRUE(PLOT_ONLY)) {
  if (!file.exists(SUMMARY_CSV)) {
    stop("Plot-only mode requested, but summary CSV was not found: ", SUMMARY_CSV)
  }
  summary_dt <- fread(SUMMARY_CSV)
  message(sprintf("Plot-only mode: loaded cached summary CSV: %s", normalizePath(SUMMARY_CSV, mustWork = FALSE)))
} else {
baseline_truth_by_class <- setNames(vector("list", length(PASQUILL_CLASSES_TO_RUN)), PASQUILL_CLASSES_TO_RUN)
for (pc in PASQUILL_CLASSES_TO_RUN) {
  baseline_met_cfg <- build_one_step_met_cfg(
    pasquill_class = pc,
    u_ref_m_s = U_REF_M_S,
    wd_base_deg = WD_BASE_DEG
  )
  baseline_met_ts <- generate_gaussian_met_timeseries(
    duration_days = baseline_met_cfg$duration_days,
    dt_minutes = baseline_met_cfg$dt_minutes,
    start_hour = baseline_met_cfg$start_hour,
    U_ref_m_s = baseline_met_cfg$U_ref_m_s,
    wd_base_deg = baseline_met_cfg$wd_base_deg,
    wd_noise_sd_deg = baseline_met_cfg$wd_noise_sd_deg,
    pasquill_class = baseline_met_cfg$pasquill_class,
    seed = SEED_BASE
  )
  baseline_truth_dt <- quietly(run_forward_model_gaussian(
    emission_result = emission_result,
    met_ts = baseline_met_ts,
    receptors = receptors,
    scenario_label = sprintf("nearby_contam_clean_%s", pc),
    save_output = FALSE,
    z_ref_m = baseline_met_cfg$z_ref_m,
    nearby_lake_config = NULL
  ))
  baseline_truth_dt[, ch4_bg_ppb := 0]
  baseline_truth_dt[, ch4_total_ppb := ch4_enhancement_ppb]
  baseline_truth_by_class[[pc]] <- baseline_truth_dt
}

run_pb <- make_progress(nrow(scenario_grid) * N_MC, label = "Nearby-lake sweep")
run_rows <- vector("list", nrow(scenario_grid) * N_MC)
row_idx <- 0L


# -----------------------------------------------------------------------------
# 4. Run sweep
# -----------------------------------------------------------------------------
for (grid_idx in seq_len(nrow(scenario_grid))) {
  grid_row <- scenario_grid[grid_idx]
  met_cfg <- build_one_step_met_cfg(
    pasquill_class = grid_row$pasquill_class,
    u_ref_m_s = U_REF_M_S,
    wd_base_deg = WD_BASE_DEG
  )

  nearby_lake_config <- resolve_nearby_lake_position(
    nearby_lake_cfg = list(
      diameter_m = NEARBY_DIAMETER_M,
      distance_m = grid_row$distance_m,
      angle_deg = WD_BASE_DEG,
      flux_mg_m2_d = NEARBY_FLUX_MG_D
    ),
    main_radius_m = lk$diameter_m / 2,
    ring_radius_m = lk$diameter_m / 2 + 10,
    wd_base_deg = WD_BASE_DEG,
    main_base_flux = BASE_FLUX_MG_D
  )

  met_ts <- generate_gaussian_met_timeseries(
    duration_days = met_cfg$duration_days,
    dt_minutes = met_cfg$dt_minutes,
    start_hour = met_cfg$start_hour,
    U_ref_m_s = met_cfg$U_ref_m_s,
    wd_base_deg = met_cfg$wd_base_deg,
    wd_noise_sd_deg = met_cfg$wd_noise_sd_deg,
    pasquill_class = met_cfg$pasquill_class,
    seed = SEED_BASE + grid_idx
  )

  scenario_label_base <- sprintf(
    "nearby_contam_%s_%04dm",
    grid_row$pasquill_class,
    as.integer(round(grid_row$distance_m))
  )

  truth_dt <- quietly(run_forward_model_gaussian(
    emission_result = emission_result,
    met_ts = met_ts,
    receptors = receptors,
    scenario_label = scenario_label_base,
    save_output = FALSE,
    z_ref_m = met_cfg$z_ref_m,
    nearby_lake_config = nearby_lake_config
  ))

  # Keep the analysis focused on plume enhancement and nearby-lake contamination.
  truth_dt[, ch4_bg_ppb := 0]
  truth_dt[, ch4_total_ppb := ch4_enhancement_ppb]

  for (mc_idx in seq_len(N_MC)) {
    scenario_label <- sprintf("%s_mc%03d", scenario_label_base, mc_idx)
    baseline_truth_dt <- baseline_truth_by_class[[grid_row$pasquill_class]]
    sensor_seed <- SEED_BASE + grid_idx * 1000L + mc_idx

    meas_dt <- apply_sensor_model(
      truth_dt = truth_dt,
      sensor_cfg = sensor_cfg,
      design_cfg = design_cfg,
      seed = sensor_seed
    )
    meas_dt <- prepare_measurements_for_inverse(meas_dt, truth_dt, receptors)

    baseline_meas_dt <- apply_sensor_model(
      truth_dt = baseline_truth_dt,
      sensor_cfg = sensor_cfg,
      design_cfg = design_cfg,
      seed = sensor_seed
    )
    baseline_meas_dt <- prepare_measurements_for_inverse(baseline_meas_dt, baseline_truth_dt, receptors)

    flux_dt <- quietly(retrieve_flux(
      meas_dt = meas_dt,
      lake_diameter = lk$diameter_m,
      distance_from_edge_m = 10,
      n_azimuth = design_cfg$n_azimuth,
      wd_uncertainty_deg = ret_cfg$wd_uncertainty_deg,
      ws_uncertainty_m_s = ret_cfg$ws_uncertainty_m_s,
      use_observed_ws = ret_cfg$use_observed_ws,
      ws_observation_error_m_s = ret_cfg$ws_observation_error_m_s,
      use_height_avg_wind = ret_cfg$use_height_avg_wind,
      scenario_label = scenario_label,
      save_output = FALSE
    ))

    baseline_flux_dt <- quietly(retrieve_flux(
      meas_dt = baseline_meas_dt,
      lake_diameter = lk$diameter_m,
      distance_from_edge_m = 10,
      n_azimuth = design_cfg$n_azimuth,
      wd_uncertainty_deg = ret_cfg$wd_uncertainty_deg,
      ws_uncertainty_m_s = ret_cfg$ws_uncertainty_m_s,
      use_observed_ws = ret_cfg$use_observed_ws,
      ws_observation_error_m_s = ret_cfg$ws_observation_error_m_s,
      use_height_avg_wind = ret_cfg$use_height_avg_wind,
      scenario_label = scenario_label,
      save_output = FALSE
    ))

    contam_metrics <- compute_contamination_metrics(
      truth_dt = truth_dt,
      receptors = receptors,
      scenario_label = scenario_label,
      ret_cfg = ret_cfg
    )
    contamination_bias_pct <- if (is.finite(flux_dt$flux_retrieved[1]) &&
                                  is.finite(baseline_flux_dt$flux_retrieved[1])) {
      100 * (flux_dt$flux_retrieved[1] - baseline_flux_dt$flux_retrieved[1]) / BASE_FLUX_MG_D
    } else {
      NA_real_
    }

    row_idx <- row_idx + 1L
    run_rows[[row_idx]] <- data.table(
      pasquill_class = grid_row$pasquill_class,
      distance_m = grid_row$distance_m,
      nearby_diameter_m = NEARBY_DIAMETER_M,
      nearby_flux_mg_m2_d = NEARBY_FLUX_MG_D,
      mc_id = mc_idx,
      flux_retrieved_mg_m2_d = flux_dt$flux_retrieved[1],
      flux_true_mg_m2_d = flux_dt$flux_true[1],
      qc_flag = flux_dt$qc_flag[1],
      mean_enhancement_ppb = flux_dt$mean_enhancement_ppb[1],
      ws_at_z_mean_m_s = flux_dt$ws_at_z_mean_m_s[1],
      dwell_s = unique(meas_dt$dwell_s)[1],
      baseline_flux_retrieved_mg_m2_d = baseline_flux_dt$flux_retrieved[1],
      contamination_bias_pct = contamination_bias_pct,
      flux_vector_contamination_pct = contam_metrics$flux_vector_contamination_pct[1],
      lateral_exit_fraction_pct = contam_metrics$lateral_exit_fraction_pct[1]
    )

    run_pb$tick()
  }
}

run_dt <- rbindlist(run_rows, fill = TRUE)
summary_dt <- compute_summary(run_dt, true_flux_mg_m2_d = BASE_FLUX_MG_D)
baseline_dt <- build_nearby_baseline_table(run_dt, true_flux_mg_m2_d = BASE_FLUX_MG_D)
key_dt <- build_nearby_key_table(summary_dt, baseline_dt)


# -----------------------------------------------------------------------------
# 5. Save tables
# -----------------------------------------------------------------------------
fwrite(run_dt, RUNS_CSV)
fwrite(summary_dt, SUMMARY_CSV)
fwrite(key_dt, KEY_CSV)
fwrite(baseline_dt, BASELINE_CSV)
message(sprintf("Saved run-level CSV: %s", normalizePath(RUNS_CSV, mustWork = FALSE)))
message(sprintf("Saved summary CSV:   %s", normalizePath(SUMMARY_CSV, mustWork = FALSE)))
message(sprintf("Saved key table:     %s", normalizePath(KEY_CSV, mustWork = FALSE)))
message(sprintf("Saved baseline table: %s", normalizePath(BASELINE_CSV, mustWork = FALSE)))
if (isTRUE(maybe_write_full_analysis_synthesis(OUTPUT_DIR))) {
  message("Updated cross-analysis synthesis table.")
}
}

write_nearby_paper_stats(summary_dt, PAPER_STATS_CSV)


# -----------------------------------------------------------------------------
# 6. Build figure
# -----------------------------------------------------------------------------
plot_dt <- copy(summary_dt)
plot_dt[, pasquill_class := factor(pasquill_class, levels = PASQUILL_CLASSES_TO_RUN)]
plot_dt_1km <- plot_dt[distance_m <= 1000]
plot_dt_10km <- plot_dt[distance_m <= 10000]

rep_pasquill_class <- if ("D" %in% PASQUILL_CLASSES_TO_RUN) "D" else PASQUILL_CLASSES_TO_RUN[1]
rep_distance_m <- min(DISTANCES_M, na.rm = TRUE)
rep_lowest_height_m <- min(design_cfg$heights_m, na.rm = TRUE)
ring_radius_m <- lk$diameter_m / 2 + 10

rep_met_cfg <- build_one_step_met_cfg(
  pasquill_class = rep_pasquill_class,
  u_ref_m_s = U_REF_M_S,
  wd_base_deg = WD_BASE_DEG
)
rep_nearby_lake_config <- resolve_nearby_lake_position(
  nearby_lake_cfg = list(
    diameter_m = NEARBY_DIAMETER_M,
    distance_m = rep_distance_m,
    angle_deg = WD_BASE_DEG,
    flux_mg_m2_d = NEARBY_FLUX_MG_D
  ),
  main_radius_m = lk$diameter_m / 2,
  ring_radius_m = ring_radius_m,
  wd_base_deg = WD_BASE_DEG,
  main_base_flux = BASE_FLUX_MG_D
)
rep_met_ts <- generate_gaussian_met_timeseries(
  duration_days = rep_met_cfg$duration_days,
  dt_minutes = rep_met_cfg$dt_minutes,
  start_hour = rep_met_cfg$start_hour,
  U_ref_m_s = rep_met_cfg$U_ref_m_s,
  wd_base_deg = rep_met_cfg$wd_base_deg,
  wd_noise_sd_deg = rep_met_cfg$wd_noise_sd_deg,
  pasquill_class = rep_met_cfg$pasquill_class,
  seed = SEED_BASE
)
rep_truth_dt <- quietly(run_forward_model_gaussian(
  emission_result = emission_result,
  met_ts = rep_met_ts,
  receptors = receptors,
  scenario_label = sprintf("nearby_contam_rep_%s_%04dm", rep_pasquill_class, as.integer(round(rep_distance_m))),
  save_output = FALSE,
  z_ref_m = rep_met_cfg$z_ref_m,
  nearby_lake_config = rep_nearby_lake_config
))
rep_truth_dt[, `:=`(ch4_bg_ppb = 0, ch4_total_ppb = ch4_enhancement_ppb)]

rep_subtitle <- sprintf(
  "Pasquill %s | %dm separation",
  rep_pasquill_class,
  as.integer(round(rep_distance_m))
)

main_grid_1km <- build_composite_figure(
  plot_dt = plot_dt_1km,
  rep_truth_dt = rep_truth_dt,
  rep_nearby_lake_config = rep_nearby_lake_config,
  rep_met_cfg = rep_met_cfg,
  rep_subtitle = rep_subtitle,
  distance_label = "distance sweep to 1 km"
)
main_grid_10km <- build_distance_response_figure(
  plot_dt = plot_dt_10km,
  distance_label = "distance sweep to 10 km"
)

ggsave_osse(FIG_PNG, main_grid_1km, width = 14, height = 12.5, dpi = 220, bg = "white")
message(sprintf("Saved figure:        %s", normalizePath(FIG_PNG, mustWork = FALSE)))
ggsave_osse(FIG_10KM_PNG, main_grid_10km, width = 14, height = 6.5, dpi = 220, bg = "white")
message(sprintf("Saved 10 km response figure: %s", normalizePath(FIG_10KM_PNG, mustWork = FALSE)))

message("Analysis complete.")
