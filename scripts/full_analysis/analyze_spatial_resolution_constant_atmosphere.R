#!/usr/bin/env Rscript
###############################################################################
# analyze_spatial_resolution_constant_atmosphere.R
#
# Full-analysis sweep: how does measurement spatial resolution affect Conley
# inversion accuracy under a constant atmosphere?
#
# Experiment configuration:
#   - 100 m lake
#   - constant whole-lake flux = 50 mg CH4 m^-2 d^-1
#   - one Gaussian-plume time step
#   - Pasquill classes: B, D, F
#   - LI-7810 sensor model
#   - Conley inverse with wind speed + direction error and height-averaged wind
#   - 100 Monte Carlo repeats per sampling design
#   - each repeat randomizes plume direction uniformly over 0-360 degrees
#   - vertical domain tiled from 0 to 100 m
#
# Outputs (under output/):
#   - spatial_resolution_constant_atmosphere_runs.csv
#   - spatial_resolution_constant_atmosphere_summary.csv
#   - spatial_resolution_constant_atmosphere_composite.png
#   - spatial_resolution_constant_atmosphere_precision.png
#
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

parse_chr_vec <- function(x) {
  if (is.null(x) || !nzchar(x)) return(character(0))
  trimws(strsplit(x, ",", fixed = TRUE)[[1]])
}

parse_bool <- function(x) {
  if (is.null(x) || !nzchar(x)) return(FALSE)
  tolower(trimws(x)) %in% c("true", "t", "1", "yes", "y")
}

cli <- parse_named_args(commandArgs(trailingOnly = TRUE))

LAKE_NAME       <- "small"
SENSOR_NAME     <- "li7810"
RETRIEVAL_NAME  <- "height_avg_wind"
BASE_FLUX_MG_D  <- if (!is.null(cli$flux_mg_m2_d)) as.numeric(cli$flux_mg_m2_d) else 50
TOP_HEIGHT_M    <- if (!is.null(cli$top_height_m)) as.numeric(cli$top_height_m) else 100
CYCLE_MIN       <- if (!is.null(cli$cycle_min)) as.numeric(cli$cycle_min) else 30
N_MC            <- if (!is.null(cli$n_mc)) as.integer(cli$n_mc) else 100L
AZIMUTH_COUNTS  <- {
  vals <- parse_num_vec(cli$azimuths)
  if (length(vals)) as.integer(vals) else c(8L, 16L, 24L, 36L, 48L, 72L, 96L)
}
HEIGHT_COUNTS   <- {
  vals <- parse_num_vec(cli$height_counts)
  if (length(vals)) as.integer(vals) else c(5L, 10L, 20L, 25L, 40L, 50L)
}
PASQUILL_CLASSES_TO_RUN <- {
  vals <- parse_chr_vec(cli$classes)
  if (length(vals)) vals else c("B", "D", "F")
}
U_REF_M_S       <- 2
WD_BASE_DEG     <- if (!is.null(cli$wd_base_deg)) as.numeric(cli$wd_base_deg) else 225
SEED_BASE       <- if (!is.null(cli$seed)) as.integer(cli$seed) else 20260414L
PAPER_STATS_ONLY <- parse_bool(cli$paper_stats_only)
OUTPUT_STEM     <- if (!is.null(cli$output_stem) && nzchar(cli$output_stem)) {
  cli$output_stem
} else {
  "spatial_resolution_constant_atmosphere"
}

stopifnot(LAKE_NAME %in% names(LAKE_CONFIGS))
stopifnot(SENSOR_NAME %in% names(SENSOR_CONFIGS))
stopifnot(RETRIEVAL_NAME %in% names(RETRIEVAL_CONFIGS))
stopifnot(all(PASQUILL_CLASSES_TO_RUN %in% names(PASQUILL_CLASSES)))

OUTPUT_DIR <- proj_path("output")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

RUNS_CSV    <- file.path(OUTPUT_DIR, sprintf("%s_runs.csv", OUTPUT_STEM))
SUMMARY_CSV <- file.path(OUTPUT_DIR, sprintf("%s_summary.csv", OUTPUT_STEM))
FIG_PNG     <- file.path(OUTPUT_DIR, sprintf("%s_composite.png", OUTPUT_STEM))
PRECISION_FIG_PNG <- file.path(OUTPUT_DIR, sprintf("%s_precision.png", OUTPUT_STEM))
KEY_CSV     <- file.path(OUTPUT_DIR, sprintf("%s_key_table.csv", OUTPUT_STEM))
MIN_AZ_CSV  <- file.path(OUTPUT_DIR, sprintf("%s_min_azimuth_by_height.csv", OUTPUT_STEM))
MIN_H_CSV   <- file.path(OUTPUT_DIR, sprintf("%s_min_height_by_azimuth.csv", OUTPUT_STEM))
PAPER_STATS_CSV <- file.path(OUTPUT_DIR, sprintf("%s_paper_stats.csv", OUTPUT_STEM))

FIG_BASE_SIZE <- 12
FIG_POINT_SIZE <- 2.8
FIG_LINEWIDTH <- 0.9
FIG_LABEL_SIZE <- 13
REPRESENTATIVE_PASQUILL_CLASS <- if ("B" %in% PASQUILL_CLASSES_TO_RUN) "B" else PASQUILL_CLASSES_TO_RUN[1]


# -----------------------------------------------------------------------------
# 2. Helpers
# -----------------------------------------------------------------------------
build_height_midpoints <- function(n_levels, top_height_m) {
  step_m <- top_height_m / n_levels
  seq(step_m / 2, top_height_m - step_m / 2, length.out = n_levels)
}

build_design_cfg <- function(n_azimuth, n_height_levels, top_height_m, cycle_min) {
  list(
    name = sprintf("res_%03daz_%03dh_%03dm", n_azimuth, n_height_levels, as.integer(round(top_height_m))),
    n_azimuth = as.integer(n_azimuth),
    heights_m = build_height_midpoints(as.integer(n_height_levels), top_height_m),
    cycle_min = cycle_min,
    inner_ring_fracs = NULL,
    background_heights_m = NULL
  )
}

geometry_metrics <- function(lake_diameter, design_cfg, distance_from_edge_m = 10) {
  lake_radius <- lake_diameter / 2
  ring_radius <- lake_radius + distance_from_edge_m
  height_bands <- compute_height_bands(design_cfg$heights_m)
  vertical_spacing_m <- mean(height_bands$dz)
  arc_spacing_m <- 2 * pi * ring_radius / design_cfg$n_azimuth
  effective_cell_area_m2 <- arc_spacing_m * vertical_spacing_m
  cylinder_surface_area_m2 <- 2 * pi * ring_radius * sum(height_bands$dz)

  data.table(
    n_azimuth = as.integer(design_cfg$n_azimuth),
    n_height_levels = as.integer(length(design_cfg$heights_m)),
    top_height_m = max(height_bands$z + height_bands$dz / 2),
    vertical_spacing_m = vertical_spacing_m,
    arc_spacing_m = arc_spacing_m,
    effective_cell_area_m2 = effective_cell_area_m2,
    n_receptors = as.integer(design_cfg$n_azimuth * length(design_cfg$heights_m)),
    cylinder_surface_area_m2 = cylinder_surface_area_m2
  )
}

build_one_step_met_cfg <- function(pasquill_class, u_ref_m_s, wd_base_deg) {
  start_hour <- if (pasquill_class %in% c("E", "F")) 0 else 12
  list(
    name = sprintf("gp_class_%s_one_step", pasquill_class),
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

sample_plume_direction_deg <- function(seed, offset_deg = 0) {
  (offset_deg + with_preserved_seed(runif(1, min = 0, max = 360), seed = seed)) %% 360
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

compute_summary <- function(run_dt, true_flux_mg_m2_d) {
  run_dt[, valid := qc_flag == "ok" & is.finite(flux_retrieved_mg_m2_d)]
  run_dt[, rel_error_pct := 100 * (flux_retrieved_mg_m2_d - true_flux_mg_m2_d) / true_flux_mg_m2_d]

  run_dt[, {
    valid_flux <- flux_retrieved_mg_m2_d[valid]
    valid_bias <- rel_error_pct[valid]
    scenario_seed <- string_to_seed(sprintf(
      "spatial|%s|%d|%d|%.4f",
      pasquill_class[1], n_azimuth[1], n_height_levels[1], effective_cell_area_m2[1]
    ))
    flux_stats <- distribution_detail_cols(
      valid_flux,
      base_name = "flux",
      suffix = "_mg_m2_d",
      seed = scenario_seed + 11L,
      include_abs_mean = FALSE
    )
    bias_stats <- distribution_detail_cols(
      valid_bias,
      base_name = "bias",
      suffix = "_pct",
      seed = scenario_seed + 29L,
      include_abs_mean = TRUE
    )

    c(list(
      flux_true_mg_m2_d = true_flux_mg_m2_d,
      mc_runs = .N,
      n_valid = sum(valid),
      valid_fraction = mean(valid),
      flux_mean_mg_m2_d = if (length(valid_flux)) mean(valid_flux, na.rm = TRUE) else NA_real_,
      flux_sd_mg_m2_d = safe_sd(valid_flux),
      bias_mean_pct = if (length(valid_bias)) mean(valid_bias, na.rm = TRUE) else NA_real_,
      bias_sd_pct = safe_sd(valid_bias),
      rmse_pct = if (length(valid_bias)) sqrt(mean(valid_bias^2, na.rm = TRUE)) else NA_real_,
      bias_direction = bias_direction_label(if (length(valid_bias)) mean(valid_bias, na.rm = TRUE) else NA_real_),
      rmse_target_met = if (length(valid_bias)) sqrt(mean(valid_bias^2, na.rm = TRUE)) < RESULTS_RMSE_TARGET_PCT else FALSE,
      mean_dwell_s = mean(dwell_s, na.rm = TRUE)
    ), flux_stats, bias_stats)
  }, by = .(
    pasquill_class,
    n_azimuth,
    n_height_levels,
    top_height_m,
    vertical_spacing_m,
    arc_spacing_m,
    effective_cell_area_m2,
    n_receptors,
    cylinder_surface_area_m2
  )][order(pasquill_class, effective_cell_area_m2, n_azimuth, n_height_levels)]
}

build_spatial_paper_stats <- function(summary_dt) {
  dt <- copy(summary_dt)
  primary_classes <- c("B", "D", "F")
  rows <- list()

  add_rows <- function(x) {
    if (nrow(x)) rows[[length(rows) + 1L]] <<- x
  }

  precision_dt <- dt[
    pasquill_class %in% primary_classes &
      n_height_levels == 50L &
      n_azimuth %in% c(8L, 96L)
  ]
  add_rows(precision_dt[, {
    precision_pct_true <- 100 * flux_sd_mg_m2_d / flux_true_mg_m2_d
    data.table(
      stat_id = ifelse(
        n_azimuth == 8L,
        "spatial_01_precision_lowest_azimuth_50_heights",
        "spatial_02_precision_highest_azimuth_50_heights"
      ),
      stat_description = ifelse(
        n_azimuth == 8L,
        "Precision at lowest azimuthal resolution: +/- 1 std dev of retrieved flux as percent of true flux at 8 azimuth cells and 50 vertical levels.",
        "Precision at highest azimuthal resolution: +/- 1 std dev of retrieved flux as percent of true flux at 96 azimuth cells and 50 vertical levels."
      ),
      mc_runs = mc_runs,
      flux_true_mg_m2_d = flux_true_mg_m2_d,
      flux_mean_mg_m2_d = flux_mean_mg_m2_d,
      flux_sd_mg_m2_d = flux_sd_mg_m2_d,
      precision_1sd_pct_true = precision_pct_true,
      bias_mean_pct = bias_mean_pct,
      bias_sd_pct = bias_sd_pct,
      value_text = sprintf("+/- %.2f%% of true flux", precision_pct_true),
      note = "precision_1sd_pct_true = 100 * flux_sd_mg_m2_d / flux_true_mg_m2_d"
    )
  }, by = .(pasquill_class, n_azimuth, n_height_levels)])

  bias_by_az_dt <- dt[
    pasquill_class %in% primary_classes &
      n_height_levels == 50L
  ]
  add_rows(bias_by_az_dt[, {
    mean_bias <- mean(bias_mean_pct, na.rm = TRUE)
    sd_bias <- stats::sd(bias_mean_pct, na.rm = TRUE)
    data.table(
      stat_id = "spatial_03_systematic_bias_across_azimuths_50_heights",
      stat_description = "Systematic bias across azimuthal resolutions: mean bias averaged across azimuthal configurations at 50 vertical levels.",
      n_height_levels = 50L,
      n_azimuth_configurations = .N,
      bias_mean_across_azimuths_pct = mean_bias,
      bias_sd_across_azimuths_pct = sd_bias,
      value_text = sprintf("%.2f +/- %.2f%%", mean_bias, sd_bias),
      note = "mean +/- std dev of bias_mean_pct across all azimuthal resolutions at 50 vertical levels"
    )
  }, by = pasquill_class])

  height_bias_dt <- dt[
    pasquill_class %in% c("F", "D") &
      n_azimuth == 96L &
      n_height_levels %in% c(40L, 50L)
  ]
  add_rows(height_bias_dt[, {
    bias_50 <- bias_mean_pct[n_height_levels == 50L][1]
    bias_40 <- bias_mean_pct[n_height_levels == 40L][1]
    data.table(
      stat_id = sprintf(
        "spatial_0%d_bias_change_50_to_40_heights_%s_96az",
        ifelse(pasquill_class[1] == "F", 4L, 5L),
        pasquill_class[1]
      ),
      stat_description = sprintf(
        "Bias change from 50 to 40 vertical levels for Pasquill %s at 96 azimuthal cells.",
        pasquill_class[1]
      ),
      n_azimuth = 96L,
      bias_50_height_levels_pct = bias_50,
      bias_40_height_levels_pct = bias_40,
      bias_change_40_minus_50_pct = bias_40 - bias_50,
      value_text = sprintf(
        "50 levels: %.2f%%; 40 levels: %.2f%%; change: %.2f percentage points",
        bias_50,
        bias_40,
        bias_40 - bias_50
      ),
      note = "bias_change_40_minus_50_pct = bias at 40 vertical levels minus bias at 50 vertical levels"
    )
  }, by = pasquill_class])

  rbindlist(rows, fill = TRUE)
}

write_spatial_paper_stats <- function(summary_dt, output_csv = PAPER_STATS_CSV) {
  paper_stats_dt <- build_spatial_paper_stats(summary_dt)
  fwrite(paper_stats_dt, output_csv)
  message(sprintf("Saved paper stats CSV: %s", normalizePath(output_csv, mustWork = FALSE)))
  invisible(paper_stats_dt)
}

if (isTRUE(PAPER_STATS_ONLY)) {
  if (!file.exists(SUMMARY_CSV)) {
    stop("Paper-stats-only mode requested, but summary CSV was not found: ", SUMMARY_CSV)
  }
  summary_dt <- fread(SUMMARY_CSV)
  write_spatial_paper_stats(summary_dt, PAPER_STATS_CSV)
  quit(save = "no", status = 0)
}

build_accuracy_plot <- function(summary_dt, true_flux_mg_m2_d) {
  plot_dt <- copy(summary_dt)
  plot_dt[, class_label := factor(
    paste("Pasquill", pasquill_class),
    levels = paste("Pasquill", PASQUILL_CLASSES_TO_RUN)
  )]
  plot_dt[, azimuth_label := factor(
    paste0(n_azimuth, " az"),
    levels = paste0(sort(unique(n_azimuth)), " az")
  )]
  plot_dt[, height_label := factor(
    paste0(n_height_levels, " heights"),
    levels = paste0(sort(unique(n_height_levels)), " heights")
  )]

  ggplot(
    plot_dt,
    aes(
      x = effective_cell_area_m2,
      y = flux_mean_mg_m2_d,
      colour = azimuth_label,
      fill = azimuth_label,
      shape = height_label,
      group = azimuth_label
    )
    ) +
    geom_ribbon(
      aes(
        ymin = flux_mean_mg_m2_d - flux_sd_mg_m2_d,
        ymax = flux_mean_mg_m2_d + flux_sd_mg_m2_d
      ),
      alpha = 0.14,
      linewidth = 0,
      show.legend = FALSE
    ) +
    geom_line(linewidth = FIG_LINEWIDTH, alpha = 0.9, show.legend = TRUE) +
    geom_point(size = FIG_POINT_SIZE, stroke = 0.3, show.legend = c(colour = FALSE, fill = FALSE, shape = TRUE)) +
    geom_hline(
      data = data.frame(true_flux_mg_m2_d = true_flux_mg_m2_d),
      aes(yintercept = true_flux_mg_m2_d, linetype = "true_flux"),
      inherit.aes = FALSE,
      colour = "black",
      linewidth = 0.75
    ) +
    facet_grid(class_label ~ ., scales = "free_y", switch = "y") +
    scale_x_log10() +
    scale_colour_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2") +
    labs(
      title = "Retrieval accuracy vs measurement spatial resolution",
      subtitle = sprintf("Ribbon = +/-1 stdev across %d replicates", N_MC),
      x = expression("Receptor cell area ("*m^2*")"),
      y = expression("Retrieved flux ("*mg~CH[4]~m^{-2}~d^{-1}*")"),
      colour = "Azimuth count",
      shape = "Vertical levels"
    ) +
    scale_linetype_manual(
      values = c(true_flux = "dashed"),
      labels = expression("True CH"[4]~"flux"),
      name = NULL
    ) +
    theme_osse(base_size = FIG_BASE_SIZE) +
    theme(
      legend.position = "right",
      panel.grid.minor = element_blank(),
      strip.placement = "outside",
      strip.background = element_rect(fill = "grey95", colour = NA)
    ) +
    guides(
      colour = guide_legend(order = 1),
      fill = "none",
      shape = guide_legend(order = 2),
      linetype = guide_legend(order = 3, override.aes = list(colour = "black", linewidth = 0.75))
    )
}

build_precision_plot <- function(summary_dt) {
  point_dt <- copy(summary_dt)
  point_dt[, class_label := factor(
    paste("Pasquill", pasquill_class),
    levels = paste("Pasquill", PASQUILL_CLASSES_TO_RUN)
  )]

  line_dt <- point_dt[, .(
    retrieval_precision_sd_mg_m2_d = mean(flux_sd_mg_m2_d, na.rm = TRUE)
  ), by = .(pasquill_class, class_label, effective_cell_area_m2)][
    order(pasquill_class, effective_cell_area_m2)
  ]

  class_colors <- c("Pasquill B" = "#B26A00", "Pasquill D" = "#2C7FB8", "Pasquill F" = "#238B45")
  class_colors <- class_colors[names(class_colors) %in% levels(point_dt$class_label)]

  ggplot() +
    geom_line(
      data = line_dt,
      aes(
        x = effective_cell_area_m2,
        y = retrieval_precision_sd_mg_m2_d,
        colour = class_label,
        group = class_label
      ),
      linewidth = 1.0
    ) +
    geom_point(
      data = point_dt,
      aes(
        x = effective_cell_area_m2,
        y = flux_sd_mg_m2_d,
        colour = class_label
      ),
      size = FIG_POINT_SIZE,
      alpha = 0.85
    ) +
    scale_x_log10() +
    scale_colour_manual(values = class_colors) +
    labs(
      title = "Retrieval precision vs measurement spatial resolution",
      x = expression("Receptor cell area ("*m^2*")"),
      y = expression("Retrieval precision (stdev; "*mg~CH[4]~m^{-2}~d^{-1}*")"),
      colour = "Stability class"
    ) +
    theme_osse(base_size = FIG_BASE_SIZE) +
    theme(
      legend.position = "right",
      panel.grid.minor = element_blank()
    )
}


# -----------------------------------------------------------------------------
# 3. Experiment setup
# -----------------------------------------------------------------------------
lk <- LAKE_CONFIGS[[LAKE_NAME]]
ret_cfg <- RETRIEVAL_CONFIGS[[RETRIEVAL_NAME]]
sensor_cfg <- SENSOR_CONFIGS[[SENSOR_NAME]]
sensor_cfg$drift_ppb_hr <- 0
sensor_cfg$cal_offset_ppb <- 0
lake_grid_resolution_m <- 1

message("========================================")
message(" Spatial Resolution Full Analysis")
message("========================================")
message(sprintf(" Lake:        %s (%d m diameter)", LAKE_NAME, lk$diameter_m))
message(sprintf(" Sensor:      %s", sensor_cfg$name))
message(sprintf(" Wind speed:  %.1f m/s (fixed true wind)", U_REF_M_S))
message(" Sensor err:  drift = 0 ppb/hr, calibration offset = 0 ppb")
message(sprintf(" Retrieval:   %s", ret_cfg$name))
message(sprintf(" Flux:        %.1f mg CH4 m^-2 d^-1", BASE_FLUX_MG_D))
message(sprintf(" Lake grid:   %.2f m", lake_grid_resolution_m))
message(sprintf(" Height top:  %.1f m", TOP_HEIGHT_M))
message(sprintf(" Cycle time:  %.1f min", CYCLE_MIN))
message(sprintf(" Monte Carlo: %d repeats", N_MC))
message(sprintf(" Classes:     %s", paste(PASQUILL_CLASSES_TO_RUN, collapse = ", ")))
message(sprintf(" Azimuths:    %s", paste(AZIMUTH_COUNTS, collapse = ", ")))
message(sprintf(" Heights:     %s", paste(HEIGHT_COUNTS, collapse = ", ")))
message("========================================")

emission_result <- quietly(generate_single_lake_emissions(
  lake_name = LAKE_NAME,
  save_output = FALSE,
  base_flux_override = BASE_FLUX_MG_D
))

# Use a fixed 1 m lake mesh in this experiment so the emission map and plotted
# lake pixels resolve the 100 m lake at high spatial detail.
emission_times_hour <- if ("hour" %in% names(emission_result$flux_field)) {
  sort(unique(emission_result$flux_field$hour))
} else {
  0:23
}
emission_result$lake_geom <- define_lake_geometry(
  diameter_m = lk$diameter_m,
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

design_grid <- CJ(
  pasquill_class = PASQUILL_CLASSES_TO_RUN,
  n_azimuth = sort(unique(AZIMUTH_COUNTS)),
  n_height_levels = sort(unique(HEIGHT_COUNTS))
)

run_pb <- make_progress(nrow(design_grid) * N_MC, label = "Resolution sweep")
run_rows <- vector("list", nrow(design_grid) * N_MC)
row_idx <- 0L

representative <- NULL
best_rep_area <- Inf


# -----------------------------------------------------------------------------
# 4. Run sweep
# -----------------------------------------------------------------------------
for (grid_idx in seq_len(nrow(design_grid))) {
  grid_row <- design_grid[grid_idx]
  design_cfg <- build_design_cfg(
    n_azimuth = grid_row$n_azimuth,
    n_height_levels = grid_row$n_height_levels,
    top_height_m = TOP_HEIGHT_M,
    cycle_min = CYCLE_MIN
  )
  geom_dt <- geometry_metrics(lk$diameter_m, design_cfg)
  for (mc_idx in seq_len(N_MC)) {
    plume_wd_deg <- sample_plume_direction_deg(
      seed = SEED_BASE + grid_idx * 1000L + mc_idx * 17L,
      offset_deg = WD_BASE_DEG
    )
    met_cfg <- build_one_step_met_cfg(
      pasquill_class = grid_row$pasquill_class,
      u_ref_m_s = U_REF_M_S,
      wd_base_deg = plume_wd_deg
    )
    scenario_label <- sprintf(
      "spatial_resolution_%s_%03daz_%03dh_wd%03d_mc%03d",
      grid_row$pasquill_class,
      grid_row$n_azimuth,
      grid_row$n_height_levels,
      as.integer(round(plume_wd_deg)),
      mc_idx
    )
    forward <- quietly(run_forward_scenario_gaussian(
      lake_name = LAKE_NAME,
      design_cfg = design_cfg,
      met_cfg = met_cfg,
      sensor_name = SENSOR_NAME,
      emission_result = emission_result,
      save_output = FALSE,
      met_seed = SEED_BASE + grid_idx * 1000L + mc_idx
    ))

    truth_dt <- copy(forward$truth_dt)
    receptors <- copy(forward$receptors)

    # This analysis studies enhancement sampling under a constant atmosphere.
    # Set the ambient background to zero so the inverse reflects plume
    # structure, sensor noise, wind uncertainty, and randomized plume
    # orientation only.
    truth_dt[, ch4_bg_ppb := 0]
    truth_dt[, ch4_total_ppb := ch4_enhancement_ppb]

    meas_dt <- apply_sensor_model(
      truth_dt = truth_dt,
      sensor_cfg = sensor_cfg,
      design_cfg = design_cfg,
      seed = SEED_BASE + grid_idx * 1000L + mc_idx
    )
    meas_dt <- prepare_measurements_for_inverse(meas_dt, truth_dt, receptors)

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

    row_idx <- row_idx + 1L
    run_rows[[row_idx]] <- cbind(
      data.table(
        pasquill_class = grid_row$pasquill_class,
        mc_id = mc_idx,
        plume_wd_deg = plume_wd_deg,
        flux_retrieved_mg_m2_d = flux_dt$flux_retrieved[1],
        flux_true_mg_m2_d = flux_dt$flux_true[1],
        qc_flag = flux_dt$qc_flag[1],
        mean_enhancement_ppb = flux_dt$mean_enhancement_ppb[1],
        ws_at_z_mean_m_s = flux_dt$ws_at_z_mean_m_s[1],
        dwell_s = unique(meas_dt$dwell_s)[1]
      ),
      geom_dt
    )

    candidate_is_representative_class <- grid_row$pasquill_class == REPRESENTATIVE_PASQUILL_CLASS
    current_is_representative_class <- !is.null(representative) &&
      representative$pasquill_class == REPRESENTATIVE_PASQUILL_CLASS
    if (is.null(representative) ||
        (candidate_is_representative_class &&
         mc_idx == 1L &&
         (!current_is_representative_class ||
          geom_dt$effective_cell_area_m2[1] < best_rep_area))) {
      representative <- list(
        measurements = copy(meas_dt),
        truth = copy(truth_dt),
        receptors = copy(receptors),
        scenario_label = scenario_label,
        pasquill_class = grid_row$pasquill_class,
        design_cfg = design_cfg,
        plume_wd_deg = plume_wd_deg,
        geometry = copy(geom_dt)
      )
      best_rep_area <- geom_dt$effective_cell_area_m2[1]
    }

    run_pb$tick()
  }
}

run_dt <- rbindlist(run_rows, fill = TRUE)
summary_dt <- compute_summary(run_dt, true_flux_mg_m2_d = BASE_FLUX_MG_D)
key_dt <- build_spatial_key_table(summary_dt)
min_az_dt <- build_spatial_min_az_table(summary_dt)
min_height_dt <- build_spatial_min_height_table(summary_dt)


# -----------------------------------------------------------------------------
# 5. Save tables
# -----------------------------------------------------------------------------
fwrite(run_dt, RUNS_CSV)
fwrite(summary_dt, SUMMARY_CSV)
fwrite(key_dt, KEY_CSV)
fwrite(min_az_dt, MIN_AZ_CSV)
fwrite(min_height_dt, MIN_H_CSV)
write_spatial_paper_stats(summary_dt, PAPER_STATS_CSV)
message(sprintf("Saved run-level CSV: %s", normalizePath(RUNS_CSV, mustWork = FALSE)))
message(sprintf("Saved summary CSV:   %s", normalizePath(SUMMARY_CSV, mustWork = FALSE)))
message(sprintf("Saved key table:     %s", normalizePath(KEY_CSV, mustWork = FALSE)))
message(sprintf("Saved azimuth table: %s", normalizePath(MIN_AZ_CSV, mustWork = FALSE)))
message(sprintf("Saved height table:  %s", normalizePath(MIN_H_CSV, mustWork = FALSE)))
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
  ws_observation_error_m_s = ret_cfg$ws_observation_error_m_s
)

left_plot <- plot_measured_isometric_over_emission_map(
  measurements_dt = representative$measurements,
  receptors = representative$receptors,
  flux_field = emission_result$flux_field,
  lake_diameter = lk$diameter_m,
  time_min = min(representative$measurements$time_min),
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
      "Pasquill %s | %d azimuths x %d heights | cell area = %.1f m² | receptor height = %.0f m",
      representative$pasquill_class,
      representative$design_cfg$n_azimuth,
      length(representative$design_cfg$heights_m),
      representative$geometry$effective_cell_area_m2[1],
      representative$geometry$top_height_m[1]
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

right_plot <- build_accuracy_plot(summary_dt, BASE_FLUX_MG_D)
precision_plot <- build_precision_plot(summary_dt)

left_legend <- cowplot::get_legend(left_plot)
right_legend <- cowplot::get_legend(right_plot)

left_main <- left_plot + theme(legend.position = "none")
right_main <- right_plot + theme(legend.position = "none")

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

right_column <- if (!is.null(right_legend)) {
  cowplot::plot_grid(
    right_main,
    right_legend,
    ncol = 2,
    rel_widths = c(1, 0.34)
  )
} else {
  right_main
}

composite_plot <- cowplot::plot_grid(
  left_column,
  right_column,
  ncol = 2,
  rel_widths = c(1.05, 1.25),
  labels = c("", "b"),
  label_size = FIG_LABEL_SIZE,
  label_x = c(0.02, 0.005),
  label_y = c(0.995, 0.995),
  hjust = 0,
  vjust = 1
)

ggsave_osse(FIG_PNG, composite_plot, width = 15, height = 11, dpi = 220)
message(sprintf("Saved composite figure: %s", normalizePath(FIG_PNG, mustWork = FALSE)))

ggsave_osse(PRECISION_FIG_PNG, precision_plot, width = 9.5, height = 6.5, dpi = 220)
message(sprintf("Saved precision figure: %s", normalizePath(PRECISION_FIG_PNG, mustWork = FALSE)))

message("Analysis complete.")
