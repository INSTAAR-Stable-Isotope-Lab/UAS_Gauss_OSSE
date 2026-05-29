suppressPackageStartupMessages({
  library(data.table)
})

RESULTS_RMSE_TARGET_PCT <- 25
RESULTS_ABS_BIAS_TARGET_PCT <- 10
RESULTS_CONTAM_BIAS_THRESHOLDS_PCT <- c(10, 5)
RESULTS_FV_THRESHOLDS_PCT <- c(10)
RESULTS_BOOTSTRAP_REPS <- 1000L

safe_sd <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) >= 2L) stats::sd(x) else 0
}

safe_quantile <- function(x, probs) {
  x <- x[is.finite(x)]
  if (!length(x)) return(rep(NA_real_, length(probs)))
  as.numeric(stats::quantile(x, probs = probs, na.rm = TRUE, names = FALSE, type = 7))
}

safe_min_value <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  min(x)
}

safe_max_value <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  max(x)
}

safe_span <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  max(x) - min(x)
}

safe_which_min <- function(x) {
  idx <- which.min(replace(x, !is.finite(x), Inf))
  if (!length(idx) || !is.finite(x[idx])) return(NA_integer_)
  idx
}

safe_which_max <- function(x) {
  idx <- which.max(replace(x, !is.finite(x), -Inf))
  if (!length(idx) || !is.finite(x[idx])) return(NA_integer_)
  idx
}

bootstrap_mean_ci <- function(x,
                              seed = 1L,
                              n_boot = RESULTS_BOOTSTRAP_REPS,
                              conf = 0.95) {
  x <- x[is.finite(x)]
  if (!length(x)) return(c(NA_real_, NA_real_))
  if (length(x) == 1L) return(c(x[1], x[1]))

  set.seed(seed)
  n <- length(x)
  boot_means <- vapply(seq_len(n_boot), function(i) {
    mean(sample(x, size = n, replace = TRUE))
  }, numeric(1))
  alpha <- (1 - conf) / 2
  safe_quantile(boot_means, c(alpha, 1 - alpha))
}

distribution_detail_cols <- function(x,
                                     base_name,
                                     suffix = "",
                                     seed = 1L,
                                     include_abs_mean = TRUE,
                                     include_boot_ci = TRUE) {
  x <- x[is.finite(x)]
  out <- list()

  median_name <- paste0(base_name, "_median", suffix)
  abs_mean_name <- paste0(base_name, "_abs_mean", suffix)
  p025_name <- paste0(base_name, "_p025", suffix)
  p05_name <- paste0(base_name, "_p05", suffix)
  p95_name <- paste0(base_name, "_p95", suffix)
  p975_name <- paste0(base_name, "_p975", suffix)
  boot_lo_name <- paste0(base_name, "_mean_boot_lo95", suffix)
  boot_hi_name <- paste0(base_name, "_mean_boot_hi95", suffix)

  if (!length(x)) {
    out[[median_name]] <- NA_real_
    if (include_abs_mean) out[[abs_mean_name]] <- NA_real_
    out[[p025_name]] <- NA_real_
    out[[p05_name]] <- NA_real_
    out[[p95_name]] <- NA_real_
    out[[p975_name]] <- NA_real_
    if (include_boot_ci) {
      out[[boot_lo_name]] <- NA_real_
      out[[boot_hi_name]] <- NA_real_
    }
    return(out)
  }

  qs <- safe_quantile(x, c(0.025, 0.05, 0.95, 0.975))
  out[[median_name]] <- stats::median(x)
  if (include_abs_mean) out[[abs_mean_name]] <- mean(abs(x))
  out[[p025_name]] <- qs[1]
  out[[p05_name]] <- qs[2]
  out[[p95_name]] <- qs[3]
  out[[p975_name]] <- qs[4]

  if (include_boot_ci) {
    ci <- bootstrap_mean_ci(x, seed = seed)
    out[[boot_lo_name]] <- ci[1]
    out[[boot_hi_name]] <- ci[2]
  }

  out
}

bias_direction_label <- function(x, near_zero_tol_pct = 1) {
  if (!is.finite(x)) return(NA_character_)
  if (abs(x) <= near_zero_tol_pct) return("near-unbiased")
  if (x > 0) "overestimate" else "underestimate"
}

nearest_tested_value <- function(x) {
  x <- sort(unique(x[is.finite(x)]))
  if (!length(x)) return(NA_real_)
  target <- stats::median(x)
  x[which.min(abs(x - target))]
}

pick_reference_row <- function(dt, x_col, reference_value, metric_col = NULL) {
  d <- copy(dt)[is.finite(get(x_col))]
  if (!nrow(d) || !is.finite(reference_value)) return(d[0])

  d[, .dist_from_ref := abs(get(x_col) - reference_value)]
  if (!is.null(metric_col) && metric_col %in% names(d)) {
    ord <- order(d$.dist_from_ref, d[[metric_col]])
  } else {
    ord <- order(d$.dist_from_ref)
  }
  out <- d[ord][1]
  out[, .dist_from_ref := NULL]
  out
}

pick_threshold_row <- function(dt,
                               x_col,
                               metric_col,
                               threshold,
                               direction = c("largest", "smallest"),
                               abs_metric = FALSE,
                               bias_col = NULL) {
  direction <- match.arg(direction)
  d <- copy(dt)[is.finite(get(x_col)) & is.finite(get(metric_col))]
  if (!nrow(d)) return(d[0])

  if (abs_metric) {
    d <- d[abs(get(metric_col)) < threshold]
  } else {
    d <- d[get(metric_col) < threshold]
  }
  if (!nrow(d)) return(d)

  if (is.null(bias_col)) {
    bias_candidates <- c("bias_mean_pct", "retrieval_bias_mean_pct", "contamination_bias_mean_pct")
    bias_col <- bias_candidates[bias_candidates %in% names(d)][1]
  }

  d[, .sort_x := get(x_col)]
  d[, .sort_valid := if ("valid_fraction" %in% names(d)) valid_fraction else 0]
  d[, .sort_bias := if (!is.na(bias_col) && bias_col %in% names(d)) abs(get(bias_col)) else 0]

  ord <- if (direction == "largest") {
    order(-d$.sort_x, -d$.sort_valid, d$.sort_bias, d[[metric_col]])
  } else {
    order(d$.sort_x, -d$.sort_valid, d$.sort_bias, d[[metric_col]])
  }

  out <- d[ord][1]
  out[, c(".sort_x", ".sort_valid", ".sort_bias") := NULL]
  out
}

endpoint_slope <- function(x, y, log_x = FALSE) {
  dt <- data.table(x = x, y = y)[is.finite(x) & is.finite(y)]
  if (!nrow(dt)) return(NA_real_)

  dt <- dt[, .(y = mean(y, na.rm = TRUE)), by = x][order(x)]
  if (nrow(dt) < 2L) return(NA_real_)
  n_rows <- nrow(dt)

  x0 <- dt$x[1]
  x1 <- dt$x[n_rows]
  if (log_x) {
    if (x0 <= 0 || x1 <= 0) return(NA_real_)
    x0 <- log10(x0)
    x1 <- log10(x1)
  }

  denom <- x1 - x0
  if (!is.finite(denom) || denom == 0) return(NA_real_)
  (dt$y[n_rows] - dt$y[1]) / denom
}

half_step_change_from_lm <- function(x, y, log_x = FALSE) {
  dt <- data.table(x = x, y = y)[is.finite(x) & is.finite(y)]
  if (nrow(dt) < 2L) return(NA_real_)
  if (log_x) {
    dt <- dt[x > 0]
    if (nrow(dt) < 2L || uniqueN(dt$x) < 2L) return(NA_real_)
    fit <- stats::lm(y ~ log10(x), data = dt)
    coef(fit)[["log10(x)"]] * log10(0.5)
  } else {
    if (uniqueN(dt$x) < 2L) return(NA_real_)
    fit <- stats::lm(y ~ x, data = dt)
    coef(fit)[["x"]] * (-0.5)
  }
}

percent_reduction <- function(start, end, absolute = TRUE) {
  if (!is.finite(start) || !is.finite(end)) return(NA_real_)
  denom <- if (absolute) abs(start) else start
  numer <- if (absolute) abs(start) - abs(end) else start - end
  if (!is.finite(denom) || denom == 0) return(NA_real_)
  100 * numer / denom
}

fmt_num <- function(x, digits = 2) {
  if (!is.finite(x)) return("NA")
  formatC(x, format = "f", digits = digits, big.mark = ",")
}

fmt_pct <- function(x, digits = 1) {
  if (!is.finite(x)) return("NA")
  paste0(fmt_num(x, digits = digits), "%")
}

fmt_hours <- function(x, digits = 1) {
  if (!is.finite(x)) return("NA")
  paste0(fmt_num(x, digits = digits), " h")
}

fmt_distance <- function(x, digits = 0) {
  if (!is.finite(x)) return("NA")
  paste0(fmt_num(x, digits = digits), " m")
}

fmt_design <- function(n_azimuth, n_heights) {
  if (!is.finite(n_azimuth) || !is.finite(n_heights)) return("NA")
  sprintf("%d azimuths x %d heights", as.integer(round(n_azimuth)), as.integer(round(n_heights)))
}

build_spatial_key_table <- function(summary_dt) {
  class_levels <- sort(unique(summary_dt$pasquill_class))

  rbindlist(lapply(class_levels, function(pc) {
    dt_class <- copy(summary_dt[pasquill_class == pc])[order(effective_cell_area_m2, n_azimuth, n_height_levels)]
    recommended <- pick_threshold_row(
      dt_class,
      x_col = "effective_cell_area_m2",
      metric_col = "rmse_pct",
      threshold = RESULTS_RMSE_TARGET_PCT,
      direction = "largest",
      abs_metric = FALSE,
      bias_col = "bias_mean_pct"
    )
    if (!nrow(recommended)) {
      best_idx <- safe_which_min(dt_class$rmse_pct)
      recommended <- if (is.finite(best_idx)) dt_class[best_idx] else dt_class[1]
    }

    reference_area <- nearest_tested_value(dt_class$effective_cell_area_m2)
    reference_row <- pick_reference_row(dt_class, "effective_cell_area_m2", reference_area, metric_col = "rmse_pct")
    finest_row <- dt_class[order(effective_cell_area_m2)][1]
    best_rmse_idx <- safe_which_min(dt_class$rmse_pct)
    best_rmse_row <- if (is.finite(best_rmse_idx)) dt_class[best_rmse_idx] else dt_class[1]

    data.table(
      pasquill_class = pc,
      rmse_target_pct = RESULTS_RMSE_TARGET_PCT,
      recommended_n_azimuth = recommended$n_azimuth[1],
      recommended_n_height_levels = recommended$n_height_levels[1],
      recommended_effective_cell_area_m2 = recommended$effective_cell_area_m2[1],
      recommended_flux_mean_mg_m2_d = recommended$flux_mean_mg_m2_d[1],
      recommended_bias_mean_pct = recommended$bias_mean_pct[1],
      recommended_bias_median_pct = recommended$bias_median_pct[1],
      recommended_rmse_pct = recommended$rmse_pct[1],
      recommended_valid_fraction = recommended$valid_fraction[1],
      recommended_mean_dwell_s = recommended$mean_dwell_s[1],
      best_rmse_pct = best_rmse_row$rmse_pct[1],
      best_bias_abs_mean_pct = safe_min_value(dt_class$bias_abs_mean_pct),
      worst_rmse_pct = safe_max_value(dt_class$rmse_pct),
      reference_effective_cell_area_m2 = reference_row$effective_cell_area_m2[1],
      bias_at_reference_pct = reference_row$bias_mean_pct[1],
      rmse_at_reference_pct = reference_row$rmse_pct[1],
      finest_effective_cell_area_m2 = finest_row$effective_cell_area_m2[1],
      finest_rmse_pct = finest_row$rmse_pct[1],
      finest_bias_mean_pct = finest_row$bias_mean_pct[1],
      rmse_improvement_to_finest_pct_points = recommended$rmse_pct[1] - finest_row$rmse_pct[1],
      bias_change_to_finest_pct_points = recommended$bias_mean_pct[1] - finest_row$bias_mean_pct[1],
      rmse_sensitivity_pct_per_log10_area = endpoint_slope(dt_class$effective_cell_area_m2, dt_class$rmse_pct, log_x = TRUE),
      bias_sensitivity_pct_per_log10_area = endpoint_slope(dt_class$effective_cell_area_m2, dt_class$bias_mean_pct, log_x = TRUE),
      rmse_change_halving_area_pct_points = half_step_change_from_lm(dt_class$effective_cell_area_m2, dt_class$rmse_pct, log_x = TRUE),
      bias_change_halving_area_pct_points = half_step_change_from_lm(dt_class$effective_cell_area_m2, dt_class$bias_mean_pct, log_x = TRUE)
    )
  }), fill = TRUE)
}

build_spatial_min_az_table <- function(summary_dt) {
  height_grid <- CJ(
    pasquill_class = sort(unique(summary_dt$pasquill_class)),
    n_height_levels = sort(unique(summary_dt$n_height_levels))
  )
  chosen <- copy(summary_dt[rmse_pct < RESULTS_RMSE_TARGET_PCT])[
    order(pasquill_class, n_height_levels, n_azimuth, -valid_fraction, abs(bias_mean_pct), rmse_pct)
  ][, .SD[1], by = .(pasquill_class, n_height_levels)]

  out <- merge(height_grid, chosen, by = c("pasquill_class", "n_height_levels"), all.x = TRUE)
  out[, .(
    pasquill_class,
    n_height_levels,
    min_azimuth_for_rmse25 = n_azimuth,
    selected_effective_cell_area_m2 = effective_cell_area_m2,
    selected_rmse_pct = rmse_pct,
    selected_bias_mean_pct = bias_mean_pct,
    selected_valid_fraction = valid_fraction
  )][order(pasquill_class, n_height_levels)]
}

build_spatial_min_height_table <- function(summary_dt) {
  az_grid <- CJ(
    pasquill_class = sort(unique(summary_dt$pasquill_class)),
    n_azimuth = sort(unique(summary_dt$n_azimuth))
  )
  chosen <- copy(summary_dt[rmse_pct < RESULTS_RMSE_TARGET_PCT])[
    order(pasquill_class, n_azimuth, n_height_levels, -valid_fraction, abs(bias_mean_pct), rmse_pct)
  ][, .SD[1], by = .(pasquill_class, n_azimuth)]

  out <- merge(az_grid, chosen, by = c("pasquill_class", "n_azimuth"), all.x = TRUE)
  out[, .(
    pasquill_class,
    n_azimuth,
    min_height_levels_for_rmse25 = n_height_levels,
    selected_effective_cell_area_m2 = effective_cell_area_m2,
    selected_rmse_pct = rmse_pct,
    selected_bias_mean_pct = bias_mean_pct,
    selected_valid_fraction = valid_fraction
  )][order(pasquill_class, n_azimuth)]
}

build_nearby_baseline_table <- function(run_dt, true_flux_mg_m2_d) {
  base_dt <- copy(run_dt)
  base_dt[, baseline_bias_pct := 100 * (baseline_flux_retrieved_mg_m2_d - true_flux_mg_m2_d) / true_flux_mg_m2_d]

  base_dt[, {
    baseline_flux <- baseline_flux_retrieved_mg_m2_d[is.finite(baseline_flux_retrieved_mg_m2_d)]
    baseline_bias <- baseline_bias_pct[is.finite(baseline_bias_pct)]
    scenario_seed <- string_to_seed(sprintf("nearby_baseline|%s", pasquill_class[1]))

    flux_stats <- distribution_detail_cols(
      baseline_flux,
      base_name = "clean_baseline_flux",
      suffix = "_mg_m2_d",
      seed = scenario_seed + 5L,
      include_abs_mean = FALSE
    )
    bias_stats <- distribution_detail_cols(
      baseline_bias,
      base_name = "clean_baseline_bias",
      suffix = "_pct",
      seed = scenario_seed + 17L,
      include_abs_mean = TRUE
    )

    c(list(
      clean_baseline_flux_mean_mg_m2_d = if (length(baseline_flux)) mean(baseline_flux, na.rm = TRUE) else NA_real_,
      clean_baseline_flux_sd_mg_m2_d = safe_sd(baseline_flux),
      clean_baseline_bias_mean_pct = if (length(baseline_bias)) mean(baseline_bias, na.rm = TRUE) else NA_real_,
      clean_baseline_bias_sd_pct = safe_sd(baseline_bias)
    ), flux_stats, bias_stats)
  }, by = pasquill_class][order(pasquill_class)]
}

build_nearby_key_table <- function(summary_dt, baseline_dt) {
  class_levels <- sort(unique(summary_dt$pasquill_class))

  rbindlist(lapply(class_levels, function(pc) {
    dt_class <- copy(summary_dt[pasquill_class == pc])[order(distance_m)]
    baseline_row <- baseline_dt[pasquill_class == pc][1]
    nearest <- dt_class[1]
    farthest <- dt_class[.N]
    ref_distance <- nearest_tested_value(dt_class$distance_m)
    ref_row <- pick_reference_row(dt_class, "distance_m", ref_distance, metric_col = "contamination_bias_abs_mean_pct")
    bias10_row <- pick_threshold_row(
      dt_class,
      "distance_m",
      "contamination_bias_mean_pct",
      RESULTS_CONTAM_BIAS_THRESHOLDS_PCT[1],
      direction = "smallest",
      abs_metric = TRUE,
      bias_col = "contamination_bias_mean_pct"
    )
    bias5_row <- pick_threshold_row(
      dt_class,
      "distance_m",
      "contamination_bias_mean_pct",
      RESULTS_CONTAM_BIAS_THRESHOLDS_PCT[2],
      direction = "smallest",
      abs_metric = TRUE,
      bias_col = "contamination_bias_mean_pct"
    )
    fv10_row <- pick_threshold_row(
      dt_class,
      "distance_m",
      "flux_vector_contamination_mean_pct",
      RESULTS_FV_THRESHOLDS_PCT[1],
      direction = "smallest",
      abs_metric = TRUE,
      bias_col = "contamination_bias_mean_pct"
    )

    data.table(
      pasquill_class = pc,
      clean_baseline_flux_mean_mg_m2_d = baseline_row$clean_baseline_flux_mean_mg_m2_d[1],
      clean_baseline_bias_mean_pct = baseline_row$clean_baseline_bias_mean_pct[1],
      nearest_distance_m = nearest$distance_m[1],
      farthest_distance_m = farthest$distance_m[1],
      reference_distance_m = ref_row$distance_m[1],
      retrieved_flux_nearest_mg_m2_d = nearest$retrieved_flux_mean_mg_m2_d[1],
      retrieved_flux_farthest_mg_m2_d = farthest$retrieved_flux_mean_mg_m2_d[1],
      retrieval_rmse_nearest_pct = nearest$retrieval_rmse_pct[1],
      retrieval_rmse_farthest_pct = farthest$retrieval_rmse_pct[1],
      valid_fraction_nearest = nearest$valid_fraction[1],
      valid_fraction_farthest = farthest$valid_fraction[1],
      contamination_bias_nearest_pct = nearest$contamination_bias_mean_pct[1],
      contamination_bias_farthest_pct = farthest$contamination_bias_mean_pct[1],
      contamination_bias_reference_pct = ref_row$contamination_bias_mean_pct[1],
      flux_vector_contamination_reference_pct = ref_row$flux_vector_contamination_mean_pct[1],
      lateral_exit_fraction_reference_pct = ref_row$lateral_exit_fraction_mean_pct[1],
      distance_threshold_bias10_m = bias10_row$distance_m[1],
      distance_threshold_bias5_m = bias5_row$distance_m[1],
      distance_threshold_fv10_m = fv10_row$distance_m[1],
      excess_bias_nearest_pct = nearest$contamination_bias_mean_pct[1],
      contamination_bias_reduction_pct = percent_reduction(
        nearest$contamination_bias_mean_pct[1],
        farthest$contamination_bias_mean_pct[1],
        absolute = TRUE
      ),
      distance_sensitivity_pct_per_m = endpoint_slope(dt_class$distance_m, dt_class$contamination_bias_mean_pct, log_x = FALSE)
    )
  }), fill = TRUE)
}

maybe_write_full_analysis_synthesis <- function(output_dir) {
  paths <- list(
    spatial_key = file.path(output_dir, "spatial_resolution_constant_atmosphere_key_table.csv"),
    temporal_key = file.path(output_dir, "temporal_measurement_changing_atmosphere_key_table.csv"),
    nearby_key = file.path(output_dir, "nearby_lake_contamination_distance_key_table.csv"),
    spatial_summary = file.path(output_dir, "spatial_resolution_constant_atmosphere_summary.csv"),
    temporal_summary = file.path(output_dir, "temporal_measurement_changing_atmosphere_summary.csv"),
    nearby_summary = file.path(output_dir, "nearby_lake_contamination_distance_summary.csv"),
    synthesis_table = file.path(output_dir, "full_analysis_results_synthesis_table.csv")
  )

  required_inputs <- unlist(paths[c("spatial_key", "temporal_key", "nearby_key",
                                    "spatial_summary", "temporal_summary", "nearby_summary")])
  if (!all(file.exists(required_inputs))) {
    return(invisible(FALSE))
  }

  spatial_key <- fread(paths$spatial_key)
  temporal_key <- fread(paths$temporal_key)
  nearby_key <- fread(paths$nearby_key)
  spatial_summary <- fread(paths$spatial_summary)
  temporal_summary <- fread(paths$temporal_summary)
  nearby_summary <- fread(paths$nearby_summary)

  rows <- list()
  row_idx <- 0L

  for (i in seq_len(nrow(spatial_key))) {
    row_idx <- row_idx + 1L
    rows[[row_idx]] <- data.table(
      analysis = "spatial",
      group = spatial_key$pasquill_class[i],
      recommendation = fmt_design(spatial_key$recommended_n_azimuth[i], spatial_key$recommended_n_height_levels[i]),
      threshold_metric = "coarsest design meeting RMSE < 25%",
      threshold_value = spatial_key$recommended_effective_cell_area_m2[i],
      threshold_unit = "m^2",
      mean_bias_pct = spatial_key$recommended_bias_mean_pct[i],
      rmse_pct = spatial_key$recommended_rmse_pct[i],
      valid_fraction = spatial_key$recommended_valid_fraction[i]
    )
  }

  for (i in seq_len(nrow(temporal_key))) {
    temporal_group <- if ("experiment_label" %in% names(temporal_key)) {
      sprintf("%s | %s", temporal_key$transition[i], temporal_key$experiment_label[i])
    } else {
      temporal_key$transition[i]
    }
    temporal_metric <- if ("experiment_family" %in% names(temporal_key) &&
                           identical(temporal_key$experiment_family[i], "repeat_profile")) {
      "maximum full-profile duration meeting RMSE < 25%"
    } else {
      "maximum mission duration meeting RMSE < 25%"
    }
    row_idx <- row_idx + 1L
    rows[[row_idx]] <- data.table(
      analysis = "temporal",
      group = temporal_group,
      recommendation = ifelse(is.finite(temporal_key$max_duration_rmse25_h[i]),
                              fmt_hours(temporal_key$max_duration_rmse25_h[i]),
                              "No tested duration met RMSE < 25%"),
      threshold_metric = temporal_metric,
      threshold_value = temporal_key$max_duration_rmse25_h[i],
      threshold_unit = "h",
      mean_bias_pct = temporal_key$bias_at_longest_pct[i],
      rmse_pct = temporal_key$rmse_at_longest_pct[i],
      valid_fraction = temporal_key$valid_fraction_at_longest[i]
    )
  }

  for (i in seq_len(nrow(nearby_key))) {
    row_idx <- row_idx + 1L
    rows[[row_idx]] <- data.table(
      analysis = "nearby",
      group = nearby_key$pasquill_class[i],
      recommendation = ifelse(is.finite(nearby_key$distance_threshold_bias10_m[i]),
                              fmt_distance(nearby_key$distance_threshold_bias10_m[i]),
                              "No tested distance met |contamination bias| < 10%"),
      threshold_metric = "minimum distance meeting |contamination bias| < 10%",
      threshold_value = nearby_key$distance_threshold_bias10_m[i],
      threshold_unit = "m",
      mean_bias_pct = nearby_key$contamination_bias_nearest_pct[i],
      rmse_pct = nearby_key$retrieval_rmse_nearest_pct[i],
      valid_fraction = nearby_key$valid_fraction_nearest[i]
    )
  }

  spatial_abs_bias <- abs(spatial_summary$bias_mean_pct)
  temporal_abs_bias <- abs(temporal_summary$retrieval_bias_mean_pct)
  nearby_retrieval_abs_bias <- abs(nearby_summary$retrieval_bias_mean_pct)

  retrieval_rmse_all <- c(spatial_summary$rmse_pct,
                          temporal_summary$retrieval_rmse_pct,
                          nearby_summary$retrieval_rmse_pct)
  valid_all <- c(spatial_summary$valid_fraction,
                 temporal_summary$valid_fraction,
                 nearby_summary$valid_fraction)

  row_idx <- row_idx + 1L
  rows[[row_idx]] <- data.table(
    analysis = "overall",
    group = "all",
    recommendation = "Combined capability envelope",
    threshold_metric = "summary",
    threshold_value = NA_real_,
    threshold_unit = "",
    mean_bias_pct = safe_max_value(c(spatial_abs_bias, temporal_abs_bias, nearby_retrieval_abs_bias)),
    rmse_pct = safe_max_value(retrieval_rmse_all),
    valid_fraction = safe_min_value(valid_all)
  )

  synthesis_dt <- rbindlist(rows, fill = TRUE)
  fwrite(synthesis_dt, paths$synthesis_table)
  invisible(TRUE)
}

