#!/usr/bin/env Rscript
###############################################################################
# analyze_dimensionless_framework.R
#
# Full-analysis sweep for the dimensionless UAS cylindrical mass-balance
# framework. The experiment maps transitional Gaussian-plume retrieval
# performance into Pi_S / Pi_G space while holding the temporal-analysis
# receptor layout fixed.
#
# Experiment v1:
#   - Main transitions: B-D, C-D, and C-F
#   - Supplemental reverse-direction transitions: F-D and D-B
#   - Wind speed is swept; wind direction can follow the temporal AmeriFlux
#     variability model and is stochastic by default
#   - Receptors: design_96az_many_heights_30min
#   - Single start-of-transition flight durations: 30, 60, 90 min
#   - Default cylinder top height: 100 m, using the temporal-analysis
#     azimuth count and number of height levels
#   - One wind/noise realization per dimensional scenario; no Monte Carlo repeats
#   - Ring radius is always outside the lake
#   - Pi_S is theoretical, using the 1 Hz sensor precision without the
#     practical precision floor
#   - Primary binning uses signal-weighted transitional Pi values by default
#
# Outputs (under output/):
#   - dimensionless_framework_runs.csv
#   - dimensionless_framework_base_timings.csv
#   - dimensionless_framework_binned_summary.csv
#   - dimensionless_framework_composite.png
#   - dimensionless_framework_additional_transitions_2x2.png
#   - dimensionless_framework_report.md
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
suppressMessages(source(file.path(proj_root, "R", "04_conley_inverse_flux.R")))
suppressMessages(source(file.path(proj_root, "R", "06_visualization.R")))

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

log_grid <- function(min_val, max_val, n) {
  if (n <= 1L) return(min_val)
  exp(seq(log(min_val), log(max_val), length.out = n))
}

format_num_vec <- function(x) {
  paste(format(x, trim = TRUE, scientific = FALSE), collapse = ", ")
}

cli <- parse_named_args(commandArgs(trailingOnly = TRUE))

PILOT <- parse_bool(cli$pilot)
SENSOR_NAME <- if (!is.null(cli$sensor)) cli$sensor else "li7810"
DESIGN_NAME <- if (!is.null(cli$design)) cli$design else "design_96az_many_heights_30min"
TRANSITION_DURATION_MIN <- if (!is.null(cli$transition_duration_min)) as.numeric(cli$transition_duration_min) else 720
if (!is.null(cli$observation_window_min)) {
  warning("--observation_window_min is ignored; each mission now uses one start-of-transition retrieval window.")
}
MISSION_DURATIONS_MIN <- {
  vals <- parse_num_vec(cli$mission_durations_min)
  if (length(vals)) vals else c(30, 60, 90)
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
FLUX_GRID_MG_M2_D <- {
  vals <- parse_num_vec(cli$flux_grid_mg_m2_d)
  if (length(vals)) vals else signif(log_grid(5, 500, 15), 6)
}
LAKE_RADII_M <- {
  vals <- parse_num_vec(cli$lake_radii_m)
  if (length(vals)) vals else signif(log_grid(10, 1000, 15), 6)
}
RING_RADIUS_RATIOS <- {
  vals <- parse_num_vec(cli$ring_radius_ratios)
  if (length(vals)) vals else c(1.01, 1.03, 1.05, 1.075, 1.1, 1.15, 1.25, 1.5, 2, 3)
}
TOP_HEIGHT_M <- if (!is.null(cli$top_height_m)) as.numeric(cli$top_height_m) else 100
HEIGHT_MODE <- if (!is.null(cli$height_mode)) cli$height_mode else "uniform_top"
U_REF_GRID_M_S <- {
  vals <- parse_num_vec(cli$u_ref_m_s)
  if (length(vals)) vals else seq(2, 5, by = 0.5)
}
WD_DEG <- if (!is.null(cli$wd_deg)) as.numeric(cli$wd_deg) else 225
STOCHASTIC_WIND <- if (!is.null(cli$stochastic_wind)) parse_bool(cli$stochastic_wind) else TRUE
MET_DT_MIN <- if (!is.null(cli$met_dt_min)) as.numeric(cli$met_dt_min) else 30
SAMPLE_RATE_HZ <- if (!is.null(cli$sample_rate_hz)) as.numeric(cli$sample_rate_hz) else 1
SEED_BASE <- if (!is.null(cli$seed)) as.integer(cli$seed) else 20260414L
N_PI_S_BINS <- if (!is.null(cli$n_pi_s_bins)) as.integer(cli$n_pi_s_bins) else 20L
N_PI_G_BINS <- if (!is.null(cli$n_pi_g_bins)) as.integer(cli$n_pi_g_bins) else 20L
MIN_BIN_N_FOR_SPREAD_PLOT <- if (!is.null(cli$min_bin_n_for_spread)) {
  as.integer(cli$min_bin_n_for_spread)
} else if (!is.null(cli$min_bin_n)) {
  as.integer(cli$min_bin_n)
} else {
  5L
}
PI_SUMMARY <- if (!is.null(cli$pi_summary)) cli$pi_summary else "weighted"
N_WORKERS_REQUESTED <- if (!is.null(cli$n_workers)) as.integer(cli$n_workers) else NA_integer_
BATCH_SIZE_REQUESTED <- if (!is.null(cli$batch_size)) as.integer(cli$batch_size) else NA_integer_
WORKER_SPARE_CORES <- if (!is.null(cli$worker_spare_cores)) as.integer(cli$worker_spare_cores) else 2L
OUTPUT_STEM <- if (!is.null(cli$output_stem) && nzchar(cli$output_stem)) {
  cli$output_stem
} else if (isTRUE(PILOT)) {
  "dimensionless_framework_pilot"
} else {
  "dimensionless_framework"
}

stopifnot(SENSOR_NAME %in% names(SENSOR_CONFIGS))
stopifnot(DESIGN_NAME %in% names(SAMPLING_CONFIGS))
stopifnot(all(vapply(TRANSITIONS, function(tr) {
  length(tr) == 2L && all(tr %in% names(PASQUILL_CLASSES))
}, logical(1))))
stopifnot(all(is.finite(FLUX_GRID_MG_M2_D)), all(FLUX_GRID_MG_M2_D > 0))
stopifnot(all(is.finite(LAKE_RADII_M)), all(LAKE_RADII_M > 0))
stopifnot(all(is.finite(RING_RADIUS_RATIOS)), all(RING_RADIUS_RATIOS > 1))
stopifnot(all(is.finite(MISSION_DURATIONS_MIN)), all(MISSION_DURATIONS_MIN > 0))
stopifnot(is.finite(TOP_HEIGHT_M), TOP_HEIGHT_M > 0)
stopifnot(HEIGHT_MODE %in% c("uniform_top", "native", "scaled_native"))
stopifnot(is.finite(TRANSITION_DURATION_MIN), TRANSITION_DURATION_MIN > 0)
stopifnot(all(MISSION_DURATIONS_MIN <= TRANSITION_DURATION_MIN + 1e-9))
stopifnot(all(is.finite(U_REF_GRID_M_S)), all(U_REF_GRID_M_S > 0))
stopifnot(is.finite(SAMPLE_RATE_HZ), SAMPLE_RATE_HZ > 0)
stopifnot(PI_SUMMARY %in% c("mean", "weighted", "start", "mid", "end"))

OUTPUT_DIR <- proj_path("output")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

BESTFIT_PATH <- file.path(OUTPUT_DIR, "ameriflux_wind_variability_bestfit.csv")
if (isTRUE(STOCHASTIC_WIND) && !file.exists(BESTFIT_PATH)) {
  stop("AmeriFlux best-fit CSV not found: ", BESTFIT_PATH,
       "\nRun process_ameriflux_wind_variability.R first or pass --stochastic_wind=FALSE.")
}

RUNS_CSV <- file.path(OUTPUT_DIR, sprintf("%s_runs.csv", OUTPUT_STEM))
TIMINGS_CSV <- file.path(OUTPUT_DIR, sprintf("%s_base_timings.csv", OUTPUT_STEM))
BINNED_CSV <- file.path(OUTPUT_DIR, sprintf("%s_binned_summary.csv", OUTPUT_STEM))
PI_DIAGNOSTIC_CSV <- file.path(OUTPUT_DIR, sprintf("%s_pi_diagnostic_summary.csv", OUTPUT_STEM))
FIG_PNG <- file.path(OUTPUT_DIR, sprintf("%s_composite.png", OUTPUT_STEM))
ADDITIONAL_FIG_PNG <- file.path(OUTPUT_DIR, sprintf("%s_additional_transitions_2x2.png", OUTPUT_STEM))
REPORT_MD <- file.path(OUTPUT_DIR, sprintf("%s_report.md", OUTPUT_STEM))


# -----------------------------------------------------------------------------
# 2. Helpers
# -----------------------------------------------------------------------------
PASQUILL_ZETA_ANCHORS <- c(A = -2.0, B = -1.0, C = -0.5,
                           D =  0.0, E =  0.5, F =  2.0)

safe_sd <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) <= 1L) return(NA_real_)
  stats::sd(x)
}

safe_mean <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  mean(x)
}

safe_min <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  min(x)
}

safe_max <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  max(x)
}

safe_quantile <- function(x, probs) {
  x <- x[is.finite(x)]
  if (!length(x)) return(rep(NA_real_, length(probs)))
  as.numeric(stats::quantile(x, probs = probs, na.rm = TRUE, names = FALSE))
}

safe_fraction <- function(x) {
  x <- x[!is.na(x)]
  if (!length(x)) return(NA_real_)
  mean(x)
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
                                   u_mean_m_s, wd_mean_deg,
                                   bestfit_interp, class_start, class_end,
                                   seed) {
  set.seed(seed)
  boundaries <- seq(0, transition_duration_min, by = met_dt_min)
  if (tail(boundaries, 1) < transition_duration_min) {
    boundaries <- c(boundaries, transition_duration_min)
  }

  rows <- lapply(boundaries, function(t_now) {
    zeta_now <- zeta_at_time(t_now, transition_duration_min, class_start, class_end)
    sigma_u_now <- bestfit_interp$sigma_u(zeta_now)
    sigma_theta_now <- bestfit_interp$sigma_theta(zeta_now)
    data.table(
      time_min = t_now,
      zeta = zeta_now,
      sigma_u = sigma_u_now,
      sigma_theta = sigma_theta_now,
      # Wind speed is intentionally fixed here; sigma_u is retained as metadata.
      ws = u_mean_m_s,
      wd = (wd_mean_deg + rnorm(1, 0, sigma_theta_now)) %% 360
    )
  })

  rbindlist(rows)
}

angle_diff_deg <- function(to_deg, from_deg) {
  ((to_deg - from_deg + 180) %% 360) - 180
}

circular_mean_deg <- function(deg, weights = NULL) {
  keep <- is.finite(deg)
  deg <- deg[keep]
  if (!length(deg)) return(NA_real_)
  if (is.null(weights)) {
    weights <- rep(1, length(deg))
  } else {
    weights <- weights[keep]
    if (length(weights) != length(deg)) return(NA_real_)
    weights[!is.finite(weights) | weights < 0] <- 0
    if (sum(weights) <= 0) weights <- rep(1, length(deg))
  }
  rad <- deg * pi / 180
  mean_rad <- atan2(sum(weights * sin(rad)), sum(weights * cos(rad)))
  (mean_rad * 180 / pi) %% 360
}

circular_sd_deg <- function(deg) {
  deg <- deg[is.finite(deg)]
  if (length(deg) <= 1L) return(NA_real_)
  rad <- deg * pi / 180
  r_bar <- sqrt(mean(cos(rad))^2 + mean(sin(rad))^2)
  if (!is.finite(r_bar) || r_bar <= 0) return(NA_real_)
  sqrt(-2 * log(pmin(r_bar, 1))) * 180 / pi
}

wd_at_time <- function(t_min, wind_tl) {
  wd <- wind_tl$wd
  if (length(wd) <= 1L) return(rep(wd[1], length(t_min)))
  wd_unwrapped <- c(wd[1], wd[1] + cumsum(angle_diff_deg(wd[-1], wd[-length(wd)])))
  approx(wind_tl$time_min, wd_unwrapped, xout = t_min, rule = 2)$y %% 360
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

interpolate_briggs_table <- function(t_min, transition_duration_min, class_start, class_end) {
  w <- pmin(pmax(t_min / transition_duration_min, 0), 1)
  cs <- PASQUILL_CLASSES[[class_start]]
  ce <- PASQUILL_CLASSES[[class_end]]
  data.table(
    sy_a = cs$sy_a + w * (ce$sy_a - cs$sy_a),
    sy_b = cs$sy_b + w * (ce$sy_b - cs$sy_b),
    sz_a = cs$sz_a + w * (ce$sz_a - cs$sz_a),
    sz_b = cs$sz_b + w * (ce$sz_b - cs$sz_b),
    sz_p = cs$sz_p + w * (ce$sz_p - cs$sz_p),
    wind_exp = cs$wind_exp + w * (ce$wind_exp - cs$wind_exp),
    zeta = zeta_at_time(t_min, transition_duration_min, class_start, class_end)
  )
}

build_disk_quadrature <- function(source_radius, flux_ug_m2_s,
                                  n_rings = 20L, n_src_az = 36L) {
  r_edges <- seq(0, source_radius, length.out = n_rings + 1L)
  r_mids <- (r_edges[-1L] + r_edges[-(n_rings + 1L)]) / 2
  theta_mids <- seq(0, 2 * pi * (1 - 1 / n_src_az), length.out = n_src_az)
  d_theta <- 2 * pi / n_src_az

  rows <- vector("list", n_rings)
  for (kr in seq_len(n_rings)) {
    r_lo <- r_edges[kr]
    r_hi <- r_edges[kr + 1L]
    r_mid <- r_mids[kr]
    dA_cell <- (r_hi^2 - r_lo^2) / 2 * d_theta
    rows[[kr]] <- data.table(
      sx = r_mid * cos(theta_mids),
      sy = r_mid * sin(theta_mids),
      q_cell = flux_ug_m2_s * dA_cell
    )
  }
  rbindlist(rows)
}

gaussian_plume_disk_timevarying <- function(rec_x, rec_y, rec_z,
                                            source_radius, flux_ug_m2_s,
                                            u_ref, wd_deg, coefs_dt,
                                            n_rings = 20L, n_src_az = 36L,
                                            chunk_size = 512L) {
  n_rec <- length(rec_x)
  if (length(wd_deg) == 1L) {
    wd_deg <- rep(wd_deg, n_rec)
  }
  if (length(wd_deg) != n_rec) {
    stop("wd_deg must be scalar or have one value per receptor.")
  }

  conc <- numeric(n_rec)
  quad <- build_disk_quadrature(
    source_radius = source_radius,
    flux_ug_m2_s = flux_ug_m2_s,
    n_rings = n_rings,
    n_src_az = n_src_az
  )

  wd_to_rad <- ((wd_deg + 180) %% 360) * pi / 180
  wind_dx <- sin(wd_to_rad)
  wind_dy <- cos(wd_to_rad)
  src_n <- nrow(quad)

  starts <- seq(1L, n_rec, by = chunk_size)
  for (chunk_start in starts) {
    idx <- chunk_start:min(chunk_start + chunk_size - 1L, n_rec)
    nr <- length(idx)

    dx <- outer(rec_x[idx], quad$sx, "-")
    dy <- outer(rec_y[idx], quad$sy, "-")
    x_down <- sweep(dx, 1L, wind_dx[idx], "*") + sweep(dy, 1L, wind_dy[idx], "*")
    y_cross <- sweep(dx, 1L, wind_dy[idx], "*") - sweep(dy, 1L, wind_dx[idx], "*")

    valid <- x_down > 1e-3
    xv <- x_down
    xv[!valid] <- 1e-3

    sy_mat <- pmax(coefs_dt$sy_a[idx] * xv / sqrt(1 + coefs_dt$sy_b[idx] * xv), 0.01)
    sz_mat <- pmax(coefs_dt$sz_a[idx] * xv * (1 + coefs_dt$sz_b[idx] * xv)^coefs_dt$sz_p[idx], 0.01)

    z_mat <- matrix(rec_z[idx], nrow = nr, ncol = src_n)
    term <- (1 / (pi * sy_mat * sz_mat * u_ref)) *
      exp(-y_cross^2 / (2 * sy_mat^2)) *
      exp(-z_mat^2 / (2 * sz_mat^2))
    term[!valid] <- 0
    term <- sweep(term, 2L, quad$q_cell, "*")
    conc[idx] <- rowSums(term)
  }

  conc
}

build_receptors_for_ring <- function(lake_radius_m, ring_radius_m, design_cfg) {
  define_receptor_config(
    lake_diameter = 2 * lake_radius_m,
    n_azimuth = design_cfg$n_azimuth,
    heights_m = design_cfg$heights_m,
    distance_from_edge_m = ring_radius_m - lake_radius_m,
    inner_ring_fracs = NULL,
    background_heights_m = NULL
  )
}

compute_pi_at_times <- function(sample_times_min,
                                flux_mg_m2_d,
                                lake_radius_m,
                                ring_radius_m,
                                top_height_m,
                                mission_duration_min,
                                transition_duration_min,
                                class_start,
                                class_end,
                                u_ref_m_s,
                                sensor_cfg,
                                sample_rate_hz,
                                sample_weights = NULL) {
  sigma1_ug_m3 <- ppb_to_ug_m3(sensor_cfg$precision_1s_ppb)
  flux_ug_m2_s <- flux_mg_d_to_ug_s(flux_mg_m2_d)
  tau_flight_s <- mission_duration_min * 60

  one_time <- function(t_min) {
    coefs <- interpolate_briggs(t_min, transition_duration_min, class_start, class_end)
    sigma_y <- pg_sigma_y_interp(ring_radius_m, coefs)
    sigma_z <- pg_sigma_z_interp(ring_radius_m, coefs)
    pi_s <- (flux_ug_m2_s * lake_radius_m^2 * sqrt(tau_flight_s * sample_rate_hz)) /
      (2 * sqrt(2 * pi) * ring_radius_m * u_ref_m_s * sigma_z * sigma1_ug_m3)
    pi_g <- (lake_radius_m + sigma_y) / ring_radius_m
    data.table(
      sigma_y_m = sigma_y,
      sigma_z_m = sigma_z,
      pi_s = pi_s,
      pi_g = pi_g
    )
  }

  key_times <- c(start = 0, mid = mission_duration_min / 2, end = mission_duration_min)
  key_dt <- rbindlist(lapply(key_times, one_time), idcol = "time_label")
  key_dt[, time_label := names(key_times)]

  sample_dt <- rbindlist(lapply(sample_times_min, one_time))
  sample_mean <- sample_dt[, .(
    sigma_y_m = mean(sigma_y_m, na.rm = TRUE),
    sigma_z_m = mean(sigma_z_m, na.rm = TRUE),
    pi_s = mean(pi_s, na.rm = TRUE),
    pi_g = mean(pi_g, na.rm = TRUE)
  )]
  if (is.null(sample_weights)) {
    sample_weights <- rep(1, nrow(sample_dt))
  } else {
    sample_weights <- as.numeric(sample_weights)
    sample_weights[!is.finite(sample_weights) | sample_weights < 0] <- 0
    if (length(sample_weights) != nrow(sample_dt)) {
      stop("sample_weights must match sample_times_min length.")
    }
    if (sum(sample_weights) <= 0) sample_weights <- rep(1, nrow(sample_dt))
  }
  sample_weighted <- data.table(
    sigma_y_m = stats::weighted.mean(sample_dt$sigma_y_m, sample_weights, na.rm = TRUE),
    sigma_z_m = stats::weighted.mean(sample_dt$sigma_z_m, sample_weights, na.rm = TRUE),
    pi_s = stats::weighted.mean(sample_dt$pi_s, sample_weights, na.rm = TRUE),
    pi_g = stats::weighted.mean(sample_dt$pi_g, sample_weights, na.rm = TRUE)
  )

  out <- data.table(
    pi_integration_min = mission_duration_min,
    pi_s_start = key_dt[time_label == "start", pi_s],
    pi_g_start = key_dt[time_label == "start", pi_g],
    pi_h_start = top_height_m / key_dt[time_label == "start", sigma_z_m],
    sigma_y_start_m = key_dt[time_label == "start", sigma_y_m],
    sigma_z_start_m = key_dt[time_label == "start", sigma_z_m],
    pi_s_mid = key_dt[time_label == "mid", pi_s],
    pi_g_mid = key_dt[time_label == "mid", pi_g],
    pi_h_mid = top_height_m / key_dt[time_label == "mid", sigma_z_m],
    sigma_y_mid_m = key_dt[time_label == "mid", sigma_y_m],
    sigma_z_mid_m = key_dt[time_label == "mid", sigma_z_m],
    pi_s_end = key_dt[time_label == "end", pi_s],
    pi_g_end = key_dt[time_label == "end", pi_g],
    pi_h_end = top_height_m / key_dt[time_label == "end", sigma_z_m],
    sigma_y_end_m = key_dt[time_label == "end", sigma_y_m],
    sigma_z_end_m = key_dt[time_label == "end", sigma_z_m],
    pi_s_mean = sample_mean$pi_s,
    pi_g_mean = sample_mean$pi_g,
    pi_h_mean = top_height_m / sample_mean$sigma_z_m,
    sigma_y_mean_m = sample_mean$sigma_y_m,
    sigma_z_mean_m = sample_mean$sigma_z_m,
    pi_s_weighted = sample_weighted$pi_s,
    pi_g_weighted = sample_weighted$pi_g,
    pi_h_weighted = top_height_m / sample_weighted$sigma_z_m,
    sigma_y_weighted_m = sample_weighted$sigma_y_m,
    sigma_z_weighted_m = sample_weighted$sigma_z_m
  )
  out[]
}

build_base_profile_truth <- function(mission_duration_min,
                                     design_cfg,
                                     lake_radius_m,
                                     ring_radius_m,
                                     transition_duration_min,
                                     class_start,
                                     class_end,
                                     u_ref_m_s,
                                     wd_deg,
                                     wind_tl = NULL) {
  receptors <- build_receptors_for_ring(lake_radius_m, ring_radius_m, design_cfg)
  setorder(receptors, z, azimuth_deg)

  n_samples <- nrow(receptors)
  dwell_min <- mission_duration_min / n_samples
  dwell_s <- dwell_min * 60
  unit_flux_ug_s <- flux_mg_d_to_ug_s(1)
  sample_times_min <- (seq_len(n_samples) - 1L) * dwell_min

  wd_vec <- if (is.null(wind_tl)) {
    rep(wd_deg, length(sample_times_min))
  } else {
    wd_at_time(sample_times_min, wind_tl)
  }
  coefs_dt <- interpolate_briggs_table(
    t_min = sample_times_min,
    transition_duration_min = transition_duration_min,
    class_start = class_start,
    class_end = class_end
  )
  enh_ug_m3 <- gaussian_plume_disk_timevarying(
    rec_x = receptors$x,
    rec_y = receptors$y,
    rec_z = receptors$z,
    source_radius = lake_radius_m,
    flux_ug_m2_s = unit_flux_ug_s,
    u_ref = u_ref_m_s,
    wd_deg = wd_vec,
    coefs_dt = coefs_dt
  )

  base_dt <- data.table(
    receptor_id = receptors$receptor_id,
    x = receptors$x,
    y = receptors$y,
    z = receptors$z,
    azimuth_deg = receptors$azimuth_deg,
    ring_type = receptors$ring_type,
    time_min = 0,
    sample_time_min = sample_times_min,
    ch4_enhancement_ppb_per_mg_m2_d = ug_m3_to_ppb(enh_ug_m3),
    dwell_s = dwell_s,
    ws_m_s = u_ref_m_s,
    ws_at_z_m_s = u_ref_m_s,
    wd_deg = wd_vec,
    ustar_m_s = NA_real_,
    L_m = NA_real_,
    z0_m = 0.001,
    zeta = coefs_dt$zeta,
    class_start = class_start,
    class_end = class_end,
    stability_class = "transitional"
  )

  list(
    base_dt = base_dt,
    receptors = receptors,
    n_samples_per_profile = n_samples,
    total_profile_time_min = mission_duration_min
  )
}

build_measurements_for_flux <- function(base_dt, flux_mg_m2_d, sensor_cfg, seed) {
  set.seed(seed)
  dt <- copy(base_dt)
  dwell_s <- unique(dt$dwell_s)[1]
  precision_at_dwell <- sensor_cfg$precision_1s_ppb / sqrt(max(dwell_s, 1))
  precision_at_dwell <- max(precision_at_dwell, sensor_cfg$precision_opt_ppb)
  noise_ppb <- rnorm(nrow(dt), 0, precision_at_dwell)

  enh_ppb <- dt$ch4_enhancement_ppb_per_mg_m2_d * flux_mg_m2_d
  dt[, `:=`(
    ch4_bg_ppb = 0,
    ch4_enhancement_ppb = round(enh_ppb, 4),
    ch4_contamination_ppb = 0,
    ch4_total_ppb = round(enh_ppb, 4),
    ch4_measured_ppb = round(enh_ppb + noise_ppb, 4),
    noise_ppb = round(noise_ppb, 4),
    drift_ppb = 0,
    offset_ppb = 0,
    precision_ppb = round(precision_at_dwell, 4),
    flux_true_mg_m2_d = flux_mg_m2_d
  )]
  dt[]
}

run_flux_retrieval <- function(base_dt,
                               flux_mg_m2_d,
                               sensor_cfg,
                               seed,
                               lake_radius_m,
                               ring_radius_m,
                               n_azimuth,
                               scenario_label) {
  meas_dt <- build_measurements_for_flux(
    base_dt = base_dt,
    flux_mg_m2_d = flux_mg_m2_d,
    sensor_cfg = sensor_cfg,
    seed = seed
  )

  flux_dt <- quietly(retrieve_flux(
    meas_dt = meas_dt,
    lake_diameter = 2 * lake_radius_m,
    distance_from_edge_m = ring_radius_m - lake_radius_m,
    n_azimuth = n_azimuth,
    use_observed_ws = TRUE,
    ws_observation_error_m_s = 0,
    wd_uncertainty_deg = 0,
    use_height_avg_wind = TRUE,
    scenario_label = scenario_label,
    save_output = FALSE
  ))

  list(measurements = meas_dt, flux = flux_dt)
}

summarise_retrieval_windows <- function(flux_dt) {
  valid_idx <- flux_dt$qc_flag == "ok" & is.finite(flux_dt$flux_retrieved)
  valid_flux <- flux_dt$flux_retrieved[valid_idx]
  valid_enh <- flux_dt$mean_enhancement_ppb[valid_idx]
  valid_ws <- flux_dt$ws_at_z_mean_m_s[valid_idx]

  qc_summary <- if (any(valid_idx)) {
    "ok"
  } else if (nrow(flux_dt) && all(flux_dt$qc_flag == "low_ws")) {
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

make_full_base_grid <- function() {
  transition_dt <- rbindlist(lapply(TRANSITIONS, function(tr) {
    data.table(
      class_start = tr[1],
      class_end = tr[2],
      transition = paste(tr, collapse = "-")
    )
  }))

  CJ(
    transition = transition_dt$transition,
    lake_radius_m = sort(unique(LAKE_RADII_M)),
    ring_radius_ratio = sort(unique(RING_RADIUS_RATIOS)),
    mission_duration_min = sort(unique(MISSION_DURATIONS_MIN)),
    u_ref_m_s = sort(unique(U_REF_GRID_M_S))
  )[transition_dt, on = "transition"][order(transition, lake_radius_m, ring_radius_ratio,
                                            mission_duration_min, u_ref_m_s)]
}

make_pilot_base_grid <- function(full_base_grid) {
  candidates <- rbindlist(lapply(TRANSITION_LABELS, function(tr) {
    rbindlist(list(
      full_base_grid[
        transition == tr &
          lake_radius_m == min(LAKE_RADII_M) &
          ring_radius_ratio == min(RING_RADIUS_RATIOS) &
          mission_duration_min == min(MISSION_DURATIONS_MIN) &
          u_ref_m_s == min(U_REF_GRID_M_S)
      ],
      full_base_grid[
        transition == tr &
          lake_radius_m == max(LAKE_RADII_M) &
          ring_radius_ratio == max(RING_RADIUS_RATIOS) &
          mission_duration_min == max(MISSION_DURATIONS_MIN) &
          u_ref_m_s == max(U_REF_GRID_M_S)
      ]
    ), fill = TRUE)
  }), fill = TRUE)
  unique(candidates)
}

build_bins <- function(x, n_bins, log10_scale = FALSE) {
  vals <- if (log10_scale) log10(x[is.finite(x) & x > 0]) else x[is.finite(x)]
  if (!length(vals)) stop("Cannot build bins: no finite values.")
  rng <- range(vals, na.rm = TRUE)
  span <- diff(rng)
  scale <- max(1, max(abs(rng), na.rm = TRUE))
  if (!is.finite(span) || span <= sqrt(.Machine$double.eps) * scale) {
    pad <- 0.05 * scale
    rng <- rng + c(-pad, pad)
  }
  breaks <- seq(rng[1], rng[2], length.out = n_bins + 1L)
  if (anyDuplicated(breaks)) {
    center <- mean(rng)
    pad <- 0.1 * scale
    breaks <- seq(center - pad, center + pad, length.out = n_bins + 1L)
  }
  breaks
}

pi_column_names <- function(pi_summary) {
  list(
    pi_s = paste0("pi_s_", pi_summary),
    pi_g = paste0("pi_g_", pi_summary),
    pi_h = paste0("pi_h_", pi_summary)
  )
}

summarise_binned_results <- function(run_dt, n_pi_s_bins, n_pi_g_bins, pi_summary = "weighted") {
  dt <- copy(run_dt)
  pi_cols <- pi_column_names(pi_summary)
  if (!all(unlist(pi_cols) %in% names(dt))) {
    stop("Pi summary columns not found for pi_summary = ", pi_summary)
  }
  dt[, pi_s_for_bin := get(pi_cols$pi_s)]
  dt[, pi_g_for_bin := get(pi_cols$pi_g)]
  dt[, pi_h_for_summary := get(pi_cols$pi_h)]
  dt <- dt[is.finite(pi_s_for_bin) & pi_s_for_bin > 0 & is.finite(pi_g_for_bin)]
  pi_s_breaks <- build_bins(dt$pi_s_for_bin, n_pi_s_bins, log10_scale = TRUE)
  pi_g_breaks <- build_bins(dt$pi_g_for_bin, n_pi_g_bins, log10_scale = FALSE)

  dt[, pi_s_bin := cut(log10(pi_s_for_bin), breaks = pi_s_breaks, include.lowest = TRUE, labels = FALSE)]
  dt[, pi_g_bin := cut(pi_g_for_bin, breaks = pi_g_breaks, include.lowest = TRUE, labels = FALSE)]
  dt <- dt[!is.na(pi_s_bin) & !is.na(pi_g_bin)]

  out <- dt[, {
    valid_vals <- qc_flag == "ok" & is.finite(flux_retrieved_mg_m2_d)
    bias_vals <- retrieval_bias_pct[valid_vals]
    abs_error_vals <- abs(bias_vals)
    abs_error_q <- safe_quantile(abs_error_vals, c(0.25, 0.5, 0.75, 0.90, 0.95))
    success_lt10 <- safe_fraction(abs_error_vals < 10)
    success_lt25 <- safe_fraction(abs_error_vals < 25)
    success_lt50 <- safe_fraction(abs_error_vals < 50)
    physical_error_vals <- flux_error_mg_m2_d[valid_vals]
    data.table(
      bin_n = .N,
      n_valid = sum(valid_vals),
      valid_fraction = mean(valid_vals),
      pi_s_center = 10^mean(range(pi_s_breaks[pi_s_bin[1]:(pi_s_bin[1] + 1L)])),
      pi_g_center = mean(range(pi_g_breaks[pi_g_bin[1]:(pi_g_bin[1] + 1L)])),
      pi_s_min = 10^pi_s_breaks[pi_s_bin[1]],
      pi_s_max = 10^pi_s_breaks[pi_s_bin[1] + 1L],
      pi_g_min = pi_g_breaks[pi_g_bin[1]],
      pi_g_max = pi_g_breaks[pi_g_bin[1] + 1L],
      bias_mean_pct = safe_mean(bias_vals),
      bias_median_pct = safe_quantile(bias_vals, 0.5),
      bias_sd_pct = safe_sd(bias_vals),
      rmse_pct = if (length(bias_vals)) sqrt(mean(bias_vals^2, na.rm = TRUE)) else NA_real_,
      abs_error_mean_pct = safe_mean(abs_error_vals),
      abs_error_median_pct = abs_error_q[2],
      abs_error_p25_pct = abs_error_q[1],
      abs_error_p75_pct = abs_error_q[3],
      abs_error_p90_pct = abs_error_q[4],
      abs_error_p95_pct = abs_error_q[5],
      abs_error_iqr_pct = abs_error_q[3] - abs_error_q[1],
      abs_error_sd_pct = safe_sd(abs_error_vals),
      success_abs_error_lt10_frac = success_lt10,
      success_abs_error_lt10_pct = 100 * success_lt10,
      success_abs_error_lt25_frac = success_lt25,
      success_abs_error_lt25_pct = 100 * success_lt25,
      success_abs_error_lt50_frac = success_lt50,
      success_abs_error_lt50_pct = 100 * success_lt50,
      rmse_mg_m2_d = if (length(physical_error_vals)) sqrt(mean(physical_error_vals^2, na.rm = TRUE)) else NA_real_,
      abs_error_sd_mg_m2_d = safe_sd(abs(physical_error_vals)),
      flux_true_min_mg_m2_d = min(flux_true_mg_m2_d, na.rm = TRUE),
      flux_true_max_mg_m2_d = max(flux_true_mg_m2_d, na.rm = TRUE),
      pi_h_median = median(pi_h_for_summary, na.rm = TRUE),
      pi_h_min = min(pi_h_for_summary, na.rm = TRUE),
      pi_h_max = max(pi_h_for_summary, na.rm = TRUE),
      n_retrieval_windows = sum(n_retrieval_windows, na.rm = TRUE),
      n_valid_retrieval_windows = sum(n_valid_retrieval_windows, na.rm = TRUE)
    )
  }, by = .(transition, pi_s_bin, pi_g_bin)][order(transition, pi_s_bin, pi_g_bin)]

  out[, pi_summary := pi_summary]

  attr(out, "pi_s_breaks_log10") <- pi_s_breaks
  attr(out, "pi_g_breaks") <- pi_g_breaks
  out[]
}

summarise_performance_stats <- function(dt, pi_s_col, pi_g_col, pi_h_col) {
  valid_vals <- dt$qc_flag == "ok" & is.finite(dt$flux_retrieved_mg_m2_d) &
    is.finite(dt$retrieval_bias_pct)
  bias_vals <- dt$retrieval_bias_pct[valid_vals]
  abs_error_vals <- abs(bias_vals)
  physical_error_vals <- dt$flux_error_mg_m2_d[valid_vals]
  abs_error_q <- safe_quantile(abs_error_vals, c(0.25, 0.5, 0.75, 0.90, 0.95))
  bias_q <- safe_quantile(bias_vals, c(0.25, 0.5, 0.75))
  success_lt10 <- safe_fraction(abs_error_vals < 10)
  success_lt25 <- safe_fraction(abs_error_vals < 25)
  success_lt50 <- safe_fraction(abs_error_vals < 50)

  pi_s_vals <- dt[[pi_s_col]]
  pi_g_vals <- dt[[pi_g_col]]
  pi_h_vals <- dt[[pi_h_col]]
  u_ref_vals <- if ("u_ref_m_s" %in% names(dt)) dt$u_ref_m_s else NA_real_
  mission_vals <- if ("mission_duration_min" %in% names(dt)) dt$mission_duration_min else NA_real_

  data.table(
    n_scenarios = nrow(dt),
    n_valid = sum(valid_vals, na.rm = TRUE),
    valid_fraction = safe_fraction(valid_vals),
    rmse_pct = if (length(bias_vals)) sqrt(mean(bias_vals^2, na.rm = TRUE)) else NA_real_,
    bias_mean_pct = safe_mean(bias_vals),
    bias_median_pct = bias_q[2],
    bias_p25_pct = bias_q[1],
    bias_p75_pct = bias_q[3],
    bias_sd_pct = safe_sd(bias_vals),
    abs_error_mean_pct = safe_mean(abs_error_vals),
    abs_error_median_pct = abs_error_q[2],
    abs_error_p25_pct = abs_error_q[1],
    abs_error_p75_pct = abs_error_q[3],
    abs_error_p90_pct = abs_error_q[4],
    abs_error_p95_pct = abs_error_q[5],
    abs_error_iqr_pct = abs_error_q[3] - abs_error_q[1],
    abs_error_sd_pct = safe_sd(abs_error_vals),
    success_abs_error_lt10_frac = success_lt10,
    success_abs_error_lt10_pct = 100 * success_lt10,
    success_abs_error_lt25_frac = success_lt25,
    success_abs_error_lt25_pct = 100 * success_lt25,
    success_abs_error_lt50_frac = success_lt50,
    success_abs_error_lt50_pct = 100 * success_lt50,
    rmse_mg_m2_d = if (length(physical_error_vals)) sqrt(mean(physical_error_vals^2, na.rm = TRUE)) else NA_real_,
    abs_error_sd_mg_m2_d = safe_sd(abs(physical_error_vals)),
    flux_true_min_mg_m2_d = safe_min(dt$flux_true_mg_m2_d),
    flux_true_max_mg_m2_d = safe_max(dt$flux_true_mg_m2_d),
    pi_s_median = safe_quantile(pi_s_vals, 0.5),
    pi_s_min = safe_min(pi_s_vals),
    pi_s_max = safe_max(pi_s_vals),
    pi_g_median = safe_quantile(pi_g_vals, 0.5),
    pi_g_min = safe_min(pi_g_vals),
    pi_g_max = safe_max(pi_g_vals),
    pi_h_median = safe_quantile(pi_h_vals, 0.5),
    pi_h_min = safe_min(pi_h_vals),
    pi_h_max = safe_max(pi_h_vals),
    mission_duration_min_min = safe_min(mission_vals),
    mission_duration_min_max = safe_max(mission_vals),
    u_ref_m_s_min = safe_min(u_ref_vals),
    u_ref_m_s_max = safe_max(u_ref_vals)
  )
}

build_pi_diagnostic_summary <- function(run_dt, pi_summary = "weighted") {
  dt <- copy(run_dt)
  pi_cols <- pi_column_names(pi_summary)
  if (!all(unlist(pi_cols) %in% names(dt))) {
    stop("Pi summary columns not found for pi_summary = ", pi_summary)
  }
  dt[, `:=`(
    pi_s_for_summary = get(pi_cols$pi_s),
    pi_g_for_summary = get(pi_cols$pi_g),
    pi_h_for_summary = get(pi_cols$pi_h)
  )]
  dt <- dt[is.finite(pi_s_for_summary) & is.finite(pi_g_for_summary) &
             is.finite(pi_h_for_summary)]

  dt[, pi_s_class := fcase(
    pi_s_for_summary < 1, "pi_s_lt_1",
    pi_s_for_summary < 5, "pi_s_1_to_5",
    pi_s_for_summary < 10, "pi_s_5_to_10",
    pi_s_for_summary < 50, "pi_s_10_to_50",
    default = "pi_s_gte_50"
  )]
  dt[, pi_g_class := fcase(
    pi_g_for_summary < 0.75, "pi_g_lt_0.75",
    pi_g_for_summary < 1, "pi_g_0.75_to_1",
    pi_g_for_summary < 1.25, "pi_g_1_to_1.25",
    default = "pi_g_gte_1.25"
  )]
  dt[, pi_h_class := fcase(
    pi_h_for_summary < 1, "pi_h_lt_1",
    pi_h_for_summary < 2, "pi_h_1_to_2",
    pi_h_for_summary < 5, "pi_h_2_to_5",
    default = "pi_h_gte_5"
  )]

  rows <- list()
  add_slice <- function(dt_sub,
                        summary_group,
                        condition_label,
                        condition_description,
                        transition_label = "all",
                        pi_s_class_label = "all",
                        pi_g_class_label = "all",
                        pi_h_class_label = "all") {
    if (!nrow(dt_sub)) return(invisible(NULL))
    stats <- summarise_performance_stats(
      dt_sub,
      pi_s_col = "pi_s_for_summary",
      pi_g_col = "pi_g_for_summary",
      pi_h_col = "pi_h_for_summary"
    )
    rows[[length(rows) + 1L]] <<- cbind(
      data.table(
        summary_group = summary_group,
        condition_label = condition_label,
        condition_description = condition_description,
        transition = transition_label,
        pi_s_class = pi_s_class_label,
        pi_g_class = pi_g_class_label,
        pi_h_class = pi_h_class_label,
        pi_summary = pi_summary
      ),
      stats
    )
    invisible(NULL)
  }

  transition_vals <- sort(unique(dt$transition))
  threshold_conditions <- list(
    list(label = "all_scenarios", desc = "All finite Pi scenarios", idx = rep(TRUE, nrow(dt))),
    list(label = "pi_s_gt_5", desc = "Pi_S > 5", idx = dt$pi_s_for_summary > 5),
    list(label = "pi_s_gt_10", desc = "Pi_S > 10", idx = dt$pi_s_for_summary > 10),
    list(label = "pi_g_lt_1", desc = "Pi_G < 1", idx = dt$pi_g_for_summary < 1),
    list(label = "pi_h_gt_1", desc = "H / sigma_z > 1", idx = dt$pi_h_for_summary > 1),
    list(label = "pi_h_gt_2", desc = "H / sigma_z > 2", idx = dt$pi_h_for_summary > 2),
    list(
      label = "good_pi_s_gt_5_pi_g_lt_1_pi_h_gt_1",
      desc = "Pi_S > 5, Pi_G < 1, and H / sigma_z > 1",
      idx = dt$pi_s_for_summary > 5 & dt$pi_g_for_summary < 1 & dt$pi_h_for_summary > 1
    ),
    list(
      label = "good_pi_s_gt_10_pi_g_lt_1_pi_h_gt_1",
      desc = "Pi_S > 10, Pi_G < 1, and H / sigma_z > 1",
      idx = dt$pi_s_for_summary > 10 & dt$pi_g_for_summary < 1 & dt$pi_h_for_summary > 1
    ),
    list(
      label = "good_pi_s_gt_10_pi_g_lt_1_pi_h_gt_2",
      desc = "Pi_S > 10, Pi_G < 1, and H / sigma_z > 2",
      idx = dt$pi_s_for_summary > 10 & dt$pi_g_for_summary < 1 & dt$pi_h_for_summary > 2
    )
  )
  for (cond in threshold_conditions) {
    add_slice(
      dt[cond$idx],
      summary_group = "threshold_success_rates",
      condition_label = cond$label,
      condition_description = cond$desc
    )
    for (tr in transition_vals) {
      add_slice(
        dt[cond$idx & transition == tr],
        summary_group = "threshold_success_rates_by_transition",
        condition_label = cond$label,
        condition_description = cond$desc,
        transition_label = tr
      )
    }
  }

  pi_s_levels <- c("pi_s_lt_1", "pi_s_1_to_5", "pi_s_5_to_10",
                   "pi_s_10_to_50", "pi_s_gte_50")
  for (lev in pi_s_levels) {
    add_slice(
      dt[pi_s_class == lev],
      summary_group = "conditional_pi_s_class",
      condition_label = lev,
      condition_description = "Pi_S class, all geometries",
      pi_s_class_label = lev
    )
    add_slice(
      dt[pi_s_class == lev & pi_g_for_summary < 1 & pi_h_for_summary > 1],
      summary_group = "conditional_pi_s_class_closed_geometry",
      condition_label = paste0(lev, "_pi_g_lt_1_pi_h_gt_1"),
      condition_description = "Pi_S class after filtering to Pi_G < 1 and H / sigma_z > 1",
      pi_s_class_label = lev,
      pi_g_class_label = "pi_g_lt_1",
      pi_h_class_label = "pi_h_gt_1"
    )
    for (tr in transition_vals) {
      add_slice(
        dt[transition == tr & pi_s_class == lev],
        summary_group = "conditional_pi_s_class_by_transition",
        condition_label = lev,
        condition_description = "Pi_S class within transition",
        transition_label = tr,
        pi_s_class_label = lev
      )
    }
  }

  pi_g_levels <- c("pi_g_lt_0.75", "pi_g_0.75_to_1",
                   "pi_g_1_to_1.25", "pi_g_gte_1.25")
  for (lev in pi_g_levels) {
    add_slice(
      dt[pi_s_for_summary > 10 & pi_g_class == lev],
      summary_group = "conditional_pi_g_class_pi_s_gt_10",
      condition_label = paste0(lev, "_pi_s_gt_10"),
      condition_description = "Pi_G class after filtering to Pi_S > 10",
      pi_s_class_label = "pi_s_gt_10",
      pi_g_class_label = lev
    )
    add_slice(
      dt[pi_s_for_summary > 10 & pi_h_for_summary > 2 & pi_g_class == lev],
      summary_group = "conditional_pi_g_class_pi_s_gt_10_pi_h_gt_2",
      condition_label = paste0(lev, "_pi_s_gt_10_pi_h_gt_2"),
      condition_description = "Pi_G class after filtering to Pi_S > 10 and H / sigma_z > 2",
      pi_s_class_label = "pi_s_gt_10",
      pi_g_class_label = lev,
      pi_h_class_label = "pi_h_gt_2"
    )
    for (tr in transition_vals) {
      add_slice(
        dt[transition == tr & pi_s_for_summary > 10 & pi_g_class == lev],
        summary_group = "conditional_pi_g_class_pi_s_gt_10_by_transition",
        condition_label = paste0(lev, "_pi_s_gt_10"),
        condition_description = "Pi_G class within transition after filtering to Pi_S > 10",
        transition_label = tr,
        pi_s_class_label = "pi_s_gt_10",
        pi_g_class_label = lev
      )
    }
  }

  for (tr in transition_vals) {
    add_slice(
      dt[transition == tr],
      summary_group = "transition_pi_performance",
      condition_label = "transition_all_scenarios",
      condition_description = "All finite Pi scenarios within transition",
      transition_label = tr
    )
    add_slice(
      dt[transition == tr & pi_s_for_summary > 10 & pi_g_for_summary < 1 & pi_h_for_summary > 1],
      summary_group = "transition_pi_performance_good_region_h_gt_1",
      condition_label = "pi_s_gt_10_pi_g_lt_1_pi_h_gt_1",
      condition_description = "Transition performance in Pi_S > 10, Pi_G < 1, and H / sigma_z > 1 region",
      transition_label = tr,
      pi_s_class_label = "pi_s_gt_10",
      pi_g_class_label = "pi_g_lt_1",
      pi_h_class_label = "pi_h_gt_1"
    )
    add_slice(
      dt[transition == tr & pi_s_for_summary > 10 & pi_g_for_summary < 1 & pi_h_for_summary > 2],
      summary_group = "transition_pi_performance_good_region_h_gt_2",
      condition_label = "pi_s_gt_10_pi_g_lt_1_pi_h_gt_2",
      condition_description = "Transition performance in Pi_S > 10, Pi_G < 1, and H / sigma_z > 2 region",
      transition_label = tr,
      pi_s_class_label = "pi_s_gt_10",
      pi_g_class_label = "pi_g_lt_1",
      pi_h_class_label = "pi_h_gt_2"
    )
  }

  if ("u_ref_m_s" %in% names(dt)) {
    wind_speed_vals <- sort(unique(dt$u_ref_m_s[is.finite(dt$u_ref_m_s)]))
    for (u_now in wind_speed_vals) {
      u_label <- sprintf("u_ref_m_s_%s", gsub("\\.", "p", format(u_now, trim = TRUE, scientific = FALSE)))
      add_slice(
        dt[u_ref_m_s == u_now],
        summary_group = "wind_speed_performance",
        condition_label = u_label,
        condition_description = sprintf("All finite Pi scenarios at U = %s m/s",
                                        format(u_now, trim = TRUE, scientific = FALSE))
      )
      add_slice(
        dt[u_ref_m_s == u_now & pi_s_for_summary > 10 &
             pi_g_for_summary < 1 & pi_h_for_summary > 1],
        summary_group = "wind_speed_performance_good_region_h_gt_1",
        condition_label = paste0(u_label, "_pi_s_gt_10_pi_g_lt_1_pi_h_gt_1"),
        condition_description = sprintf(
          "U = %s m/s within Pi_S > 10, Pi_G < 1, and H / sigma_z > 1 region",
          format(u_now, trim = TRUE, scientific = FALSE)
        ),
        pi_s_class_label = "pi_s_gt_10",
        pi_g_class_label = "pi_g_lt_1",
        pi_h_class_label = "pi_h_gt_1"
      )
      add_slice(
        dt[u_ref_m_s == u_now & pi_s_for_summary > 10 &
             pi_g_for_summary < 1 & pi_h_for_summary > 2],
        summary_group = "wind_speed_performance_good_region_h_gt_2",
        condition_label = paste0(u_label, "_pi_s_gt_10_pi_g_lt_1_pi_h_gt_2"),
        condition_description = sprintf(
          "U = %s m/s within Pi_S > 10, Pi_G < 1, and H / sigma_z > 2 region",
          format(u_now, trim = TRUE, scientific = FALSE)
        ),
        pi_s_class_label = "pi_s_gt_10",
        pi_g_class_label = "pi_g_lt_1",
        pi_h_class_label = "pi_h_gt_2"
      )
    }
  }

  rbindlist(rows, fill = TRUE)
}

build_dimensionless_plot <- function(binned_dt, transition_labels = TRANSITION_LABELS) {
  transition_labels <- intersect(transition_labels, unique(binned_dt$transition))
  if (!length(transition_labels)) {
    return(NULL)
  }

  plot_dt <- copy(binned_dt[transition %in% transition_labels])
  plot_dt[, transition_facet := factor(transition, levels = transition_labels)]
  plot_dt[, rmse_log_for_plot := fifelse(is.finite(rmse_pct) & rmse_pct > 0,
                                         rmse_pct, NA_real_)]
  spread_dt <- plot_dt[bin_n >= MIN_BIN_N_FOR_SPREAD_PLOT]
  spread_dt[, abs_error_sd_for_plot := abs_error_sd_pct]
  spread_dt[, abs_error_sd_log_for_plot := fifelse(is.finite(abs_error_sd_for_plot) & abs_error_sd_for_plot > 0,
                                                   abs_error_sd_for_plot, NA_real_)]

  base_theme <- theme_osse(base_size = 11) +
    theme(
      legend.position = "right",
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold")
    )

  rmse_plot <- ggplot(plot_dt) +
    geom_rect(aes(xmin = pi_s_min, xmax = pi_s_max,
                  ymin = pi_g_min, ymax = pi_g_max,
                  fill = rmse_log_for_plot)) +
    geom_vline(xintercept = 5, linetype = "dashed", linewidth = 0.55, colour = "grey20") +
    geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.55, colour = "grey20") +
    facet_wrap(~transition_facet, ncol = 1) +
    scale_x_log10() +
    scale_fill_gradient(
      low = "#F7FCF0",
      high = "#084081",
      trans = "log10",
      breaks = c(1, 3, 10, 30, 100, 300, 1000),
      na.value = "grey88"
    ) +
    labs(
      title = "a) Retrieval RMSE",
      x = expression(Pi[S]),
      y = expression(Pi[G]),
      fill = "RMSE (%)"
    ) +
    base_theme

  rmse_sd_plot <- ggplot(plot_dt) +
    geom_rect(data = spread_dt,
              aes(xmin = pi_s_min, xmax = pi_s_max,
                  ymin = pi_g_min, ymax = pi_g_max,
                  fill = abs_error_sd_log_for_plot)) +
    geom_vline(xintercept = 5, linetype = "dashed", linewidth = 0.55, colour = "grey20") +
    geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.55, colour = "grey20") +
    facet_wrap(~transition_facet, ncol = 1) +
    scale_x_log10() +
    scale_fill_gradient(
      low = "#F7FCF0",
      high = "#084081",
      trans = "log10",
      breaks = c(1, 3, 10, 30, 100, 300, 1000),
      na.value = "grey88"
    ) +
    labs(
      title = "b) Within-bin absolute-error spread",
      x = expression(Pi[S]),
      y = expression(Pi[G]),
      fill = "SD (%)"
    ) +
    base_theme

  cowplot::plot_grid(rmse_plot, rmse_sd_plot, ncol = 2, rel_widths = c(1, 1))
}

format_duration <- function(seconds) {
  if (!is.finite(seconds)) return("NA")
  if (seconds < 60) return(sprintf("%.1f sec", seconds))
  if (seconds < 3600) return(sprintf("%.1f min", seconds / 60))
  sprintf("%.2f h", seconds / 3600)
}

build_uniform_height_midpoints <- function(n_levels, top_height_m) {
  step_m <- top_height_m / n_levels
  seq(step_m / 2, top_height_m - step_m / 2, length.out = n_levels)
}

resolve_dimensionless_design_cfg <- function(base_design_cfg, top_height_m, height_mode) {
  heights_native <- base_design_cfg$heights_m
  heights_m <- switch(
    height_mode,
    uniform_top = build_uniform_height_midpoints(length(heights_native), top_height_m),
    native = heights_native,
    scaled_native = heights_native / max(heights_native) * top_height_m
  )
  out <- base_design_cfg
  out$name <- sprintf("%s_%s_top%gm", base_design_cfg$name, height_mode, top_height_m)
  out$heights_m <- heights_m
  out
}

build_report_lines <- function(run_dt, timing_dt, binned_dt, pi_diagnostic_dt,
                               full_counts, is_pilot, n_workers = 1L,
                               wall_time_s = NA_real_) {
  valid_dt <- run_dt[qc_flag == "ok" & is.finite(flux_retrieved_mg_m2_d)]
  median_base_s <- stats::median(timing_dt$base_profile_seconds, na.rm = TRUE)
  median_flux_s <- stats::median(timing_dt$flux_loop_seconds / pmax(timing_dt$n_fluxes, 1), na.rm = TRUE)
  est_full_serial_s <- median_base_s * full_counts$n_base + median_flux_s * full_counts$n_scenarios
  est_full_parallel_s <- est_full_serial_s / max(1L, as.integer(n_workers))
  main_transition_labels <- intersect(PRIMARY_TRANSITION_LABELS, TRANSITION_LABELS)
  if (!length(main_transition_labels)) {
    main_transition_labels <- TRANSITION_LABELS
  }
  additional_transition_labels <- intersect(ADDITIONAL_TRANSITION_LABELS, TRANSITION_LABELS)
  supplemental_label_text <- if (length(additional_transition_labels)) {
    paste(additional_transition_labels, collapse = ", ")
  } else {
    "none"
  }

  c(
    "# Dimensionless framework analysis",
    "",
    sprintf("- Mode: %s", ifelse(isTRUE(is_pilot), "pilot timing subset", "full sweep")),
    sprintf("- Transitions: %s", paste(TRANSITION_LABELS, collapse = ", ")),
    sprintf("- Main figure transitions: %s",
            paste(main_transition_labels, collapse = ", ")),
    sprintf("- Supplemental figure transitions: %s",
            supplemental_label_text),
    sprintf("- Receptor design: %s (%d azimuths x %d heights = %d receptors)",
            DESIGN_NAME, SAMPLING_CONFIGS[[DESIGN_NAME]]$n_azimuth,
            length(SAMPLING_CONFIGS[[DESIGN_NAME]]$heights_m),
            SAMPLING_CONFIGS[[DESIGN_NAME]]$n_azimuth * length(SAMPLING_CONFIGS[[DESIGN_NAME]]$heights_m)),
    sprintf("- Height mode: %s; cylinder top: %.1f m", HEIGHT_MODE, TOP_HEIGHT_M),
    sprintf("- Mission durations: %s min", paste(MISSION_DURATIONS_MIN, collapse = ", ")),
    "- Note: each mission is one start-of-transition retrieval window; duration changes dwell time and transition exposure, not sampling density.",
    sprintf("- Wind speed sweep: %s m/s; mean direction %.1f deg; stochastic direction = %s",
            format_num_vec(U_REF_GRID_M_S), WD_DEG, STOCHASTIC_WIND),
    sprintf("- Pi binning summary: %s", PI_SUMMARY),
    sprintf("- Run rows completed: %d", nrow(run_dt)),
    sprintf("- Retrieval windows completed: %d", sum(run_dt$n_retrieval_windows, na.rm = TRUE)),
    sprintf("- Valid retrieval fraction: %.3f", nrow(valid_dt) / max(nrow(run_dt), 1)),
    sprintf("- Run-level bias range: %.2f to %.2f%%",
            min(run_dt$retrieval_bias_pct, na.rm = TRUE),
            max(run_dt$retrieval_bias_pct, na.rm = TRUE)),
    sprintf("- Binned RMSE range: %.2f to %.2f%%",
            min(binned_dt$rmse_pct, na.rm = TRUE),
            max(binned_dt$rmse_pct, na.rm = TRUE)),
    sprintf("- Binned rows: %d", nrow(binned_dt)),
    sprintf("- Pi diagnostic summary rows: %d", nrow(pi_diagnostic_dt)),
    "",
    "## Runtime",
    "",
    sprintf("- Observed wall-clock time for this run: %s", format_duration(wall_time_s)),
    sprintf("- Worker count: %d", as.integer(n_workers)),
    sprintf("- Median base plume time: %s", format_duration(median_base_s)),
    sprintf("- Median per-flux retrieval time: %s", format_duration(median_flux_s)),
    sprintf("- Estimated full-grid serial-equivalent runtime from these timings: %s", format_duration(est_full_serial_s)),
    sprintf("- Estimated full-grid parallel runtime from these timings: %s", format_duration(est_full_parallel_s)),
    sprintf("- Full-grid base cases: %d", full_counts$n_base),
    sprintf("- Full-grid scenario rows: %d", full_counts$n_scenarios),
    "",
    "## Next steps",
    "",
    "- Compare weighted transitional Pi binning against start, mid, end, and unweighted mean Pi binning.",
    "- Decide whether stochastic wind-direction variability should remain the default or become a sensitivity mode.",
    "- Add a Monte Carlo mode if uncertainty is needed rather than single-run bin RMSE.",
    "- Add a second experiment that varies receptor count if sampling-density effects need explicit testing.",
    "- Keep real-lake overlay points deferred until literature values and assumed mission designs are agreed."
  )
}


# -----------------------------------------------------------------------------
# 3. Experiment setup
# -----------------------------------------------------------------------------
base_design_cfg <- SAMPLING_CONFIGS[[DESIGN_NAME]]
design_cfg <- resolve_dimensionless_design_cfg(base_design_cfg, TOP_HEIGHT_M, HEIGHT_MODE)
sensor_cfg <- SENSOR_CONFIGS[[SENSOR_NAME]]
sensor_cfg$drift_ppb_hr <- 0
sensor_cfg$cal_offset_ppb <- 0
bestfit_interp <- if (isTRUE(STOCHASTIC_WIND)) {
  load_bestfit_interpolators(BESTFIT_PATH)
} else {
  NULL
}

n_receptors <- as.integer(design_cfg$n_azimuth * length(design_cfg$heights_m))
full_base_grid <- make_full_base_grid()
full_n_base <- nrow(full_base_grid)
full_n_scenarios <- full_n_base * length(FLUX_GRID_MG_M2_D)
base_grid <- if (isTRUE(PILOT)) make_pilot_base_grid(full_base_grid) else full_base_grid
flux_grid_to_run <- if (isTRUE(PILOT)) {
  sort(unique(c(min(FLUX_GRID_MG_M2_D), max(FLUX_GRID_MG_M2_D))))
} else {
  sort(unique(FLUX_GRID_MG_M2_D))
}
available_cores <- parallel::detectCores(logical = TRUE)
if (!is.finite(available_cores) || available_cores < 1L) available_cores <- 1L
WORKER_SPARE_CORES <- max(0L, min(as.integer(WORKER_SPARE_CORES), max(0L, available_cores - 1L)))
N_WORKERS <- if (is.finite(N_WORKERS_REQUESTED) && !is.na(N_WORKERS_REQUESTED)) {
  max(1L, as.integer(N_WORKERS_REQUESTED))
} else if (isTRUE(PILOT)) {
  1L
} else {
  max(1L, available_cores - WORKER_SPARE_CORES)
}
N_WORKERS <- min(N_WORKERS, nrow(base_grid))
BATCH_SIZE <- if (is.finite(BATCH_SIZE_REQUESTED) && !is.na(BATCH_SIZE_REQUESTED)) {
  max(1L, as.integer(BATCH_SIZE_REQUESTED))
} else {
  max(1L, N_WORKERS)
}

message("========================================")
message(" Dimensionless Framework Analysis")
message("========================================")
message(sprintf(" Mode:        %s", ifelse(PILOT, "pilot", "full")))
message(sprintf(" Sensor:      %s", sensor_cfg$name))
message(sprintf(" Design:      %s", DESIGN_NAME))
message(sprintf(" Height mode: %s | top height %.1f m", HEIGHT_MODE, TOP_HEIGHT_M))
message(sprintf(" Sampling:    %d azimuths x %d heights = %d receptors",
                design_cfg$n_azimuth, length(design_cfg$heights_m), n_receptors))
message(sprintf(" Missions:    %s min, one retrieval window starting at t = 0",
                paste(MISSION_DURATIONS_MIN, collapse = ", ")))
message(sprintf(" Wind speeds: %s m/s | mean direction %.1f deg | stochastic direction: %s",
                format_num_vec(U_REF_GRID_M_S), WD_DEG, STOCHASTIC_WIND))
message(sprintf(" Pi summary:  %s", PI_SUMMARY))
message(sprintf(" Transition:  %.1f h", TRANSITION_DURATION_MIN / 60))
message(sprintf(" Cases:       %d base plume(s), %d scenario row(s)",
                nrow(base_grid), nrow(base_grid) * length(flux_grid_to_run)))
message(sprintf(" Full grid:   %d base plume(s), %d scenario row(s)",
                full_n_base, full_n_scenarios))
message(sprintf(" Parallel:    %d worker%s on %d detected core%s | batch size %d",
                N_WORKERS, ifelse(N_WORKERS == 1L, "", "s"),
                available_cores, ifelse(available_cores == 1L, "", "s"),
                BATCH_SIZE))
if (!is.finite(N_WORKERS_REQUESTED) || is.na(N_WORKERS_REQUESTED)) {
  message(sprintf(" Spare cores: %d (auto workers = detected cores - spare cores)", WORKER_SPARE_CORES))
}
message("========================================")

run_base_case <- function(base_spec) {
  base_idx <- base_spec$base_idx
  base_row <- as.data.table(base_spec$base_row)
  lake_radius_m <- base_row$lake_radius_m
  ring_radius_m <- base_row$ring_radius_ratio * lake_radius_m
  distance_from_edge_m <- ring_radius_m - lake_radius_m
  u_ref_m_s <- base_row$u_ref_m_s
  wind_seed <- SEED_BASE + base_idx * 9973L
  wind_tl <- if (isTRUE(STOCHASTIC_WIND)) {
    generate_wind_timeline(
      transition_duration_min = TRANSITION_DURATION_MIN,
      met_dt_min = MET_DT_MIN,
      u_mean_m_s = u_ref_m_s,
      wd_mean_deg = WD_DEG,
      bestfit_interp = bestfit_interp,
      class_start = base_row$class_start,
      class_end = base_row$class_end,
      seed = wind_seed
    )
  } else {
    NULL
  }

  base_start <- Sys.time()
  base_profile <- build_base_profile_truth(
    mission_duration_min = base_row$mission_duration_min,
    design_cfg = design_cfg,
    lake_radius_m = lake_radius_m,
    ring_radius_m = ring_radius_m,
    transition_duration_min = TRANSITION_DURATION_MIN,
    class_start = base_row$class_start,
    class_end = base_row$class_end,
    u_ref_m_s = u_ref_m_s,
    wd_deg = WD_DEG,
    wind_tl = wind_tl
  )
  base_seconds <- as.numeric(difftime(Sys.time(), base_start, units = "secs"))
  wd_mean_case <- circular_mean_deg(base_profile$base_dt$wd_deg)
  wd_sd_case <- circular_sd_deg(base_profile$base_dt$wd_deg)

  local_rows <- vector("list", length(flux_grid_to_run))
  flux_start <- Sys.time()
  for (flux_idx in seq_along(flux_grid_to_run)) {
    flux_now <- flux_grid_to_run[flux_idx]
    scenario_label <- sprintf(
      "dimensionless_%s_r%07.2f_R%04.2f_t%03d_u%04.1f_q%08.3f",
      base_row$transition,
      lake_radius_m,
      base_row$ring_radius_ratio,
      as.integer(round(base_row$mission_duration_min)),
      u_ref_m_s,
      flux_now
    )
    scenario_seed <- SEED_BASE + base_idx * 10000L + flux_idx

    retrieval <- run_flux_retrieval(
      base_dt = base_profile$base_dt,
      flux_mg_m2_d = flux_now,
      sensor_cfg = sensor_cfg,
      seed = scenario_seed,
      lake_radius_m = lake_radius_m,
      ring_radius_m = ring_radius_m,
      n_azimuth = design_cfg$n_azimuth,
      scenario_label = scenario_label
    )
    flux_dt <- retrieval$flux
    retrieval_summary <- summarise_retrieval_windows(flux_dt)

    pi_dt <- compute_pi_at_times(
      sample_times_min = base_profile$base_dt$sample_time_min,
      flux_mg_m2_d = flux_now,
      lake_radius_m = lake_radius_m,
      ring_radius_m = ring_radius_m,
      top_height_m = TOP_HEIGHT_M,
      mission_duration_min = base_row$mission_duration_min,
      transition_duration_min = TRANSITION_DURATION_MIN,
      class_start = base_row$class_start,
      class_end = base_row$class_end,
      u_ref_m_s = u_ref_m_s,
      sensor_cfg = sensor_cfg,
      sample_rate_hz = SAMPLE_RATE_HZ,
      sample_weights = pmax(base_profile$base_dt$ch4_enhancement_ppb_per_mg_m2_d, 0)
    )
    pi_cols <- pi_column_names(PI_SUMMARY)
    pi_dt[, `:=`(
      pi_s_selected = get(pi_cols$pi_s),
      pi_g_selected = get(pi_cols$pi_g),
      pi_h_selected = get(pi_cols$pi_h)
    )]

    flux_retrieved <- retrieval_summary$flux_retrieved_mg_m2_d
    flux_error <- flux_retrieved - flux_now
    bias_pct <- 100 * flux_error / flux_now
    dwell_s <- unique(base_profile$base_dt$dwell_s)[1]

    local_rows[[flux_idx]] <- cbind(
      data.table(
        scenario_label = scenario_label,
        transition = base_row$transition,
        class_start = base_row$class_start,
        class_end = base_row$class_end,
        flux_true_mg_m2_d = flux_now,
        flux_retrieved_mg_m2_d = flux_retrieved,
        flux_error_mg_m2_d = flux_error,
        retrieval_bias_pct = bias_pct,
        squared_rel_error = (flux_error / flux_now)^2,
        lake_radius_m = lake_radius_m,
        lake_diameter_m = 2 * lake_radius_m,
        ring_radius_m = ring_radius_m,
        ring_radius_ratio = base_row$ring_radius_ratio,
        distance_from_edge_m = distance_from_edge_m,
        mission_duration_min = base_row$mission_duration_min,
        total_profile_time_min = base_profile$total_profile_time_min,
        transition_duration_min = TRANSITION_DURATION_MIN,
        sample_rate_hz = SAMPLE_RATE_HZ,
        u_ref_m_s = u_ref_m_s,
        wd_mean_deg = WD_DEG,
        wd_sample_mean_deg = wd_mean_case,
        wd_sample_sd_deg = wd_sd_case,
        stochastic_wind = STOCHASTIC_WIND,
        wind_seed = wind_seed,
        n_azimuth = design_cfg$n_azimuth,
        n_height_levels = length(design_cfg$heights_m),
        n_receptors = n_receptors,
        n_samples_per_profile = base_profile$n_samples_per_profile,
        n_measurements_total = nrow(base_profile$base_dt),
        n_retrieval_windows = retrieval_summary$n_retrieval_windows,
        n_valid_retrieval_windows = retrieval_summary$n_valid_retrieval_windows,
        retrieval_window_valid_fraction = retrieval_summary$retrieval_window_valid_fraction,
        min_height_m = min(design_cfg$heights_m),
        max_height_m = max(design_cfg$heights_m),
        top_height_m = TOP_HEIGHT_M,
        height_mode = HEIGHT_MODE,
        dwell_s = dwell_s,
        mean_enhancement_ppb = retrieval_summary$mean_enhancement_ppb,
        ws_at_z_mean_m_s = retrieval_summary$ws_at_z_mean_m_s,
        within_run_flux_sd_mg_m2_d = retrieval_summary$within_run_flux_sd_mg_m2_d,
        qc_flag = retrieval_summary$qc_flag,
        seed = scenario_seed
      ),
      pi_dt
    )
  }
  flux_loop_seconds <- as.numeric(difftime(Sys.time(), flux_start, units = "secs"))

  timing <- data.table(
    base_idx = base_idx,
    transition = base_row$transition,
    lake_radius_m = lake_radius_m,
    ring_radius_ratio = base_row$ring_radius_ratio,
    ring_radius_m = ring_radius_m,
    mission_duration_min = base_row$mission_duration_min,
    u_ref_m_s = u_ref_m_s,
    total_profile_time_min = base_profile$total_profile_time_min,
    n_fluxes = length(flux_grid_to_run),
    n_receptors = n_receptors,
    n_measurements_total = nrow(base_profile$base_dt),
    stochastic_wind = STOCHASTIC_WIND,
    wind_seed = wind_seed,
    wd_sample_mean_deg = wd_mean_case,
    wd_sample_sd_deg = wd_sd_case,
    base_profile_seconds = base_seconds,
    flux_loop_seconds = flux_loop_seconds,
    total_seconds = base_seconds + flux_loop_seconds
  )

  list(
    run_rows = rbindlist(local_rows, fill = TRUE),
    timing = timing
  )
}


# -----------------------------------------------------------------------------
# 4. Run sweep
# -----------------------------------------------------------------------------
pb <- make_progress(nrow(base_grid), label = "Dimensionless base plumes")
sweep_start <- Sys.time()
base_specs <- lapply(seq_len(nrow(base_grid)), function(idx) {
  list(base_idx = idx, base_row = base_grid[idx])
})
base_chunks <- split(base_specs, ceiling(seq_along(base_specs) / BATCH_SIZE))

cl <- NULL
if (N_WORKERS > 1L) {
  worker_outfile <- tempfile(pattern = "dimensionless_parallel_workers_", fileext = ".log")
  cl <- parallel::makePSOCKcluster(
    N_WORKERS,
    outfile = worker_outfile,
    rscript_args = c("--vanilla")
  )
  on.exit(parallel::stopCluster(cl), add = TRUE)
  on.exit(unlink(worker_outfile, force = TRUE), add = TRUE)
  parallel::clusterExport(cl, varlist = c("proj_root"), envir = environment())
  parallel::clusterEvalQ(cl, {
    library(data.table)
    options(uas_gauss_osse_root = proj_root)
    suppressMessages(source(file.path(proj_root, "R", "utils.R")))
    suppressMessages(source(file.path(proj_root, "R", "04_conley_inverse_flux.R")))
    NULL
  })
  parallel::clusterExport(
    cl,
    varlist = c(
      "safe_sd",
      "quietly",
      "interpolate_briggs",
      "zeta_at_time",
      "generate_wind_timeline",
      "angle_diff_deg",
      "circular_mean_deg",
      "circular_sd_deg",
      "wd_at_time",
      "pg_sigma_y_interp",
      "pg_sigma_z_interp",
      "interpolate_briggs_table",
      "build_disk_quadrature",
      "gaussian_plume_disk_timevarying",
      "build_receptors_for_ring",
      "compute_pi_at_times",
      "build_base_profile_truth",
      "build_measurements_for_flux",
      "run_flux_retrieval",
      "summarise_retrieval_windows",
      "pi_column_names",
      "run_base_case",
      "PASQUILL_ZETA_ANCHORS",
      "design_cfg",
      "sensor_cfg",
      "bestfit_interp",
      "n_receptors",
      "flux_grid_to_run",
      "TRANSITION_DURATION_MIN",
      "MET_DT_MIN",
      "U_REF_GRID_M_S",
      "WD_DEG",
      "STOCHASTIC_WIND",
      "SEED_BASE",
      "SAMPLE_RATE_HZ",
      "PI_SUMMARY",
      "TOP_HEIGHT_M",
      "HEIGHT_MODE"
    ),
    envir = environment()
  )
}

result_chunks <- vector("list", length(base_chunks))
for (chunk_idx in seq_along(base_chunks)) {
  chunk <- base_chunks[[chunk_idx]]
  chunk_results <- if (!is.null(cl)) {
    parallel::parLapplyLB(cl, chunk, run_base_case)
  } else {
    lapply(chunk, run_base_case)
  }
  result_chunks[[chunk_idx]] <- chunk_results
  for (tick_idx in seq_along(chunk)) pb$tick()
}

flat_results <- unlist(result_chunks, recursive = FALSE)
run_dt <- rbindlist(lapply(flat_results, `[[`, "run_rows"), fill = TRUE)
timing_dt <- rbindlist(lapply(flat_results, `[[`, "timing"), fill = TRUE)
binned_dt <- summarise_binned_results(run_dt, N_PI_S_BINS, N_PI_G_BINS, PI_SUMMARY)
pi_diagnostic_dt <- build_pi_diagnostic_summary(run_dt, PI_SUMMARY)
sweep_wall_time_s <- as.numeric(difftime(Sys.time(), sweep_start, units = "secs"))


# -----------------------------------------------------------------------------
# 5. Save tables and report
# -----------------------------------------------------------------------------
fwrite(run_dt, RUNS_CSV)
fwrite(timing_dt, TIMINGS_CSV)
fwrite(binned_dt, BINNED_CSV)
fwrite(pi_diagnostic_dt, PI_DIAGNOSTIC_CSV)
report_lines <- build_report_lines(
  run_dt = run_dt,
  timing_dt = timing_dt,
  binned_dt = binned_dt,
  pi_diagnostic_dt = pi_diagnostic_dt,
  full_counts = list(n_base = full_n_base, n_scenarios = full_n_scenarios),
  is_pilot = PILOT,
  n_workers = N_WORKERS,
  wall_time_s = sweep_wall_time_s
)
writeLines(report_lines, REPORT_MD)

message(sprintf("Saved run-level CSV:      %s", normalizePath(RUNS_CSV, mustWork = FALSE)))
message(sprintf("Saved timing CSV:         %s", normalizePath(TIMINGS_CSV, mustWork = FALSE)))
message(sprintf("Saved binned summary CSV: %s", normalizePath(BINNED_CSV, mustWork = FALSE)))
message(sprintf("Saved Pi diagnostics CSV: %s", normalizePath(PI_DIAGNOSTIC_CSV, mustWork = FALSE)))
message(sprintf("Saved report:             %s", normalizePath(REPORT_MD, mustWork = FALSE)))


# -----------------------------------------------------------------------------
# 6. Build composite figure
# -----------------------------------------------------------------------------
main_transition_labels <- intersect(PRIMARY_TRANSITION_LABELS, TRANSITION_LABELS)
if (!length(main_transition_labels)) {
  main_transition_labels <- TRANSITION_LABELS
}
additional_transition_labels <- intersect(ADDITIONAL_TRANSITION_LABELS, TRANSITION_LABELS)

composite_plot <- build_dimensionless_plot(
  binned_dt = binned_dt,
  transition_labels = main_transition_labels
)
if (!is.null(composite_plot)) {
  ggsave_osse(FIG_PNG, composite_plot, width = 14, height = 10.5, dpi = 220)
  message(sprintf("Saved composite figure:   %s", normalizePath(FIG_PNG, mustWork = FALSE)))
}

additional_plot <- build_dimensionless_plot(
  binned_dt = binned_dt,
  transition_labels = additional_transition_labels
)
if (!is.null(additional_plot)) {
  ggsave_osse(ADDITIONAL_FIG_PNG, additional_plot, width = 14, height = 8.5, dpi = 220)
  message(sprintf("Saved supplemental figure: %s", normalizePath(ADDITIONAL_FIG_PNG, mustWork = FALSE)))
}

message("Dimensionless framework analysis complete.")

