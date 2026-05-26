#!/usr/bin/env Rscript
###############################################################################
# process_ameriflux_wind_variability.R
#
# Unzip AmeriFlux BASE ZIP archives, load HH/HR CSVs, compute wind variability
# and stability metrics, and write summary CSVs + best-fit figure.
#
# Usage (from repository root):
#   Rscript --vanilla scripts/process_ameriflux_wind_variability.R
#
# Use --vanilla to avoid loading a project .Rprofile that sources utils.R.
#
# Requires: data.table, ggplot2 (install.packages(c("data.table", "ggplot2")))
#
# Notes:
#   - Recursively finds all .zip under INPUT_DIR (e.g. data/fluxData).
#   - HH/HR CSVs must include U_SIGMA, V_SIGMA, and MO_LENGTH (or L) for rows
###############################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# -----------------------------------------------------------------------------
# CONFIGURATION
# -----------------------------------------------------------------------------
INPUT_DIR  <- normalizePath("data/fluxData", mustWork = FALSE)
OUTPUT_DIR <- normalizePath("output", mustWork = FALSE)

WS_MIN           <- 1.0   # m/s; minimum wind speed
DEFAULT_Z        <- 10.0  # m; fallback sonic height if not parsed from metadata
QC_MAX           <- 1L    # keep QC in {0, 1}; exclude if QC > QC_MAX when column exists
ZETA_CLIP        <- 5.0   # clip |ZETA| to this value
MIN_OBS_PER_BIN  <- 30L   # minimum rows per ZETA bin for fitting binned medians

# TI = U_SIGMA/WS and sigma_theta ~ V_SIGMA/WS blow up when WS is at WS_MIN but
# sigmas are large; drop non-physical / extreme ratio rows.
TI_MAX_RETAIN          <- 10.0   # dimensionless; max turbulence intensity kept
SIGMA_THETA_ABS_MAX    <- 90.0   # degrees; small-angle validity + plot sanity

# Figure ribbons: hard clip displayed quantiles (defense in depth vs. outliers)
RIBBON_U_SIGMA_MAX     <- 15.0   # m/s
RIBBON_TI_MAX          <- 10.0
RIBBON_SIGMA_THETA_ABS <- 90.0   # degrees
RIBBON_MIN_N           <- 10L    # skip ribbon in sparse bins (unstable quantiles)

# Optional: set to a single site folder name to inspect one site only (NULL = all)
DEBUG_SITE_ID <- NULL

# Optional CLI: Rscript script.R [INPUT_DIR]
# If the first argument is an existing directory, it overrides INPUT_DIR.
.local_args <- commandArgs(trailingOnly = TRUE)
if (length(.local_args) >= 1L && dir.exists(.local_args[1L])) {
  INPUT_DIR <- normalizePath(.local_args[1L], mustWork = TRUE)
}

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------
parse_site_id_from_zip <- function(zip_path) {
  bn <- basename(zip_path)
  m  <- regexec("^AMF_([A-Za-z0-9-]+)_BASE", bn)
  r  <- regmatches(bn, m)[[1]]
  if (length(r) < 2L) {
    warning("Could not parse site ID from: ", bn)
    return(NA_character_)
  }
  r[2]
}

find_zip_files <- function(root_dir) {
  if (!dir.exists(root_dir)) return(character(0))
  sort(list.files(root_dir, pattern = "\\.zip$", recursive = TRUE, full.names = TRUE,
                  ignore.case = TRUE))
}

unzip_to_site_dir <- function(zip_path, input_root, site_id, overwrite = FALSE) {
  dest <- file.path(dirname(zip_path), site_id)
  if (!dir.exists(dest)) dir.create(dest, recursive = TRUE, showWarnings = FALSE)

  # Already unzipped: look for BASE HH/HR CSV
  existing <- list.files(dest, pattern = "BASE_(HH|HR)_", full.names = TRUE)
  if (length(existing) > 0L && !overwrite) {
    message("  [skip unzip] CSV already present in ", basename(dest))
    return(dest)
  }

  message("  Unzipping -> ", dest)
  tryCatch(
    unzip(zip_path, exdir = dest, overwrite = TRUE),
    error = function(e) {
      warning("Unzip failed for ", zip_path, ": ", e$message)
    }
  )
  dest
}

#' Count leading # comment lines for fread skip=
count_leading_hash_lines <- function(path) {
  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
  n <- 0L
  for (i in seq_along(lines)) {
    if (grepl("^\\s*#", lines[i])) n <- n + 1L else break
  }
  n
}

#' Parse sensor height (m) from AmeriFlux # metadata lines
parse_sensor_height_from_metadata <- function(path, ws_var) {
  lines <- readLines(path, warn = FALSE, n = 500L, encoding = "UTF-8")
  txt <- paste(lines[grepl("^\\s*#", lines)], collapse = " ")
  esc <- gsub("([][{}()^$|*+?.\\\\])", "\\\\\\1", ws_var, perl = TRUE)
  m <- regexpr(paste0(esc, "[^0-9]{0,40}([0-9]+\\.?[0-9]*)\\s*m"), txt, perl = TRUE, ignore.case = TRUE)
  if (m[1] > 0) {
    cap <- regmatches(txt, regexec(paste0(esc, "[^0-9]{0,40}([0-9]+\\.?[0-9]*)\\s*m"),
                                  txt, perl = TRUE, ignore.case = TRUE))[[1]]
    if (length(cap) >= 2) return(as.numeric(cap[2]))
  }
  m2 <- regexpr("([0-9]+\\.?[0-9]*)\\s*m\\b", txt, perl = TRUE)
  if (m2[1] > 0) {
    cap2 <- regmatches(txt, regexec("([0-9]+\\.?[0-9]*)\\s*m\\b", txt, perl = TRUE))[[1]]
    if (length(cap2) >= 2) {
      z <- as.numeric(cap2[2])
      if (is.finite(z) && z > 0.5 && z < 200) return(z)
    }
  }
  NA_real_
}

replace_missing_flags <- function(dt) {
  for (nm in names(dt)) {
    if (!is.numeric(dt[[nm]])) next
    dt[[nm]][dt[[nm]] <= -9000] <- NA_real_
  }
  invisible(dt)
}

#' WS columns: WS_H_V_R (exclude gap-filled _F, PI variants _PI)
find_ws_candidates <- function(names) {
  ws <- grep("^WS_[0-9]+_[0-9]+_[0-9]+$", names, value = TRUE)
  ws <- ws[!grepl("_F$|_PI$", ws, ignore.case = TRUE)]
  unique(ws)
}

parse_hvr <- function(var) {
  p <- strsplit(var, "_", fixed = FALSE)[[1]]
  if (length(p) < 4) return(NULL)
  tryCatch(
    list(H = as.integer(p[2]), V = as.integer(p[3]), R = as.integer(p[4])),
    error = function(...) NULL
  )
}

#' Select WS with minimum V (highest physical level), tie-break R == 1
select_ws_variable <- function(ws_cols) {
  if (length(ws_cols) == 0L) return(NULL)
  rows <- lapply(ws_cols, function(w) {
    p <- parse_hvr(w)
    if (is.null(p)) return(NULL)
    data.table(ws = w, V = p$V, R = p$R, H = p$H)
  })
  rows <- rbindlist(rows[!vapply(rows, is.null, logical(1))])
  if (nrow(rows) == 0L) return(NULL)
  min_v <- min(rows$V, na.rm = TRUE)
  rows <- rows[V == min_v]
  if (any(rows$R == 1L)) rows[R == 1L][1]$ws else rows[which.min(R)]$ws
}

col_at_hvr <- function(prefix, H, V, R, names_dt) {
  exact <- sprintf("%s_%d_%d_%d", prefix, H, V, R)
  if (exact %in% names_dt) return(exact)
  # Alternate naming US_SIGMA vs U_SIGMA
  alt <- if (prefix == "U_SIGMA") {
    c(sprintf("US_SIGMA_%d_%d_%d", H, V, R))
  } else character(0)
  for (a in alt) if (a %in% names_dt) return(a)
  NA_character_
}

find_first_matching <- function(pattern, names_dt) {
  hit <- grep(pattern, names_dt, value = TRUE)[1]
  if (is.na(hit)) NA_character_ else hit
}

#' Prefer MO_LENGTH / USTAR at same H,V as WS; else any MO_LENGTH / USTAR
resolve_mo_length <- function(dt, H, V, ws_name) {
  nm <- names(dt)
  cand <- c(
    sprintf("MO_LENGTH_%d_%d_%d", H, V, 1L),
    sprintf("MO_LENGTH_%d_%d_%d", H, V, 2L),
    "MO_LENGTH"
  )
  for (c in cand) if (c %in% nm) return(c)
  m <- grep("^MO_LENGTH(_|$)", nm, value = TRUE)
  if (length(m)) return(m[1])
  if ("L" %in% nm) return("L")
  NA_character_
}

resolve_ustar <- function(dt, H, V) {
  nm <- names(dt)
  cand <- c(
    sprintf("USTAR_%d_%d_%d", H, V, 1L),
    sprintf("USTAR_%d_%d_%d", H, V, 2L),
    "USTAR"
  )
  for (c in cand) if (c %in% nm) return(c)
  m <- grep("^USTAR(_|$)", nm, value = TRUE)
  if (length(m)) m[1] else NA_character_
}

resolve_wd <- function(dt, H, V) {
  nm <- names(dt)
  cand <- c(
    sprintf("WD_%d_%d_%d", H, V, 1L),
    "WD"
  )
  for (c in cand) if (c %in% nm) return(c)
  NA_character_
}

find_qc_column <- function(ws_col, names_dt) {
  cand <- paste0(ws_col, "_QC")
  if (cand %in% names_dt) return(cand)
  NA_character_
}

pick_preferred_csv <- function(site_dir) {
  f <- list.files(site_dir, pattern = "BASE_(HH|HR)_.*\\.csv$", full.names = TRUE,
                  ignore.case = TRUE)
  if (length(f) == 0L) return(NA_character_)
  hh <- f[grepl("_HH_", basename(f), ignore.case = TRUE)]
  if (length(hh)) return(hh[1])
  hr <- f[grepl("_HR_", basename(f), ignore.case = TRUE)]
  if (length(hr)) return(hr[1])
  f[1]
}

read_ameriflux_csv <- function(path) {
  skip <- count_leading_hash_lines(path)
  dt <- tryCatch(
    data.table::fread(path, skip = skip, na.strings = c("NA", ""), showProgress = FALSE),
    error = function(e) {
      warning("fread failed: ", path, " - ", e$message)
      NULL
    }
  )
  if (is.null(dt)) return(NULL)
  setnames(dt, trimws(names(dt)))
  if (!"TIMESTAMP_START" %in% names(dt)) {
    warning("TIMESTAMP_START missing in ", path)
    return(NULL)
  }
  replace_missing_flags(dt)
  dt
}

parse_timestamp <- function(x) {
  x <- as.character(x)
  # %Y%m%d%H%M: 12 chars typical; some files may use 10 (hourly?)
  as.POSIXct(x, format = "%Y%m%d%H%M", tz = "UTC")
}

stability_class_from_zeta <- function(zeta) {
  fcase(
    is.na(zeta), NA_character_,
    zeta < -1.0, "strongly_unstable",
    zeta < -0.1, "unstable",
    zeta <  0.1, "near_neutral",
    zeta <  1.0, "stable",
    default = "strongly_stable"
  )
}

#' Piecewise cubic on z < 0 and z >= 0, fitted to (z, y) binned medians
fit_piecewise_poly <- function(z, y, degree = 3L) {
  ok <- is.finite(z) & is.finite(y)
  z <- z[ok]; y <- y[ok]
  fit_side <- function(idx) {
    if (sum(idx) < degree + 1L) return(NULL)
    zs <- z[idx]; ys <- y[idx]
    tryCatch(
      lm(ys ~ poly(zs, degree = degree, raw = TRUE)),
      error = function(...) NULL
    )
  }
  list(
    neg = fit_side(z < 0),
    pos = fit_side(z >= 0)
  )
}

predict_piecewise_poly <- function(z, fits, degree = 3L) {
  out <- rep(NA_real_, length(z))
  if (is.null(fits)) return(out)
  pred1 <- function(fit, zz) {
    if (is.null(fit)) return(rep(NA_real_, length(zz)))
    tryCatch({
      df <- data.frame(zs = zz)
      predict(fit, newdata = df)
    }, error = function(...) rep(NA_real_, length(zz)))
  }
  neg <- z < 0 & is.finite(z)
  pos <- z >= 0 & is.finite(z)
  if (any(neg)) out[neg] <- pred1(fits$neg, z[neg])
  if (any(pos)) out[pos] <- pred1(fits$pos, z[pos])
  out
}

fit_quality <- function(y, yhat) {
  ok <- is.finite(y) & is.finite(yhat)
  y <- y[ok]
  yhat <- yhat[ok]
  if (length(y) == 0L) {
    return(list(n = 0L, rmse = NA_real_, r2 = NA_real_))
  }
  resid <- y - yhat
  rmse <- sqrt(mean(resid^2))
  tss <- sum((y - mean(y))^2)
  r2 <- if (is.finite(tss) && tss > 0) 1 - sum(resid^2) / tss else NA_real_
  list(n = length(y), rmse = rmse, r2 = r2)
}

# -----------------------------------------------------------------------------
# Per-site processing
# -----------------------------------------------------------------------------
process_one_site <- function(site_id, csv_path, meta_path_hint = csv_path) {
  dt <- read_ameriflux_csv(csv_path)
  if (is.null(dt)) return(NULL)

  n_raw <- nrow(dt)

  dt[, TIMESTAMP_START := parse_timestamp(TIMESTAMP_START)]

  nm <- names(dt)
  ws_col <- NULL
  H <- V <- R <- NA_integer_

  ws_cand <- find_ws_candidates(nm)
  if (length(ws_cand)) {
    ws_col <- select_ws_variable(ws_cand)
    p <- parse_hvr(ws_col)
    H <- p$H; V <- p$V; R <- p$R
  } else if ("WS" %in% nm) {
    ws_col <- "WS"
    H <- 1L; V <- 1L; R <- 1L
  } else {
    warning("Site ", site_id, ": no WS_H_V_R or WS column - skipped.")
    return(NULL)
  }

  u_sig <- col_at_hvr("U_SIGMA", H, V, R, nm)
  if (is.na(u_sig)) u_sig <- find_first_matching("^U_?SIGMA", nm)
  v_sig <- col_at_hvr("V_SIGMA", H, V, R, nm)
  if (is.na(v_sig)) v_sig <- find_first_matching("^V_?SIGMA", nm)

  mo_col <- resolve_mo_length(dt, H, V, ws_col)
  wd_col <- resolve_wd(dt, H, V)

  qc_col <- find_qc_column(ws_col, nm)

  # Sensor height
  z_m <- parse_sensor_height_from_metadata(meta_path_hint, ws_col)
  if (!is.finite(z_m)) z_m <- DEFAULT_Z

  # Pull vectors
  WS <- dt[[ws_col]]
  U_SIGMA <- if (!is.na(u_sig)) dt[[u_sig]] else rep(NA_real_, n_raw)
  V_SIGMA <- if (!is.na(v_sig)) dt[[v_sig]] else rep(NA_real_, n_raw)
  MO <- if (!is.na(mo_col)) dt[[mo_col]] else rep(NA_real_, n_raw)

  if (!is.na(qc_col)) {
    qc <- dt[[qc_col]]
    qc_ok <- is.na(qc) | (qc <= QC_MAX)
  } else {
    qc_ok <- rep(TRUE, n_raw)
  }

  # ZETA
  ZETA <- z_m / MO
  ZETA[!is.finite(MO) | MO == 0] <- NA_real_
  ZETA <- pmax(-ZETA_CLIP, pmin(ZETA_CLIP, ZETA))

  # Derived (TI and sigma_theta scale as 1/WS; exclude ratio blow-ups at low WS)
  TI <- U_SIGMA / WS
  SIGMA_THETA_DEG <- (V_SIGMA / WS) * (180 / pi)

  STAB <- stability_class_from_zeta(ZETA)

  keep <- is.finite(WS) & WS >= WS_MIN &
    is.finite(U_SIGMA) & is.finite(V_SIGMA) &
    is.finite(MO) & MO != 0 &
    qc_ok &
    is.finite(TI) & TI >= 0 & TI <= TI_MAX_RETAIN &
    is.finite(SIGMA_THETA_DEG) & abs(SIGMA_THETA_DEG) <= SIGMA_THETA_ABS_MAX

  out <- data.table(
    SITE_ID           = site_id,
    TIMESTAMP_START   = dt$TIMESTAMP_START,
    WS                = WS,
    WS_VAR_NAME       = ws_col,
    SENSOR_HEIGHT_M   = z_m,
    U_SIGMA           = U_SIGMA,
    V_SIGMA           = V_SIGMA,
    TI                = TI,
    SIGMA_THETA_DEG   = SIGMA_THETA_DEG,
    MO_LENGTH         = MO,
    ZETA              = ZETA,
    STABILITY_CLASS   = STAB
  )

  out <- out[keep]

  n_keep <- nrow(out)

  mo_warn <- is.na(mo_col) ||
    (!is.na(mo_col) && all(is.na(dt[[mo_col]])))
  if (mo_warn && !is.na(mo_col)) {
    warning("Site ", site_id, ": MO_LENGTH column present but all NA - 0 rows retained.")
  } else if (mo_warn && is.na(mo_col)) {
    warning("Site ", site_id, ": no MO_LENGTH / L column - 0 rows retained.")
  }

  list(
    data = out,
    log = list(
      site_id = site_id,
      ws_var = ws_col,
      sensor_height_m = z_m,
      n_raw = n_raw,
      n_keep = n_keep,
      mo_warn = mo_warn
    )
  )
}

# -----------------------------------------------------------------------------
# Aggregation & fits
# -----------------------------------------------------------------------------
global_summary <- function(all_dt) {
  all_dt[, .(
    N_OBS = .N,
    ZETA_MEDIAN = median(ZETA, na.rm = TRUE),
    WS_MEAN = mean(WS, na.rm = TRUE),
    WS_STD = sd(WS, na.rm = TRUE),
    U_SIGMA_MEDIAN = median(U_SIGMA, na.rm = TRUE),
    U_SIGMA_MEAN = mean(U_SIGMA, na.rm = TRUE),
    U_SIGMA_STD = sd(U_SIGMA, na.rm = TRUE),
    U_SIGMA_P25 = quantile(U_SIGMA, 0.25, na.rm = TRUE),
    U_SIGMA_P75 = quantile(U_SIGMA, 0.75, na.rm = TRUE),
    U_SIGMA_P5 = quantile(U_SIGMA, 0.05, na.rm = TRUE),
    U_SIGMA_P95 = quantile(U_SIGMA, 0.95, na.rm = TRUE),
    TI_MEDIAN = median(TI, na.rm = TRUE),
    TI_MEAN = mean(TI, na.rm = TRUE),
    TI_STD = sd(TI, na.rm = TRUE),
    TI_P25 = quantile(TI, 0.25, na.rm = TRUE),
    TI_P75 = quantile(TI, 0.75, na.rm = TRUE),
    TI_P5 = quantile(TI, 0.05, na.rm = TRUE),
    TI_P95 = quantile(TI, 0.95, na.rm = TRUE),
    SIGMA_THETA_MEDIAN = median(SIGMA_THETA_DEG, na.rm = TRUE),
    SIGMA_THETA_MEAN = mean(SIGMA_THETA_DEG, na.rm = TRUE),
    SIGMA_THETA_STD = sd(SIGMA_THETA_DEG, na.rm = TRUE),
    SIGMA_THETA_P25 = quantile(SIGMA_THETA_DEG, 0.25, na.rm = TRUE),
    SIGMA_THETA_P75 = quantile(SIGMA_THETA_DEG, 0.75, na.rm = TRUE),
    SIGMA_THETA_P5 = quantile(SIGMA_THETA_DEG, 0.05, na.rm = TRUE),
    SIGMA_THETA_P95 = quantile(SIGMA_THETA_DEG, 0.95, na.rm = TRUE)
  ), by = STABILITY_CLASS]
}

bin_zeta_metrics <- function(all_dt, n_bin = 40L, zmin = -5, zmax = 5) {
  edges <- seq(zmin, zmax, length.out = n_bin + 1L)
  centers <- (edges[-1L] + edges[-length(edges)]) / 2
  d <- copy(all_dt)
  d[, bin_id := findInterval(ZETA, edges, all.inside = TRUE)]
  agg <- d[is.finite(ZETA) & bin_id >= 1L & bin_id <= n_bin, .(
    n = .N,
    u_med = median(U_SIGMA, na.rm = TRUE),
    ti_med = median(TI, na.rm = TRUE),
    th_med = median(SIGMA_THETA_DEG, na.rm = TRUE),
    u_p5 = quantile(U_SIGMA, 0.05, na.rm = TRUE),
    u_p25 = quantile(U_SIGMA, 0.25, na.rm = TRUE),
    u_p75 = quantile(U_SIGMA, 0.75, na.rm = TRUE),
    u_p95 = quantile(U_SIGMA, 0.95, na.rm = TRUE),
    ti_p5 = quantile(TI, 0.05, na.rm = TRUE),
    ti_p25 = quantile(TI, 0.25, na.rm = TRUE),
    ti_p75 = quantile(TI, 0.75, na.rm = TRUE),
    ti_p95 = quantile(TI, 0.95, na.rm = TRUE),
    th_p5 = quantile(SIGMA_THETA_DEG, 0.05, na.rm = TRUE),
    th_p25 = quantile(SIGMA_THETA_DEG, 0.25, na.rm = TRUE),
    th_p75 = quantile(SIGMA_THETA_DEG, 0.75, na.rm = TRUE),
    th_p95 = quantile(SIGMA_THETA_DEG, 0.95, na.rm = TRUE)
  ), by = bin_id]
  agg[, ZETA_BIN_CENTER := centers[bin_id]]
  agg
}

build_bestfit_table <- function(bin_dt, min_obs = MIN_OBS_PER_BIN) {
  fit_metric <- function(zc, ymed, n) {
    use <- n >= min_obs & is.finite(zc) & is.finite(ymed)
    zfit <- zc[use]
    yfit <- ymed[use]
    pp <- fit_piecewise_poly(zfit, yfit, 3L)
    poly_pred <- predict_piecewise_poly(zc, pp, 3L)
    poly_fit <- predict_piecewise_poly(zfit, pp, 3L)
    list(poly = pp, poly_y = poly_pred, stats = fit_quality(yfit, poly_fit))
  }

  u  <- fit_metric(bin_dt$ZETA_BIN_CENTER, bin_dt$u_med, bin_dt$n)
  ti <- fit_metric(bin_dt$ZETA_BIN_CENTER, bin_dt$ti_med, bin_dt$n)
  th <- fit_metric(bin_dt$ZETA_BIN_CENTER, bin_dt$th_med, bin_dt$n)

  list(
    table = data.table(
      ZETA_BIN_CENTER = bin_dt$ZETA_BIN_CENTER,
      U_SIGMA_MEDIAN = bin_dt$u_med,
      TI_MEDIAN = bin_dt$ti_med,
      SIGMA_THETA_MEDIAN = bin_dt$th_med,
      U_SIGMA_POLY_FIT = u$poly_y,
      TI_POLY_FIT = ti$poly_y,
      SIGMA_THETA_POLY_FIT = th$poly_y
    ),
    fits = list(U_SIGMA = u, TI = ti, SIGMA_THETA_DEG = th),
    fit_stats = rbindlist(list(
      data.table(METRIC = "U_SIGMA", FIT_N_BINS = u$stats$n,
                 FIT_RMSE = u$stats$rmse, FIT_R2 = u$stats$r2),
      data.table(METRIC = "TI", FIT_N_BINS = ti$stats$n,
                 FIT_RMSE = ti$stats$rmse, FIT_R2 = ti$stats$r2),
      data.table(METRIC = "SIGMA_THETA_DEG", FIT_N_BINS = th$stats$n,
                 FIT_RMSE = th$stats$rmse, FIT_R2 = th$stats$r2)
    ))
  )
}

lm_coef_string <- function(fit) {
  if (is.null(fit)) return(NA_character_)
  paste(format(coef(fit), scientific = TRUE, digits = 8), collapse = ";")
}

build_coefficients_table <- function(bf_out) {
  rows <- list()
  for (nm in names(bf_out$fits)) {
    f <- bf_out$fits[[nm]]
    st <- bf_out$fit_stats[METRIC == nm]
    rows[[length(rows) + 1L]] <- data.table(
      METRIC = nm,
      MODEL = "piecewise_cubic_best_fit",
      fit_n_bins = st$FIT_N_BINS,
      fit_rmse = st$FIT_RMSE,
      fit_r2 = st$FIT_R2,
      poly_neg_coef = lm_coef_string(f$poly$neg),
      poly_pos_coef = lm_coef_string(f$poly$pos)
    )
  }
  rbindlist(rows, fill = TRUE)
}

# -----------------------------------------------------------------------------
# Plot
# -----------------------------------------------------------------------------
#' Clip ribbon bounds so ymin <= ymax, enforce physical ranges, and NA-out thin bins.
clip_metric_ribbons <- function(lo, hi, lo2, hi2, lo_lim, hi_lim, n, n_min) {
  thin <- n < n_min
  clip2 <- function(a, b) {
    a <- pmax(a, lo_lim)
    b <- pmin(b, hi_lim)
    b <- pmax(b, a)
    list(a = a, b = b)
  }
  inner <- clip2(lo, hi)
  outer <- clip2(lo2, hi2)
  inner$a[thin] <- inner$b[thin] <- NA_real_
  outer$a[thin] <- outer$b[thin] <- NA_real_
  list(
    ymin = inner$a, ymax = inner$b,
    y5 = outer$a, y95 = outer$b
  )
}

make_bestfit_figure <- function(bin_dt, bestfit_dt, fit_stats_dt, path_out) {
  zlim <- c(-3, 3)
  bin_dt <- bin_dt[is.finite(ZETA_BIN_CENTER) &
                     ZETA_BIN_CENTER >= zlim[1] & ZETA_BIN_CENTER <= zlim[2]]

  u <- clip_metric_ribbons(
    bin_dt$u_p25, bin_dt$u_p75, bin_dt$u_p5, bin_dt$u_p95,
    0, RIBBON_U_SIGMA_MAX, bin_dt$n, RIBBON_MIN_N
  )
  th <- clip_metric_ribbons(
    bin_dt$th_p25, bin_dt$th_p75, bin_dt$th_p5, bin_dt$th_p95,
    -RIBBON_SIGMA_THETA_ABS, RIBBON_SIGMA_THETA_ABS, bin_dt$n, RIBBON_MIN_N
  )

  long_rib <- rbindlist(list(
    data.table(metric = "sigma_u", z = bin_dt$ZETA_BIN_CENTER,
               ymin = u$ymin, ymax = u$ymax, y5 = u$y5, y95 = u$y95),
    data.table(metric = "sigma_theta", z = bin_dt$ZETA_BIN_CENTER,
               ymin = th$ymin, ymax = th$ymax, y5 = th$y5, y95 = th$y95)
  ))
  long_rib <- long_rib[order(metric, z)]

  pts <- rbindlist(list(
    data.table(metric = "sigma_u", z = bestfit_dt$ZETA_BIN_CENTER,
               y = bestfit_dt$U_SIGMA_MEDIAN,
               poly = bestfit_dt$U_SIGMA_POLY_FIT),
    data.table(metric = "sigma_theta", z = bestfit_dt$ZETA_BIN_CENTER,
               y = bestfit_dt$SIGMA_THETA_MEDIAN,
               poly = bestfit_dt$SIGMA_THETA_POLY_FIT)
  ))
  pts <- pts[z >= zlim[1] & z <= zlim[2]]
  pts <- pts[order(metric, z)]

  fit_label <- function(metric_name) {
    st <- fit_stats_dt[METRIC == metric_name]
    if (nrow(st) == 0L) return(expression(atop("Piecewise cubic best fit", R^2)))
    substitute(
      atop("Piecewise cubic best fit", R^2 == r2),
      list(r2 = round(st$FIT_R2, 2))
    )
  }

  make_panel <- function(metric_key, y_label, fit_metric_name, unit_label,
                         show_x_title = TRUE) {
    rib <- long_rib[metric == metric_key]
    pt <- pts[metric == metric_key]
    p <- ggplot() +
      geom_ribbon(data = rib, aes(x = z, ymin = y5, ymax = y95),
                  fill = "grey90", alpha = 0.9, na.rm = TRUE) +
      geom_ribbon(data = rib, aes(x = z, ymin = ymin, ymax = ymax),
                  fill = "grey75", alpha = 0.85, na.rm = TRUE) +
      geom_point(data = pt, aes(x = z, y = y), size = 1.6, color = "black") +
      geom_line(data = pt, aes(x = z, y = poly),
                color = "#2166AC", linewidth = 0.8) +
      annotate(
        "text",
        x = zlim[2] - 0.03 * diff(zlim), y = Inf,
        label = as.character(as.expression(fit_label(fit_metric_name))),
        parse = TRUE,
        hjust = 1, vjust = 1.25, size = 3.2
      ) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
      coord_cartesian(xlim = zlim, clip = "off") +
      labs(
        x = if (show_x_title) expression(zeta == z / L) else NULL,
        y = y_label
      ) +
      theme_bw(base_size = 11) +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.title.x = element_text(size = 11),
        axis.title.y = element_text(size = 12),
        panel.border = element_rect(color = "black"),
        plot.margin = margin(6, 8, 4, 8)
      )
    if (!show_x_title) {
      p <- p + theme(axis.title.x = element_blank())
    }
    p
  }

  p_sigma_u <- make_panel(
    "sigma_u",
    expression(sigma[u]~"(m"~s^{-1}*")"),
    "U_SIGMA",
    "m/s",
    show_x_title = FALSE
  )
  p_sigma_theta <- make_panel(
    "sigma_theta",
    expression(sigma[theta]~"(degrees)"),
    "SIGMA_THETA_DEG",
    "deg",
    show_x_title = TRUE
  )

  grDevices::png(path_out, width = 8, height = 7.2, units = "in", res = 150)
  on.exit(grDevices::dev.off(), add = TRUE)
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(
    layout = grid::grid.layout(
      nrow = 2, ncol = 1,
      heights = grid::unit(c(1, 1.08), "null")
    )
  ))
  print(p_sigma_u, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(p_sigma_theta, vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 1))
  invisible(path_out)
}

# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------
main <- function() {
  if (!dir.exists(INPUT_DIR)) {
    stop("INPUT_DIR does not exist: ", INPUT_DIR,
         "\nCreate it or set INPUT_DIR in the script.")
  }
  dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

  zips <- find_zip_files(INPUT_DIR)
  if (length(zips) == 0L) {
    message("No .zip files under ", INPUT_DIR)
  }

  site_logs <- list()
  all_list <- list()

  for (zp in zips) {
    site_id <- parse_site_id_from_zip(zp)
    if (is.na(site_id)) next
    if (!is.null(DEBUG_SITE_ID) && site_id != DEBUG_SITE_ID) next

    dest <- unzip_to_site_dir(zp, INPUT_DIR, site_id)
    csv_path <- pick_preferred_csv(dest)
    if (is.na(csv_path)) {
      warning("No BASE HH/HR CSV in ", dest)
      next
    }

    message("Processing site ", site_id, " - ", basename(csv_path))
    res <- process_one_site(site_id, csv_path, meta_path_hint = csv_path)
    if (is.null(res)) next

    site_logs[[length(site_logs) + 1L]] <- res$log
    if (nrow(res$data) > 0L) all_list[[length(all_list) + 1L]] <- res$data

    lg <- res$log
    message(sprintf(
      "  WS=%s, z=%.2f m, raw=%d, kept=%d%s",
      lg$ws_var, lg$sensor_height_m, lg$n_raw, lg$n_keep,
      if (lg$mo_warn) " [WARN: MO_LENGTH]" else ""
    ))
  }

  if (length(all_list) == 0L) {
    warning("No site data combined - check CSV columns (U_SIGMA, V_SIGMA, MO_LENGTH).")
    return(invisible(NULL))
  }

  all_dt <- rbindlist(all_list, use.names = TRUE, fill = TRUE)

  fwrite(all_dt, file.path(OUTPUT_DIR, "ameriflux_wind_variability_all.csv"))

  gs <- global_summary(all_dt)
  fwrite(gs, file.path(OUTPUT_DIR, "ameriflux_wind_variability_global_summary.csv"))

  binned <- bin_zeta_metrics(all_dt, 40L, -ZETA_CLIP, ZETA_CLIP)
  bf_out <- build_bestfit_table(binned, MIN_OBS_PER_BIN)
  fwrite(bf_out$table, file.path(OUTPUT_DIR, "ameriflux_wind_variability_bestfit.csv"))
  fwrite(build_coefficients_table(bf_out),
         file.path(OUTPUT_DIR, "ameriflux_wind_variability_bestfit_coefficients.csv"))

  n_sites <- uniqueN(all_dt$SITE_ID)
  n_obs <- nrow(all_dt)

  make_bestfit_figure(
    binned, bf_out$table, bf_out$fit_stats,
    file.path(OUTPUT_DIR, "ameriflux_wind_variability_bestfit.png")
  )

  message(sprintf("Combined sites = %d; retained observations = %s",
                  n_sites, format(n_obs, big.mark = ",")))
  for (metric_name in c("U_SIGMA", "SIGMA_THETA_DEG")) {
    st <- bf_out$fit_stats[METRIC == metric_name]
    message(sprintf(
      "Piecewise cubic best fit (%s): R^2 = %.3f, RMSE = %.3f, bins = %d",
      metric_name, st$FIT_R2, st$FIT_RMSE, st$FIT_N_BINS
    ))
  }

  message("Done. Outputs in ", OUTPUT_DIR)
  invisible(TRUE)
}

if (sys.nframe() == 0L) {
  main()
}

