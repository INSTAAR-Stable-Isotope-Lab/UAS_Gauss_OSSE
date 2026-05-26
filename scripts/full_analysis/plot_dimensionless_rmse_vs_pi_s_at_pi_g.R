#!/usr/bin/env Rscript
###############################################################################
# plot_dimensionless_rmse_vs_pi_s_at_pi_g.R
#
# Build a line plot from the already computed
# dimensionless-framework binned summary. The plot shows retrieval RMSE across
# Pi_S at one selected Pi_G bin, with within-bin absolute-error standard
# deviation shown as error bars.
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

if (!exists("theme_osse", mode = "function")) {
  theme_osse <- function(base_size = 11) {
    theme_minimal(base_size = base_size) +
      theme(
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey92", colour = NA),
        legend.position = "bottom",
        plot.title = element_text(face = "bold", size = base_size + 1)
      )
  }
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

parse_bool <- function(x) {
  if (is.null(x) || !nzchar(x)) return(FALSE)
  tolower(trimws(x)) %in% c("true", "t", "1", "yes", "y")
}

parse_transition_vec <- function(x) {
  if (is.null(x) || !nzchar(x)) return(character(0))
  trimws(strsplit(x, ",", fixed = TRUE)[[1]])
}

format_target <- function(x) {
  gsub("\\.", "p", format(round(x, 3), trim = TRUE, scientific = FALSE))
}

cli <- parse_named_args(commandArgs(trailingOnly = TRUE))

INPUT_CSV <- if (!is.null(cli$input) && nzchar(cli$input)) {
  cli$input
} else {
  file.path(proj_root, "output", "dimensionless_framework_binned_summary.csv")
}
if (!file.exists(INPUT_CSV)) {
  stop("Binned summary CSV not found: ", INPUT_CSV,
       "\nRun analyze_dimensionless_framework.R first.")
}

PI_G_TARGET <- if (!is.null(cli$pi_g_target)) as.numeric(cli$pi_g_target) else 0.99
MIN_BIN_N <- if (!is.null(cli$min_bin_n)) as.integer(cli$min_bin_n) else 5L
TRANSITIONS <- parse_transition_vec(cli$transitions)
LOG_Y <- if (!is.null(cli$log_y)) parse_bool(cli$log_y) else FALSE
ERROR_SD_COL <- if (!is.null(cli$error_sd_col) && nzchar(cli$error_sd_col)) {
  cli$error_sd_col
} else {
  "abs_error_sd_pct"
}

default_stem <- sprintf("dimensionless_framework_rmse_vs_pi_s_pi_g_%s",
                        format_target(PI_G_TARGET))
OUTPUT_PNG <- if (!is.null(cli$output) && nzchar(cli$output)) {
  cli$output
} else {
  file.path(proj_root, "output", paste0(default_stem, ".png"))
}
OUTPUT_CSV <- if (!is.null(cli$output_csv) && nzchar(cli$output_csv)) {
  cli$output_csv
} else {
  file.path(proj_root, "output", paste0(default_stem, ".csv"))
}

stopifnot(is.finite(PI_G_TARGET))
stopifnot(is.finite(MIN_BIN_N), MIN_BIN_N >= 1L)


# -----------------------------------------------------------------------------
# 2. Select the Pi_G bin and build plot data
# -----------------------------------------------------------------------------
dt <- fread(INPUT_CSV)
required_cols <- c("transition", "pi_s_center", "pi_g_bin", "pi_g_min",
                   "pi_g_center", "pi_g_max", "rmse_pct", "bin_n",
                   ERROR_SD_COL)
missing_cols <- setdiff(required_cols, names(dt))
if (length(missing_cols)) {
  stop("Input CSV is missing required column(s): ",
       paste(missing_cols, collapse = ", "))
}

numeric_cols <- c("pi_s_center", "pi_g_bin", "pi_g_min", "pi_g_center",
                  "pi_g_max", "rmse_pct", "bin_n", ERROR_SD_COL)
for (col in numeric_cols) {
  set(dt, j = col, value = as.numeric(dt[[col]]))
}

if (length(TRANSITIONS)) {
  dt <- dt[transition %in% TRANSITIONS]
}
if (!nrow(dt)) stop("No rows remain after transition filtering.")

g_bins <- unique(dt[, .(pi_g_bin, pi_g_min, pi_g_center, pi_g_max)])
g_bins <- g_bins[is.finite(pi_g_bin) & is.finite(pi_g_center)]
g_bins[, contains_target := pi_g_min <= PI_G_TARGET & PI_G_TARGET <= pi_g_max]
g_bins[, center_distance := abs(pi_g_center - PI_G_TARGET)]

if (any(g_bins$contains_target)) {
  selected_bin <- g_bins[contains_target == TRUE][order(center_distance)][1]
  selection_mode <- "containing"
} else {
  selected_bin <- g_bins[order(center_distance)][1]
  selection_mode <- "nearest"
}

plot_dt <- dt[
  pi_g_bin == selected_bin$pi_g_bin &
    is.finite(pi_s_center) &
    is.finite(rmse_pct) &
    is.finite(get(ERROR_SD_COL)) &
    bin_n >= MIN_BIN_N
]
if (!nrow(plot_dt)) {
  stop("No rows remain for selected Pi_G bin and min_bin_n = ", MIN_BIN_N)
}

setorder(plot_dt, transition, pi_s_center)
plot_dt[, rmse_sd_pct := get(ERROR_SD_COL)]
plot_dt[, rmse_ymin_pct := pmax(0, rmse_pct - rmse_sd_pct)]
plot_dt[, rmse_ymax_pct := rmse_pct + rmse_sd_pct]
plot_dt[, rmse_ymin_plot_pct := rmse_ymin_pct]
if (isTRUE(LOG_Y)) {
  min_positive_y <- min(plot_dt$rmse_pct[is.finite(plot_dt$rmse_pct) & plot_dt$rmse_pct > 0],
                        na.rm = TRUE)
  plot_dt[, rmse_ymin_plot_pct := pmax(rmse_ymin_pct, min_positive_y / 10)]
}
plot_dt[, pi_g_target := PI_G_TARGET]
plot_dt[, pi_g_selection_mode := selection_mode]
plot_dt[, selected_pi_g_bin := selected_bin$pi_g_bin]
plot_dt[, selected_pi_g_min := selected_bin$pi_g_min]
plot_dt[, selected_pi_g_center := selected_bin$pi_g_center]
plot_dt[, selected_pi_g_max := selected_bin$pi_g_max]
plot_dt[, error_bar_definition := ERROR_SD_COL]


# -----------------------------------------------------------------------------
# 3. Plot and write outputs
# -----------------------------------------------------------------------------
base_theme <- theme_osse(base_size = 11) +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

p <- ggplot(plot_dt, aes(x = pi_s_center, y = rmse_pct,
                         colour = transition, group = transition)) +
  geom_errorbar(aes(ymin = rmse_ymin_plot_pct, ymax = rmse_ymax_pct),
                width = 0, linewidth = 0.35, alpha = 0.65) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2.1) +
  scale_x_log10() +
  labs(
    x = expression(Pi[S]),
    y = "Retrieval RMSE (%)",
    colour = "Stability Class\nTransition"
  ) +
  base_theme

if (isTRUE(LOG_Y)) {
  p <- p + scale_y_log10()
}

dir.create(dirname(OUTPUT_PNG), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(OUTPUT_CSV), recursive = TRUE, showWarnings = FALSE)
ggsave_osse(OUTPUT_PNG, p, width = 8.5, height = 5.2, dpi = 220)
fwrite(plot_dt, OUTPUT_CSV)

message("Saved RMSE line plot: ", normalizePath(OUTPUT_PNG, mustWork = FALSE))
message("Saved plot data CSV:  ", normalizePath(OUTPUT_CSV, mustWork = FALSE))
message(sprintf(
  "Selected Pi_G bin: target %.3f, %s bin %.0f [%.3f, %.3f], center %.3f",
  PI_G_TARGET, selection_mode, selected_bin$pi_g_bin,
  selected_bin$pi_g_min, selected_bin$pi_g_max, selected_bin$pi_g_center
))
message(sprintf(
  "Rows plotted: %d | min bin_n: %d | error bars: +/- %s",
  nrow(plot_dt), MIN_BIN_N, ERROR_SD_COL
))

