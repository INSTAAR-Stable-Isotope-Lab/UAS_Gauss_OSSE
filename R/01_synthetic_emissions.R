###############################################################################
# 01_synthetic_emissions.R - Build lake-averaged methane emission fields
#
# Scientific notes:
#   - Ebullition is treated as part of an effective, lake-averaged surface
#     flux (mg CH4 m-2 d-1).  We do NOT resolve hotspots or fine-scale
#     spatial structure; any ebullition is assumed homogenised in space and
#     time at the lake scale.
#   - The OSSE is a publication companion for Gaussian-forward /
#     Conley-inverse manuscript experiments.
#
# Inputs :
#   - Lake size configurations (from utils.R: LAKE_CONFIGS).
#   - Emission defaults (from utils.R: EMISSION_DEFAULTS).
#
# Optional outputs when save_output = TRUE:
#   - data/raw/emissions_<size>.csv
#   - data/raw/lake_mean_flux_<size>.csv
###############################################################################

if (!exists("project_root", mode = "function")) {
  proj_root <- if (file.exists("R/utils.R")) "." else ".."
  source(file.path(proj_root, "R", "utils.R"))
}

# -----------------------------------------------------------------------------
# Main function: generate and save emission fields for specified lake configs
# -----------------------------------------------------------------------------
#' Generate emission fields for one or more lake configurations.
#'
#' @param lake_configs  List of lake config lists, each with $name, $diameter_m.
#'                      Default: LAKE_CONFIGS.
#' @param emission_params  List of emission parameters (see EMISSION_DEFAULTS).
#' @param save_output  Logical; write generated intermediate CSV files to data/raw?
#' @param base_flux_override  Optional numeric: if provided and a single lake
#'   is requested, use this as base_flux (mg CH4 m-2 d-1) instead of
#'   emission_params$base_flux_default * size_factor (e.g. for sensitivity).
#' @return A named list of lists, each containing $flux_field, $lake_mean,
#'         $lake_geom, $base_flux, $diameter_m.
generate_all_emissions <- function(lake_configs = LAKE_CONFIGS,
                                   emission_params = EMISSION_DEFAULTS,
                                   save_output = TRUE,
                                   base_flux_override = NULL) {

  results <- list()

  # times_hour must be part of the emissions-catalog key because it controls
  # the evaluated diurnal sampling (and therefore the returned flux/lake_mean).
  times_hour <- seq(0, 23.5, by = 0.5)

  cache_root <- proj_path("cache", "emissions")
  dir.create(cache_root, showWarnings = FALSE, recursive = TRUE)

  hash_emissions_config <- function(full_config_for_hash) {
    # Use a deterministic stable hash; prefer digest when available.
    old_digits <- getOption("digits")
    options(digits = 17)
    on.exit(options(digits = old_digits), add = TRUE)

    dput_str <- paste(capture.output(
      dput(full_config_for_hash, control = "all")
    ), collapse = "\n")

    if (requireNamespace("digest", quietly = TRUE)) {
      digest::digest(dput_str, algo = "md5")
    } else {
      # Fallback: deterministic but not cryptographic.
      sprintf("%08d", string_to_seed(dput_str))
    }
  }

  for (lk_name in names(lake_configs)) {
    lk <- lake_configs[[lk_name]]
    message(sprintf("  Emissions: %s lake (d=%d m)...",
                    lk$name, lk$diameter_m))

    # Build lake grid
    res_m <- if (lk$diameter_m <= 150) 5 else if (lk$diameter_m <= 400) 10 else 20
    lake_geom <- define_lake_geometry(lk$diameter_m, resolution_m = res_m)

    # Base flux: override for single-lake sensitivity, else size-scaled default
    base_flux_override_present <- !is.null(base_flux_override)
    base_flux_override_numeric <- if (is.null(base_flux_override)) NULL else as.numeric(base_flux_override)
    override_applied <- base_flux_override_present &&
      length(lake_configs) == 1L &&
      is.finite(base_flux_override_numeric)

    if (override_applied) {
      base_flux <- base_flux_override_numeric
    } else {
      size_factor <- 1 - 0.3 * (lk$diameter_m - 100) / 700
      base_flux <- emission_params$base_flux_default * size_factor
    }

    # Full config key for caching.
    # Includes:
    #   - lake_name, diameter_m
    #   - derived grid resolution (res_m) from the diameter thresholds
    #   - times_hour sampling grid
    #   - emission_params fields used by generate_effective_flux_field()
    #   - base_flux_override (NULL vs numeric) + override_applied
    #   - resolved base_flux actually used in this run
    #   - EMISSIONS_CATALOG_SOFTWARE_VERSION
    full_config_for_hash <- list(
      emissions_catalog_software_version = EMISSIONS_CATALOG_SOFTWARE_VERSION,
      lake_name    = lk$name,
      diameter_m   = lk$diameter_m,
      res_m         = res_m,
      times_hour   = times_hour,
      emission_params = list(
        littoral_factor        = emission_params$littoral_factor,
        spatial_corr_length_m = emission_params$spatial_corr_length_m,
        diel_amplitude        = emission_params$diel_amplitude,
        diel_phase_hour       = emission_params$diel_phase_hour
      ),
      base_flux_override         = base_flux_override_numeric, # NULL vs numeric preserved
      base_flux_override_applied = override_applied,
      resolved_base_flux         = base_flux
    )

    cfg_hash <- hash_emissions_config(full_config_for_hash)

    lake_tag <- gsub("[^A-Za-z0-9._-]", "_", lk$name)
    # Naming scheme mirrors the TD-catalog style:
    #   emcat__lake<name>__res<res_m>__t<n_times>__osseV<ver>__h<hash>/
    #
    # IMPORTANT: emissions caching does NOT depend on sensor config, so
    # sensor details are intentionally excluded from both the name and the
    # reuse check criteria.
    cache_subdir <- sprintf(
      "emcat__lake%s__res%d__t%d__osseV%s__h%s",
      lake_tag,
      res_m,
      length(times_hour),
      EMISSIONS_CATALOG_SOFTWARE_VERSION,
      cfg_hash
    )
    cache_dir <- file.path(cache_root, cache_subdir)

    flux_field <- NULL
    lake_mean <- NULL
    cache_hit <- FALSE
    cache_hit_source <- NULL

    cache_cfg_rds <- file.path(cache_dir, "emissions_config.rds")
    cache_flux_rds <- file.path(cache_dir, "flux_field.rds")
    cache_lake_mean_rds <- file.path(cache_dir, "lake_mean.rds")

    # Backward-compatible: attempt to load from the legacy naming scheme
    # (created by earlier implementations).
    legacy_cache_subdir <- sprintf(
      "emcat__%s__osseV%s__h%s",
      lake_tag,
      EMISSIONS_CATALOG_SOFTWARE_VERSION,
      cfg_hash
    )
    legacy_cache_dir <- file.path(cache_root, legacy_cache_subdir)
    legacy_cache_cfg_rds <- file.path(legacy_cache_dir, "emissions_config.rds")
    legacy_cache_flux_rds <- file.path(legacy_cache_dir, "flux_field.rds")
    legacy_cache_lake_mean_rds <- file.path(legacy_cache_dir, "lake_mean.rds")

    if (file.exists(cache_cfg_rds) &&
        file.exists(cache_flux_rds) &&
        file.exists(cache_lake_mean_rds)) {
      cached_cfg <- tryCatch(readRDS(cache_cfg_rds), error = function(e) NULL)

      # Cache match criteria:
      #  - version gate
      #  - config hash match (config includes NULL vs numeric override,
      #    resolved base_flux, derived res_m, times_hour, emission_params)
      if (!is.null(cached_cfg) &&
          is.list(cached_cfg) &&
          identical(cached_cfg$emissions_catalog_software_version, EMISSIONS_CATALOG_SOFTWARE_VERSION) &&
          identical(cached_cfg$full_config_hash, cfg_hash)) {
        flux_field <- tryCatch(readRDS(cache_flux_rds), error = function(e) NULL)
        lake_mean <- tryCatch(readRDS(cache_lake_mean_rds), error = function(e) NULL)
        cache_hit <- !is.null(flux_field) && !is.null(lake_mean)
        if (cache_hit) cache_hit_source <- "new"
      }
    }

    if (!cache_hit &&
        file.exists(legacy_cache_cfg_rds) &&
        file.exists(legacy_cache_flux_rds) &&
        file.exists(legacy_cache_lake_mean_rds)) {
      cached_cfg <- tryCatch(readRDS(legacy_cache_cfg_rds), error = function(e) NULL)
      if (!is.null(cached_cfg) &&
          is.list(cached_cfg) &&
          identical(cached_cfg$emissions_catalog_software_version, EMISSIONS_CATALOG_SOFTWARE_VERSION) &&
          identical(cached_cfg$full_config_hash, cfg_hash)) {
        flux_field <- tryCatch(readRDS(legacy_cache_flux_rds), error = function(e) NULL)
        lake_mean <- tryCatch(readRDS(legacy_cache_lake_mean_rds), error = function(e) NULL)
        cache_hit <- !is.null(flux_field) && !is.null(lake_mean)
        if (cache_hit) cache_hit_source <- "legacy"
      }
    }

    if (cache_hit) {
      message(sprintf("    Using cached emissions catalog (h=%s).", cfg_hash))
      # Promote legacy cache directories into the new naming scheme so
      # future runs are faster to interpret (no recomputation).
      if (identical(cache_hit_source, "legacy") && !dir.exists(cache_dir)) {
        dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
        cache_cfg <- list(
          emissions_catalog_software_version = EMISSIONS_CATALOG_SOFTWARE_VERSION,
          full_config_hash                   = cfg_hash,
          full_config                        = full_config_for_hash, # authoritative
          created_at                         = as.character(Sys.time())
        )
        saveRDS(cache_cfg, cache_cfg_rds)
        saveRDS(flux_field, cache_flux_rds)
        saveRDS(lake_mean, cache_lake_mean_rds)

        cache_cfg_txt <- file.path(cache_dir, "emissions_config.txt")
        txt_lines <- c(
          sprintf("emissions_catalog_software_version: %s", EMISSIONS_CATALOG_SOFTWARE_VERSION),
          sprintf("full_config_hash: %s", cfg_hash),
          sprintf("created_at: %s", as.character(cache_cfg$created_at)),
          "",
          "full_config (dput):",
          paste(capture.output(dput(full_config_for_hash, control = "all")), collapse = "\n")
        )
        writeLines(txt_lines, cache_cfg_txt)
      }
    } else {
      if (dir.exists(cache_dir)) {
        message(sprintf("    Cache miss or mismatch; regenerating emissions (h=%s).", cfg_hash))
      } else {
        message(sprintf("    Cache miss; generating emissions (h=%s).", cfg_hash))
      }

      # Generate effective flux field
      flux_field <- generate_effective_flux_field(
        lake_geom             = lake_geom,
        base_flux             = base_flux,
        littoral_factor       = emission_params$littoral_factor,
        spatial_corr_length_m = emission_params$spatial_corr_length_m,
        diel_amplitude        = emission_params$diel_amplitude,
        diel_phase_hour       = emission_params$diel_phase_hour,
        times_hour            = times_hour
      )

      # Compute lake-mean time series
      lake_mean <- flux_field[, .(
        lake_mean_flux = mean(flux_mg_m2_d),
        flux_sd        = sd(flux_mg_m2_d),
        flux_cv        = sd(flux_mg_m2_d) / mean(flux_mg_m2_d),
        n_cells        = .N
      ), by = hour]

      lake_mean[, `:=`(
        lake_name   = lk$name,
        diameter_m  = lk$diameter_m,
        base_flux   = base_flux
      )]

      # Write cache (authoritative)
      dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
      cache_cfg <- list(
        emissions_catalog_software_version = EMISSIONS_CATALOG_SOFTWARE_VERSION,
        full_config_hash                   = cfg_hash,
        full_config                        = full_config_for_hash, # saved for human inspection + integrity
        created_at                         = as.character(Sys.time())
      )

      saveRDS(cache_cfg, cache_cfg_rds)
      saveRDS(flux_field, cache_flux_rds)
      saveRDS(lake_mean, cache_lake_mean_rds)

      cache_cfg_txt <- file.path(cache_dir, "emissions_config.txt")
      txt_lines <- c(
        sprintf("emissions_catalog_software_version: %s", EMISSIONS_CATALOG_SOFTWARE_VERSION),
        sprintf("full_config_hash: %s", cfg_hash),
        sprintf("created_at: %s", as.character(cache_cfg$created_at)),
        "",
        "full_config (dput):",
        paste(capture.output(dput(full_config_for_hash, control = "all")), collapse = "\n")
      )
      writeLines(txt_lines, cache_cfg_txt)
    }

    if (save_output) {
      outdir <- proj_path("data", "raw")
      dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
      fwrite(flux_field,
             file.path(outdir, sprintf("emissions_%s.csv", lk$name)))
      fwrite(lake_mean,
             file.path(outdir, sprintf("lake_mean_flux_%s.csv", lk$name)))
    }

    results[[lk_name]] <- list(
      flux_field = flux_field,
      lake_mean  = lake_mean,
      lake_geom  = lake_geom,
      base_flux  = base_flux,
      diameter_m = lk$diameter_m
    )
  }

  message("  Emission fields generated for all lakes.")
  invisible(results)
}


# -----------------------------------------------------------------------------
# Generate emissions for a single lake (convenience wrapper)
# -----------------------------------------------------------------------------
#' Generate emission field for a single lake size.
#'
#' @param lake_name  One of "small", "medium", "large".
#' @param save_output  Logical.
#' @param base_flux_override  Optional numeric (mg CH4 m-2 d-1). When set, used
#'   as the lake-mean flux instead of the size-scaled default (e.g. for sensitivity).
#' @return List with $flux_field, $lake_mean, $lake_geom, $base_flux, $diameter_m.
generate_single_lake_emissions <- function(lake_name = "small",
                                           save_output = TRUE,
                                           base_flux_override = NULL) {
  if (!lake_name %in% names(LAKE_CONFIGS)) {
    stop("Unknown lake_name: ", lake_name,
         ". Choose from: ", paste(names(LAKE_CONFIGS), collapse = ", "))
  }
  lk_list <- LAKE_CONFIGS[lake_name]
  generate_all_emissions(lake_configs = lk_list, save_output = save_output,
                          base_flux_override = base_flux_override)[[1]]
}


# -----------------------------------------------------------------------------
# If run as a script, generate all emission fields
# -----------------------------------------------------------------------------
if (sys.nframe() == 0) {
  proj_root <- project_root()
  message("=== 01: Generating synthetic emission fields ===")
  all_emissions <- generate_all_emissions()
  message("Done.")
}

