#!/usr/bin/env Rscript
###############################################################################
# run_publication_analyses.R
#
# Run the manuscript companion analyses in sequence. Command-line arguments are
# forwarded to each analysis script, so lightweight settings such as
# --paper_stats_only or --output_stem=<name> can still be used where supported.
###############################################################################

proj_root <- if (file.exists("R/utils.R")) {
  normalizePath(".")
} else if (file.exists("../R/utils.R")) {
  normalizePath("..")
} else {
  stop("Cannot find project root. Run from the repository root.")
}

scripts_to_run <- file.path(
  proj_root,
  "scripts",
  "full_analysis",
  c(
    "analyze_spatial_resolution_constant_atmosphere.R",
    "analyze_nearby_lake_contamination_distance.R",
    "analyze_temporal_measurement_changing_atmosphere.R",
    "analyze_dimensionless_framework.R",
    "make_conley_methods_diagram_assets.R",
    "make_conley_face_diffusion_loss_schematic.R",
    "make_gaussian_plume_stability_comparison.R"
  )
)

missing_scripts <- scripts_to_run[!file.exists(scripts_to_run)]
if (length(missing_scripts)) {
  stop("Missing analysis script(s):\n", paste(missing_scripts, collapse = "\n"))
}

all_args <- commandArgs(trailingOnly = TRUE)
dry_run <- "--dry_run" %in% all_args || "--dry-run" %in% all_args
forwarded_args <- setdiff(all_args, c("--dry_run", "--dry-run"))
rscript_bin <- file.path(
  R.home("bin"),
  if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript"
)

if (dry_run) {
  message("Dry run: found publication scripts in execution order:")
  for (script_path in scripts_to_run) {
    message("  - ", normalizePath(script_path, winslash = "/", mustWork = TRUE))
  }
  quit(status = 0)
}

for (script_path in scripts_to_run) {
  message("========================================")
  message("Running: ", basename(script_path))
  if (length(forwarded_args)) {
    message("Args:    ", paste(forwarded_args, collapse = " "))
  }
  message("========================================")

  status <- system2(
    command = rscript_bin,
    args = c(script_path, forwarded_args),
    stdout = "",
    stderr = ""
  )

  if (!identical(status, 0L)) {
    stop("Script failed with exit status ", status, ": ", basename(script_path))
  }
}

message("Publication companion analyses completed successfully.")

