# UAS Gauss OSSE Publication Companion

This repository accompanies the manuscript **"Evaluation of area-source methane-emission quantification by Uncrewed Aerial Systems using Gauss's law"**.

The code implements an observing system simulation experiment (OSSE) for UAS-based methane flux retrieval over an idealized circular lake. A Gaussian-plume forward model generates synthetic methane enhancement fields, simulated UAS receptors sample a cylindrical control surface, and a cylindrical mass-balance inverse retrieves the lake-integrated flux.

## Repository Scope

1. Spatial sampling resolution under a constant atmosphere.
2. Nearby upwind lake contamination.
3. Sequential UAS sampling during changing atmospheric stability.
4. Retrieval-error regimes in nondimensional space.

## Layout

```text
R/
  utils.R
  01_synthetic_emissions.R
  02b_gaussian_forward.R
  03_drone_and_sensor_sim.R
  04_conley_inverse_flux.R
  06_visualization.R

scripts/
  run_publication_analyses.R
  process_ameriflux_wind_variability.R
  full_analysis/

output/
  Selected manuscript tables, reports, and figures.

docs/
  dimensionless_framework.md
```

## Requirements

Install R packages from CRAN:

```r
install.packages(c(
  "data.table",
  "ggplot2",
  "cowplot",
  "viridis",
  "scales",
  "progressr"
))
```

## Run Analyses

From the repository root:

```powershell
Rscript scripts/run_publication_analyses.R
```

Individual manuscript analyses can also be run directly:

```powershell
Rscript scripts/full_analysis/analyze_spatial_resolution_constant_atmosphere.R
Rscript scripts/full_analysis/analyze_nearby_lake_contamination_distance.R
Rscript scripts/full_analysis/analyze_temporal_measurement_changing_atmosphere.R
Rscript scripts/full_analysis/analyze_dimensionless_framework.R
```

Most scripts write CSV summaries, report markdown, and PNG figures into `output/`.

## License

Code is licensed under GPL-3.0-or-later. See `LICENSE`.
