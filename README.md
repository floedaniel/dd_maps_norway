# DD-maps Norway

Degree-day pest phenology maps for Norway at 1 km resolution, with a Shiny app.

For each pest species (base temperature T0, degree-day threshold, optional photoperiod/temperature
diapause) and each year, the model computes, for development starting on the 1st of every month:

- **adult emergence date** of the first generation (up to 4 years ahead), and
- **number of generations per year** (with a fractional last generation),

plus 1991–2020 medians and anomalies. Daily degree-days are computed from Tmin/Tmax with the
single-sine (default), single-triangle or average method.

The maps show where development is **thermally possible**. The model has no upper temperature
threshold, mortality, winter survival or host availability, and is not yet validated against
Norwegian observations.

## App

```r
shiny::runApp("4_scripts/ddmaps_app")   # from the project root
```

- **Maps**: species × year × start month; emergence date or generations; year, 30-yr median or anomaly; click a 1 km cell for details.
- **Point model**: live run for one location with preset or custom parameters.
- **Species**: parameter table with references; add your own species and compute maps for them in the background.

## Pipeline

| Script | Purpose |
|---|---|
| `4_scripts/extract_senorge.R` | Inspect seNorge2018 NetCDF; extract years to an int16 binary store (and optional GeoTIFFs) |
| `4_scripts/prep_norway_mask.R` | Norway mask and borders (GADM 4.1, downloaded via `geodata`) |
| `4_scripts/opram_core.R` | Model: degree-days, diapause, emergence, generations, normals |
| `4_scripts/run_norway.R` | Parallel run over all land cells; writes GeoTIFFs to `5_outputs/opram/` |
| `4_scripts/run_1991_2020.R` | Full batch: extract 1991–2021, all species 1991–2020, 30-yr normals |
| `4_scripts/species_db.R` | Species table (presets + user-added species) |
| `4_scripts/tests/` | `testthat` tests: model vs. a literal port of the original algorithm, app workflow |

Setup:

1. Download seNorge2018 daily TN and TX NetCDF (MET Norway) and set their folder in a `.Renviron`
   file in the project root: `SENORGE_DIR="path/to/seNorge"`.
2. `Rscript 4_scripts/extract_senorge.R 2010 2011 2012 2013` (or run `4_scripts/run_1991_2020.R`).
3. `Rscript 4_scripts/prep_norway_mask.R`
4. `Rscript 4_scripts/run_norway.R 2010`, then start the app.

Data and outputs are not in the repository (tens of GB). R packages: terra, ncdf4, future.apply,
matrixStats, ggplot2, shiny, bslib, callr, geodata, testthat.

## Data

- Climate: seNorge_2018 daily minimum/maximum temperature, MET Norway
  (Lussana et al. 2019, *Earth Syst. Sci. Data* 11:1531), CC BY 4.0 / NLOD.
- Boundaries: GADM 4.1 (gadm.org), downloaded at run time; not redistributed here.

## Acknowledgements and license

The model is an R port of **OPRAM** by Jon Yearsley
([github.com/DrJonYearsley/OPRAM](https://github.com/DrJonYearsley/OPRAM), GPL-3), the model behind
[pestmaps.ie](https://www.pestmaps.ie). The species presets in `3_metadata/species_parameters_OPRAM.csv`
come from that repository. Differences: the port matches OPRAM with `save1year = false`; days with
Tmax < Tmin are swapped.

This project is licensed under the GNU General Public License v3.0 (see `LICENSE`).
