# FinnPRIO: Degree-day pest phenology maps for Norway, 1 km

Degree-day pest phenology maps for Norway at 1 km resolution, with a Shiny app. Part of the FinnPRIO
project (repository: `dd_maps_norway`).

For each pest species and each year, the model computes, for development starting on the 1st of every month:

- **adult emergence date** of the first generation (up to 4 years ahead), and
- **number of generations per year** (with a fractional last generation),

plus, for two reference periods (1991–2020 and 2008–2022), **medians**, **anomalies** and
**probabilities** (share of years with ≥1 and ≥2 generations).

Species parameters: lower threshold T0, degree-day sum, degree-day method (single sine by default,
single triangle or simple average), optional upper threshold (horizontal cutoff), optional
photoperiod/temperature diapause, and an *emergence only* flag for parameter sets that cover part of
the life cycle (their generation counts are not shown).

New species parameters are only used after the values are confirmed in the source paper
(quote and page recorded); DOIs alone only prove that a source exists.

The maps show where development is **thermally possible**. The model has no mortality, winter
survival, host availability or microclimate, and is not yet validated against Norwegian observations.

## App

```r
shiny::runApp("4_scripts/ddmaps_app")   # from the project root
```

- **Maps**: species × year × start month; emergence date or generations; selected year, median or anomaly
  for a reference period, or probability of ≥1/≥2 generations. Click a 1 km cell for its monthly values
  and a 1991–2022 trend (with change per decade).
  The side table per start month: **Adults** (emergence date, "+1y" = following year), **Gen.** (generations
  that year), **±days** (days earlier (−) or later (+) than the reference-period median), **Years** (share of
  reference-period years with at least one generation by 31 December). Column names have hover help.
- **Point model**: live run for one place (its 1 km cell; latitude/longitude under the species, or click the
  map) with preset or custom parameters (method, upper threshold, diapause).
- **Species**: parameter table with references; add your own species and compute maps for them in the
  background (one job per species at a time); delete user species, optionally with their maps.

## Pipeline

| Script | Purpose |
|---|---|
| `4_scripts/download_senorge.R` | Download yearly seNorge v23.09 files from MET Norway's THREDDS archive (resumable) |
| `4_scripts/extract_senorge.R` | Inspect seNorge NetCDF; extract years to an int16 binary store (and optional GeoTIFFs) |
| `4_scripts/prep_norway_mask.R` | Norway mask and borders (GADM 4.1, downloaded via `geodata`) |
| `4_scripts/opram_core.R` | Model: degree-days (with upper-threshold cutoff), diapause, emergence, generations, normals |
| `4_scripts/run_norway.R` | Parallel run over all land cells → GeoTIFFs in `5_outputs/opram/`; normals and probabilities per reference period |
| `4_scripts/run_v2309.R` | Full batch: extract 1991–2025, all species 1991–2022, normals + probabilities for both periods |
| `4_scripts/run_1991_2020.R` | Earlier batch on the long single-variable files (superseded by `run_v2309.R`) |
| `4_scripts/species_db.R` | Species table: OPRAM presets + user species (`3_metadata/species_parameters_user.csv`) |
| `4_scripts/check_references_doi.R` | Find and verify DOIs for every reference in a degree-day parameter spreadsheet (Crossref/DataCite; author + year + species-in-title) |
| `4_scripts/extract_pdf_text.R` | Extract text from downloaded source papers (`1_raw_data/dd_pdf`) for checking parameter values |

Source papers are kept as `1_raw_data/dd_pdf/Author_Year_short_title.pdf` (not in the repository); `_renamed.csv`
in that folder lists the original file names with DOI or source.
| `4_scripts/tests/` | `testthat` tests: model vs. a literal port of the original algorithm, cutoff vs. numerical integration, app workflow |

Setup:

1. `Rscript 4_scripts/download_senorge.R` (yearly seNorge files, ~75 GB for 1991–2025). Optional `.Renviron`
   in the project root: `SENORGE_ARCHIVE_DIR="path/to/files"`.
2. `Rscript 4_scripts/extract_senorge.R 2010 2011 2012 2013`, or `4_scripts/run_v2309.R` for the full pipeline.
3. `Rscript 4_scripts/prep_norway_mask.R`
4. `Rscript 4_scripts/run_norway.R 2010`, then start the app.

Data and outputs are not in the repository (tens of GB). R packages: terra, ncdf4, future.apply,
matrixStats, ggplot2, shiny, bslib, callr, geodata, testthat.

## Data

- Climate: seNorge_2018 v23.09 daily minimum/maximum temperature, MET Norway
  (Lussana et al. 2019, *Earth Syst. Sci. Data* 11:1531), CC BY 4.0 / NLOD. Yearly files from the
  [THREDDS archive](https://thredds.met.no/thredds/catalog/senorge/seNorge_2018/Archive/catalog.html);
  download with `4_scripts/download_senorge.R`. Weather 1991–2025, maps 1991–2022.
- Boundaries: GADM 4.1 (gadm.org), downloaded at run time; not redistributed here.
- User species added so far: *Popillia japonica* (Ebbenga et al. 2022, doi:10.3389/finsc.2022.1075807) and
  *Agrilus planipennis* (Barker et al. 2023, doi:10.3389/finsc.2023.1239173), both EPPO A2.

## Acknowledgements and license

The model is an R port of **OPRAM** by Jon Yearsley
([github.com/DrJonYearsley/OPRAM](https://github.com/DrJonYearsley/OPRAM), GPL-3), the model behind
[pestmaps.ie](https://www.pestmaps.ie). The species presets in `3_metadata/species_parameters_OPRAM.csv`
come from that repository. Differences: the port matches OPRAM with `save1year = false`; any day with
Tmax < Tmin is swapped (303 cell-days in 1991–2025); additions: upper threshold, per-species method,
emergence-only flag, recent reference period, probability maps. The app uses the shared VKM / FinnPRIO theme.

This project is licensed under the GNU General Public License v3.0 (see `LICENSE`).
