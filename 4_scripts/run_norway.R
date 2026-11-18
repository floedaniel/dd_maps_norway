# run_norway.R -- run the OPRAM port over Norway (seNorge2018, 1 km) and write GeoTIFFs.
# Port of OPRAM (Jon Yearsley, github.com/DrJonYearsley/OPRAM), GPL-3. Run from project root.
#
# Output per parameter set and simulation year (5_outputs/opram/<name>/):
#   <name>_<year>_emerge.tif  INT2S, 12 layers (start 1st of month): first-generation emergence day
#                             (continuous day index from 1 Jan of <year>; >365 = following year)
#   <name>_<year>_ngen.tif    FLT4S, 12 layers: generations per year (fractional last generation)
#
# Speed: input is the int16 binary store (extract_senorge.R). Each worker takes a block of land
# cells, reads every needed met year ONCE, and runs all simulation years and parameter sets.
# Windows: workers are separate R processes (future multisession); RAM ~2-3 GB per worker.
# Resumable: finished (parameter set, year) GeoTIFFs are skipped unless overwrite = TRUE.

suppressPackageStartupMessages({ library(terra); library(future.apply); library(ggplot2) })
source("4_scripts/opram_core.R")
source("4_scripts/extract_senorge.R")            # sn_read_bin(), ndays(), sn_bin_file()

sn_proc   <- "2_processed_data/senorge"
opram_out <- "5_outputs/opram"
last_meteo <- 2025L                              # last full year (seNorge v23.09 yearly files)

slug <- function(x) gsub("[^a-z0-9]+", "_", tolower(trimws(x)))
set_name <- function(p) slug(if (is.na(p$name)) sprintf("T%s_DD%s", p$T0, p$threshold) else p$name)
out_file <- function(nm, y, v) file.path(opram_out, nm, sprintf("%s_%d_%s.tif", nm, y, v))

# one block: land cells i (contiguous); returns nothing, writes <tmp>/<set>__<block>.rds
run_block <- function(b, i, met, todo, params, nms, lat, tmp) {
  n <- length(i); rd <- function(v) setNames(lapply(met, function(y) sn_read_bin(v, y, i[1], n, out = sn_proc)), met)
  tn <- rd("tn"); tx <- rd("tx")                               # integer, 0.01 degC, read once
  res <- lapply(todo, function(yrs) list(em = array(NA_integer_, c(n, 12L, length(yrs))),
                                         ng = array(NA_real_, c(n, 12L, length(yrs)))))
  mY <- max(vapply(params, `[[`, 1L, "maxYears"))
  for (y in sort(unique(unlist(todo)))) {                      # window built once per year, shared by all sets
    w <- as.character(opram_years(y, mY, last_meteo)); nd <- vapply(as.integer(w), ndays, 1L)
    TN <- do.call(rbind, tn[w]) / 100; TX <- do.call(rbind, tx[w]) / 100
    lo <- pmin(TN, TX); TX <- pmax(TN, TX); TN <- lo; rm(lo)   # Tmin <= Tmax once, not per set
    dates <- do.call(c, lapply(as.integer(w), function(z) as.Date(sprintf("%d-01-01", z)) + seq_len(ndays(z)) - 1L))
    for (j in which(vapply(todo, function(t) y %in% t, TRUE))) {
      p <- params[[j]]; r <- seq_len(sum(nd[seq_len(p$maxYears)])); k <- match(y, todo[[j]])
      o <- if (length(r) == nrow(TN)) opram_run(TN, TX, dates, lat, p, swap = FALSE)   # no copy
           else opram_run(TN[r, , drop = FALSE], TX[r, , drop = FALSE], dates[r], lat, p, swap = FALSE)
      res[[j]]$em[, , k] <- o$emerge; res[[j]]$ng[, , k] <- o$ngen
    }
  }
  for (j in which(lengths(todo) > 0))
    saveRDS(res[[j]], file.path(tmp, sprintf("%s__%04d.rds", nms[j], b)), compress = FALSE)
  invisible(NULL)
}

# assemble one parameter set from block files -> one pair of GeoTIFFs per year
assemble_set <- function(nm, yrs, nb, gr, tmp, p) {
  tmpl <- rast(nrows = gr$nrow, ncols = gr$ncol, extent = ext(gr$ext), crs = gr$crs)
  em <- array(NA_integer_, c(length(gr$land), 12L, length(yrs))); ng <- array(NA_real_, dim(em)); i0 <- 0L
  for (b in seq_len(nb)) {
    z <- readRDS(file.path(tmp, sprintf("%s__%04d.rds", nm, b))); n <- dim(z$em)[1]
    em[i0 + seq_len(n), , ] <- z$em; ng[i0 + seq_len(n), , ] <- z$ng; i0 <- i0 + n
  }
  stopifnot(i0 == length(gr$land))
  dir.create(file.path(opram_out, nm), FALSE, TRUE); saveRDS(p, file.path(opram_out, nm, "params.rds"))
  for (k in seq_along(yrs)) for (v in c("emerge", "ngen")) {
    r <- rast(tmpl, nlyrs = 12); names(r) <- paste0("start_", month.abb)
    vals <- matrix(NA_real_, ncell(tmpl), 12); vals[gr$land, ] <- if (v == "emerge") em[, , k] else ng[, , k]
    values(r) <- vals
    writeRaster(r, out_file(nm, yrs[k], v), overwrite = TRUE, gdal = "COMPRESS=DEFLATE",
                datatype = if (v == "emerge") "INT2S" else "FLT4S")
  }
  unlink(file.path(tmp, sprintf("%s__%04d.rds", nm, seq_len(nb))))
  invisible(nm)
}

# run simulation years x parameter sets over Norway
opram_norway <- function(years, params, workers = 10, block_cells = 10000, overwrite = FALSE) {
  if (!is.null(params$T0)) params <- list(params)
  gr <- readRDS(file.path(sn_proc, "grid.rds")); nms <- vapply(params, set_name, "")
  if (anyDuplicated(nms)) stop("duplicate parameter-set names: ", paste(nms[duplicated(nms)], collapse = ", "))
  # parameters changed since the maps were made (e.g. a corrected user species)? -> drop old maps, redo all
  key <- function(p) p[c("T0", "threshold", "diapause_daylength", "diapause_temperature", "method", "maxYears")]
  for (j in seq_along(params)) {
    pf <- file.path(opram_out, nms[j], "params.rds")
    if (file.exists(pf) && !identical(key(readRDS(pf)), key(params[[j]]))) {
      old <- list.files(file.path(opram_out, nms[j]), "\\.tif$", full.names = TRUE)
      message(sprintf("%s: parameters changed -> removing %d old map files and recomputing", nms[j], length(old)))
      unlink(c(old, pf))
    }
  }
  todo <- lapply(nms, function(nm) years[overwrite | !file.exists(out_file(nm, years, "emerge")) |
                                           !file.exists(out_file(nm, years, "ngen"))])
  if (!any(lengths(todo))) { message("all outputs exist"); return(invisible(NULL)) }
  met <- sort(unique(unlist(lapply(seq_along(params), function(j)
    lapply(todo[[j]], opram_years, params[[j]]$maxYears, last_meteo)))))
  need <- c(sn_bin_file("tn", met, sn_proc), sn_bin_file("tx", met, sn_proc))
  if (any(miss <- !file.exists(need))) stop("extract these years first (4_scripts/extract_senorge.R):\n  ",
                                            paste(need[miss], collapse = "\n  "))
  blocks <- split(seq_along(gr$land), ceiling(seq_along(gr$land) / block_cells))
  tmp <- file.path(opram_out, paste0("_tmp_", Sys.getpid()))   # per process: concurrent jobs never collide
  dir.create(tmp, FALSE, TRUE)
  old <- plan(multisession, workers = min(workers, length(blocks))); on.exit(plan(old))
  t0 <- Sys.time()
  message(sprintf("%d sets x %d years, %d blocks, %d workers, met years %d-%d",
                  sum(lengths(todo) > 0), length(unique(unlist(todo))), length(blocks), nbrOfWorkers(), min(met), max(met)))
  future_lapply(seq_along(blocks), function(b)
    run_block(b, blocks[[b]], met, todo, params, nms, gr$lat[blocks[[b]]], tmp),
    future.seed = TRUE, future.packages = "matrixStats")
  t1 <- Sys.time(); message(sprintf("model: %.1f min", difftime(t1, t0, units = "mins")))
  j <- which(lengths(todo) > 0)
  plan(multisession, workers = min(workers, 4L, length(j)))         # assembly: ~5 GB RAM per set (32 yrs); 6 workers starved a 63 GB machine
  future_lapply(j, function(jj) assemble_set(nms[jj], todo[[jj]], length(blocks), gr, tmp, params[[jj]]),
                future.seed = TRUE, future.packages = "terra")
  message(sprintf("write GeoTIFFs: %.1f min", difftime(Sys.time(), t1, units = "mins")))
  if (!length(list.files(tmp))) unlink(tmp, recursive = TRUE)
  invisible(nms)
}

# multi-year normal (OPRAM_calculate_average.jl): per-cell median across years for each start month;
# emergence median > 366*maxYears or all-missing -> NA and generations -> 0. Anomaly = year - normal.
opram_normals <- function(name, years = 1991:2020, maxYears = 4L) {
  nm <- slug(name); d <- file.path(opram_out, nm)
  f <- function(v) out_file(nm, years, v)
  if (any(miss <- !file.exists(c(f("emerge"), f("ngen"))))) stop("missing years: ", sum(miss), " files")
  med <- function(v) rast(lapply(1:12, function(m) median(rast(lapply(f(v), function(x) rast(x)[[m]])), na.rm = TRUE)))
  em <- med("emerge"); ng <- med("ngen")
  bad <- is.na(em) | em > 366 * maxYears
  em <- mask(em, bad, maskvalue = 1); ng <- ifel(bad & !is.na(ng), 0, ng)
  names(em) <- names(ng) <- paste0("start_", month.abb)
  tag <- sprintf("%s_normal_%d_%d", nm, min(years), max(years))
  writeRaster(em, file.path(d, paste0(tag, "_emerge.tif")), overwrite = TRUE, datatype = "FLT4S", gdal = "COMPRESS=DEFLATE")
  writeRaster(ng, file.path(d, paste0(tag, "_ngen.tif")), overwrite = TRUE, datatype = "FLT4S", gdal = "COMPRESS=DEFLATE")
  invisible(tag)
}

# probability maps: share of years in which >= k generations complete within the year, per start month.
# For k = 1 this equals P(first adults emerge within the same year) -> also valid for "emergence only" rows.
# Writes <name>_prob_<y1>_<y2>_ngen_ge<k>.tif (FLT4S, 12 layers, values 0-1).
opram_probability <- function(name, years = 1991:2020, k = c(1, 2)) {
  nm <- slug(name); f <- out_file(nm, years, "ngen")
  if (any(miss <- !file.exists(f))) stop(nm, ": missing ", sum(miss), " year files for ", min(years), "-", max(years))
  for (kk in k) {
    p <- rast(lapply(1:12, function(m) {
      s <- rast(lapply(f, function(x) rast(x)[[m]]))
      sum(s >= kk - 1e-9) / length(f)                  # ngen is float32; tolerance for exactly k
    }))
    names(p) <- paste0("start_", month.abb)
    writeRaster(p, file.path(opram_out, nm, sprintf("%s_prob_%d_%d_ngen_ge%d.tif", nm, min(years), max(years), kk)),
                overwrite = TRUE, datatype = "FLT4S", gdal = "COMPRESS=DEFLATE")
  }
  invisible(nm)
}

# compute one species from the species table (OPRAM presets + 3_metadata/species_parameters_user.csv);
# used by the app's "Compute" button (runs in a background R process via callr).
opram_compute_species <- function(name, years, workers = 6L, block_cells = 5000L) {
  source("4_scripts/species_db.R", local = TRUE)
  df <- sp_read_all("."); df <- df[df$species == name, , drop = FALSE]
  if (nrow(df) != 1) stop("species not found: ", name)
  p <- sp_params(df)[[1]]
  message(sprintf("STAGE model %s %d-%d", name, min(years), max(years)))
  opram_norway(years, p, workers = workers, block_cells = block_cells)
  if (all(1991:2020 %in% years)) {
    message("STAGE normals"); opram_normals(name, 1991:2020, p$maxYears)
    message("STAGE probability"); opram_probability(name, 1991:2020)
  }
  message("STAGE done")
  invisible(set_name(p))
}

# quick map of one layer (ggplot2 + viridis)
plot_opram <- function(file, layer = "start_Jan", title = basename(file)) {
  em <- grepl("emerge", file)
  df <- as.data.frame(rast(file)[[layer]], xy = TRUE, na.rm = TRUE); names(df)[3] <- "v"
  ggplot(df, aes(x, y, fill = v)) + geom_raster() + coord_equal() +
    scale_fill_viridis_c(name = if (em) "Emergence\n(day)" else "Generations",
                         option = if (em) "magma" else "viridis", direction = if (em) -1 else 1) +
    labs(title = title, subtitle = layer, x = NULL, y = NULL) + theme_minimal(base_size = 10) +
    theme(axis.text = element_blank(), panel.grid = element_blank())
}

if (sys.nframe() == 0) {
  a <- as.integer(commandArgs(TRUE)); yrs <- if (length(a)) a else 2010L
  ps <- list(opram_params(10, 300, name = "test_T10_DD300"),
             opram_params(8.3, 557, 14.5, 14.5, name = "test_T8.3_DD557_diapause"))
  opram_norway(yrs, ps)
}
