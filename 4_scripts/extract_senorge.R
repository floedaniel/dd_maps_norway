# extract_senorge.R -- inspect seNorge2018 TN/TX NetCDF and extract single years.
# Data: MET Norway seNorge_2018 (Lussana et al. 2019, ESSD 11:1531), CC BY 4.0 / NLOD.
# Part of an R port of OPRAM (github.com/DrJonYearsley/OPRAM). License: GPL-3. Run from project root.
#
# Outputs per year and variable (2_processed_data/senorge/):
#   bin/<v>_<year>.bin  model input: int16 (0.01 degC), land cells only, cell-major [days x land]
#                       -> any contiguous range of land cells is one seek + one readBin. NA = -32768.
#   <v>/<v>_<year>.tif  optional float32 GeoTIFF (one layer per day) for GIS use
#   grid.tif / grid.rds lat, lon, land mask / land cell index
#
# NetCDF: ~40 GB, one chunk per day [X=1195, Y=1550, 1] -> always read by day. Time is
# "hours since 1900-01-01", standard calendar (leap days present), stamped 18 UTC; day D covers
# 18 UTC D-1 .. 18 UTC D. ncdf4 returns [X, Y] with Y north->south = terra cell order.

suppressPackageStartupMessages({ library(ncdf4); library(terra) })

# folder with the seNorge2018 NetCDF files; set SENORGE_DIR in a (git-ignored) .Renviron in the project root
sn_dir  <- Sys.getenv("SENORGE_DIR", "1_raw_data/senorge")
sn_out  <- "2_processed_data/senorge"   # relative to project root (phenology_model.Rproj)
sn_files <- c(tn = file.path(sn_dir, "seNorge2018_TN_1957_2022.nc"),
              tx = file.path(sn_dir, "seNorge2018_TX_1957_2021.nc"))
sn_crs  <- "EPSG:25833"   # file proj4 says WGS84/UTM33 (EPSG:32633); <1 m apart
sn_na16 <- -32768L

ndays <- function(y) 365L + (y %% 4 == 0 & (y %% 100 != 0 | y %% 400 == 0))

# dates from a CF time axis "<unit> since <origin>"
sn_dates <- function(nc) {
  u <- strsplit(ncatt_get(nc, "time", "units")$value, " since ")[[1]]
  s <- c(seconds = 1, hours = 3600, days = 86400)[[u[1]]]
  as.Date(as.POSIXct(u[2], tz = "UTC") + ncvar_get(nc, "time") * s)
}

sn_var <- function(nc) intersect(c("tn", "tx"), names(nc$var))

sn_template <- function(nc) {
  x <- ncvar_get(nc, "X"); y <- ncvar_get(nc, "Y")
  rast(nrows = length(y), ncols = length(x), crs = sn_crs,
       extent = ext(min(x) - 500, max(x) + 500, min(y) - 500, max(y) + 500))
}

# print header summary + consistency checks
sn_inspect <- function(files = sn_files) {
  ncs <- lapply(files, nc_open); on.exit(lapply(ncs, nc_close))
  for (k in names(ncs)) {
    nc <- ncs[[k]]; v <- sn_var(nc); d <- sn_dates(nc); a <- nc$var[[v]]
    cat(sprintf("\n[%s] %s\n  var=%s dims=%s units=%s fill=%s version=%s\n  time: %s | calendar=%s | %s .. %s (n=%d, contiguous=%s, Feb29=%d)\n  lat/lon vars: %s\n",
      k, basename(files[[k]]), v, paste(sapply(a$dim, `[[`, "name"), sapply(a$dim, `[[`, "len"), sep = "=", collapse = " "),
      a$units, a$missval, ncatt_get(nc, 0, "title")$value,
      ncatt_get(nc, "time", "units")$value, ncatt_get(nc, "time", "calendar")$value,
      min(d), max(d), length(d), all(diff(d) == 1), sum(format(d, "%m-%d") == "02-29"),
      paste(grep("^lat|^lon", names(nc$var), value = TRUE), collapse = ",")))
  }
  if (length(ncs) == 2) {
    same <- function(v) identical(ncvar_get(ncs[[1]], v), ncvar_get(ncs[[2]], v))
    d1 <- sn_dates(ncs[[1]]); d2 <- sn_dates(ncs[[2]]); n <- min(length(d1), length(d2))
    cat(sprintf("\nTN/TX: X identical=%s, Y identical=%s, shared dates identical=%s, common span %s .. %s\n",
                same("X"), same("Y"), identical(d1[1:n], d2[1:n]), max(min(d1), min(d2)), d1[n]))
  }
  invisible(NULL)
}

# grid file: lat, lon, land mask (+ grid.rds with land cell index and coords)
sn_grid <- function(out = sn_out, file = sn_files[["tn"]]) {
  f <- file.path(out, "grid.tif"); fr <- file.path(out, "grid.rds")
  if (file.exists(f) && file.exists(fr)) return(readRDS(fr))
  dir.create(out, FALSE, TRUE)
  nc <- nc_open(file); on.exit(nc_close(nc))
  lv <- grep("^lat", names(nc$var), value = TRUE)[1]; ov <- grep("^lon", names(nc$var), value = TRUE)[1]
  r <- sn_template(nc)
  lat <- as.vector(ncvar_get(nc, lv)); lon <- as.vector(ncvar_get(nc, ov))
  land <- which(!is.na(as.vector(ncvar_get(nc, sn_var(nc), c(1, 1, 1), c(-1, -1, 1)))))
  g <- rast(r, nlyrs = 3); names(g) <- c("lat", "lon", "land")
  values(g) <- cbind(lat, lon, replace(rep(NA_real_, ncell(r)), land, 1))
  writeRaster(g, f, overwrite = TRUE, datatype = "FLT4S", gdal = "COMPRESS=DEFLATE")
  gr <- list(land = land, lat = lat[land], lon = lon[land], nrow = nrow(r), ncol = ncol(r),
             ext = as.vector(ext(r)), crs = sn_crs)
  saveRDS(gr, fr); gr
}

# yearly THREDDS files (tn+tx+tg+rr, one version; see download_senorge.R) take precedence over the
# long single-variable files when present for that year
sn_archive_dir <- Sys.getenv("SENORGE_ARCHIVE_DIR", "1_raw_data/senorge_archive")
sn_source <- function(v, y, files = sn_files) {          # local archive first, then the OneDrive copy
  a <- file.path(c(sn_archive_dir, Sys.getenv("ARCHIVE_DEST", "")), sprintf("seNorge2018_%d.nc", y))
  a <- a[file.exists(a)]
  if (length(a)) a[1] else files[[v]]
}

# read one year of one variable -> matrix [land cells x days]
sn_read_year <- function(file, year, land, v = NULL) {
  nc <- nc_open(file); on.exit(nc_close(nc))
  if (is.null(v)) v <- sn_var(nc)[1]
  if (!v %in% names(nc$var)) stop(v, " not in ", basename(file))
  d <- sn_dates(nc); i <- which(as.integer(format(d, "%Y")) == year)
  if (!length(i)) stop(sprintf("%d not in %s (%s..%s)", year, basename(file), min(d), max(d)))
  stopifnot(length(i) == ndays(year), all(diff(i) == 1))
  m <- matrix(NA_real_, length(land), length(i))
  for (k in seq_along(i)) {
    s <- as.vector(ncvar_get(nc, v, c(1, 1, i[k]), c(-1, -1, 1)))
    if (sum(!is.na(s)) != sum(!is.na(s[land]))) warning(v, " ", d[i[k]], ": values outside land mask dropped")
    m[, k] <- s[land]
  }
  attr(m, "dates") <- d[i]; m
}

# ---- binary store (model input) ----------------------------------------------------------
sn_bin_file <- function(v, y, out = sn_out) file.path(out, "bin", sprintf("%s_%d.bin", v, y))

# m [land x days] (degC) -> int16 0.01 degC, written cell-major [days x land]
sn_write_bin <- function(m, f) {
  stopifnot(all(abs(m) < 300, na.rm = TRUE))
  x <- as.integer(round(t(m) * 100)); x[is.na(x)] <- sn_na16
  dir.create(dirname(f), FALSE, TRUE); tmp <- paste0(f, ".part")
  con <- file(tmp, "wb"); writeBin(x, con, size = 2L); close(con)
  file.rename(tmp, f)                    # atomic-ish on NTFS: no half-written .bin after a crash
}

# land cells i1..i1+n-1 of year y -> integer [days x n] (0.01 degC). Windows-safe binary seek.
sn_read_bin <- function(v, y, i1 = 1L, n = NULL, out = sn_out) {
  f <- sn_bin_file(v, y, out); nd <- ndays(y)
  if (is.null(n)) n <- file.size(f) / 2 / nd - (i1 - 1)
  con <- file(f, "rb"); on.exit(close(con))
  seek(con, (i1 - 1) * nd * 2, rw = "read")
  x <- readBin(con, "integer", nd * n, size = 2L, signed = TRUE)
  if (length(x) != nd * n) stop("short read: ", f)
  x[x == sn_na16] <- NA_integer_; dim(x) <- c(nd, n); x
}

# write [land x days] matrix as multi-layer float32 GeoTIFF, row block by row block
sn_write_tif <- function(m, gr, f) {
  r <- rast(nrows = gr$nrow, ncols = gr$ncol, extent = ext(gr$ext), crs = gr$crs, nlyrs = ncol(m))
  d <- attr(m, "dates"); names(r) <- format(d); time(r) <- d
  b <- writeStart(r, f, overwrite = TRUE, datatype = "FLT4S", names = format(d),
                  gdal = c("COMPRESS=DEFLATE", "PREDICTOR=3", "TILED=YES"))
  for (j in seq_len(b$n)) {
    c0 <- (b$row[j] - 1) * gr$ncol; nc <- b$nrows[j] * gr$ncol
    k <- which(gr$land > c0 & gr$land <= c0 + nc)
    v <- matrix(NA_real_, nc, ncol(m)); v[gr$land[k] - c0, ] <- m[k, ]
    writeValues(r, v, b$row[j], b$nrows[j])
  }
  writeStop(r); invisible(f)
}

# extract one year (both vars when available), QC tx >= tn; .bin always, .tif optional
sn_extract_year <- function(year, vars = c("tn", "tx"), files = sn_files, out = sn_out,
                            tif = TRUE, overwrite = FALSE) {
  gr <- sn_grid(out, files[["tn"]]); res <- list()
  for (v in vars) {
    fb <- sn_bin_file(v, year, out); ft <- file.path(out, v, sprintf("%s_%d.tif", v, year))
    need_b <- overwrite || !file.exists(fb); need_t <- tif && (overwrite || !file.exists(ft))
    if (!need_b && !need_t) next
    m <-tryCatch(sn_read_year(sn_source(v, year, files), year, gr$land, v), error = function(e) { warning(conditionMessage(e)); NULL })
    if (is.null(m)) next
    if (need_b) sn_write_bin(m, fb)
    if (need_t) { dir.create(dirname(ft), FALSE, TRUE); sn_write_tif(m, gr, ft) }
    res[[v]] <- m
  }
  if (all(c("tn", "tx") %in% names(res))) {
    bad <- res$tx < res$tn
    data.frame(year = year, tx_lt_tn = sum(bad, na.rm = TRUE), max_gap = max(c(0, (res$tn - res$tx)[which(bad)])))
  }
}

# extract many years in parallel (Windows: separate R processes; ~6 GB RAM per worker)
sn_extract <- function(years, workers = 4, ...) {
  sn_grid()                                              # create grid once, before workers start
  if (workers > 1) {
    old <- future::plan(future::multisession, workers = min(workers, length(years))); on.exit(future::plan(old))
    r <- future.apply::future_lapply(years, function(y) sn_extract_year(y, ...),
                                     future.packages = c("ncdf4", "terra"), future.seed = TRUE)
  } else r <- lapply(years, sn_extract_year, ...)
  do.call(rbind, r)
}

# read an extracted year back as SpatRaster (days as layers)
sn_load <- function(var, year, out = sn_out) rast(file.path(out, var, sprintf("%s_%d.tif", var, year)))

if (sys.nframe() == 0) {
  sn_inspect()
  a <- as.integer(commandArgs(TRUE))
  t0 <- Sys.time(); print(sn_extract(if (length(a)) a else 2010:2013))
  message("done in ", round(difftime(Sys.time(), t0, units = "mins"), 1), " min")
}
