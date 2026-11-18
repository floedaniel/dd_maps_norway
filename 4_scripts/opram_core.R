# opram_core.R -- R port of the OPRAM degree-day pest phenology model.
#
# Original model and code: Jon Yearsley, OPRAM (https://github.com/DrJonYearsley/OPRAM),
#   julia/OPRAM_main_functions.jl, julia/OPRAM_processing_functions.jl,
#   julia/OPRAM_calculate_average.jl. Copyright (C) 2026 Jon Yearsley, GPL-3.
# This port: GPL-3 (or later). Vectorised re-implementation; equivalence with the Julia
# algorithm is checked in 4_scripts/tests/test_opram_core.R.
#
# Conventions (as OPRAM):
#  * temperature matrices are [days x cells]; days run from 1 Jan of the simulation year
#    over maxYears consecutive years (day index is continuous, "emergeDOY" uses it too)
#  * a generation is "within the year" if emergence day <= 365 (leap day ignored)
#  * development starts on the 1st of each month (leap-aware DOY) unless start_doy is given

# ---- parameters -------------------------------------------------------------------------
nullna <- function(x) if (is.null(x) || !length(x) || is.na(x[1])) NULL else x

# upper: optional upper development threshold (horizontal cutoff: no extra development above it).
# emergence_only: the degree-day sum covers only part of the life cycle (e.g. overwintered larvae
#   to adult); generation counts are then not meaningful (used by the app, not by the model).
opram_params <- function(T0, threshold, diapause_daylength = NULL, diapause_temperature = NULL,
                         method = c("sine", "triangle", "average"), maxYears = 4L, name = NA_character_,
                         upper = NULL, emergence_only = FALSE) {
  method <- if (is.null(nullna(method))) "sine" else match.arg(method, c("sine", "triangle", "average"))
  stopifnot(is.numeric(T0), is.numeric(threshold), threshold > 0, maxYears >= 1)
  dl <- nullna(diapause_daylength); dt <- nullna(diapause_temperature); up <- nullna(upper)
  if (!is.null(dt) && is.null(dl)) stop("diapause_temperature requires diapause_daylength")
  if (!is.null(up) && up <= T0) stop("upper threshold must be above T0")
  list(name = name, T0 = T0, threshold = threshold, diapause_daylength = dl,
       diapause_temperature = dt, method = method, maxYears = as.integer(maxYears),
       upper = up, emergence_only = isTRUE(as.logical(emergence_only)))
}

# species presets from OPRAM data/species_parameters.csv -> list of opram_params
opram_species <- function(file, method = "sine", maxYears = 4L) {
  s <- read.csv(file, check.names = FALSE)
  setNames(lapply(seq_len(nrow(s)), function(i) with(s[i, ], opram_params(
    baseline_temperature, threshold, diapause_daylength, diapause_temperature,
    method, maxYears, species))), s$species)        # OPRAM presets: no upper threshold, full life cycle
}

# ---- degree-day methods (OPRAM: average_/triangle_/sine_degreeday) ------------------------
# Valid only on "update" days (sine/triangle: Tmax > T0; average: Tmin + Tmax > 2*T0).
dd_average <- function(tn, tx, T0) pmax((tn + tx) / 2 - T0, 0)

dd_triangle <- function(tn, tx, T0) {
  dd <- (tn + tx) / 2 - T0; i <- which(tn < T0)
  dd[i] <- (tx[i] - T0)^2 / (2 * (tx[i] - tn[i])); dd
}

dd_sine <- function(tn, tx, T0) {
  dd <- (tn + tx) / 2 - T0; i <- which(tn < T0)
  dT <- (tx[i] - tn[i]) / 2; d0 <- dd[i]
  dd[i] <- (sqrt(dT^2 - d0^2) + acos(-d0 / dT) * d0) / pi; dd
}

# ---- photoperiod / diapause (OPRAM: photoperiod) ------------------------------------------
# Forsythe et al. (1995) daylength [h], as geosphere::daylength. Returns [length(lat) x length(doy)].
daylength <- function(lat, doy) {
  P <- asin(0.39795 * cos(0.2163108 + 2 * atan(0.9671396 * tan(0.0086 * (doy - 186)))))
  r <- pi / 180
  a <- (sin(0.8333 * r) + outer(sin(lat * r), sin(P))) / outer(cos(lat * r), cos(P))
  24 - 24 / pi * acos(pmin(pmax(a, -1), 1))
}

# first DOY >= min_doy with daylength < threshold; 367 = never (no diapause)
diapause_doy <- function(lat, threshold, min_doy = 150L) {
  d <- min_doy:366; below <- daylength(lat, d) < threshold
  out <- d[max.col(below, "first")]; out[rowSums(below) == 0] <- 367L; out
}

# ---- daily degree days with diapause (OPRAM: calculate_GDD) ------------------------------
# tn, tx [days x cells]; doy day-of-year per row; lat per cell. Returns DD (0 on non-update days).
degree_days <- function(tn, tx, p, doy, lat = NULL, swap = TRUE) {
  if (swap) { lo <- pmin(tn, tx); tx <- pmax(tn, tx); tn <- lo }   # enforce Tmin <= Tmax (seNorge TN/TX versions differ)
  upd <- if (p$method == "average") tn + tx > 2 * p$T0 else tx > p$T0
  if (!is.null(p$diapause_daylength)) {
    off <- outer(doy, diapause_doy(lat, p$diapause_daylength), ">=")
    if (!is.null(p$diapause_temperature)) off <- off & tn + tx <= 2 * p$diapause_temperature
    upd <- upd & !off
  }
  upd[is.na(upd)] <- FALSE
  f <- switch(p$method, sine = dd_sine, triangle = dd_triangle, average = dd_average)
  dd <- array(0, dim(tn)); i <- which(upd); dd[i] <- f(tn[i], tx[i], p$T0)
  if (!is.null(p$upper)) {   # horizontal cutoff: subtract the degree-days above the upper threshold
    above <- upd & (if (p$method == "average") tn + tx > 2 * p$upper else tx > p$upper)
    j <- which(above); dd[j] <- dd[j] - f(tn[j], tx[j], p$upper)
  }
  dd
}

# ---- emergence (OPRAM: emergence / location_loop) -----------------------------------------
# Column-wise cumulative sums stacked into one non-decreasing vector by adding a per-column
# offset, so any number of "first day C[t, j] >= target" queries is a single findInterval().
cum_index <- function(dd) {
  C <- if (requireNamespace("matrixStats", quietly = TRUE)) matrixStats::colCumsums(dd)   # ~6x faster
       else apply(dd, 2, cumsum)
  if (is.null(dim(C))) C <- matrix(C, 1)
  n <- nrow(C); big <- max(C[n, ], 0) + 1
  list(C = C, n = n, big = big, g = as.vector(C) + rep((seq_len(ncol(C)) - 1) * big, each = n))
}

# first day t with C[t, j] >= C[start - 1, j] + threshold, for cells j (NA = no emergence)
emerge_day <- function(ci, j, start, threshold) {
  base <- ifelse(start > 1, ci$C[cbind(pmax(start - 1, 1), j)], 0)
  tgt <- base + threshold
  o <- (j - 1) * ci$big; tgt <- pmin(tgt, ci$big - 0.5)   # targets past column total -> no emergence
  # g is sorted and NA-free by construction: skip findInterval's full-vector checks (~2x faster)
  t <- findInterval(tgt + o, ci$g, left.open = TRUE, checkSorted = FALSE, checkNA = FALSE) - (j - 1) * ci$n + 1L
  t[t > ci$n | base + threshold > ci$C[ci$n, j]] <- NA_integer_; as.integer(t)
}

# ---- generations from one start day (OPRAM: extract_results) -----------------------------
generations <- function(ci, start, threshold) {
  k <- ncol(ci$C); s <- rep(as.integer(start), k)
  ngen <- numeric(k); em1 <- rep(NA_integer_, k); first <- TRUE; act <- seq_len(k)
  while (length(act)) {
    e <- emerge_day(ci, act, s[act], threshold)
    ok <- !is.na(e); full <- ok & e <= 365; part <- ok & !full
    ngen[act[full]] <- ngen[act[full]] + 1
    ngen[act[part]] <- ngen[act[part]] + (365 - s[act[part]]) / (e[part] - s[act[part]])
    if (first) { em1[act] <- e; first <- FALSE }
    ns <- integer(length(act)); ns[full] <- e[full] + 1L; ns[ns >= 365] <- 0L
    s[act] <- ns; act <- act[ns > 0]
  }
  ngen[is.na(em1)] <- 0
  list(emerge = em1, ngen = ngen)
}

# ---- one simulation year for a block of cells --------------------------------------------
# tn, tx [days x cells] for p$maxYears years from 1 Jan; dates per row; lat per cell.
# Returns list(emerge, ngen): [cells x starts] (emerge = continuous day index from 1 Jan).
# swap = FALSE if the caller already enforced Tmin <= Tmax (saves time when looping parameter sets).
opram_run <- function(tn, tx, dates, lat, p, start_doy = NULL, swap = TRUE) {
  stopifnot(nrow(tn) == length(dates), identical(dim(tn), dim(tx)))
  y1 <- as.integer(format(dates[1], "%Y")); doy <- as.integer(format(dates, "%j"))
  if (is.null(start_doy)) start_doy <- month_starts(y1)
  ci <- cum_index(degree_days(tn, tx, p, doy, lat, swap))
  r <- lapply(start_doy, function(s) generations(ci, s, p$threshold))
  list(emerge = sapply(r, `[[`, "emerge"), ngen = sapply(r, `[[`, "ngen"), start_doy = start_doy)
}

month_starts <- function(year) as.integer(format(as.Date(sprintf("%d-%02d-01", year, 1:12)), "%j"))

# years of meteo used for simulation year y (OPRAM read_meteo: capped at last available year)
opram_years <- function(y, maxYears, last_year) pmin(y:(y + maxYears - 1L), last_year)

# multi-year normal (OPRAM_calculate_average): median across years; emerge NA if beyond 366*maxYears
opram_normal <- function(emerge, ngen, maxYears) {   # [cells x years] each
  med <- function(m) apply(m, 1, median, na.rm = TRUE)
  em <- med(emerge); ng <- med(ngen)
  bad <- is.na(em) | em > 366 * maxYears; em[bad] <- NA; ng[bad] <- 0
  list(emerge = em, ngen = pmax(ng, 0))
}
