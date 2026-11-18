# Tests for opram_core.R. Run from project root:
#   Rscript -e "testthat::test_file('4_scripts/tests/test_opram_core.R')"
library(testthat)
src <- c("4_scripts/opram_core.R", "../opram_core.R", "../../4_scripts/opram_core.R")
source(src[file.exists(src)][1])

# ---- literal scalar port of the Julia pipeline for one location (reference) --------------
# emergence() -> cleanup_results(save1year = TRUE) -> extract_results()
jl_emergence <- function(cum, cumdoy, thr) {
  e <- integer(length(cum)); i <- 1; em <- which(cum >= thr)[1]
  while (!is.na(em)) { e[i] <- cumdoy[em]; i <- i + 1; em <- which(cum >= cum[i - 1] + thr)[1] }
  e
}
jl_location <- function(tn, tx, doy, lat, p, start, save1year = FALSE) {
  lo <- pmin(tn, tx); tx <- pmax(tn, tx); tn <- lo
  upd <- if (p$method == "average") tn + tx > 2 * p$T0 else tx > p$T0
  if (!is.null(p$diapause_daylength)) {
    dd0 <- diapause_doy(lat, p$diapause_daylength)
    ow <- doy >= dd0
    if (!is.null(p$diapause_temperature)) ow <- ow & (tn + tx <= 2 * p$diapause_temperature)
    upd <- upd & !ow
  }
  f <- switch(p$method, sine = dd_sine, triangle = dd_triangle, average = dd_average)
  rDOY <- which(upd); cum <- cumsum(f(tn[upd], tx[upd], p$T0))
  rEm <- jl_emergence(cum, rDOY, p$threshold)
  n <- length(rDOY); if (!n) return(c(emerge = NA, ngen = 0))
  keep <- rEm > 0 & c(rEm[-n] != rEm[-1], TRUE)
  if (save1year) keep <- keep & c(TRUE, rDOY[-1] <= 365 | (rDOY[-n] < 365 & rDOY[-1] >= 365))
  D <- rDOY[keep]; E <- rEm[keep]
  s <- start; ngen <- 0; em1 <- NA; first <- TRUE
  while (s > 0) {
    i3 <- which(D >= s)[1]; done <- !is.na(i3)
    within <- done && E[i3] <= 365
    if (done && within) ngen <- ngen + 1
    if (done && !within) ngen <- ngen + (365 - s) / (E[i3] - s)
    if (first) { if (done) em1 <- E[i3]; first <- FALSE }
    s <- if (done && within) E[i3] + 1 else 0
    if (s >= 365) s <- 0
  }
  if (is.na(em1)) ngen <- 0
  c(emerge = em1, ngen = ngen)
}

# synthetic seasonal climate [days x cells]
synth <- function(y1, years, ncell, warm = 0, seed = 1) {
  set.seed(seed)
  dates <- seq(as.Date(sprintf("%d-01-01", y1)), as.Date(sprintf("%d-12-31", y1 + years - 1)), "day")
  doy <- as.integer(format(dates, "%j")); nd <- length(dates)
  mu <- outer(-cos(2 * pi * (doy - 15) / 365.25) * 10, runif(ncell, 0, 6) + warm, `+`) + 4
  tn <- mu + matrix(rnorm(nd * ncell, 0, 3), nd)
  tx <- tn + matrix(runif(nd * ncell, -0.3, 12), nd)        # a few Tmax < Tmin on purpose
  list(tn = tn, tx = tx, dates = dates, lat = runif(ncell, 58, 71.2))
}

# ---- tests --------------------------------------------------------------------------------
test_that("degree-day methods match known values", {
  expect_equal(dd_sine(5, 15, 10), 5 / pi)                   # Tmin = T0 - a, Tmax = T0 + a: a/pi
  expect_equal(dd_sine(12, 20, 10), 6)                       # Tmin >= T0: mean - T0
  expect_equal(dd_triangle(5, 15, 10), 25 / 20)
  expect_equal(dd_triangle(12, 20, 10), 6)
  expect_equal(dd_average(c(5, 2), c(15, 10), 8), c(2, 0))
  # sine method equals the integral of max(T(t) - T0, 0) over a sinusoidal day
  for (v in list(c(2, 14, 10), c(-3, 11, 10), c(8, 30, 10), c(9.9, 10.1, 10))) {
    m <- mean(v[1:2]); a <- diff(v[1:2]) / 2
    ref <- integrate(function(t) pmax(m + a * sin(2 * pi * t) - v[3], 0), 0, 1,
                     subdivisions = 5000L, rel.tol = 1e-12)$value
    expect_equal(dd_sine(v[1], v[2], v[3]), ref, tolerance = 1e-6)
  }
})

test_that("upper threshold = horizontal cutoff, checked by numerical integration", {
  t <- (seq_len(20000) - 0.5) / 20000
  sine_T <- function(tn, tx) (tn + tx) / 2 + (tx - tn) / 2 * sin(2 * pi * t)
  tri_T  <- function(tn, tx) tn + (tx - tn) * (1 - abs(2 * t - 1))
  for (v in list(c(5, 30), c(14, 26), c(18, 24), c(2, 12), c(-3, 9))) {
    for (m in c("sine", "triangle", "average")) {
      p <- opram_params(10, 100, method = m, upper = 21.7)
      got <- degree_days(matrix(v[1]), matrix(v[2]), p, doy = 1)[1, 1]
      ref <- switch(m,
        sine     = mean(pmin(pmax(sine_T(v[1], v[2]) - 10, 0), 11.7)),
        triangle = mean(pmin(pmax(tri_T(v[1], v[2]) - 10, 0), 11.7)),
        average  = min(max(mean(v) - 10, 0), 11.7))
      expect_equal(got, ref, tolerance = 1e-4, info = paste(m, v[1], v[2]))
    }
  }
  # no upper threshold -> identical to before; upper never reached -> no change
  p0 <- opram_params(10, 100); p1 <- opram_params(10, 100, upper = 40)
  tn <- matrix(c(5, 8, 12), 3); tx <- matrix(c(20, 25, 30), 3)
  expect_equal(degree_days(tn, tx, p0, 1:3), degree_days(tn, tx, p1, 1:3))
  expect_error(opram_params(10, 100, upper = 9), "above T0")
  expect_true(opram_params(10, 100, emergence_only = TRUE)$emergence_only)
  expect_equal(opram_params(10, 100, method = NA)$method, "sine")
})

test_that("daylength matches geosphere::daylength, incl. polar day/night", {
  skip_if_not_installed("geosphere")
  # DOY 366 excluded: geosphere wraps it to DOY 1, OPRAM (Julia) uses the raw formula
  lat <- c(0, 45, 58, 63.4, 66.6, 69.6, 71.2); doy <- c(1, 80, 150, 172, 200, 266, 355, 365)
  ref <- t(sapply(lat, function(l) geosphere::daylength(l, doy)))
  expect_equal(daylength(lat, doy), ref, tolerance = 1e-10)
  expect_equal(daylength(71, 172)[1], 24); expect_equal(daylength(71, 355)[1], 0)
  expect_equal(diapause_doy(c(60, 71), 0), c(367L, 367L))    # never below 0 h -> no diapause
  expect_true(all(diapause_doy(c(60, 71), 14.5) %in% 150:366))
})

test_that("emerge_day equals brute-force search on random data", {
  set.seed(42); nd <- 400; k <- 60
  dd <- matrix(rexp(nd * k) * rbinom(nd * k, 1, 0.6), nd)
  ci <- cum_index(dd); C <- ci$C
  for (thr in c(5, 50, 200)) {
    st <- sample(1:nd, k, TRUE)
    bf <- sapply(1:k, function(j) { b <- if (st[j] > 1) C[st[j] - 1, j] else 0; which(C[, j] >= b + thr)[1] })
    expect_identical(emerge_day(ci, 1:k, st, thr), as.integer(bf))
  }
})

test_that("vectorised model equals literal Julia-pipeline port", {
  ps <- list(opram_params(10, 300), opram_params(8.3, 557, 14.5, 14.5),
             opram_params(10, 300, 12, method = "triangle"), opram_params(6.6, 81, method = "average"),
             opram_params(12.9, 625, 15, 20, maxYears = 2))
  for (cfg in list(list(y = 2011, w = 2), list(y = 2012, w = 8))) {        # 2012 = leap
    for (p in ps) {
      d <- synth(cfg$y, p$maxYears, 25, cfg$w, seed = cfg$y)
      r <- opram_run(d$tn, d$tx, d$dates, d$lat, p)
      doy <- as.integer(format(d$dates, "%j"))
      for (m in seq_along(r$start_doy)) {
        ref <- t(sapply(1:25, function(j) jl_location(d$tn[, j], d$tx[, j], doy, d$lat[j], p, r$start_doy[m])))
        expect_identical(r$emerge[, m], as.integer(ref[, "emerge"]), info = paste(p$method, cfg$y, m))
        expect_equal(r$ngen[, m], ref[, "ngen"], tolerance = 1e-12, info = paste(p$method, cfg$y, m))
      }
    }
  }
})

test_that("Julia save1year=TRUE only differs where emergence run crosses day 365 (documented quirk)", {
  p <- opram_params(10, 300); d <- synth(2011, 4, 25, 2, seed = 2011)
  r <- opram_run(d$tn, d$tx, d$dates, d$lat, p); doy <- as.integer(format(d$dates, "%j"))
  ref <- sapply(1:25, function(j) jl_location(d$tn[, j], d$tx[, j], doy, d$lat[j], p, r$start_doy[11], TRUE))
  diff <- which(!(r$emerge[, 11] %in% ref["emerge", ]) | is.na(ref["emerge", ]) != is.na(r$emerge[, 11]))
  cat("\n save1year=TRUE: Nov-1 start, cells differing:", length(diff), "of 25\n")
  expect_true(all(is.na(ref["emerge", diff])))                 # quirk only ever drops results
})

test_that("cold cells give no emergence and zero generations", {
  d <- synth(2013, 4, 5, warm = -20)
  r <- opram_run(d$tn, d$tx, d$dates, d$lat, opram_params(10, 300))
  expect_true(all(is.na(r$emerge))); expect_true(all(r$ngen == 0))
})

test_that("helpers", {
  expect_equal(month_starts(2012)[3], 61L); expect_equal(month_starts(2013)[3], 60L)
  expect_equal(opram_years(2020, 4, 2021), c(2020, 2021, 2021, 2021))
  expect_error(opram_params(10, 300, diapause_temperature = 14))
  n <- opram_normal(cbind(c(100, NA, 2000), c(120, NA, 2000), c(110, 50, 2000)), cbind(c(2, 0, 1), c(1, 0, 1), c(3, 1, 1)), 4)
  expect_equal(n$emerge, c(110, 50, NA)); expect_equal(n$ngen, c(2, 0, 0))
})
