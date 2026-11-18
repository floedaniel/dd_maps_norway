# Headless smoke test of the Shiny server. Run from project root:
#   Rscript -e "testthat::test_file('4_scripts/tests/test_app.R')"
library(testthat); library(shiny)
f_ips <- normalizePath(Filter(file.exists, file.path(c(".", "../.."), "5_outputs/opram/ips_typographus/ips_typographus_2010_emerge.tif"))[1])
app <- shinyAppDir(if (file.exists("4_scripts/ddmaps_app")) "4_scripts/ddmaps_app" else "../ddmaps_app")

test_that("maps, reference periods, trend and point model render", {
  testServer(app, {
    session$setInputs(set = "ips_typographus", year = "2010", month = "start_Jan", var = "emerge",
                      view = "year", period = "1991_2020")
    expect_gt(nrow(mapdf()), 1e5)
    expect_match(output$map$src, "^data:image/png")                              # plot rendered
    session$setInputs(map_click = list(x = 262000, y = 6650000))                 # near Oslo
    expect_match(output$cellinfo, "Gen")
    session$setInputs(var = "ngen", view = "normal")
    expect_gt(sum(!is.na(mapdf()$v)), 1e5)                                       # 1991-2020 median available
    expect_match(output$maptitle, "median 1991.2020")
    session$setInputs(view = "p1")                                               # probability maps (0-100 %)
    v <- mapdf()$v; expect_true(all(v >= 0 & v <= 100, na.rm = TRUE)); expect_gt(max(v, na.rm = TRUE), 50)
    expect_match(output$maptitle, "probability of ≥1 generation")
    expect_match(output$cellinfo, ">Years<"); expect_match(output$cellinfo, "title=")
    session$setInputs(period = "2008_2022")                                      # recent reference period
    expect_match(output$maptitle, "2008"); expect_match(output$maptitle, "2022"); expect_gt(sum(!is.na(mapdf()$v)), 1e5)
    expect_match(output$cellnote$html, "2008"); expect_match(output$cellnote$html, "median")
    session$setInputs(view = "p2"); expect_match(output$map$src, "^data:image/png")
    session$setInputs(set = "test_t10_dd300", view = "p1")
    expect_error(layer(), "No probability maps")                                 # set without a full period
    session$setInputs(set = "ips_typographus", view = "year", period = "1991_2020")
    z <- trend_data(); expect_gte(nrow(z), 30); expect_true(all(z$ngen >= 0))    # trend 1991-2022 at the cell
    expect_match(output$trend$src, "^data:image/png"); expect_match(output$trendhead, "^Trend 1991")
    session$setInputs(preset = "Ips typographus", T0 = 8.3, thr = 557, dia = TRUE, dial = 14.5, diat_on = TRUE,
                      diat = 14.5, method = "sine", pupper = NA, pemonly = FALSE, maxY = 4, pyear = "2010", pmonth = "1",
                      lat = 59.91, lon = 10.75, go = 1)
    z <- prun()
    expect_equal(length(z$C), 1461)
    pt <- terra::project(terra::vect(cbind(10.75, 59.91), crs = "EPSG:4326"), "EPSG:25833")
    pre <- unlist(terra::extract(terra::rast(f_ips), pt)[, -1])
    expect_equal(as.numeric(z$r$emerge), as.numeric(pre))                         # live = precomputed
    expect_match(output$ptable, "Generations"); expect_match(output$cumplot$src, "^data:image/png")
    session$setInputs(pupper = 15, go = 2)                                        # upper threshold slows development
    expect_true(all(prun()$r$emerge >= pre, na.rm = TRUE))
  })
})

test_that("emergence-only species hide generation counts", {
  skip_if_not(file.exists(gsub("ips_typographus", "popillia_japonica", f_ips)), "Popillia maps not computed")
  testServer(app, {
    session$setInputs(set = "popillia_japonica", year = "2010", month = "start_Jan", var = "emerge",
                      view = "year", period = "1991_2020", map_click = list(x = 262000, y = 6650000))
    expect_true(emonly()); expect_gt(sum(!is.na(mapdf()$v)), 1e4)
    expect_false(grepl(">Gen.<", output$cellinfo, fixed = TRUE))              # no generation column
    session$setInputs(var = "ngen"); expect_error(layer(), "emergence only")
    session$setInputs(var = "emerge", view = "p2"); expect_error(layer(), "emergence only")
  })
})

test_that("user can add a species, compute it in the background, see it in Maps, delete it", {
  nm <- "Testus testus"; sl <- "testus_testus"
  rt <- normalizePath(dirname(dirname(dirname(f_ips))), "/")           # 5_outputs/opram -> project root
  rt <- normalizePath(file.path(rt, ".."), "/"); od <- file.path(rt, "5_outputs/opram")
  source(file.path(rt, "4_scripts/species_db.R"), local = TRUE)
  on.exit({ sp_delete(rt, nm); unlink(file.path(od, sl), recursive = TRUE)
            unlink(file.path(od, "_jobs", paste0(sl, c(".log", ".lock")))) }, add = TRUE)
  testServer(app, {
    n0 <- nrow(sp_df())
    # invalid: temperature without daylength, missing degree-days, upper below T0 -> rejected
    session$setInputs(sp_name = nm, sp_T0 = 9, sp_thr = NA, sp_dl = NA, sp_dt = 12, sp_upper = 5, sp_add = 1)
    expect_equal(nrow(sp_df()), n0)
    # valid, with method, upper threshold and emergence-only flag
    session$setInputs(sp_name = nm, sp_common = "Test beetle", sp_event = "Egg to adult", sp_T0 = 9, sp_thr = 400,
                      sp_dtype = "facultative", sp_dl = 14, sp_dt = NA, sp_ref1 = "https://doi.org/10.1000/xyz",
                      sp_ref2 = "Smith 2020, field data", sp_ref3 = "", sp_method = "triangle", sp_upper = 30,
                      sp_emonly = TRUE, sp_add = 2)
    expect_equal(nrow(sp_df()), n0 + 1)
    row <- sp_read_all(rt); row <- row[row$species == nm, ]
    expect_equal(row$method, "triangle"); expect_equal(row$upper_threshold, 30); expect_true(row$emergence_only)
    expect_match(output$sptable, "Testus testus"); expect_match(output$sptable, "Smith 2020, field data")
    expect_match(output$sptable, "emergence only")
    expect_equal(species()[[nm]]$upper, 30)
    # duplicate name rejected
    session$setInputs(sp_add = 3); expect_equal(nrow(sp_df()), n0 + 1)
    # compute 2010 in the background
    session$setInputs(cp_species = nm, cp_years = c(2010, 2010), cp_workers = 2, cp_go = 1)
    j <- job(); expect_true(j$proc$is_alive())
    t0 <- Sys.time(); while (j$proc$is_alive() && Sys.time() - t0 < 600) Sys.sleep(2)
    cat("\n  compute job:", round(as.numeric(difftime(Sys.time(), t0, units = "secs"))), "s\n")
    expect_equal(j$proc$get_exit_status(), 0)
    expect_false(file.exists(file.path(od, "_jobs", paste0(sl, ".lock"))))       # lock removed after the job
    session$elapse(4000)                                    # let invalidateLater() refresh the status panel
    expect_match(output$cp_status$html, "done")
    expect_true(all(file.exists(file.path(od, sl, paste0(sl, "_2010_", c("emerge", "ngen"), ".tif")))))
    expect_equal(readRDS(file.path(od, sl, "params.rds"))$upper, 30)
    session$flushReact()
    expect_true(sl %in% list.dirs(od, FALSE, FALSE))
    session$setInputs(set = sl, year = "2010", month = "start_Jan", var = "emerge", view = "year", period = "1991_2020")
    expect_gt(sum(!is.na(mapdf()$v)), 1e4)
    # delete row and maps
    session$setInputs(del_species = nm, del_maps = TRUE, del_go = 1)
    expect_false(nm %in% sp_df()$species); expect_false(dir.exists(file.path(od, sl)))
  })
})
