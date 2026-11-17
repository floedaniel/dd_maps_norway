# Headless smoke test of the Shiny server. Run from project root:
#   Rscript -e "testthat::test_file('4_scripts/tests/test_app.R')"
library(testthat); library(shiny)
f_ips <- normalizePath(Filter(file.exists, file.path(c(".", "../.."), "5_outputs/opram/ips_typographus/ips_typographus_2010_emerge.tif"))[1])
app <- shinyAppDir(if (file.exists("4_scripts/ddmaps_app")) "4_scripts/ddmaps_app" else "../ddmaps_app")

test_that("maps and point model render", {
  testServer(app, {
    session$setInputs(set = "ips_typographus", year = "2010", month = "start_Jan", var = "emerge",
                      view = "year")
    expect_gt(nrow(mapdf()), 1e5)
    expect_match(output$map$src, "^data:image/png")                              # plot rendered
    session$setInputs(map_click = list(x = 262000, y = 6650000))                 # near Oslo
    expect_match(output$cellinfo, "Gen")
    session$setInputs(var = "ngen", view = "normal")
    expect_error(layer())                                                        # no normal yet -> validate()
    session$setInputs(preset = "Ips typographus", T0 = 8.3, thr = 557, dia = TRUE, dial = 14.5, diat_on = TRUE,
                      diat = 14.5, method = "sine", maxY = 4, pyear = "2010", pmonth = "1", lat = 59.91, lon = 10.75, go = 1)
    z <- prun()
    expect_equal(length(z$C), 1461)
    pt <- terra::project(terra::vect(cbind(10.75, 59.91), crs = "EPSG:4326"), "EPSG:25833")
    pre <- unlist(terra::extract(terra::rast(f_ips), pt)[, -1])
    expect_equal(as.numeric(z$r$emerge), as.numeric(pre))                         # live = precomputed
    expect_match(output$ptable, "generations"); expect_match(output$cumplot$src, "^data:image/png")
  })
})

test_that("user can add a species, compute it in the background, see it in Maps, delete it", {
  nm <- "Testus testus"; sl <- "testus_testus"
  rt <- normalizePath(dirname(dirname(dirname(f_ips))), "/")           # 5_outputs/opram -> project root
  rt <- normalizePath(file.path(rt, ".."), "/"); od <- file.path(rt, "5_outputs/opram")
  source(file.path(rt, "4_scripts/species_db.R"), local = TRUE)
  on.exit({ sp_delete(rt, nm); unlink(file.path(od, sl), recursive = TRUE)
            unlink(file.path(od, "_jobs", paste0(sl, ".log"))) }, add = TRUE)
  testServer(app, {
    n0 <- nrow(sp_df())
    # invalid: temperature without daylength, missing threshold -> rejected
    session$setInputs(sp_name = nm, sp_T0 = 9, sp_thr = NA, sp_dl = NA, sp_dt = 12, sp_add = 1)
    expect_equal(nrow(sp_df()), n0)
    # valid
    session$setInputs(sp_name = nm, sp_common = "Test beetle", sp_event = "Egg to adult", sp_T0 = 9, sp_thr = 400,
                      sp_dtype = "facultative", sp_dl = 14, sp_dt = NA, sp_ref1 = "https://doi.org/10.1000/xyz",
                      sp_ref2 = "Smith 2020, field data", sp_ref3 = "", sp_add = 2)
    expect_equal(nrow(sp_df()), n0 + 1)
    expect_true(nm %in% sp_read_all(rt)$species)
    expect_match(output$sptable, "Testus testus"); expect_match(output$sptable, "Smith 2020, field data")
    expect_true(nm %in% names(species()))
    # duplicate name rejected
    session$setInputs(sp_add = 3); expect_equal(nrow(sp_df()), n0 + 1)
    # compute 2010 in the background
    session$setInputs(cp_species = nm, cp_years = c(2010, 2010), cp_workers = 2, cp_go = 1)
    j <- job(); expect_true(j$proc$is_alive())
    t0 <- Sys.time(); while (j$proc$is_alive() && Sys.time() - t0 < 600) Sys.sleep(2)
    cat("\n  compute job:", round(as.numeric(difftime(Sys.time(), t0, units = "secs"))), "s\n")
    expect_equal(j$proc$get_exit_status(), 0)
    session$elapse(4000)                                    # let invalidateLater() refresh the status panel
    expect_match(output$cp_status$html, "done")
    expect_true(all(file.exists(file.path(od, sl, paste0(sl, "_2010_", c("emerge", "ngen"), ".tif")))))
    session$flushReact()
    expect_true(sl %in% list.dirs(od, FALSE, FALSE))
    session$setInputs(set = sl, year = "2010", month = "start_Jan", var = "ngen", view = "year")
    expect_gt(sum(!is.na(mapdf()$v)), 1e5)
    # delete
    session$setInputs(del_species = nm, del_go = 1)
    expect_false(nm %in% sp_df()$species)
  })
})
