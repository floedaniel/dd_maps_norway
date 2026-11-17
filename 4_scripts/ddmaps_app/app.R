# DD-maps Norway -- Shiny app: degree-day pest phenology maps for Norway (1 km, seNorge2018).
# Model: R port of OPRAM (Jon Yearsley, github.com/DrJonYearsley/OPRAM), GPL-3; app idea after pestmaps.ie.
# Climate: MET Norway seNorge_2018 (CC BY 4.0). This code: GPL-3.
# Run from project root:  shiny::runApp("4_scripts/ddmaps_app")

suppressPackageStartupMessages({ library(shiny); library(bslib); library(terra); library(ggplot2) })

root <- normalizePath(if (file.exists("4_scripts/opram_core.R")) "." else "../..", "/")
source(file.path(root, "4_scripts/opram_core.R"))
source(file.path(root, "4_scripts/extract_senorge.R"))           # sn_read_bin(), ndays()
proc <- file.path(root, "2_processed_data/senorge"); outd <- file.path(root, "5_outputs/opram")
gr <- readRDS(file.path(proc, "grid.rds"))
tmpl <- rast(nrows = gr$nrow, ncols = gr$ncol, extent = ext(gr$ext), crs = gr$crs)
nor <- readRDS(file.path(proc, "norway_mask.rds"))                # from 4_scripts/prep_norway_mask.R
ncell_nor <- gr$land[nor$keep]                                    # map cells: Norway only
xy_nor <- xyFromCell(tmpl, ncell_nor)
# geom_raster needs every column/row present with equal spacing: pad empty ones with transparent cells
pad_nor <- local({
  gx <- setdiff(seq(min(xy_nor[, 1]), max(xy_nor[, 1]), 1000), xy_nor[, 1])
  gy <- setdiff(seq(min(xy_nor[, 2]), max(xy_nor[, 2]), 1000), xy_nor[, 2])
  data.frame(x = c(gx, rep(min(xy_nor[, 1]), length(gy))), y = c(rep(min(xy_nor[, 2]), length(gx)), gy), v = NA_real_)
})
# grey land underlay: shows through where a value is NA (= no emergence within maxYears)
base_nor <- rbind(data.frame(x = xy_nor[, 1], y = xy_nor[, 2], a = 1), data.frame(pad_nor[, 1:2], a = 0))
landmap <- match(seq_len(ncell(tmpl)), gr$land)                   # cell -> land index
landmap[gr$land[!nor$keep]] <- NA                                 # clicks outside Norway ignored
source(file.path(root, "4_scripts/species_db.R"))                # OPRAM presets + user species
n_blocks <- ceiling(length(gr$land) / 5000)                       # progress denominator for compute jobs

# species overview table (inputs + references) as HTML-ready data frame
species_tbl <- function(s) {
  fix <- function(u) {                                   # repair malformed links in the upstream CSV
    u <- trimws(u); u <- sub("^https://doi\\..*?org/", "https://doi.org/", u)   # pasted local path inside URL
    sub("doi\\.org/doi:", "doi.org/", u)
  }
  esc <- function(x) htmltools::htmlEscape(ifelse(is.na(x), "", x))
  link <- function(u) ifelse(is.na(u) | u == "none", "",
    ifelse(grepl("^https?://", trimws(u)),
      sprintf("<a href='%s' target='_blank'>%s</a>", esc(fix(u)), esc(sub("^https?://(dx\\.)?(doi\\.org/)?", "", fix(u)))),
      esc(u)))                                           # free-text reference (user input): shown as text
  refs <- apply(matrix(sapply(c("Reference1", "Reference2", "Reference3"), function(k) link(s[[k]])), nrow(s)), 1,
                function(r) paste(r[nzchar(r)], collapse = "<br>"))
  num <- function(x, d = 1) ifelse(is.na(x), "–", formatC(x, format = "f", digits = d, drop0trailing = TRUE))
  data.frame(Species = sprintf("<i>%s</i>", esc(s$species)), `Common name` = esc(s$commonName),
             Event = esc(s$event), `T0 (°C)` = num(s$baseline_temperature, 2),
             `Threshold (DD)` = num(s$threshold, 2), Diapause = ifelse(is.na(s$diapause), "–", esc(s$diapause)),
             `Daylength (h)` = num(s$diapause_daylength), `Temp. (°C)` = num(s$diapause_temperature),
             References = refs, Source = ifelse(s$source == "user", "<span class='badge bg-info'>user</span>", "OPRAM"),
             check.names = FALSE)
}
last_meteo <- 2021L
slug <- function(x) gsub("[^a-z0-9]+", "_", tolower(trimws(x)))
by <- function(v) sub("^(tn|tx)_(\\d+)\\.bin$", "\\2", list.files(file.path(proc, "bin"), sprintf("^%s_\\d+\\.bin$", v)))
met_years <- sort(as.integer(intersect(by("tn"), by("tx"))))
sim_years <- function(mY) met_years[vapply(met_years, function(y) all(opram_years(y, mY, last_meteo) %in% met_years), TRUE)]

result_sets <- function() {
  d <- list.dirs(outd, FALSE, FALSE)
  d[vapply(d, function(x) length(list.files(file.path(outd, x), "_emerge\\.tif$")) > 0, TRUE)]
}
set_label <- function(d, nms) { n <- nms[slug(nms) == d]; if (length(n)) n[1] else d }
set_years <- function(d) sort(unique(sub(".*_(\\d{4})_emerge\\.tif$", "\\1",
                                         list.files(file.path(outd, d), "_\\d{4}_emerge\\.tif$"))))
set_normal <- function(d) sub("_emerge\\.tif$", "", list.files(file.path(outd, d), "_normal_\\d+_\\d+_emerge\\.tif$"))[1]

day_lab <- function(y, fmt = "%d %b %Y") function(v) ifelse(is.na(v), "", format(as.Date(sprintf("%d-01-01", y)) + v - 1, fmt))
short_lab <- function(y) function(v) {                             # "07 Jul", "12 Jun +1y"
  d <- as.Date(sprintf("%d-01-01", y)) + v - 1; k <- as.integer(format(d, "%Y")) - y
  ifelse(is.na(v), "–", paste0(format(d, "%d %b"), ifelse(k > 0, sprintf(" +%dy", k), "")))
}
xy_of <- function(lon, lat) crds(project(vect(cbind(lon, lat), crs = "EPSG:4326"), gr$crs))

# ---------------------------------------------------------------------------------------- UI
ui <- page_navbar(
  title = tags$span("DD-maps Norway", tags$small(class = "ms-2 opacity-75", "degree-day pest phenology")),
  window_title = "DD-maps Norway", theme = bs_theme(version = 5, bootswatch = "flatly"),
  nav_panel("Maps", layout_sidebar(
    sidebar = sidebar(width = 300,
      selectInput("set", "Species / parameter set", NULL),
      selectInput("year", "Year", NULL),
      selectInput("month", "Development starts", setNames(paste0("start_", month.abb), paste("1", month.name))),
      radioButtons("var", "Show", c("Adult emergence date" = "emerge", "Generations per year" = "ngen")),
      radioButtons("view", "Values", c("Selected year" = "year", "30-yr median" = "normal", "Anomaly (year - median)" = "anom")),
      helpText("Click the map to inspect a 1 km cell; the cell is also used for the point model."),
      tags$small("Degree-day model ported from OPRAM (J. Yearsley, GPL-3). Climate: seNorge_2018, MET Norway. Borders: GADM 4.1.")),
    layout_columns(col_widths = c(9, 3), fill = TRUE,
      card(full_screen = TRUE, card_header(textOutput("maptitle")),
           plotOutput("map", height = "780px", click = "map_click")),
      card(card_header(textOutput("cellhead")), class = "small",
           tableOutput("cellinfo"))))),
  nav_panel("Point model", layout_sidebar(
    sidebar = sidebar(width = 320,
      selectInput("preset", "Preset", c("Custom", sp_read_all(root)$species), "Ips typographus"),
      numericInput("T0", "Base temperature T0 (°C)", 8.3, step = 0.1),
      numericInput("thr", "Threshold (degree-days)", 557, min = 1),
      checkboxInput("dia", "Diapause (photoperiod)", TRUE),
      conditionalPanel("input.dia", numericInput("dial", "Diapause daylength (h)", 14.5, step = 0.1),
                       checkboxInput("diat_on", "Also require mean temp ≤", TRUE),
                       conditionalPanel("input.diat_on", numericInput("diat", "Diapause temperature (°C)", 14.5, step = 0.1))),
      selectInput("method", "Degree-day method", c("sine", "triangle", "average")),
      sliderInput("maxY", "maxYears", 1, 4, 4, step = 1),
      selectInput("pyear", "Year", NULL),
      selectInput("pmonth", "Development starts", setNames(1:12, paste("1", month.name))),
      fluidRow(column(6, numericInput("lat", "Lat", 59.91, step = 0.01)), column(6, numericInput("lon", "Lon", 10.75, step = 0.01))),
      actionButton("go", "Run model", class = "btn-primary")),
    card(card_header(textOutput("ptitle")), plotOutput("cumplot", height = "330px"), plotOutput("tplot", height = "260px")),
    card(card_header("Results for all start months"), tableOutput("ptable")))),
  nav_panel("Species", layout_sidebar(
    sidebar = sidebar(width = 340, open = "always",
      accordion(id = "sp_acc", open = "Add a species",
        accordion_panel("Add a species", icon = icon("plus"),
          textInput("sp_name", "Species (scientific name) *"),
          textInput("sp_common", "Common name"),
          textInput("sp_event", "Life stage / event", placeholder = "e.g. Egg to adult"),
          fluidRow(column(6, numericInput("sp_T0", "T0 (°C) *", NA, step = 0.1)),
                   column(6, numericInput("sp_thr", "Threshold (DD) *", NA, min = 1))),
          selectInput("sp_dtype", "Diapause type (label only)", c("", "facultative", "obligate", "quiescence")),
          fluidRow(column(6, numericInput("sp_dl", "Daylength (h)", NA, min = 0, max = 24, step = 0.1)),
                   column(6, numericInput("sp_dt", "Temp. (°C)", NA, step = 0.1))),
          helpText("Leave daylength empty for no diapause. A temperature requires a daylength."),
          textInput("sp_ref1", "Reference 1", placeholder = "https://doi.org/..."),
          textInput("sp_ref2", "Reference 2"), textInput("sp_ref3", "Reference 3"),
          actionButton("sp_add", "Save species", class = "btn-primary", icon = icon("floppy-disk"))),
        accordion_panel("Compute maps", icon = icon("play"),
          selectInput("cp_species", "Species", NULL),
          sliderInput("cp_years", "Years", 1991, 2020, c(2018, 2020), step = 1, sep = ""),
          numericInput("cp_workers", "Parallel workers", 6, min = 1, max = 14),
          helpText("All of Norway at 1 km: ~10-30 s per year with 6 workers. Runs in the background;",
                   "finished years are skipped. 30-yr normals are added when 1991-2020 is covered."),
          actionButton("cp_go", "Compute", class = "btn-success", icon = icon("play")),
          uiOutput("cp_status")),
        accordion_panel("Delete a user species", icon = icon("trash"),
          selectInput("del_species", "User species", NULL),
          helpText("Removes the row from the table; computed maps stay in 5_outputs/opram."),
          actionButton("del_go", "Delete", class = "btn-outline-danger")))),
    card(card_header("Species parameters and references"),
      tableOutput("sptable"),
      tags$small(class = "text-muted",
        "T0: base temperature for development; Threshold: degree-days from start of development to adult emergence. ",
        "Diapause: development stops after the first day (from DOY 150) when daylength falls below the given hours; ",
        "if a temperature is given, only on days with mean temperature at or below it. ",
        "OPRAM rows: data/species_parameters.csv (J. Yearsley, GPL-3; one malformed DOI link repaired for display). ",
        "User rows: 3_metadata/species_parameters_user.csv."))))
)

# ------------------------------------------------------------------------------------ server
server <- function(input, output, session) {
  # ---- species table (OPRAM + user), add / delete ----
  sp_df <- reactiveVal(sp_read_all(root))
  species <- reactive(sp_params(sp_df()))
  sets_tick <- reactiveVal(0)                                     # bump to refresh the Maps set list
  output$sptable <- renderTable(species_tbl(sp_df()), sanitize.text.function = identity, striped = TRUE,
                                hover = TRUE, spacing = "s", width = "100%", align = "lllrrlrrll")
  observe({
    s <- sp_df(); u <- s$species[s$source == "user"]
    updateSelectInput(session, "cp_species", choices = s$species, selected = isolate(input$cp_species) %||% u[length(u)])
    updateSelectInput(session, "del_species", choices = u)
    updateSelectInput(session, "preset", choices = c("Custom", s$species), selected = isolate(input$preset))
  })
  observeEvent(input$sp_add, {
    num <- function(x) if (is.null(x) || is.na(x)) NA_real_ else as.numeric(x)
    txt <- function(x) if (is.null(x) || !nzchar(trimws(x))) NA_character_ else trimws(x)
    r <- list(species = trimws(input$sp_name), baseline_temperature = num(input$sp_T0), threshold = num(input$sp_thr),
              diapause = txt(input$sp_dtype), diapause_daylength = num(input$sp_dl), diapause_temperature = num(input$sp_dt),
              Reference1 = txt(input$sp_ref1), Reference2 = txt(input$sp_ref2), Reference3 = txt(input$sp_ref3),
              event = txt(input$sp_event), commonName = txt(input$sp_common))
    p <- sp_add(root, r)
    if (length(p)) { showNotification(HTML(paste(p, collapse = "<br>")), type = "error", duration = 8); return() }
    sp_df(sp_read_all(root))
    updateSelectInput(session, "cp_species", selected = r$species)
    for (k in c("sp_name", "sp_common", "sp_event", "sp_ref1", "sp_ref2", "sp_ref3")) updateTextInput(session, k, value = "")
    accordion_panel_open("sp_acc", "Compute maps")
    showNotification(sprintf("Saved %s. Choose years and press Compute to make maps.", r$species), type = "message")
  })
  observeEvent(input$del_go, {
    req(input$del_species); sp_delete(root, input$del_species); sp_df(sp_read_all(root))
    showNotification(sprintf("Deleted %s from the table.", input$del_species), type = "warning")
  })
  observe({ y <- sim_years(4); req(length(y))
    updateSliderInput(session, "cp_years", min = min(y), max = max(y), value = c(max(min(y), max(y) - 2), max(y))) })

  # ---- background compute job (one at a time; callr child process, Windows-safe) ----
  job <- reactiveVal(NULL)
  observeEvent(input$cp_go, {
    j <- job()
    if (!is.null(j) && j$proc$is_alive()) { showNotification("A compute job is already running.", type = "warning"); return() }
    req(input$cp_species); yrs <- seq(input$cp_years[1], input$cp_years[2])
    miss <- setdiff(unique(unlist(lapply(yrs, opram_years, 4L, last_meteo))), met_years)
    if (length(miss)) { showNotification(paste("Climate data missing for:", paste(miss, collapse = ", ")), type = "error"); return() }
    log <- file.path(outd, "_jobs", paste0(slug(input$cp_species), ".log")); dir.create(dirname(log), FALSE, TRUE)
    p <- callr::r_bg(function(root, name, years, workers) {
      setwd(root); source("4_scripts/run_norway.R"); opram_compute_species(name, years, workers)
    }, args = list(root = root, name = input$cp_species, years = yrs, workers = as.integer(input$cp_workers)),
    stdout = log, stderr = "2>&1", supervise = TRUE)                 # job stops if the app is closed (resumable)
    job(list(proc = p, name = input$cp_species, years = range(yrs), t0 = Sys.time(), log = log, done = FALSE))
  })
  output$cp_status <- renderUI({
    j <- job(); if (is.null(j)) return(NULL)
    alive <- j$proc$is_alive(); if (alive) invalidateLater(3000)
    mins <- as.numeric(difftime(Sys.time(), j$t0, units = "mins"))
    ln <- if (file.exists(j$log)) readLines(j$log, warn = FALSE) else character()
    stage <- tail(grep("STAGE|model:|write GeoTIFFs|Error|rror", ln, value = TRUE), 1)
    nb <- length(list.files(file.path(outd, paste0("_tmp_", j$proc$get_pid()))))
    if (alive) {
      tagList(tags$hr(), tags$b(sprintf("%s %d-%d", j$name, j$years[1], j$years[2])),
        tags$div(class = "progress my-2", tags$div(class = "progress-bar progress-bar-striped progress-bar-animated",
          style = sprintf("width:%d%%", round(100 * min(nb / n_blocks, 1))), sprintf("%d%%", round(100 * nb / n_blocks)))),
        tags$small(sprintf("%.1f min – %s", mins, sub("^.*STAGE ", "", stage %||% "starting"))))
    } else {
      ok <- j$proc$get_exit_status() == 0
      if (ok && !isTRUE(j$done)) { j$done <- TRUE; job(j); sets_tick(sets_tick() + 1) }
      tagList(tags$hr(), if (ok) tags$div(class = "text-success", icon("check"),
                sprintf(" %s %d-%d done in %.1f min – see Maps.", j$name, j$years[1], j$years[2], mins))
              else tags$div(class = "text-danger", icon("triangle-exclamation"), " Job failed:",
                tags$pre(style = "white-space:pre-wrap;font-size:0.75em", paste(tail(ln, 6), collapse = "\n"))))
    }
  })

  # ---- maps ----
  observe({
    sets_tick(); invalidateLater(120000)                          # also pick up results written by batch runs
    s <- result_sets(); cur <- isolate(input$set)
    updateSelectInput(session, "set", choices = setNames(s, vapply(s, set_label, "", nms = isolate(sp_df())$species)),
                      selected = if (!is.null(cur) && cur %in% s) cur else if ("ips_typographus" %in% s) "ips_typographus" else s[1])
  })
  observe({ req(input$set); sets_tick(); invalidateLater(120000); y <- set_years(input$set); cur <- isolate(input$year)
    updateSelectInput(session, "year", choices = y, selected = if (!is.null(cur) && cur %in% y) cur else tail(y, 1)) })
  cell <- reactiveVal(cellFromXY(tmpl, xy_of(10.75, 59.91)))

  layer <- reactive({
    req(input$set, input$year, input$var, input$month); d <- file.path(outd, input$set)
    fy <- file.path(d, sprintf("%s_%s_%s.tif", input$set, input$year, input$var))
    validate(need(file.exists(fy), sprintf("No %s results for %s yet.", input$year, set_label(input$set, sp_df()$species))))
    f <- function(tag) rast(file.path(d, sprintf("%s_%s.tif", tag, input$var)))[[input$month]]
    yr <- f(sprintf("%s_%s", input$set, input$year)); nt <- set_normal(input$set)
    if (input$view != "year") validate(need(!is.na(nt), "No 30-yr normal for this set yet (run opram_normals())."))
    switch(input$view, year = yr, normal = f(nt), anom = yr - f(nt))
  })
  mapdf <- reactive({                                   # Norway cells only, native 1 km
    r <- layer()          # evaluate first: req() inside an S4 argument would surface as a red error
    v <- values(r, mat = FALSE)[ncell_nor]
    rbind(data.frame(x = xy_nor[, 1], y = xy_nor[, 2], v = v), pad_nor)
  })
  output$maptitle <- renderText(sprintf("%s – %s – %s, development from 1 %s", set_label(input$set, sp_df()$species),
    c(emerge = "adult emergence", ngen = "generations")[input$var],
    c(year = input$year, normal = "30-yr median", anom = paste(input$year, "anomaly"))[input$view],
    month.name[match(sub("start_", "", input$month), month.abb)]))
  output$map <- renderPlot({
    df <- mapdf(); em <- input$var == "emerge"; an <- input$view == "anom"
    na <- "transparent"
    sc <- if (an) scale_fill_viridis_c(option = "cividis", na.value = na, name = if (em) "Days vs\nmedian" else "Generations\nvs median")
          else if (em) scale_fill_viridis_c(option = "magma", direction = -1, na.value = na, name = "Emergence",
                                            labels = if (input$view == "year") day_lab(as.integer(input$year)) else function(v) paste("day", v))
          else scale_fill_viridis_c(name = "Generations", na.value = na)
    xy <- xyFromCell(tmpl, cell())
    ggplot(df, aes(x, y)) +
      geom_raster(data = base_nor, aes(alpha = a), fill = "grey80") + scale_alpha_identity() +
      geom_raster(aes(fill = v)) + sc +
      geom_path(data = nor$county, aes(x, y, group = id), inherit.aes = FALSE, colour = "grey55", linewidth = 0.15) +
      geom_path(data = nor$border, aes(x, y, group = id), inherit.aes = FALSE, colour = "grey25", linewidth = 0.25) +
      annotate("point", x = xy[1], y = xy[2], shape = 21, size = 4, stroke = 1.5, colour = "black", fill = NA) +
      coord_equal(xlim = nor$ext[1:2], ylim = nor$ext[3:4], expand = FALSE) +
      theme_void(base_size = 13) +
      theme(legend.position = c(0.82, 0.25), legend.key.height = unit(1.4, "cm"), legend.background = element_blank())
  }, res = 110)
  observeEvent(input$map_click, {
    c1 <- cellFromXY(tmpl, cbind(input$map_click$x, input$map_click$y))
    if (!is.na(c1) && !is.na(landmap[c1])) {
      cell(c1); ll <- c(gr$lon[landmap[c1]], gr$lat[landmap[c1]])
      updateNumericInput(session, "lat", value = round(ll[2], 3)); updateNumericInput(session, "lon", value = round(ll[1], 3))
    }
  })
  output$cellhead <- renderText({ i <- landmap[cell()]; req(!is.na(i))
    sprintf("%.2f°N %.2f°E, %s", gr$lat[i], gr$lon[i], input$year) })
  output$cellinfo <- renderTable({
    req(input$set, input$year); i <- landmap[cell()]; req(!is.na(i))
    d <- file.path(outd, input$set); tag <- sprintf("%s_%s", input$set, input$year)
    e <- unlist(rast(file.path(d, paste0(tag, "_emerge.tif")))[cell()]); g <- unlist(rast(file.path(d, paste0(tag, "_ngen.tif")))[cell()])
    out <- data.frame(Start = month.abb, Emerges = short_lab(as.integer(input$year))(e), Gen = sprintf("%.2f", g))
    nt <- set_normal(input$set)
    if (!is.na(nt)) { en <- unlist(rast(file.path(d, paste0(nt, "_emerge.tif")))[cell()])
      out[["Δd"]] <- ifelse(is.na(e - en), "–", sprintf("%+.0f", e - en)) }
    out
  }, striped = TRUE, spacing = "xs", width = "100%")

  # ---- point model ----
  observeEvent(input$preset, { p <- species()[[input$preset]]; req(p)
    updateNumericInput(session, "T0", value = p$T0); updateNumericInput(session, "thr", value = p$threshold)
    updateCheckboxInput(session, "dia", value = !is.null(p$diapause_daylength))
    if (!is.null(p$diapause_daylength)) updateNumericInput(session, "dial", value = p$diapause_daylength)
    updateCheckboxInput(session, "diat_on", value = !is.null(p$diapause_temperature))
    if (!is.null(p$diapause_temperature)) updateNumericInput(session, "diat", value = p$diapause_temperature) })
  observeEvent(input$maxY, { y <- sim_years(input$maxY); updateSelectInput(session, "pyear", choices = y, selected = tail(y, 1)) })

  prun <- eventReactive(input$go, {
    req(input$pyear); c1 <- cellFromXY(tmpl, xy_of(input$lon, input$lat)); i <- landmap[c1]
    validate(need(!is.na(i), "Location is outside Norway's land cells."))
    p <- opram_params(input$T0, input$thr, if (input$dia) input$dial, if (input$dia && input$diat_on) input$diat,
                      input$method, input$maxY)
    yrs <- opram_years(as.integer(input$pyear), p$maxYears, last_meteo)
    rd <- function(v) unlist(lapply(yrs, function(y) sn_read_bin(v, y, i, 1L, out = proc))) / 100
    dates <- do.call(c, lapply(yrs, function(y) seq(as.Date(sprintf("%d-01-01", y)), as.Date(sprintf("%d-12-31", y)), "day")))
    tn <- matrix(rd("tn")); tx <- matrix(rd("tx"))
    r <- opram_run(tn, tx, dates, gr$lat[i], p)
    dd <- degree_days(tn, tx, p, as.integer(format(dates, "%j")), gr$lat[i]); ci <- cum_index(dd)
    list(p = p, y = as.integer(input$pyear), lat = gr$lat[i], lon = gr$lon[i], tn = tn[, 1], tx = tx[, 1], dd = dd[, 1],
         C = ci$C[, 1], ci = ci, r = r, dia = if (!is.null(p$diapause_daylength)) diapause_doy(gr$lat[i], p$diapause_daylength))
  })
  chain <- function(z, s) {          # generation chain from start day s (as OPRAM extract_results)
    out <- NULL
    while (s > 0) { e <- emerge_day(z$ci, 1L, s, z$p$threshold); if (is.na(e)) break
      out <- rbind(out, data.frame(start = s, emerge = e)); s <- if (e <= 365 && e + 1 < 365) e + 1L else 0L }
    out
  }
  output$ptitle <- renderText({ z <- prun(); sprintf("%.3f°N %.3f°E – %d, T0 %.1f°C, %g DD, %s%s", z$lat, z$lon, z$y,
    z$p$T0, z$p$threshold, z$p$method, if (is.null(z$dia)) "" else sprintf(", diapause from day %d", z$dia)) })
  output$cumplot <- renderPlot({
    z <- prun(); n <- length(z$C); d0 <- as.Date(sprintf("%d-01-01", z$y)); s <- z$r$start_doy[as.integer(input$pmonth)]
    ch <- chain(z, s); df <- data.frame(date = d0 + seq_len(n) - 1, C = z$C)
    g <- ggplot(df, aes(date, C)) + geom_line(linewidth = 0.8, colour = viridisLite::viridis(1)) +
      labs(x = NULL, y = "Cumulative degree-days") + theme_minimal(base_size = 13) +
      annotate("rect", xmin = d0, xmax = d0 + 364, ymin = -Inf, ymax = Inf, alpha = 0.06)
    if (!is.null(ch)) {
      ch$gen <- factor(seq_len(nrow(ch))); ch$Cs <- c(0, z$C)[ch$start]; ch$Ce <- z$C[ch$emerge]
      g <- g + geom_segment(data = ch, aes(x = d0 + start - 1, xend = d0 + emerge - 1, y = Cs, yend = Cs + z$p$threshold, colour = gen),
                            inherit.aes = FALSE, linewidth = 1, linetype = 2) +
        geom_point(data = ch, aes(d0 + emerge - 1, Ce, colour = gen), inherit.aes = FALSE, size = 3) +
        scale_colour_viridis_d(end = 0.85, name = "Generation")
    }
    g
  }, res = 96)
  output$tplot <- renderPlot({
    z <- prun(); n <- length(z$tn); df <- data.frame(date = as.Date(sprintf("%d-01-01", z$y)) + seq_len(n) - 1, tn = z$tn, tx = z$tx)
    pal <- viridisLite::viridis(3, end = 0.8)
    ggplot(df, aes(date)) + geom_ribbon(aes(ymin = tn, ymax = tx), fill = pal[2], alpha = 0.35) +
      geom_hline(yintercept = z$p$T0, colour = pal[1], linetype = 2) +
      geom_col(aes(y = z$dd), fill = pal[3], width = 1) +
      labs(x = NULL, y = "°C  /  daily DD") + theme_minimal(base_size = 13)
  }, res = 96)
  output$ptable <- renderTable({
    z <- prun(); lab <- day_lab(z$y)
    data.frame(start = paste("1", month.abb), emergence = lab(z$r$emerge), emergence_day = z$r$emerge, generations = round(z$r$ngen, 2))
  }, na = "–", striped = TRUE, spacing = "xs")
}

shinyApp(ui, server)
