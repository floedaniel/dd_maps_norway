# FinnPRIO degree-day maps -- Shiny app: degree-day pest phenology maps for Norway (1 km, seNorge2018).
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
  data.frame(Species = paste0(sprintf("<i>%s</i>", esc(s$species)),
                             ifelse(s$source == "user", " <span class='dd-badge'>user</span>", "")),
             `Common name` = esc(s$commonName),
             Event = esc(s$event), `T0 (°C)` = num(s$baseline_temperature, 2),
             `Degree-days` = num(s$threshold, 2), `Upper (°C)` = num(s$upper_threshold),
             Method = esc(s$method), Covers = ifelse(s$emergence_only, "emergence only", "full cycle"),
             Diapause = paste0(ifelse(is.na(s$diapause), "–", esc(s$diapause)),
                               ifelse(is.na(s$diapause_daylength), "",
                                      paste0(": ", num(s$diapause_daylength), " h",
                                             ifelse(is.na(s$diapause_temperature), "", paste0(", ≤ ", num(s$diapause_temperature), " °C"))))),
             References = refs, check.names = FALSE)
}
last_meteo <- 2025L
slug <- function(x) gsub("[^a-z0-9]+", "_", tolower(trimws(x)))
by <- function(v) sub("^(tn|tx)_(\\d+)\\.bin$", "\\2", list.files(file.path(proc, "bin"), sprintf("^%s_\\d+\\.bin$", v)))
met_years <- sort(as.integer(intersect(by("tn"), by("tx"))))
# simulation years whose whole window is observed (no padding with repeated weather)
sim_years <- function(mY) met_years[vapply(met_years, function(y) y + mY - 1 <= last_meteo && all((y:(y + mY - 1)) %in% met_years), TRUE)]

result_sets <- function() {
  d <- list.dirs(outd, FALSE, FALSE)
  d[vapply(d, function(x) length(list.files(file.path(outd, x), "_emerge\\.tif$")) > 0, TRUE)]
}
set_label <- function(d, nms) { n <- nms[slug(nms) == d]; if (length(n)) n[1] else d }
set_years <- function(d) sort(unique(sub(".*_(\\d{4})_emerge\\.tif$", "\\1",
                                         list.files(file.path(outd, d), "_\\d{4}_emerge\\.tif$"))))
# reference periods (computed by opram_summaries() in run_norway.R): "y1_y2"
periods <- c("1991–2020 (standard normal)" = "1991_2020", "2008–2022 (recent 15 years)" = "2008_2022")
period_lab <- function(pr) sub("_", "\u2013", pr)
set_normal <- function(d, pr) { f <- file.path(outd, d, sprintf("%s_normal_%s", d, pr))
  if (file.exists(paste0(f, "_emerge.tif"))) f else NA }
set_prob <- function(d, k, pr) { f <- file.path(outd, d, sprintf("%s_prob_%s_ngen_ge%d.tif", d, pr, k))
  if (file.exists(f)) f else NA }
# species rows flagged "emergence only": generation counts are not meaningful
set_emonly <- function(d, df) any(df$emergence_only[slug(df$species) == d])
# a compute job for this species running (also from another session)? lock file holds the PID
job_running <- function(name) {
  f <- file.path(outd, "_jobs", paste0(slug(name), ".lock")); if (!file.exists(f)) return(FALSE)
  pid <- suppressWarnings(as.integer(readLines(f, warn = FALSE)[1])); if (is.na(pid)) return(FALSE)
  out <- suppressWarnings(system2("tasklist", c("/FI", shQuote(sprintf("PID eq %d", pid)), "/NH"), stdout = TRUE))
  any(grepl(sprintf("\\b%d\\b", pid), out))
}

day_lab <- function(y, fmt = "%d %b %Y") function(v) ifelse(is.na(v), "", format(as.Date(sprintf("%d-01-01", y)) + v - 1, fmt))
short_lab <- function(y) function(v) {                             # "07 Jul", "12 Jun +1y"
  d <- as.Date(sprintf("%d-01-01", y)) + v - 1; k <- as.integer(format(d, "%Y")) - y
  ifelse(is.na(v), "–", paste0(format(d, "%d %b"), ifelse(k > 0, sprintf(" +%dy", k), "")))
}
xy_of <- function(lon, lat) crds(project(vect(cbind(lon, lat), crs = "EPSG:4326"), gr$crs))

# ---------------------------------------------------------------------------------------- UI
source(file.path(root, "4_scripts/ddmaps_app/R/theme.R"), local = TRUE)   # shared VKM / FinnPRIO look
ui <- page(
  title = "FinnPRIO degree-day maps", theme = hm_theme(),
  tags$head(tags$link(rel = "stylesheet", href = "ddmaps.css")),
  vkm_header("FinnPRIO: Degree-day pest phenology maps for Norway, 1 km"),
  div(class = "vkm-navbar", navset_bar(id = "tabs", fillable = FALSE,
  nav_panel("Maps", layout_sidebar(
    sidebar = sidebar(width = 300,
      selectInput("set", "Species / parameter set", NULL),
      selectInput("year", "Year", NULL),
      selectInput("month", "Development starts", setNames(paste0("start_", month.abb), paste("1", month.name))),
      radioButtons("var", "Show", c("Adult emergence date" = "emerge", "Generations per year" = "ngen")),
      radioButtons("view", "Values", c("Selected year" = "year", "Median of reference period" = "normal",
                                       "Anomaly (year − median)" = "anom",
                                       "Probability: ≥1 generation" = "p1", "Probability: ≥2 generations" = "p2")),
      selectInput("period", "Reference period", periods),
      conditionalPanel("input.view == 'p1' || input.view == 'p2'",
        helpText("Share of years in the reference period in which the generation(s) complete within the calendar year.")),
      uiOutput("emonly_note"),
      helpText("Click the map to inspect a 1 km cell. The same cell is used in the point model.")),
    layout_columns(col_widths = c(8, 4), fill = TRUE,
      card(full_screen = TRUE, card_header(textOutput("maptitle")),
           plotOutput("map", height = "780px", click = "map_click"),
           card_footer(class = "small text-muted",
             "Climate: seNorge_2018 daily Tmin/Tmax, 1 km, © MET Norway (CC BY 4.0). ",
             "Borders: GADM 4.1. Details under About.")),
      div(class = "dd-right",
        card(card_header(textOutput("cellhead")), class = "dd-cell",
             tableOutput("cellinfo"),
             card_footer(class = "small text-muted", textOutput("cellnote"))),
        card(card_header(textOutput("trendhead")), class = "dd-trend",
             plotOutput("trend", height = "300px")))))),
  nav_panel("Point model", layout_sidebar(
    sidebar = sidebar(width = 320,
      selectInput("preset", "Preset", c("Custom", sp_read_all(root)$species), "Ips typographus"),
      numericInput("T0", "Lower threshold T0 (°C)", 8.3, step = 0.1),
      numericInput("thr", "Degree-days required (°C·days)", 557, min = 1),
      checkboxInput("dia", "Diapause (photoperiod)", TRUE),
      conditionalPanel("input.dia", numericInput("dial", "Diapause daylength (h)", 14.5, step = 0.1),
                       checkboxInput("diat_on", "Only on days with mean temperature at or below", TRUE),
                       conditionalPanel("input.diat_on", numericInput("diat", "Diapause temperature (°C)", 14.5, step = 0.1))),
      selectInput("method", "Degree-day method", c("sine", "triangle", "average")),
      numericInput("pupper", "Upper threshold (°C), optional", NA, step = 0.1),
      checkboxInput("pemonly", "Emergence only (degree-days cover part of the life cycle)", FALSE),
      sliderInput("maxY", "Years allowed for development", 1, 4, 4, step = 1),
      selectInput("pyear", "Year", NULL),
      selectInput("pmonth", "Development starts", setNames(1:12, paste("1", month.name))),
      fluidRow(column(6, numericInput("lat", "Latitude", 59.91, step = 0.01)), column(6, numericInput("lon", "Longitude", 10.75, step = 0.01))),
      helpText("Or click a cell on the Maps tab."),
      actionButton("go", "Run model", class = "btn-primary")),
    card(card_header(textOutput("ptitle")), plotOutput("cumplot", height = "330px"), plotOutput("tplot", height = "260px")),
    card(card_header("Results for all start months"), tableOutput("ptable")))),
  nav_panel("Species", layout_sidebar(
    sidebar = sidebar(width = 340, open = "always",
      accordion(id = "sp_acc", open = "Add a species",
        accordion_panel("Add a species",
          textInput("sp_name", "Species (scientific name) *"),
          textInput("sp_common", "Common name"),
          textInput("sp_event", "Life stage / event", placeholder = "e.g. Egg to adult"),
          fluidRow(column(6, numericInput("sp_T0", "Lower threshold T0 (°C) *", NA, step = 0.1)),
                   column(6, numericInput("sp_thr", "Degree-days required (°C·days) *", NA, min = 1))),
          helpText("T0: temperature below which no development occurs (often 8-15 °C).",
                   "Degree-days: heat sum from start of development to adult emergence (typically 60-1500)."),
          fluidRow(column(6, selectInput("sp_method", "Degree-day method", c("sine", "triangle", "average"))),
                   column(6, numericInput("sp_upper", "Upper threshold (°C)", NA, step = 0.1))),
          helpText("Use the method and upper threshold from the source. Above the upper threshold development",
                   "does not speed up further (horizontal cutoff). Leave it empty for none."),
          checkboxInput("sp_emonly", "Emergence only: the degree-days cover only part of the life cycle (e.g. overwintered larvae to adult)", FALSE),
          selectInput("sp_dtype", "Diapause type (label only)", c("", "facultative", "obligate", "quiescence")),
          fluidRow(column(6, numericInput("sp_dl", "Daylength (h)", NA, min = 0, max = 24, step = 0.1)),
                   column(6, numericInput("sp_dt", "Temp. (°C)", NA, step = 0.1))),
          helpText("Leave daylength empty for no diapause. A temperature requires a daylength."),
          textInput("sp_ref1", "Reference 1", placeholder = "https://doi.org/..."),
          textInput("sp_ref2", "Reference 2"), textInput("sp_ref3", "Reference 3"),
          actionButton("sp_add", "Save species", class = "btn-primary")),
        accordion_panel("Compute maps",
          selectInput("cp_species", "Species", NULL),
          sliderInput("cp_years", "Years", 1991, 2022, c(1991, 2022), step = 1, sep = ""),
          numericInput("cp_workers", "Parallel workers", 6, min = 1, max = 14),
          helpText("All of Norway at 1 km: ~10-30 s per year with 6 workers. Runs in the background;",
                   "finished years are skipped. Medians and probabilities are added for each reference period",
                   "(1991–2020, 2008–2022) that the computed years cover."),
          actionButton("cp_go", "Compute maps", class = "btn-primary"),
          uiOutput("cp_status")),
        accordion_panel("Delete a user species",
          selectInput("del_species", "User species", NULL),
          checkboxInput("del_maps", "Also delete its computed maps", FALSE),
          actionButton("del_go", "Delete species", class = "btn-outline-secondary")))),
    card(card_header("Species parameters and references"), class = "dd-species",
      tableOutput("sptable"),
      tags$small(class = "text-muted",
        "T0: base temperature for development; Degree-days: heat sum from start of development to adult emergence; ",
        "Upper: upper threshold (horizontal cutoff); Covers: full life cycle, or emergence only (no generation counts). ",
        "Diapause: type, then the daylength below which development stops (from DOY 150) and, if given, ",
        "the mean temperature at or below which it stops; ",

        "Rows without a badge come from OPRAM data/species_parameters.csv (J. Yearsley, GPL-3; one malformed DOI link repaired for display). ",
        "Rows marked user are your own (3_metadata/species_parameters_user.csv).")))),
  nav_panel("About", layout_columns(col_widths = c(7, 5),
    card(card_header("Climate data"),
      tags$p(tags$b("seNorge_2018"), " daily minimum (TN) and maximum (TX) 2 m air temperature, ",
             "Norwegian Meteorological Institute (MET Norway). An observational gridded dataset: station data from ",
             "MET Norway, SMHI, FMI and ECA&D, spatially interpolated to a 1 km grid with automatic quality control."),
      tags$table(class = "table table-sm small",
        tags$tbody(
          tags$tr(tags$th("Variables"), tags$td("TN daily minimum, TX daily maximum temperature (°C)")),
          tags$tr(tags$th("Version"), tags$td("seNorge_2018 v23.09 for both TN and TX (yearly files from MET Norway's THREDDS archive)")),
          tags$tr(tags$th("Grid"), tags$td("1 km, UTM zone 33N (1195 × 1550 cells); Norway land cells only in this app")),
          tags$tr(tags$th("Day definition"), tags$td("18 UTC previous day to 18 UTC; labelled with the end date")),
          tags$tr(tags$th("Period used"), tags$td("Weather 1991–2025; maps for 1991–2022, so every 4-year development window is observed weather")),
          tags$tr(tags$th("License"), tags$td("CC BY 4.0 or NLOD; credit: “MET Norway”")))),
      tags$p(tags$b("Please cite:")),
      tags$ul(class = "small",
        tags$li("Lussana, C., Tveito, O. E., Dobler, A. & Tunheim, K. (2019). seNorge_2018, daily precipitation, and ",
                "temperature datasets over Norway. ", tags$i("Earth System Science Data"), " 11, 1531–1551. ",
                tags$a(href = "https://doi.org/10.5194/essd-11-1531-2019", target = "_blank", "doi:10.5194/essd-11-1531-2019")),
        tags$li("Lussana, C., Tveito, O. E. & Uboldi, F. (2018). Three-dimensional spatial interpolation of 2 m ",
                "temperature over Norway. ", tags$i("Q. J. R. Meteorol. Soc."), " ",
                tags$a(href = "https://doi.org/10.1002/qj.3208", target = "_blank", "doi:10.1002/qj.3208")),
        tags$li("Data: MET Norway THREDDS archive, ",
                tags$a(href = "https://thredds.met.no/thredds/catalog/senorge/seNorge_2018/Archive/catalog.html",
                       target = "_blank", "senorge/seNorge_2018/Archive", .noWS = "outside"), ".")),
      tags$p(tags$b("Processing in this app:"), class = "mb-1"),
      tags$ul(class = "small",
        tags$li("Values stored at 0.01 °C resolution (rounding error ≤ 0.005 °C)."),
        tags$li("Leap days are included; a generation counts for the year if it completes by day 365."))),
    card(card_header("Model, boundaries and limitations"),
      tags$p(tags$b("Model: "), "degree-day phenology model ported to R from ",
             tags$a(href = "https://github.com/DrJonYearsley/OPRAM", target = "_blank", "OPRAM", .noWS = "after"),
             " (Jon Yearsley, GPL-3; the model behind ",
             tags$a(href = "https://www.pestmaps.ie", target = "_blank", "pestmaps.ie", .noWS = "outside"), "). ",
             "Daily degree-days from TN/TX above a lower threshold T0, with the method given per species (single sine ",
             "by default, single triangle or simple average) and an optional upper threshold (horizontal cutoff: above it ",
             "development does not speed up further). Adult emergence when the degree-day sum is reached; optional ",
             "photoperiod/temperature diapause."),
      tags$p(tags$b("Summaries: "), "medians, anomalies and probabilities are computed for two reference periods: ",
             "1991–2020 (the standard climate normal) and 2008–2022 (the most recent 15 years with fully observed ",
             "4-year windows). The trend panel shows all computed years at the selected cell with a linear trend."),
      tags$p(tags$b("Borders: "), "GADM 4.1 (", tags$a(href = "https://gadm.org", target = "_blank", "gadm.org", .noWS = "outside"), ")."),
      tags$p(tags$b("Limitations: "), "maps show where development is thermally possible. No upper temperature ",
             "threshold, mortality, winter survival, host availability or microclimate (e.g. bark or soil temperature). ",
             "Not validated against Norwegian observations. Species marked 'emergence only' describe part of the life ",
             "cycle (e.g. overwintered larvae to adult): their emergence dates are valid, generation counts are hidden."),
      tags$p(class = "small text-muted", "Source code: ",
             tags$a(href = "https://github.com/floedaniel/dd_maps_norway", target = "_blank", "github.com/floedaniel/dd_maps_norway"),
             " (GPL-3).")))
  ))),
  vkm_footer("This app")
)

# ------------------------------------------------------------------------------------ server
server <- function(input, output, session) {
  # ---- species table (OPRAM + user), add / delete ----
  sp_df <- reactiveVal(sp_read_all(root))
  species <- reactive(sp_params(sp_df()))
  sets_tick <- reactiveVal(0)                                     # bump to refresh the Maps set list
  output$sptable <- renderTable(species_tbl(sp_df()), sanitize.text.function = identity, striped = TRUE,
                                hover = TRUE, spacing = "s", width = "100%", align = "lllrrrllll")
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
              event = txt(input$sp_event), commonName = txt(input$sp_common),
              method = input$sp_method, upper_threshold = num(input$sp_upper), emergence_only = isTRUE(input$sp_emonly))
    p <- sp_add(root, r)
    if (length(p)) { showNotification(HTML(paste(p, collapse = "<br>")), type = "error", duration = 8); return() }
    sp_df(sp_read_all(root))
    updateSelectInput(session, "cp_species", selected = r$species)
    for (k in c("sp_name", "sp_common", "sp_event", "sp_ref1", "sp_ref2", "sp_ref3")) updateTextInput(session, k, value = "")
    accordion_panel_open("sp_acc", "Compute maps")
    showNotification(sprintf("Saved %s. Choose years and press Compute to make maps.", r$species), type = "message")
  })
  observeEvent(input$del_go, {
    req(input$del_species); dm <- isTRUE(input$del_maps)
    sp_delete(root, input$del_species, delete_maps = dm); sp_df(sp_read_all(root)); if (dm) sets_tick(sets_tick() + 1)
    showNotification(sprintf("Deleted %s%s.", input$del_species, if (dm) " and its maps" else " from the table (maps kept)"),
                     type = "warning")
  })
  observe({ y <- sim_years(4); req(length(y))
    updateSliderInput(session, "cp_years", min = min(y), max = max(y), value = c(min(y), max(y))) })

  # ---- background compute job (one at a time; callr child process, Windows-safe) ----
  job <- reactiveVal(NULL)
  observeEvent(input$cp_go, {
    j <- job()
    if (!is.null(j) && j$proc$is_alive()) { showNotification("A compute job is already running.", type = "warning"); return() }
    req(input$cp_species); yrs <- seq(input$cp_years[1], input$cp_years[2])
    if (job_running(input$cp_species)) {
      showNotification(sprintf("Maps for %s are already being computed (another session). Try again when it has finished.",
                               input$cp_species), type = "warning", duration = 8); return() }
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

  emonly <- reactive({ req(input$set); set_emonly(input$set, sp_df()) })
  output$emonly_note <- renderUI(if (emonly()) div(class = "dd-note",
    "This species' degree-days cover only part of the life cycle, so generation counts are not shown.",
    "Use the emergence date (1 January start) or 'Probability ≥1', which here means adults emerge in the same season."))
  layer <- reactive({
    req(input$set, input$year, input$var, input$month, input$period); d <- file.path(outd, input$set)
    validate(need(!(emonly() && (input$var == "ngen" || input$view == "p2")),
                  "Not shown: generation counts are not meaningful for an 'emergence only' species."))
    if (input$view %in% c("p1", "p2")) {                 # probability maps (precomputed, values 0-1)
      fp <- set_prob(input$set, if (input$view == "p1") 1 else 2, input$period)
      validate(need(!is.na(fp), sprintf("No probability maps for %s yet (compute all years of the period).", period_lab(input$period))))
      return(rast(fp)[[input$month]] * 100)
    }
    fy <- file.path(d, sprintf("%s_%s_%s.tif", input$set, input$year, input$var))
    validate(need(file.exists(fy), sprintf("No %s results for %s yet.", input$year, set_label(input$set, sp_df()$species))))
    f <- function(path) rast(paste0(path, "_", input$var, ".tif"))[[input$month]]
    yr <- f(file.path(d, sprintf("%s_%s", input$set, input$year))); nt <- set_normal(input$set, input$period)
    if (input$view != "year") validate(need(!is.na(nt), sprintf("No %s median for this species yet (compute all years of the period).",
                                                               period_lab(input$period))))
    switch(input$view, year = yr, normal = f(nt), anom = yr - f(nt))
  })
  mapdf <- reactive({                                   # Norway cells only, native 1 km
    r <- layer()          # evaluate first: req() inside an S4 argument would surface as a red error
    v <- values(r, mat = FALSE)[ncell_nor]
    rbind(data.frame(x = xy_nor[, 1], y = xy_nor[, 2], v = v), pad_nor)
  })
  output$maptitle <- renderText({
    mon <- month.name[match(sub("start_", "", input$month), month.abb)]
    per <- period_lab(input$period)
    if (input$view %in% c("p1", "p2"))
      return(sprintf("%s – probability of %s generation%s within the year (%s), development from 1 %s",
                     set_label(input$set, sp_df()$species), if (input$view == "p1") "≥1" else "≥2",
                     if (input$view == "p1") "" else "s", per, mon))
    sprintf("%s – %s – %s, development from 1 %s", set_label(input$set, sp_df()$species),
      c(emerge = "adult emergence", ngen = "generations")[input$var],
      c(year = input$year, normal = paste("median", per), anom = sprintf("%s vs median %s", input$year, per))[input$view], mon)
  })
  output$map <- renderPlot({
    df <- mapdf(); em <- input$var == "emerge"; an <- input$view == "anom"; pr <- input$view %in% c("p1", "p2")
    na <- "transparent"
    # viridis everywhere; yellow = more pest pressure (earlier emergence / more generations / higher probability)
    pmax_ <- if (pr) max(10, ceiling(round(max(df$v, na.rm = TRUE), 3) / 10) * 10) else 100   # stretch low probabilities (float32-safe)
    sc <- if (pr) scale_fill_viridis_c(option = "viridis", limits = c(0, pmax_), na.value = na,
                                       labels = function(v) paste0(v, "%"),
                                       name = if (input$view == "p1") "Years with\n≥1 generation" else "Years with\n≥2 generations")
          else if (an) scale_fill_viridis_c(option = "viridis", direction = if (em) -1 else 1, na.value = na,
                                       name = if (em) "Days vs\nmedian" else "Generations\nvs median")
          else if (em) {                                       # first season in colour; later = darkest class
            hi <- 365; lo <- min(floor(min(df$v, na.rm = TRUE)), hi - 60)
            lab <- if (input$view == "year") day_lab(as.integer(input$year), "%d %b")
                   else function(v) format(as.Date("2001-01-01") + v - 1, "%d %b")
            scale_fill_viridis_c(option = "viridis", direction = -1, na.value = na, name = "Emergence",
                                 limits = c(lo, hi), oob = scales::squish,
                                 breaks = function(l) { b <- pretty(l, 4); c(b[b >= l[1] & b < hi - 25], hi) },
                                 labels = function(b) ifelse(b >= hi, "31 Dec or later", lab(b)))
          }
          else scale_fill_viridis_c(option = "viridis", name = "Generations", na.value = na)
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
    if (emonly()) out$Gen <- NULL
    nt <- set_normal(input$set, input$period)
    if (!is.na(nt)) { en <- unlist(rast(paste0(nt, "_emerge.tif"))[cell()])
      out[["Δd"]] <- ifelse(is.na(e - en), "–", sprintf("%+.0f", e - en)) }
    fp <- set_prob(input$set, 1, input$period)
    if (!is.na(fp)) out[["P≥1"]] <- sprintf("%.0f%%", 100 * unlist(rast(fp)[cell()]))
    out
  }, striped = TRUE, spacing = "xs", width = "100%")

  output$cellnote <- renderText({ per <- period_lab(input$period)
    sprintf("Δd: days earlier (−) or later (+) than the %s median. P≥1: share of years %s with at least one %s.",
            per, per, if (emonly()) "emergence in the same season" else "generation") })

  # ---- trend at the selected cell: all computed years, selected start month ----
  trend_data <- reactive({
    req(input$set, input$month); i <- landmap[cell()]; req(!is.na(i)); sets_tick()
    y <- as.integer(set_years(input$set)); req(length(y) > 2); d <- file.path(outd, input$set)
    m <- match(input$month, paste0("start_", month.abb))
    rd <- function(v) vapply(y, function(yy) unlist(rast(file.path(d, sprintf("%s_%d_%s.tif", input$set, yy, v)))[cell()])[[m]], 0)
    data.frame(year = y, emerge = rd("emerge"), ngen = rd("ngen"))
  })
  output$trendhead <- renderText({ z <- trend_data()
    sprintf("Trend %d–%d at this cell, development from 1 %s", min(z$year), max(z$year),
            month.name[match(sub("start_", "", input$month), month.abb)]) })
  output$trend <- renderPlot({
    z <- trend_data(); pal <- viridisLite::viridis(3, end = 0.8)
    long <- rbind(data.frame(year = z$year, v = z$emerge, what = "Adult emergence (day of year)"),
                  if (!emonly()) data.frame(year = z$year, v = z$ngen, what = "Generations per year"))
    long <- long[!is.na(long$v), ]; validate(need(nrow(long) > 2, "No emergence at this cell in the computed years."))
    slope <- function(w) { k <- long[long$what == w, ]; if (nrow(k) < 5) return(NA); 10 * coef(lm(v ~ year, k))[[2]] }
    se <- slope("Adult emergence (day of year)"); sg <- if (!emonly()) slope("Generations per year") else NA
    sub <- paste(c(if (!is.na(se)) sprintf("Emergence %+.1f days per decade", se),
                   if (!is.na(sg)) sprintf("Generations %+.2f per decade", sg)), collapse = "\n")
    ggplot(long, aes(year, v)) + geom_line(colour = pal[2], linewidth = 0.4) + geom_point(colour = pal[1], size = 1.6) +
      geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = pal[3], linewidth = 0.7, linetype = 2) +
      facet_wrap(~what, ncol = 1, scales = "free_y") + labs(x = NULL, y = NULL, subtitle = sub) +
      theme_minimal(base_size = 11) + theme(strip.text = element_text(hjust = 0, face = "bold"), panel.grid.minor = element_blank())
  }, res = 96)

  # ---- point model ----
  observeEvent(input$preset, { p <- species()[[input$preset]]; req(p)
    updateSelectInput(session, "method", selected = p$method)
    updateNumericInput(session, "pupper", value = if (is.null(p$upper)) NA else p$upper)
    updateCheckboxInput(session, "pemonly", value = isTRUE(p$emergence_only))
    updateNumericInput(session, "T0", value = p$T0); updateNumericInput(session, "thr", value = p$threshold)
    updateCheckboxInput(session, "dia", value = !is.null(p$diapause_daylength))
    if (!is.null(p$diapause_daylength)) updateNumericInput(session, "dial", value = p$diapause_daylength)
    updateCheckboxInput(session, "diat_on", value = !is.null(p$diapause_temperature))
    if (!is.null(p$diapause_temperature)) updateNumericInput(session, "diat", value = p$diapause_temperature) })
  observeEvent(input$maxY, { y <- sim_years(input$maxY); updateSelectInput(session, "pyear", choices = y, selected = tail(y, 1)) })

  # run once automatically when the year list is ready (Oslo, preset species), then on "Run model"
  auto_run <- reactiveVal(0)
  observeEvent(input$pyear, if (auto_run() == 0 && nzchar(input$pyear)) auto_run(1), ignoreInit = TRUE)
  prun <- eventReactive(list(input$go, auto_run()), {
    req(auto_run() > 0 || input$go > 0)
    req(input$pyear); c1 <- cellFromXY(tmpl, xy_of(input$lon, input$lat)); i <- landmap[c1]
    validate(need(!is.na(i), "Location is outside Norway's land cells."))
    p <- opram_params(input$T0, input$thr, if (input$dia) input$dial, if (input$dia && input$diat_on) input$diat,
                      input$method, input$maxY, upper = input$pupper, emergence_only = isTRUE(input$pemonly))
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
  output$ptitle <- renderText({ z <- prun(); sprintf("%.3f°N %.3f°E – %d, T0 %.1f°C, %g DD, %s%s%s%s", z$lat, z$lon, z$y,
    z$p$T0, z$p$threshold, z$p$method, if (is.null(z$p$upper)) "" else sprintf(", upper %.1f°C", z$p$upper),
    if (is.null(z$dia)) "" else sprintf(", diapause from day %d", z$dia), if (z$p$emergence_only) ", emergence only" else "") })
  output$cumplot <- renderPlot({
    z <- prun(); n <- length(z$C); d0 <- as.Date(sprintf("%d-01-01", z$y)); s <- z$r$start_doy[as.integer(input$pmonth)]
    ch <- chain(z, s); df <- data.frame(date = d0 + seq_len(n) - 1, C = z$C)
    if (z$p$emergence_only && !is.null(ch)) ch <- ch[1, , drop = FALSE]   # partial life cycle: first emergence only
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
    out <- data.frame(Start = paste("1", month.abb), Emergence = lab(z$r$emerge), `Day of year` = z$r$emerge,
                      Generations = round(z$r$ngen, 2), check.names = FALSE)
    if (z$p$emergence_only) out$Generations <- NULL
    out
  }, na = "–", striped = TRUE, spacing = "xs")
}

shinyApp(ui, server)
