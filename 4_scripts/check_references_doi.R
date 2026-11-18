# check_references_doi.R -- find and verify DOIs for every reference in a degree-day parameter table.
# Part of DD-maps Norway (FinnPRIO), GPL-3. Run from project root:
#   Rscript 4_scripts/check_references_doi.R <input.xlsx> [output.xlsx]
#
# Per reference (References column, split on ";"):
#   - DOI written in the citation   -> checked in Crossref (then DataCite)            -> "DOI in citation, verified"
#   - otherwise, with author + year -> Crossref bibliographic search (citation + species);
#       accepted only if first-author surname AND year (±1) match                     -> "found (author + year match)"
#       only one of them matches                                                       -> "possible - check"
#   - EPPO/USDA/web/model documents without year                                       -> "no DOI expected"
# Output: copy of the input with columns "DOI" (verified links), "DOI to check" (possible matches)
#         and "DOI check" (summary), plus a sheet "Reference check" (one row per reference).
# Only bibliographic text is sent to Crossref/DataCite; responses are cached (resumable).

suppressPackageStartupMessages({ library(jsonlite); library(openxlsx) })
a <- commandArgs(TRUE)
f_in  <- if (length(a)) a[1] else "1_raw_data/degree_day_parameters.xlsx"   # or pass the path as first argument
f_out <- if (length(a) > 1) a[2] else sub("\\.xlsx$", "_with_doi.xlsx", f_in)
cache_f <- "2_processed_data/doi_cache.rds"
cache <- if (file.exists(cache_f)) readRDS(cache_f) else list()
ua <- "dd-maps-norway reference checker (R)"

get_json <- function(u) {
  if (!is.null(cache[[u]])) return(cache[[u]])
  h <- curl::new_handle(useragent = ua, timeout = 30)
  r <- tryCatch({ x <- curl::curl_fetch_memory(u, h)
                  if (x$status_code == 200) fromJSON(rawToChar(x$content), simplifyVector = FALSE) else list(status = x$status_code) },
                error = function(e) NULL)
  Sys.sleep(0.15)
  if (!is.null(r)) { cache[[u]] <<- r; if (length(cache) %% 25 == 0) saveRDS(cache, cache_f) }
  r
}
plain <- function(x) tolower(iconv(x, "UTF-8", "ASCII//TRANSLIT"))
surname <- function(ref) {                      # first capitalised word that is not an organisation
  w <- regmatches(ref, regexpr("[A-Z][A-Za-z\u00C0-\u017F'\u2019-]{1,}", ref))
  if (!length(w) || w %in% c("EPPO", "USDA", "CABI", "EFSA", "CDFA", "APHIS", "OSU", "IPPC", "USA", "NPN", "The", "Nature",
                             "Front", "Embrapa", "Florida", "Effects", "Study", "Forest", "Journal", "Annual", "Crop", "Plant"))
    return(NA_character_)
  w
}
year_of <- function(ref) { y <- regmatches(ref, regexpr("(19|20)[0-9]{2}", ref)); if (length(y)) as.integer(y) else NA_integer_ }
doi_rx <- "10\\.[0-9]{4,9}/[^\\s;,)\\]]+"

first <- function(x) tryCatch(if (length(x)) as.character(x[[1]]) else "", error = function(e) "")   # missing fields -> ""
meta <- function(m) {                           # Crossref work -> compact fields (any field may be missing)
  fam <- tryCatch(vapply(m$author %||% list(), function(z) as.character(z$family %||% z$name %||% ""), ""),
                  error = function(e) character())
  yr <- tryCatch(as.integer(m$issued$`date-parts`[[1]][[1]]), error = function(e) NA_integer_)
  list(doi = tolower(first(m$DOI)), title = gsub("<[^>]+>|\\s+", " ", first(m$title)), authors = paste(head(fam, 4), collapse = ", "),
       families = fam, year = yr, journal = first(m$`container-title`))
}

# does a title name the species? genus, or the species epithet (handles abbreviated genus "A. planipennis")
names_species <- function(title, species) {
  w <- strsplit(trimws(species), "\\s+")[[1]]; keys <- unique(c(w[1], if (length(w) > 1) w[2]))
  any(vapply(keys, function(k) grepl(k, title, ignore.case = TRUE, fixed = FALSE), TRUE))
}

check_ref <- function(ref, species) {
  out <- list(reference = ref, doi = NA, status = NA, title = NA, authors = NA, year = NA, journal = NA)
  d <- regmatches(ref, regexpr(doi_rx, ref, perl = TRUE))
  if (length(d)) {
    d <- sub("[.]$", "", d); r <- get_json(paste0("https://api.crossref.org/works/", URLencode(d, reserved = TRUE)))
    if (!is.null(r$message)) { m <- meta(r$message)
      st <- if (names_species(m$title, species)) "DOI in citation, verified (Crossref)"
            else "possible - check (DOI in citation exists, but its title does not name the species)"
      return(modifyList(out, list(doi = m$doi, status = st, title = m$title,
                                  authors = m$authors, year = m$year, journal = m$journal))) }
    r <- get_json(paste0("https://api.datacite.org/dois/", URLencode(d, reserved = TRUE)))
    if (!is.null(r$data)) { at <- r$data$attributes
      return(modifyList(out, list(doi = tolower(d), status = "DOI in citation, verified (DataCite)",
                                  title = at$titles[[1]]$title %||% "", year = at$publicationYear %||% NA))) }
    return(modifyList(out, list(doi = tolower(d), status = "DOI in citation, NOT found in Crossref/DataCite")))
  }
  sn <- surname(ref); yr <- year_of(ref)
  if (is.na(yr) || is.na(sn)) {
    return(modifyList(out, list(status = if (grepl("EPPO|datasheet|USDA|CDFA|APHIS|CABI|model|forecast|extension|website|http|reporting service|regulat|regs",
                                                   ref, ignore.case = TRUE)) "no DOI expected (datasheet/agency/web/model)"
                                         else "cannot search (no author or year in citation)")))
  }
  q <- URLencode(paste(ref, species), reserved = TRUE)
  r <- get_json(paste0("https://api.crossref.org/works?rows=5&select=DOI,title,author,issued,container-title&query.bibliographic=", q))
  items <- r$message$items
  if (!length(items)) return(modifyList(out, list(status = "not found")))
  cand <- lapply(items, meta)
  key <- plain(sn); if (is.na(key)) key <- tolower(sn)            # transliteration can fail -> fall back
  au <- vapply(cand, function(m) isTRUE(any(grepl(key, plain(m$families), fixed = TRUE))), TRUE)
  ye <- vapply(cand, function(m) isTRUE(!is.na(m$year) && abs(m$year - yr) <= 1), TRUE)
  sp <- vapply(cand, function(m) names_species(m$title, species), TRUE)       # title names genus or species
  pick <- if (any(au & ye & sp)) which(au & ye & sp)[1] else if (any(au & ye)) which(au & ye)[1]
          else if (any((au | ye) & sp)) which((au | ye) & sp)[1] else NA
  if (is.na(pick)) return(modifyList(out, list(status = "not found (no candidate with matching author or year)")))
  m <- cand[[pick]]
  st <- if (au[pick] && ye[pick] && sp[pick]) "found (author + year + species in title)"
        else if (au[pick] && ye[pick]) "possible - check (author + year match, title does not name the species)"
        else "possible - check (only author or only year matches)"
  modifyList(out, list(doi = m$doi, status = st, title = m$title, authors = m$authors, year = m$year, journal = m$journal))
}

d <- read.xlsx(f_in, sheet = "Degree-day parameters", check.names = FALSE, sep.names = " ")
rows <- list(); t0 <- Sys.time()
for (i in seq_len(nrow(d))) {
  refs <- trimws(unlist(strsplit(d$References[i] %||% "", ";"))); refs <- refs[nzchar(refs)]
  for (rf in refs) {
    res <- tryCatch(check_ref(rf, d$Species[i]), error = function(e)
      list(reference = rf, doi = NA, status = paste("error:", conditionMessage(e)), title = NA, authors = NA, year = NA, journal = NA))
    rows[[length(rows) + 1]] <- c(list(row = i, species = d$Species[i], eppo = d$`EPPO code`[i]), res)
  }
  if (i %% 10 == 0) message(sprintf("%d/%d species, %d references, %.1f min", i, nrow(d), length(rows),
                                    as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
saveRDS(cache, cache_f)
rc <- do.call(rbind, lapply(rows, function(z) as.data.frame(lapply(z, function(v) if (is.null(v) || !length(v)) NA else v[[1]]))))
ok <- grepl("^found|verified", rc$status)

link <- function(x) ifelse(is.na(x), NA, paste0("https://doi.org/", x))
d$DOI <- vapply(seq_len(nrow(d)), function(i) paste(unique(link(rc$doi[rc$row == i & ok])), collapse = "\n"), "")
d$`DOI to check` <- vapply(seq_len(nrow(d)), function(i) paste(unique(link(rc$doi[rc$row == i & grepl("^possible", rc$status)])), collapse = "\n"), "")
d$`DOI check` <- vapply(seq_len(nrow(d)), function(i) { s <- rc$status[rc$row == i]
  paste(sprintf("%d %s", as.integer(table(s)), names(table(s))), collapse = "; ") }, "")

wb <- createWorkbook()
addWorksheet(wb, "Degree-day parameters"); writeData(wb, 1, d)
addWorksheet(wb, "Reference check")
rc_out <- rc[, c("species", "eppo", "reference", "status", "doi", "title", "authors", "year", "journal")]
rc_out$doi_link <- link(rc_out$doi); class(rc_out$doi_link) <- "hyperlink"
writeData(wb, 2, rc_out)
s <- read.xlsx(f_in, sheet = "Summary"); addWorksheet(wb, "Summary"); writeData(wb, 3, s)
wrap <- createStyle(wrapText = TRUE, valign = "top")
addStyle(wb, 1, wrap, rows = 1:(nrow(d) + 1), cols = 1:ncol(d), gridExpand = TRUE)
setColWidths(wb, 1, cols = which(names(d) %in% c("DOI", "DOI to check")), widths = 45)
setColWidths(wb, 2, cols = 1:ncol(rc_out), widths = c(24, 9, 60, 30, 28, 60, 30, 6, 30, 40))
freezePane(wb, 1, firstRow = TRUE); freezePane(wb, 2, firstRow = TRUE)
saveWorkbook(wb, f_out, overwrite = TRUE)

cat(sprintf("\nreferences: %d | verified/found: %d | possible: %d | no DOI expected: %d | not found: %d | cannot search: %d\n",
            nrow(rc), sum(ok), sum(grepl("^possible", rc$status)), sum(grepl("no DOI expected", rc$status)),
            sum(grepl("^not found", rc$status)), sum(grepl("cannot search", rc$status))))
cat("species with >=1 verified DOI:", sum(nzchar(d$DOI)), "of", nrow(d), "\nwritten:", f_out, "\n")
