# species_db.R -- species parameter table: OPRAM presets (read-only) + user-added species.
# Same columns as OPRAM data/species_parameters.csv. Part of the R port of OPRAM, GPL-3.

sp_cols <- c("species", "baseline_temperature", "threshold", "diapause", "diapause_daylength",
             "diapause_temperature", "Reference1", "Reference2", "Reference3", "event", "commonName")
sp_file_opram <- function(root) file.path(root, "3_metadata/species_parameters_OPRAM.csv")
sp_file_user  <- function(root) file.path(root, "3_metadata/species_parameters_user.csv")

sp_read <- function(f) {
  if (!file.exists(f)) return(NULL)
  s <- read.csv(f, check.names = FALSE, na.strings = c("NA", ""), colClasses = "character")
  for (k in setdiff(sp_cols, names(s))) s[[k]] <- NA_character_
  for (k in c("baseline_temperature", "threshold", "diapause_daylength", "diapause_temperature"))
    s[[k]] <- as.numeric(s[[k]])
  s[, sp_cols, drop = FALSE]
}

# all species with a source column ("OPRAM" / "user")
sp_read_all <- function(root) {
  o <- sp_read(sp_file_opram(root)); u <- sp_read(sp_file_user(root))
  rbind(cbind(o, source = "OPRAM"), if (!is.null(u) && nrow(u)) cbind(u, source = "user"))
}

# list of opram_params, named by species
sp_params <- function(df, method = "sine", maxYears = 4L)
  setNames(lapply(seq_len(nrow(df)), function(i) with(df[i, ], opram_params(
    baseline_temperature, threshold, diapause_daylength, diapause_temperature, method, maxYears, species))),
    df$species)

sp_slug <- function(x) gsub("[^a-z0-9]+", "_", tolower(trimws(x)))

# validate one new row against existing table; returns character vector of problems (empty = ok)
sp_validate <- function(r, all) {
  p <- character()
  if (!nzchar(trimws(r$species))) p <- c(p, "Species name is required.")
  if (sp_slug(r$species) %in% sp_slug(all$species)) p <- c(p, "A species with this name already exists.")
  if (!is.finite(r$baseline_temperature) || r$baseline_temperature < -10 || r$baseline_temperature > 40)
    p <- c(p, "T0 must be between -10 and 40 °C.")
  if (!is.finite(r$threshold) || r$threshold <= 0 || r$threshold > 20000) p <- c(p, "Degree-days required must be > 0.")
  else if (r$threshold < 40) p <- c(p, paste0("Degree-days required = ", r$threshold, " looks like a temperature. ",
    "Enter the degree-day SUM to complete development (OPRAM species: 61-1131 DD), ",
    "not the upper temperature threshold."))
  if (!is.na(r$diapause_daylength) && (r$diapause_daylength < 0 || r$diapause_daylength > 24))
    p <- c(p, "Diapause daylength must be 0-24 h.")
  if (!is.na(r$diapause_temperature) && is.na(r$diapause_daylength))
    p <- c(p, "A diapause temperature needs a diapause daylength.")
  p
}

# append a user species (validated); returns problems (empty on success)
sp_add <- function(root, r) {
  all <- sp_read_all(root); p <- sp_validate(r, all); if (length(p)) return(p)
  u <- sp_read(sp_file_user(root)); r <- as.data.frame(r[sp_cols], check.names = FALSE)
  write.csv(rbind(u, r), sp_file_user(root), row.names = FALSE, na = "NA")
  character()
}

# delete a user species (OPRAM presets are read-only); results on disk are kept
sp_delete <- function(root, name) {
  u <- sp_read(sp_file_user(root)); if (is.null(u)) return(invisible(FALSE))
  write.csv(u[u$species != name, , drop = FALSE], sp_file_user(root), row.names = FALSE, na = "NA")
  invisible(TRUE)
}
