# species_db.R -- species parameter table: OPRAM presets (read-only) + user-added species.
# Columns of OPRAM data/species_parameters.csv plus three optional model columns:
#   method          degree-day method: sine (default) / triangle / average
#   upper_threshold upper development threshold, horizontal cutoff (°C; empty = none)
#   emergence_only  TRUE if the degree-day sum covers only part of the life cycle
#                   (e.g. overwintered larvae -> adult): emergence dates valid, generation counts not
# Part of DD-maps Norway (R port of OPRAM), GPL-3.

sp_cols <- c("species", "baseline_temperature", "threshold", "diapause", "diapause_daylength",
             "diapause_temperature", "Reference1", "Reference2", "Reference3", "event", "commonName",
             "method", "upper_threshold", "emergence_only")
sp_file_opram <- function(root) file.path(root, "3_metadata/species_parameters_OPRAM.csv")
sp_file_user  <- function(root) file.path(root, "3_metadata/species_parameters_user.csv")

sp_read <- function(f) {
  if (!file.exists(f)) return(NULL)
  s <- read.csv(f, check.names = FALSE, na.strings = c("NA", ""), colClasses = "character")
  for (k in setdiff(sp_cols, names(s))) s[[k]] <- NA_character_
  for (k in c("baseline_temperature", "threshold", "diapause_daylength", "diapause_temperature", "upper_threshold"))
    s[[k]] <- as.numeric(s[[k]])
  s$method[is.na(s$method)] <- "sine"
  s$emergence_only <- toupper(s$emergence_only) %in% c("TRUE", "T", "1", "YES")
  s[, sp_cols, drop = FALSE]
}

# all species with a source column ("OPRAM" / "user")
sp_read_all <- function(root) {
  o <- sp_read(sp_file_opram(root)); u <- sp_read(sp_file_user(root))
  rbind(cbind(o, source = "OPRAM"), if (!is.null(u) && nrow(u)) cbind(u, source = "user"))
}

# list of opram_params, named by species (method / upper / emergence_only per row)
sp_params <- function(df, maxYears = 4L)
  setNames(lapply(seq_len(nrow(df)), function(i) with(df[i, ], opram_params(
    baseline_temperature, threshold, diapause_daylength, diapause_temperature, method, maxYears, species,
    upper = upper_threshold, emergence_only = emergence_only))), df$species)

sp_slug <- function(x) gsub("[^a-z0-9]+", "_", tolower(trimws(x)))

# validate one row against the existing table; returns problems (empty = ok).
# replace = name of the row being edited (excluded from the duplicate check).
sp_validate <- function(r, all, replace = NULL) {
  p <- character(); others <- setdiff(all$species, replace)
  if (!nzchar(trimws(r$species))) p <- c(p, "Species name is required.")
  if (sp_slug(r$species) %in% sp_slug(others)) p <- c(p, "A species with this name already exists.")
  if (!is.finite(r$baseline_temperature) || r$baseline_temperature < -10 || r$baseline_temperature > 40)
    p <- c(p, "T0 must be between -10 and 40 °C.")
  if (!is.finite(r$threshold) || r$threshold <= 0 || r$threshold > 20000) p <- c(p, "Degree-days required must be > 0.")
  else if (r$threshold < 40) p <- c(p, paste0("Degree-days required = ", r$threshold, " looks like a temperature. ",
    "Enter the degree-day SUM to complete development (OPRAM species: 61-1131 DD); ",
    "an upper temperature threshold goes in its own field."))
  if (!is.null(r$upper_threshold) && !is.na(r$upper_threshold) && is.finite(r$baseline_temperature) &&
      r$upper_threshold <= r$baseline_temperature)
    p <- c(p, "The upper threshold must be above T0.")
  if (!is.null(r$method) && !is.na(r$method) && !r$method %in% c("sine", "triangle", "average"))
    p <- c(p, "Method must be sine, triangle or average.")
  if (!is.na(r$diapause_daylength) && (r$diapause_daylength < 0 || r$diapause_daylength > 24))
    p <- c(p, "Diapause daylength must be 0-24 h.")
  if (!is.na(r$diapause_temperature) && is.na(r$diapause_daylength))
    p <- c(p, "A diapause temperature needs a diapause daylength.")
  p
}

as_row <- function(r) {
  for (k in setdiff(sp_cols, names(r))) r[[k]] <- NA
  if (is.na(r$method)) r$method <- "sine"
  r$emergence_only <- isTRUE(as.logical(r$emergence_only))
  as.data.frame(r[sp_cols], check.names = FALSE)
}

# append a user species (validated); returns problems (empty on success)
sp_add <- function(root, r) {
  all <- sp_read_all(root); p <- sp_validate(r, all); if (length(p)) return(p)
  u <- sp_read(sp_file_user(root))
  write.csv(rbind(u, as_row(r)), sp_file_user(root), row.names = FALSE, na = "NA")
  character()
}

# replace an existing user species row (same or new name); returns problems
sp_update <- function(root, name, r) {
  u <- sp_read(sp_file_user(root)); if (is.null(u) || !name %in% u$species) return("Not a user species.")
  p <- sp_validate(r, sp_read_all(root), replace = name); if (length(p)) return(p)
  u[u$species == name, ] <- as_row(r)
  write.csv(u, sp_file_user(root), row.names = FALSE, na = "NA"); character()
}

# delete a user species (OPRAM presets are read-only); maps are kept unless delete_maps = TRUE
sp_delete <- function(root, name, delete_maps = FALSE) {
  u <- sp_read(sp_file_user(root)); if (is.null(u)) return(invisible(FALSE))
  write.csv(u[u$species != name, , drop = FALSE], sp_file_user(root), row.names = FALSE, na = "NA")
  if (delete_maps) unlink(file.path(root, "5_outputs/opram", sp_slug(name)), recursive = TRUE)
  invisible(TRUE)
}
