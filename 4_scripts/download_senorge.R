# download_senorge.R -- download yearly seNorge_2018 files (tn, tx, tg, rr; v23.09) from MET Norway THREDDS.
# Data: MET Norway, CC BY 4.0 / NLOD (credit "MET Norway"). Part of DD-maps Norway, GPL-3.
#
# Usage (from project root):
#   Rscript 4_scripts/download_senorge.R              # default years 1991-2025 (~75 GB, ~1 h)
#   Rscript 4_scripts/download_senorge.R 2018 2025    # from-to
# Destination: SENORGE_ARCHIVE_DIR in .Renviron (git-ignored), else 1_raw_data/senorge_archive.
#
# - ~2.2-2.4 GB per year, ~20 MB/s -> ~2 min per year
# - resumable: complete files are skipped; partial downloads (.part) are continued (curl) or restarted
# - each file is checked: size == server Content-Length, opens with ncdf4, has tn+tx and 365/366 days

base <- "https://thredds.met.no/thredds/fileServer/senorge/seNorge_2018/Archive"
dest <- Sys.getenv("SENORGE_ARCHIVE_DIR", "1_raw_data/senorge_archive")
a <- as.integer(commandArgs(TRUE))
years <- if (length(a) == 2) a[1]:a[2] else if (length(a) == 1) a else 1991:2025   # default: full period
dir.create(dest, FALSE, TRUE)
curl <- Sys.which("curl")                                   # curl.exe ships with Windows 10/11
options(timeout = 7200)

remote_size <- function(url) {
  h <- curlGetHeaders(url, redirect = TRUE)
  as.numeric(sub(".*:\\s*", "", grep("^content-length", h, ignore.case = TRUE, value = TRUE) |> tail(1)))
}

check_nc <- function(f, y) {
  nc <- ncdf4::nc_open(f); on.exit(ncdf4::nc_close(nc))
  nd <- nc$dim$time$len; ok <- all(c("tn", "tx") %in% names(nc$var)) && nd == 365 + (y %% 4 == 0)
  if (!ok) stop(sprintf("%s: unexpected content (vars: %s; days: %d)", basename(f), paste(names(nc$var), collapse = ","), nd))
  sprintf("tn+tx, %d days, version %s", nd, ncdf4::ncatt_get(nc, 0, "version")$value)
}

for (y in years) {
  url <- sprintf("%s/seNorge2018_%d.nc", base, y)
  f <- file.path(dest, sprintf("seNorge2018_%d.nc", y)); part <- paste0(f, ".part")
  size <- remote_size(url)
  if (file.exists(f) && file.size(f) == size) { message(y, ": already complete, skipped"); next }
  message(sprintf("%d: downloading %.2f GB ...", y, size / 1e9)); t0 <- Sys.time()
  if (nzchar(curl)) {
    st <- system2(curl, c("-L", "-s", "-S", "--retry", "5", "--retry-delay", "10", "-C", "-", "-o", shQuote(part), url))
    if (st != 0) stop(y, ": curl failed (exit ", st, "); rerun the script to resume")
  } else {
    unlink(part); download.file(url, part, mode = "wb", quiet = TRUE)   # mode = "wb": required on Windows
  }
  if (file.size(part) != size) stop(sprintf("%d: size mismatch (%.0f vs %.0f bytes); rerun to resume", y, file.size(part), size))
  file.rename(part, f)
  message(sprintf("%d: done in %.1f min (%.0f MB/s) - %s", y, as.numeric(difftime(Sys.time(), t0, units = "mins")),
                  size / 1e6 / as.numeric(difftime(Sys.time(), t0, units = "secs")), check_nc(f, y)))
}
message("All requested years in ", normalizePath(dest))
