# run_v2309.R -- migrate to seNorge v23.09 (yearly THREDDS files) and recompute everything.
#  1. re-extract 1991-2025 TN/TX into the .bin store from 1_raw_data/senorge_archive (all v23.09)
#  2. move the raw yearly files to ARCHIVE_DEST (OneDrive) once every .bin is verified
#  3. recompute all species (OPRAM presets + user species) for 1991-2022 (fully observed 4-yr windows)
#  4. 1991-2020 normals + probability maps
# Resumable per step. Run detached from project root:
#   Start-Process Rscript -ArgumentList "4_scripts/run_v2309.R" -WindowStyle Hidden `
#     -RedirectStandardOutput 5_outputs/logs/v2309.out -RedirectStandardError 5_outputs/logs/v2309.err
source("4_scripts/run_norway.R"); source("4_scripts/species_db.R")
dir.create("5_outputs/logs", FALSE, TRUE); logf <- "5_outputs/logs/run_v2309.log"
say <- function(...) { m <- sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...)); message(m); cat(m, "\n", file = logf, append = TRUE) }
src <- sn_archive_dir
dest <- Sys.getenv("ARCHIVE_DEST", "")                                   # set in .Renviron
met <- 1991:2025; t0 <- Sys.time(); say("start (archive: ", src, ")")

# 1. extraction (skip years whose .bin is already from v23.09 -> marker file)
mark <- file.path(sn_out, "bin", "v2309_done.txt")
done <- if (file.exists(mark)) as.integer(readLines(mark)) else integer()
todo <- setdiff(met, done)
if (length(todo)) {
  stopifnot(all(file.exists(file.path(src, sprintf("seNorge2018_%d.nc", todo)))))
  for (ch in split(todo, ceiling(seq_along(todo) / 4))) {                # 4 years at a time (~8 GB RAM each)
    qc <- sn_extract(ch, workers = 4, tif = FALSE, overwrite = TRUE)
    cat(ch, file = mark, sep = "\n", append = TRUE)
    say("extracted ", paste(ch, collapse = ","), if (!is.null(qc)) paste0("; TX<TN cell-days: ", sum(qc$tx_lt_tn)))
  }
}
say("extraction done: ", round(difftime(Sys.time(), t0, units = "mins"), 1), " min")

# 2. verify .bin store, then move raw files to OneDrive
gr <- readRDS(file.path(sn_out, "grid.rds"))
ok <- vapply(met, function(y) all(file.size(sn_bin_file(c("tn", "tx"), y)) == 2 * ndays(y) * length(gr$land)), TRUE)
if (!all(ok)) stop("bin store incomplete for: ", paste(met[!ok], collapse = ","))
if (nzchar(dest)) {
  dir.create(dest, FALSE, TRUE)
  for (f in list.files(src, "^seNorge2018_\\d{4}\\.nc$", full.names = TRUE)) {
    to <- file.path(dest, basename(f))
    if (!file.rename(f, to)) {                                            # same volume: instant; else copy+check
      if (!file.copy(f, to, overwrite = TRUE) || file.size(to) != file.size(f)) stop("move failed: ", f)
      unlink(f)
    }
  }
  say("raw files moved to ", dest)
}

# 3. recompute all species 1991-2022
sp <- sp_params(sp_read_all("."))
t1 <- Sys.time()
opram_norway(1991:2022, sp, workers = 14, block_cells = 5000, overwrite = TRUE)
say("model done (", length(sp), " sets x 1991-2022): ", round(difftime(Sys.time(), t1, units = "mins"), 1), " min")

# 4. normals + probabilities (parallel over species)
t2 <- Sys.time(); plan(multisession, workers = 6)
r <- future_lapply(sp, function(p) tryCatch({ opram_normals(p$name, 1991:2020, p$maxYears)
  opram_probability(p$name, 1991:2020); "ok" }, error = function(e) conditionMessage(e)), future.packages = "terra")
say("normals + probabilities: ", paste(names(table(unlist(r))), table(unlist(r)), collapse = "; "), " in ",
    round(difftime(Sys.time(), t2, units = "mins"), 1), " min; total ", round(difftime(Sys.time(), t0, units = "mins"), 1), " min")
