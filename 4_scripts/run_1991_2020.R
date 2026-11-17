# run_1991_2020.R -- full batch: extract 1991-2021, run all OPRAM species for 1991-2020, 30-yr normals.
# Resumable (finished files are skipped). Run from project root, e.g. in the background on Windows:
#   Start-Process Rscript -ArgumentList "4_scripts/run_1991_2020.R" -WindowStyle Hidden `
#     -RedirectStandardOutput 5_outputs/logs/run.out -RedirectStandardError 5_outputs/logs/run.err
source("4_scripts/run_norway.R")
dir.create("5_outputs/logs", FALSE, TRUE)
logf <- "5_outputs/logs/run_1991_2020.log"
say <- function(...) { m <- sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...)); message(m); cat(m, "\n", file = logf, append = TRUE) }

years <- 1991:2020; met <- 1991:last_meteo
say("start")

t0 <- Sys.time()
qc <- sn_extract(met, workers = 4, tif = FALSE)             # ~6 GB RAM per worker
if (!is.null(qc)) write.csv(qc, "5_outputs/logs/senorge_tx_lt_tn_qc.csv", row.names = FALSE)
say("extraction done: ", round(difftime(Sys.time(), t0, units = "mins"), 1), " min")

sp <- opram_species("3_metadata/species_parameters_OPRAM.csv")
t1 <- Sys.time()
opram_norway(years, sp, workers = 14, block_cells = 5000)   # ~1.5 GB RAM per worker
say("model done: ", round(difftime(Sys.time(), t1, units = "mins"), 1), " min")

t2 <- Sys.time()
for (s in sp) opram_normals(s$name, years, s$maxYears)
say("normals done: ", round(difftime(Sys.time(), t2, units = "mins"), 1), " min; total ",
    round(difftime(Sys.time(), t0, units = "mins"), 1), " min")
