# prep_norway_mask.R -- Norway mask + outline on the seNorge 1 km grid (for maps / app).
# Boundary: GADM 4.1 (gadm.org; free for non-commercial use), NOR level 1 via geodata. Run from project root.
suppressPackageStartupMessages(library(terra))
proc <- "2_processed_data/senorge"
gr <- readRDS(file.path(proc, "grid.rds"))
tmpl <- rast(nrows = gr$nrow, ncols = gr$ncol, extent = ext(gr$ext), crs = gr$crs)

dir.create("1_raw_data/gadm", FALSE, TRUE)
nor <- geodata::gadm("NOR", level = 1, path = "1_raw_data/gadm")       # cached after first call
nor <- project(nor, gr$crs)

# land cells touched by Norway (keeps coastal/island cells whose centre lies in the sea)
inside <- !is.na(values(rasterize(nor, tmpl, touches = TRUE), mat = FALSE))
keep <- inside[gr$land]
e <- ext(tmpl, cells = gr$land[keep]); e <- ext(e[1] - 10000, e[2] + 10000, e[3] - 10000, e[4] + 10000)

# outline (counties + coast) as simplified path data frame for ggplot2
lines_df <- function(v, tol) {
  g <- geom(as.lines(simplifyGeom(v, tolerance = tol, preserveTopology = TRUE)))
  data.frame(x = g[, "x"], y = g[, "y"], id = paste(g[, "geom"], g[, "part"]))
}
saveRDS(list(keep = keep, ext = as.vector(e), border = lines_df(aggregate(nor), 300), county = lines_df(nor, 500)),
        file.path(proc, "norway_mask.rds"))
cat("Norway land cells:", sum(keep), "of", length(keep), " extent:", round(as.vector(e)), "\n")
