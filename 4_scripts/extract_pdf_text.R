# extract_pdf_text.R -- text from downloaded source papers (1_raw_data/dd_pdf) for parameter verification.
# Writes 2_processed_data/dd_pdf_text/<file>.txt (pages separated by "\f") and an index CSV with
# pages, characters, text layer yes/no and the first DOI found. Part of DD-maps Norway (FinnPRIO), GPL-3.
src <- "1_raw_data/dd_pdf"; out <- "2_processed_data/dd_pdf_text"; dir.create(out, FALSE, TRUE)
pdfs <- list.files(src, "\\.pdf$", full.names = TRUE, ignore.case = TRUE)
idx <- do.call(rbind, lapply(pdfs, function(f) {
  tx <- tryCatch(pdftools::pdf_text(f), error = function(e) character())
  info <- tryCatch(pdftools::pdf_info(f), error = function(e) list(pages = NA))
  all <- paste(tx, collapse = "\f")
  writeLines(all, file.path(out, sub("\\.pdf$", ".txt", basename(f), ignore.case = TRUE)), useBytes = TRUE)
  doi <- regmatches(all, regexpr("10\\.[0-9]{4,9}/[^\\s\"<>,;]+", all, perl = TRUE))
  data.frame(file = basename(f), pages = info$pages, chars = nchar(all),
             text_layer = nchar(gsub("\\s", "", all)) > 500 * max(1, info$pages %||% 1) / 4,
             doi = if (length(doi)) sub("[.)]+$", "", doi) else NA,
             first_line = substr(trimws(strsplit(trimws(tx[1] %||% ""), "\n")[[1]][1]), 1, 90))
}))
write.csv(idx, file.path(out, "_index.csv"), row.names = FALSE)
options(width = 220); print(idx[, c("file", "pages", "chars", "text_layer", "doi")], row.names = FALSE)
