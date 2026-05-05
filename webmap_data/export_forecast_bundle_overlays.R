#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(jsonlite)
})

root <- "webmap_data/forecast_bundle"
exporter <- "webmap_data/export_web_overlay_from_tif.R"

if (!file.exists(root)) {
  stop("Forecast bundle folder not found: ", root)
}
if (!file.exists(exporter)) {
  stop("Exporter script not found: ", exporter)
}

prob_breaks <- "0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9"
prob_colors <- "#fff7ec,#fee8c8,#fdd49e,#fdbb84,#fc8d59,#ef6548,#d7301f,#b30000,#7f0000,#4d0000"
prob_labels <- "<= 0.10,0.10-0.20,0.20-0.30,0.30-0.40,0.40-0.50,0.50-0.60,0.60-0.70,0.70-0.80,0.80-0.90,> 0.90"

risk_breaks <- "0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9"
risk_colors <- "#fee8c800,#fdd49e,#fdbb84,#fc8d59,#ef6548,#d7301f,#b30000,#7f0000,#4d0000"
risk_labels <- "<= 0.55,0.55-0.60,0.60-0.65,0.65-0.70,0.70-0.75,0.75-0.80,0.80-0.85,0.85-0.90,> 0.90"

run_export <- function(input_tif, breaks_arg, colors_arg, labels_arg, max_dim = "768") {
  output_png <- sub("\\.tif$", ".png", input_tif, ignore.case = TRUE)
  output_meta <- sub("\\.tif$", "_meta.json", input_tif, ignore.case = TRUE)

  cmd <- paste(
    c(
      "Rscript",
      shQuote(exporter),
      shQuote(input_tif),
      shQuote(output_png),
      shQuote(output_meta),
      shQuote(breaks_arg),
      shQuote(colors_arg),
      shQuote(labels_arg),
      shQuote(max_dim)
    ),
    collapse = " "
  )

  status <- system(cmd)
  if (!identical(status, 0L)) {
    stop("Export failed for ", input_tif)
  }

  meta <- fromJSON(output_meta, simplifyVector = TRUE)
  meta$value_range <- list(
    min = if (grepl("probability_ensemble_mean\\.tif$", input_tif)) 0 else NA_real_,
    max = 1
  )
  write_json(meta, output_meta, pretty = TRUE, auto_unbox = TRUE)

  cat("Prepared:", output_png, "\n")
}

tif_files <- list.files(root, pattern = "\\.tif$", recursive = TRUE, full.names = TRUE)
if (!length(tif_files)) {
  stop("No tif files found in forecast bundle.")
}

for (input_tif in tif_files) {
  if (grepl("probability_ensemble_mean\\.tif$", input_tif)) {
    run_export(input_tif, prob_breaks, prob_colors, prob_labels)
  } else if (grepl("risk_visible_t[0-9]+\\.tif$", input_tif)) {
    run_export(input_tif, risk_breaks, risk_colors, risk_labels)
  }
}
