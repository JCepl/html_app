#!/usr/bin/env Rscript

# Builds the forecast-bundle web overlays (PNG + *_meta.json) from the TIFFs.
#   - risk_visible_t*  -> thresholded "risk highlight" layers (STILL USED).
#   - probability_ensemble_mean -> classified probability tiles. SUPERSEDED:
#     the app's probability view is now rendered CONTINUOUSLY by
#     colorize_prob_continuous.R, so re-run THAT (not this) for probability.
# The RdYlBu "bold" ramp below must stay in sync with RISK_RAMP in index.html
# and ramp_colors in colorize_prob_continuous.R / export_wind_overlays.R.
# Usage: Rscript webmap_data/export_forecast_bundle_overlays.R

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

# Intuitive risk palette RdYlBu reversed (kept in sync with RISK_RAMP in index.html):
# blue = low probability -> pale yellow = mid -> red = high probability.
prob_breaks <- "0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9"
prob_colors <- "#14155f,#2a45c2,#3f7af0,#6f74ee,#9a5ee0,#c94fc0,#ec4a86,#ff7a1f,#f5331a,#8f0000"
prob_labels <- "<= 0.10,0.10-0.20,0.20-0.30,0.30-0.40,0.40-0.50,0.50-0.60,0.60-0.70,0.70-0.80,0.80-0.90,> 0.90"

risk_breaks <- "0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9"
risk_colors <- "#00000000,#fef08a,#fecc5c,#fd8d3c,#f03b20,#bd0026,#800026,#4d0000,#1a0000"
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
