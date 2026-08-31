#!/usr/bin/env Rscript

# Builds the forecast-bundle web overlays (PNG + *_meta.json) from the TIFFs.
#   - risk_visible_t*  -> thresholded "risk highlight" layers (STILL USED).
#   - probability_ensemble_mean -> classified probability tiles. SUPERSEDED:
#     the app's probability view is now rendered CONTINUOUSLY by
#     colorize_prob_continuous.R, so re-run THAT (not this) for probability.
# The ramp below is bark beetle's own scale (BARK_BEETLE_RAMP in
# color_scales.R -- see that file for how it relates to wind's WIND_RAMP).
# Usage: Rscript webmap_data/export_forecast_bundle_overlays.R

suppressPackageStartupMessages({
  library(jsonlite)
})
source("webmap_data/color_scales.R")

root <- "webmap_data/forecast_bundle"
exporter <- "webmap_data/export_web_overlay_from_tif.R"

if (!file.exists(root)) {
  stop("Forecast bundle folder not found: ", root)
}
if (!file.exists(exporter)) {
  stop("Exporter script not found: ", exporter)
}

prob_breaks <- "0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9"
# NOTE: this 10-stop resample of BARK_BEETLE_RAMP had already drifted from
# it (different intermediate hex values) before this refactor -- left as
# its pre-existing historical value rather than silently reharmonized,
# since this path is SUPERSEDED (colorize_prob_continuous.R is what
# actually runs now) and re-resampling dead code isn't worth the risk of
# a surprise visual change if it's ever revived. Fix explicitly if reviving.
prob_colors <- "#5e4fa2,#388eba,#75c8a5,#bfe5a0,#f1f9a9,#feeea2,#fdbf6f,#f67b4a,#d8434e,#9e0142"
prob_labels <- "<= 0.10,0.10-0.20,0.20-0.30,0.30-0.40,0.40-0.50,0.50-0.60,0.60-0.70,0.70-0.80,0.80-0.90,> 0.90"

risk_breaks <- ramp_to_csv(BARK_BEETLE_CLASSIFIED_BREAKS)
risk_colors <- ramp_to_csv(BARK_BEETLE_CLASSIFIED_COLORS)
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
