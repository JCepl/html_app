#!/usr/bin/env Rscript

# Build the missing RCP8.5 forecast ENSEMBLE from the available rcp85 model
# probability rasters (xgb, rf, gam) — a multi-model mean (each input is itself
# a GCM-mean probability). Mirrors the existing ensemble_rcp26/45 outputs:
#   probability_ensemble_mean.tif/.png/_meta.json
#   risk_visible_t{50..90}.tif/.png/_meta.json
#   entry.json
# Historical evaluation metrics are scenario-independent, so the manifest entry
# (added separately) reuses ensemble/metrics_confusion.json.
#
# Usage: Rscript webmap_data/build_ensemble_rcp85.R

suppressPackageStartupMessages({
  library(terra)
  library(jsonlite)
})
source("webmap_data/color_scales.R")

root        <- "webmap_data/forecast_bundle"
out_dir     <- file.path(root, "ensemble", "rcp85")
exporter    <- "webmap_data/export_web_overlay_from_tif.R"
models      <- c("xgb", "rf", "gam")   # xgb first -> reference grid (full extent)
periods     <- c("2011_2040", "2041_2070", "2071_2100")
cutoffs     <- c(0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9)
max_dim     <- 768
fade_knee   <- 0.3
fade_floor  <- 90

# Bark beetle's own scale -- see color_scales.R.
ramp_colors <- BARK_BEETLE_RAMP
color_ramp  <- colorRamp(ramp_colors, space = "rgb")

# Classified "risk highlight" scale for the thresholded view (matches
# export_forecast_bundle_overlays.R).
risk_breaks <- ramp_to_csv(BARK_BEETLE_CLASSIFIED_BREAKS)
risk_colors <- ramp_to_csv(BARK_BEETLE_CLASSIFIED_COLORS)
risk_labels <- "<= 0.55,0.55-0.60,0.60-0.65,0.65-0.70,0.70-0.75,0.75-0.80,0.80-0.85,0.85-0.90,> 0.90"

# Render a continuous probability raster to an RGBA web PNG + meta JSON,
# matching colorize_prob_continuous.R (sRGB ramp + low-risk alpha floor).
colorize_continuous <- function(r, out_png, out_meta, src_name) {
  fact <- ceiling(max(ncol(r), nrow(r)) / max_dim)
  if (fact > 1) r <- aggregate(r, fact = fact, fun = "mean", na.rm = TRUE)
  r <- project(r, "EPSG:3857", method = "bilinear")
  merc_ext <- ext(r)

  nx <- ncol(r); ny <- nrow(r)
  m  <- matrix(values(r, mat = FALSE), nrow = ny, ncol = nx, byrow = TRUE)
  finite_mask <- is.finite(m)
  norm <- pmin(pmax(m, 0), 1)

  rgba <- matrix("#00000000", nrow = ny, ncol = nx)
  if (any(finite_mask)) {
    vv  <- norm[finite_mask]
    rgb <- color_ramp(vv)
    alpha <- round(fade_floor + (255 - fade_floor) * pmin(1, vv / fade_knee))
    rgba[finite_mask] <- sprintf("#%02X%02X%02X%02X",
                                 round(rgb[, 1]), round(rgb[, 2]), round(rgb[, 3]), alpha)
  }

  png(out_png, width = nx, height = ny, bg = "transparent")
  par(mar = c(0, 0, 0, 0), xaxs = "i", yaxs = "i")
  plot.new()
  rasterImage(as.raster(rgba), 0, 0, 1, 1, interpolate = FALSE)
  dev.off()

  corner_xy <- rbind(
    c(xmin(merc_ext), ymin(merc_ext)), c(xmin(merc_ext), ymax(merc_ext)),
    c(xmax(merc_ext), ymin(merc_ext)), c(xmax(merc_ext), ymax(merc_ext))
  )
  ll <- project(corner_xy, from = "EPSG:3857", to = "EPSG:4326")
  meta <- list(
    image = basename(out_png),
    bounds = list(south = min(ll[, 2]), west = min(ll[, 1]),
                  north = max(ll[, 2]), east = max(ll[, 1])),
    dimensions = list(width = nx, height = ny),
    source = src_name,
    value_range = list(min = 0, max = 1)
  )
  write_json(meta, out_meta, pretty = TRUE, auto_unbox = TRUE, digits = 8)
}

for (per in periods) {
  paths <- character(0)
  for (mdl in models) {
    f <- file.path(root, mdl, "rcp85", per, "probability_ensemble_mean.tif")
    if (file.exists(f)) paths <- c(paths, f)
  }
  if (!length(paths)) { cat("skip", per, "(no models)\n"); next }

  ref <- rast(paths[1])
  layers <- list(ref)
  if (length(paths) > 1) {
    for (k in 2:length(paths)) {
      r <- rast(paths[k])
      if (!isTRUE(try(compareGeom(r, ref, stopOnError = FALSE), silent = TRUE))) {
        r <- resample(r, ref, method = "bilinear")
      }
      layers[[length(layers) + 1]] <- r
    }
  }
  ens <- mean(rast(layers), na.rm = TRUE)   # cell-wise multi-model mean

  pdir <- file.path(out_dir, per)
  dir.create(pdir, recursive = TRUE, showWarnings = FALSE)
  prob_tif <- file.path(pdir, "probability_ensemble_mean.tif")
  writeRaster(ens, prob_tif, overwrite = TRUE)
  colorize_continuous(ens, file.path(pdir, "probability_ensemble_mean.png"),
                      file.path(pdir, "probability_ensemble_mean_meta.json"),
                      "ensemble rcp85 multi-model mean")

  tvis <- list()
  for (ct in cutoffs) {
    tt  <- sprintf("t%02d", round(ct * 100))
    rv  <- ens
    rv[rv < ct] <- NA
    rv_tif <- file.path(pdir, paste0("risk_visible_", tt, ".tif"))
    writeRaster(rv, rv_tif, overwrite = TRUE)
    system(paste("Rscript", shQuote(exporter), shQuote(rv_tif),
                 shQuote(sub("\\.tif$", ".png", rv_tif)),
                 shQuote(sub("\\.tif$", "_meta.json", rv_tif)),
                 shQuote(risk_breaks), shQuote(risk_colors), shQuote(risk_labels),
                 shQuote("768")), ignore.stdout = TRUE)
    tvis[[as.character(ct)]] <- rv_tif
  }

  entry <- list(
    model = "ensemble", scenario = "rcp85", period = per,
    gcms = I(character(0)),
    source_probability_files = I(unname(paths)),
    ensemble_probability_tif = prob_tif,
    threshold_visible_tifs = tvis,
    forecast_confusion_by_gcm = I(list())
  )
  write_json(entry, file.path(pdir, "entry.json"), pretty = TRUE, auto_unbox = TRUE)
  cat("built ensemble rcp85", per, "from", length(paths), "models\n")
}
cat("done\n")
