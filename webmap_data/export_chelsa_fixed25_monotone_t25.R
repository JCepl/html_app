#!/usr/bin/env Rscript
# Adds a risk_visible_t25 raster to the existing chelsa_fixed25_monotone
# forecast entries -- the model's OWN CV-calibrated operating threshold
# (0.25, see metrics_confusion.json / export_chelsa_fixed25_monotone.R's
# comments) never had a matching classified raster exported; only
# risk_visible_t50 exists, which uses BARK_BEETLE_CLASSIFIED_BREAKS' own
# floor (0.55) as its alpha cutoff -- well above 0.25, so a user could never
# actually SEE the map at the threshold the model was validated at, only
# read about it in the metrics note.
#
# Builds directly from the already-exported probability_ensemble_mean.tif
# (already reprojected/downsampled) rather than re-reading the large
# CHELSA-native source rasters -- pixel-identical alignment, much cheaper.
#
# Colour scheme: reuses BARK_BEETLE_CLASSIFIED_COLORS unchanged for the
# 0.55+ bins (visual consistency with every other model's classified view),
# with ONE new pale band prepended for 0.25-0.55 -- this is intentionally
# scoped to this script, NOT added to the shared BARK_BEETLE_CLASSIFIED_*
# constants in color_scales.R, since no other model asked for a sub-0.55 bin
# and changing the shared constant would silently alter every other model's
# t50+ raster too.
#
# Usage: Rscript webmap_data/export_chelsa_fixed25_monotone_t25.R

suppressPackageStartupMessages({
  library(terra)
  library(jsonlite)
})
source("webmap_data/color_scales.R")

out_root <- "webmap_data/forecast_bundle/chelsa_fixed25_monotone"
periods <- c("1983_2010", "2071_2100")

# BARK_BEETLE_CLASSIFIED_COLORS with one new low band prepended for 0.25-0.55
# (a pale tint preceding the palette's existing lowest real color, #5e4fa2).
risk_breaks <- c(0.25, BARK_BEETLE_CLASSIFIED_BREAKS)
risk_colors <- c("#00000000", "#ccc7ec", BARK_BEETLE_CLASSIFIED_COLORS[-1])
risk_colors_rgb <- substr(risk_colors, 1, 7)

classify_color <- function(vv) {
  idx <- findInterval(vv, risk_breaks) + 1
  t(col2rgb(risk_colors_rgb[idx]))
}
classify_alpha <- function(vv) ifelse(vv < risk_breaks[1], 0, 255)

to_png <- function(r, out_png, colorfun, alpha_fun) {
  merc_ext <- ext(r)
  nx <- ncol(r); ny <- nrow(r)
  m <- matrix(values(r, mat = FALSE), nrow = ny, ncol = nx, byrow = TRUE)
  finite_mask <- is.finite(m)
  rgba <- matrix("#00000000", nrow = ny, ncol = nx)
  if (any(finite_mask)) {
    vv <- m[finite_mask]
    rgb <- colorfun(vv)
    alpha <- alpha_fun(vv)
    rgba[finite_mask] <- sprintf("#%02X%02X%02X%02X", round(rgb[,1]), round(rgb[,2]), round(rgb[,3]), alpha)
  }
  png(out_png, width = nx, height = ny, bg = "transparent")
  par(mar = c(0, 0, 0, 0), xaxs = "i", yaxs = "i")
  plot.new()
  rasterImage(as.raster(rgba), 0, 0, 1, 1, interpolate = FALSE)
  dev.off()
  merc_ext
}

bounds_from_ext <- function(merc_ext) {
  corner_xy <- rbind(c(xmin(merc_ext), ymin(merc_ext)), c(xmin(merc_ext), ymax(merc_ext)),
                      c(xmax(merc_ext), ymin(merc_ext)), c(xmax(merc_ext), ymax(merc_ext)))
  ll <- project(corner_xy, from = "EPSG:3857", to = "EPSG:4326")
  list(south = min(ll[,2]), west = min(ll[,1]), north = max(ll[,2]), east = max(ll[,1]))
}

for (key in periods) {
  out_dir <- file.path(out_root, key)
  r <- rast(file.path(out_dir, "probability_ensemble_mean.tif"))
  cat("Building risk_visible_t25 for", key, "...\n")

  risk_tif <- file.path(out_dir, "risk_visible_t25.tif")
  writeRaster(r, risk_tif, overwrite = TRUE, datatype = "FLT4S", gdal = c("COMPRESS=LZW"))
  merc_ext <- to_png(r, file.path(out_dir, "risk_visible_t25.png"),
    colorfun = classify_color, alpha_fun = classify_alpha)
  meta <- list(image = "risk_visible_t25.png", bounds = bounds_from_ext(merc_ext),
               dimensions = list(width = ncol(r), height = nrow(r)),
               source = "CHELSA-native fixed25/monotone-XGB (see PROJECT_REFERENCE_bark_beetle_model.md)",
               value_range = list(min = 0.25, max = 1),
               note = "Classified at the model's own CV-calibrated operating threshold (0.25), not the generic 0.55 floor risk_visible_t50 uses.")
  write_json(meta, file.path(out_dir, "risk_visible_t25_meta.json"), pretty = TRUE, auto_unbox = TRUE, digits = 8)
  cat("Wrote", risk_tif, "\n")
}
cat("Done.\n")
