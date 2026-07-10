#!/usr/bin/env Rscript

# Re-render all forecast `probability_ensemble_mean.tif` files as CONTINUOUS
# sRGB overlays (instead of 10 discrete classes) so the on-map colors match the
# smooth legend gradient exactly. Probability is mapped 0..1 over the shared
# RISK_RAMP, identical interpolation to the CSS legend and the JS ramp.
#
# Usage: Rscript webmap_data/colorize_prob_continuous.R

suppressPackageStartupMessages({
  library(terra)
  library(jsonlite)
})

root <- "webmap_data/forecast_bundle"
max_dim <- 768
vmin <- 0
vmax <- 1
fade_knee <- 0.3   # low-risk end fades in across the bottom fraction
fade_floor <- 90   # minimum alpha at zero risk (stays faintly visible)

ramp_colors <- c("#14155f", "#2a45c2", "#3f7af0", "#6f74ee", "#9a5ee0",
                 "#c94fc0", "#ec4a86", "#ff7a1f", "#f5331a", "#8f0000")
color_ramp  <- colorRamp(ramp_colors, space = "rgb")

tifs <- list.files(root, pattern = "probability_ensemble_mean\\.tif$",
                   recursive = TRUE, full.names = TRUE)
cat("Found", length(tifs), "probability rasters\n")

for (in_tif in tifs) {
  r <- rast(in_tif)

  fact <- ceiling(max(ncol(r), nrow(r)) / max_dim)
  if (fact > 1) r <- aggregate(r, fact = fact, fun = "mean", na.rm = TRUE)
  r <- project(r, "EPSG:3857", method = "bilinear")
  merc_ext <- ext(r)

  nx <- ncol(r); ny <- nrow(r)
  m  <- matrix(values(r, mat = FALSE), nrow = ny, ncol = nx, byrow = TRUE)
  finite_mask <- is.finite(m)

  norm <- pmin(pmax(m, vmin), vmax) / (vmax - vmin)
  rgba <- matrix("#00000000", nrow = ny, ncol = nx)
  if (any(finite_mask)) {
    vv  <- norm[finite_mask]
    rgb <- color_ramp(vv)
    # Low-risk end fades toward (but not to) transparent so ~1% risk still
    # shows faintly (see wind export).
    alpha <- round(fade_floor + (255 - fade_floor) * pmin(1, vv / fade_knee))
    rgba[finite_mask] <- sprintf("#%02X%02X%02X%02X",
                                 round(rgb[, 1]), round(rgb[, 2]), round(rgb[, 3]), alpha)
  }

  out_png <- sub("\\.tif$", ".png", in_tif)
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
    source = basename(in_tif),
    value_range = list(min = vmin, max = vmax)
  )
  write_json(meta, sub("\\.tif$", "_meta.json", in_tif),
             pretty = TRUE, auto_unbox = TRUE, digits = 8)
}
cat("Re-rendered", length(tifs), "probability PNGs (continuous sRGB)\n")
