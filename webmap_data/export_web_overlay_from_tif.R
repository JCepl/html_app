#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(terra)
  library(jsonlite)
})

# Convert a georeferenced WGS84 GeoTIFF to web overlay assets:
# 1) RGBA PNG with transparent NA
# 2) metadata JSON with exact bounds for Leaflet imageOverlay
#
# Usage:
#   Rscript webmap_data/export_web_overlay_from_tif.R \
#     webmap_data/gdd_gen_diff_1km_wgs84.tif \
#     webmap_data/gdd_gen_diff_raster.png \
#     webmap_data/gdd_gen_diff_raster_meta.json

args <- commandArgs(trailingOnly = TRUE)
input_tif <- if (length(args) >= 1) args[[1]] else "webmap_data/gdd_gen_diff_1km_wgs84.tif"
output_png <- if (length(args) >= 2) args[[2]] else "webmap_data/gdd_gen_diff_raster.png"
output_meta <- if (length(args) >= 3) args[[3]] else "webmap_data/gdd_gen_diff_raster_meta.json"

if (!file.exists(input_tif)) {
  stop("Input tif not found: ", input_tif)
}

r <- rast(input_tif)
crs_text <- crs(r)
if (is.na(crs_text) || (!grepl("WGS 84", crs_text, fixed = TRUE) && !grepl("4326", crs_text, fixed = TRUE))) {
  stop("Input raster must be WGS84/EPSG:4326.")
}

vals <- values(r, mat = FALSE)
nx <- ncol(r)
ny <- nrow(r)

if (length(vals) != nx * ny) {
  stop("Unexpected raster value length.")
}

# Classification bins (must match app legend/click classes)
# <= -0.5: blue, (-0.5,0.5]: light, (0.5,1.5]: orange, > 1.5: red
rgba <- matrix("#00000000", nrow = ny, ncol = nx)

m <- matrix(vals, nrow = ny, ncol = nx, byrow = TRUE)
rgba[is.finite(m) & m <= -0.5] <- "#1565c0ff"
rgba[is.finite(m) & m > -0.5 & m <= 0.5] <- "#f0f4f8ff"
rgba[is.finite(m) & m > 0.5 & m <= 1.5] <- "#ff9800ff"
rgba[is.finite(m) & m > 1.5] <- "#c62828ff"

# Render exact cell grid. No interpolation.
png(output_png, width = nx, height = ny, bg = "transparent")
par(mar = c(0, 0, 0, 0), xaxs = "i", yaxs = "i")
plot.new()
rasterImage(as.raster(rgba), 0, 0, 1, 1, interpolate = FALSE)
dev.off()

e <- ext(r)
meta <- list(
  image = basename(output_png),
  bounds = list(
    south = unname(ymin(e)),
    west = unname(xmin(e)),
    north = unname(ymax(e)),
    east = unname(xmax(e))
  ),
  dimensions = list(
    width = nx,
    height = ny
  ),
  source = basename(input_tif),
  legend = list(
    labels = c("<= -0.5", "-0.5 to 0.5", "0.5 to 1.5", "> 1.5"),
    colors = c("#1565c0", "#f0f4f8", "#ff9800", "#c62828")
  )
)

write_json(meta, output_meta, pretty = TRUE, auto_unbox = TRUE)

cat("Created:", output_png, "\n")
cat("Created:", output_meta, "\n")
cat("Bounds:", xmin(e), ymin(e), xmax(e), ymax(e), "\n")
cat("Dimensions:", nx, "x", ny, "\n")
