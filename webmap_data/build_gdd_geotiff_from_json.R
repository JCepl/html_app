#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(jsonlite)
  library(terra)
})

# Build georeferenced GeoTIFF from point JSON.
# Usage:
#   Rscript webmap_data/build_gdd_geotiff_from_json.R \
#     webmap_data/GDD_GEN_1km_preview.json \
#     webmap_data/gdd_gen_diff_1km_wgs84.tif

args <- commandArgs(trailingOnly = TRUE)
input_json <- if (length(args) >= 1) args[[1]] else "webmap_data/GDD_GEN_1km_preview.json"
output_tif <- if (length(args) >= 2) args[[2]] else "webmap_data/gdd_gen_diff_1km_wgs84.tif"

if (!file.exists(input_json)) {
  stop("Input JSON not found: ", input_json)
}

x <- fromJSON(input_json)
required_cols <- c("lng", "lat", "GEN_diff")
missing_cols <- setdiff(required_cols, names(x))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

x <- x[is.finite(x$lng) & is.finite(x$lat) & is.finite(x$GEN_diff), c("lng", "lat", "GEN_diff")]
if (nrow(x) == 0) {
  stop("No valid rows in input JSON")
}

v <- vect(x, geom = c("lng", "lat"), crs = "EPSG:4326")
v3035 <- project(v, "EPSG:3035")

e <- ext(v3035)
# Build an aligned 1km analysis grid in projected meters.
e_aligned <- ext(
  floor(xmin(e) / 1000) * 1000,
  ceiling(xmax(e) / 1000) * 1000,
  floor(ymin(e) / 1000) * 1000,
  ceiling(ymax(e) / 1000) * 1000
)

r_template <- rast(e_aligned, resolution = 1000, crs = "EPSG:3035")
r_3035 <- rasterize(v3035, r_template, field = "GEN_diff", fun = "mean")

# Reproject for web mapping compatibility.
r_4326 <- project(r_3035, "EPSG:4326", method = "near")

writeRaster(
  r_4326,
  output_tif,
  overwrite = TRUE,
  wopt = list(gdal = c("COMPRESS=LZW"))
)

cat("Built GeoTIFF:", output_tif, "\n")
cat("Rows x Cols:", nrow(r_4326), "x", ncol(r_4326), "\n")
cat("Extent:", as.character(ext(r_4326)), "\n")
