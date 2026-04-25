#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(terra)
  library(jsonlite)
})

# Convert a georeferenced WGS84 GeoTIFF to web overlay assets:
# 1) RGBA PNG with transparent NA
# 2) metadata JSON with exact bounds for Leaflet imageOverlay and class breaks
#
# Usage:
#   Rscript webmap_data/export_web_overlay_from_tif.R \
#     input.tif \
#     output.png \
#     output_meta.json \
#     0,0.5,1.5,3,5 \
#     "#f7fbff,#deebf7,#c6dbef,#9ecae1,#6baed6,#2171b5" \
#     "0-0.5,0.5-1.5,1.5-3,3-5,5-7,>7"

args <- commandArgs(trailingOnly = TRUE)
input_tif <- if (length(args) >= 1) args[[1]] else "webmap_data/gdd_gen_diff_1km_wgs84.tif"
output_png <- if (length(args) >= 2) args[[2]] else "webmap_data/gdd_gen_diff_raster.png"
output_meta <- if (length(args) >= 3) args[[3]] else "webmap_data/gdd_gen_diff_raster_meta.json"
breaks_arg <- if (length(args) >= 4) args[[4]] else "-0.5,0.5,1.5"
colors_arg <- if (length(args) >= 5) args[[5]] else "#1565c0,#f0f4f8,#ff9800,#c62828"
labels_arg <- if (length(args) >= 6) args[[6]] else "<= -0.5,-0.5 to 0.5,0.5 to 1.5,> 1.5"
max_dim_arg <- if (length(args) >= 7) args[[7]] else NA_character_

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

if (!is.na(max_dim_arg) && nzchar(max_dim_arg)) {
  max_dim <- as.numeric(max_dim_arg)
  if (!is.finite(max_dim) || max_dim <= 0) {
    stop("max_dim must be a positive number when provided.")
  }
  scale_factor <- max(nx, ny) / max_dim
  if (scale_factor > 1) {
    fact <- ceiling(scale_factor)
    r <- aggregate(r, fact = fact, fun = "max", na.rm = TRUE)
    vals <- values(r, mat = FALSE)
    nx <- ncol(r)
    ny <- nrow(r)
  }
}

if (length(vals) != nx * ny) {
  stop("Unexpected raster value length.")
}

breaks <- as.numeric(trimws(strsplit(breaks_arg, ",", fixed = TRUE)[[1]]))
colors <- trimws(strsplit(colors_arg, ",", fixed = TRUE)[[1]])
labels <- trimws(strsplit(labels_arg, ",", fixed = TRUE)[[1]])

if (any(!is.finite(breaks))) {
  stop("All class breaks must be numeric.")
}
if (length(colors) != length(breaks) + 1) {
  stop("Number of colors must equal number of classes: breaks + 1.")
}
if (length(labels) != length(colors)) {
  stop("Number of labels must equal number of colors/classes.")
}

to_rgba <- function(x) {
  if (grepl("^#[0-9A-Fa-f]{8}$", x)) {
    return(x)
  }
  if (grepl("^#[0-9A-Fa-f]{6}$", x)) {
    return(paste0(x, "ff"))
  }
  stop("Colors must use #RRGGBB or #RRGGBBAA format.")
}

rgba <- matrix("#00000000", nrow = ny, ncol = nx)
m <- matrix(vals, nrow = ny, ncol = nx, byrow = TRUE)
finite_mask <- is.finite(m)

if (length(breaks) == 0) {
  rgba[finite_mask] <- to_rgba(colors[[1]])
} else {
  lower <- -Inf
  for (i in seq_along(colors)) {
    upper <- if (i <= length(breaks)) breaks[[i]] else Inf
    class_mask <- finite_mask & m > lower & m <= upper
    if (i == 1) {
      class_mask <- finite_mask & m <= upper
    }
    if (i == length(colors)) {
      class_mask <- finite_mask & m > lower
    }
    rgba[class_mask] <- to_rgba(colors[[i]])
    lower <- upper
  }
}

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
    breaks = breaks,
    labels = labels,
    colors = colors
  )
)

write_json(meta, output_meta, pretty = TRUE, auto_unbox = TRUE)

cat("Created:", output_png, "\n")
cat("Created:", output_meta, "\n")
cat("Bounds:", xmin(e), ymin(e), xmax(e), ymax(e), "\n")
cat("Dimensions:", nx, "x", ny, "\n")
