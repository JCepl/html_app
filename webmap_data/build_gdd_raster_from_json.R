#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(jsonlite)
})

# Build a georeferenced PNG + bounds JSON from point JSON.
# Usage:
#   Rscript webmap_data/build_gdd_raster_from_json.R \
#     webmap_data/GDD_GEN_1km_preview.json \
#     webmap_data/gdd_gen_diff_raster.png \
#     webmap_data/gdd_gen_diff_raster_meta.json

args <- commandArgs(trailingOnly = TRUE)
input_json <- if (length(args) >= 1) args[[1]] else "webmap_data/GDD_GEN_1km_preview.json"
output_png <- if (length(args) >= 2) args[[2]] else "webmap_data/gdd_gen_diff_raster.png"
output_meta <- if (length(args) >= 3) args[[3]] else "webmap_data/gdd_gen_diff_raster_meta.json"

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

# Robust step estimation from successive points in scan order.
lat_diff <- abs(diff(x$lat))
lng_diff <- abs(diff(x$lng))

h_steps <- lng_diff[lat_diff < 2e-4 & lng_diff > 0]
v_steps <- lat_diff[lng_diff < 2e-4 & lat_diff > 0]

if (length(h_steps) == 0 || length(v_steps) == 0) {
  stop("Could not estimate grid spacing from input point order")
}

step_lng <- as.numeric(median(h_steps))
step_lat <- as.numeric(median(v_steps))

# Bounds must be cell edges, not center coordinates.
west <- min(x$lng) - step_lng / 2
south <- min(x$lat) - step_lat / 2
east <- max(x$lng) + step_lng / 2
north <- max(x$lat) + step_lat / 2

nx <- as.integer(round((east - west) / step_lng))
ny <- as.integer(round((north - south) / step_lat))

if (nx <= 0 || ny <= 0) {
  stop("Invalid raster dimensions computed: ", nx, "x", ny)
}

# Aggregate points onto the derived raster grid.
sum_mat <- matrix(0, nrow = ny, ncol = nx)
cnt_mat <- matrix(0L, nrow = ny, ncol = nx)

x_idx <- as.integer(round((x$lng - west) / step_lng) + 1)
y_idx <- as.integer(round((north - x$lat) / step_lat) + 1)

valid <- x_idx >= 1 & x_idx <= nx & y_idx >= 1 & y_idx <= ny
x_idx <- x_idx[valid]
y_idx <- y_idx[valid]
v <- x$GEN_diff[valid]

for (i in seq_along(v)) {
  r <- y_idx[i]
  c <- x_idx[i]
  sum_mat[r, c] <- sum_mat[r, c] + v[i]
  cnt_mat[r, c] <- cnt_mat[r, c] + 1L
}

avg <- matrix(NA_real_, nrow = ny, ncol = nx)
mask <- cnt_mat > 0
avg[mask] <- sum_mat[mask] / cnt_mat[mask]

# Class colors (same bins used in app click classification).
cols <- matrix("#00000000", nrow = ny, ncol = nx)
cols[is.finite(avg) & avg <= -0.5] <- "#1565c0"
cols[is.finite(avg) & avg > -0.5 & avg <= 0.5] <- "#f0f4f8"
cols[is.finite(avg) & avg > 0.5 & avg <= 1.5] <- "#ff9800"
cols[is.finite(avg) & avg > 1.5] <- "#c62828"

png(output_png, width = nx, height = ny, bg = "transparent")
par(mar = c(0, 0, 0, 0), xaxs = "i", yaxs = "i")
plot.new()
rasterImage(as.raster(cols), 0, 0, 1, 1, interpolate = FALSE)
dev.off()

meta <- list(
  image = basename(output_png),
  bounds = list(
    south = unname(south),
    west = unname(west),
    north = unname(north),
    east = unname(east)
  ),
  cellSize = list(
    lng = step_lng,
    lat = step_lat
  ),
  dimensions = list(
    width = nx,
    height = ny
  ),
  legend = list(
    labels = c("<= -0.5", "-0.5 to 0.5", "0.5 to 1.5", "> 1.5"),
    colors = c("#1565c0", "#f0f4f8", "#ff9800", "#c62828")
  )
)

write_json(meta, output_meta, auto_unbox = TRUE, pretty = TRUE)

cat("Built raster:", output_png, "\n")
cat("Built metadata:", output_meta, "\n")
cat("Grid:", nx, "x", ny, " step(lng,lat)=", step_lng, step_lat, "\n")
