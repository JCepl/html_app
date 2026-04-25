# ------------------------------
# 0) Setup
# ------------------------------
preferred_wd <- "/Users/jaroslavcepl/Library/Mobile Documents/com~apple~CloudDocs/Documents/D_MAC_2024/CE_maps_Poland"
if (dir.exists(preferred_wd)) {
  setwd(preferred_wd)
}

output_web_dir <- file.path(getwd(), "webmap_data")
if (!dir.exists(output_web_dir)) {
  dir.create(output_web_dir, recursive = TRUE)
}

library(terra)
library(sf)
library(jsonlite)
library(magrittr)

# Central Europe bbox in WGS84
bbox_wgs84 <- ext(5, 30, 42, 55)

# Keep 7 months (Apr-Oct style window in your setup)
months_idx <- 1:7
days_per_month <- c(30, 31, 30, 31, 31, 30, 31)

# Temperature thresholds
T_base <- 8.3
T_max  <- 38.9

# ------------------------------
# 1) Load monthly rasters + crop
# ------------------------------
euro_av_temp_future <- euro_av_temp_now <- vector("list", length(months_idx))

for (i in months_idx) {
  month_code <- 3 + i
  mm <- ifelse(month_code < 10, paste0("0", month_code), as.character(month_code))
  
  future_file <- paste0("CHELSA_mpi-esm1-2-hr_r1i1p1f1_w5e5_ssp370_tas_", mm, "_2041_2070_norm.tif")
  now_file    <- paste0("CHELSA_mpi-esm1-2-hr_r1i1p1f1_w5e5_ssp370_tas_", mm, "_2011_2040_norm.tif")
  
  euro_av_temp_future[[i]] <- crop(rast(future_file), bbox_wgs84)
  euro_av_temp_now[[i]]    <- crop(rast(now_file),    bbox_wgs84)
}

# ------------------------------
# 2) Convert monthly temperature to monthly GDD
# ------------------------------
for (i in months_idx) {
  # Future
  r_f <- euro_av_temp_future[[i]]
  r_f[r_f > T_max] <- T_max
  r_f <- r_f - T_base
  r_f[r_f < 0] <- 0
  euro_av_temp_future[[i]] <- r_f * days_per_month[i]
  
  # Now
  r_n <- euro_av_temp_now[[i]]
  r_n[r_n > T_max] <- T_max
  r_n <- r_n - T_base
  r_n[r_n < 0] <- 0
  euro_av_temp_now[[i]] <- r_n * days_per_month[i]
}

# Annual (7-month) GDD
GDD_now    <- Reduce(`+`, euro_av_temp_now)
GDD_future <- Reduce(`+`, euro_av_temp_future)

# ------------------------------
# 3) Build strict 1 km grid in EPSG:3035
# ------------------------------
bbox_poly_wgs84 <- as.polygons(bbox_wgs84, crs = "EPSG:4326")
bbox_poly_3035  <- project(bbox_poly_wgs84, "EPSG:3035")
e <- ext(bbox_poly_3035)

# Align extent to exact 1000 m boundaries
e_aligned <- ext(
  floor(xmin(e) / 1000) * 1000,
  ceiling(xmax(e) / 1000) * 1000,
  floor(ymin(e) / 1000) * 1000,
  ceiling(ymax(e) / 1000) * 1000
)

grid_1km <- rast(e_aligned, resolution = 1000, crs = "EPSG:3035")

# ------------------------------
# 4) Reproject/resample to aligned 1 km grid
# ------------------------------
# Bilinear for continuous variables (GDD)
GDD_now_1km    <- project(GDD_now,    grid_1km, method = "bilinear")
GDD_future_1km <- project(GDD_future, grid_1km, method = "bilinear")

# Optional exact mask to bbox polygon
GDD_now_1km    <- mask(GDD_now_1km, bbox_poly_3035)
GDD_future_1km <- mask(GDD_future_1km, bbox_poly_3035)

# ------------------------------
# 5) Derive generations on 1 km aligned grid
# ------------------------------
GEN_now_1km    <- floor(GDD_now_1km / 557)
GEN_future_1km <- floor(GDD_future_1km / 557)
GEN_diff_1km   <- GEN_future_1km - GEN_now_1km

names(GDD_now_1km)    <- "GDD_now"
names(GDD_future_1km) <- "GDD_future"
names(GEN_now_1km)    <- "GEN_now"
names(GEN_future_1km) <- "GEN_future"
names(GEN_diff_1km)   <- "GEN_diff"

# ------------------------------
# 6) Raster export (analysis-ready)
# ------------------------------
writeRaster(GDD_now_1km,    "GDD_now_1km_EPSG3035.tif",    overwrite = TRUE, wopt = list(gdal = c("COMPRESS=LZW")))
writeRaster(GDD_future_1km, "GDD_future_1km_EPSG3035.tif", overwrite = TRUE, wopt = list(gdal = c("COMPRESS=LZW")))
writeRaster(GEN_now_1km,    "GEN_now_1km_EPSG3035.tif",    overwrite = TRUE, wopt = list(gdal = c("COMPRESS=LZW")))
writeRaster(GEN_future_1km, "GEN_future_1km_EPSG3035.tif", overwrite = TRUE, wopt = list(gdal = c("COMPRESS=LZW")))
writeRaster(GEN_diff_1km,   "GEN_diff_1km_EPSG3035.tif",   overwrite = TRUE, wopt = list(gdal = c("COMPRESS=LZW")))

# Optional: WGS84 copies for quick web interoperability
GEN_diff_wgs84 <- project(GEN_diff_1km, "EPSG:4326", method = "near")
writeRaster(GEN_diff_wgs84, "GEN_diff_1km_WGS84.tif", overwrite = TRUE, wopt = list(gdal = c("COMPRESS=LZW")))
writeRaster(GEN_diff_wgs84, file.path(output_web_dir, "gdd_gen_diff_1km_wgs84.tif"), overwrite = TRUE, wopt = list(gdal = c("COMPRESS=LZW")))

# ------------------------------
# 7) JSON export from aligned raster
# ------------------------------
# Warning: full 1 km JSON over Central Europe can be very large.
# This exports point records at cell centers.

stacked <- c(GDD_now_1km, GDD_future_1km, GEN_now_1km, GEN_future_1km, GEN_diff_1km)
df_3035 <- as.data.frame(stacked, xy = TRUE, na.rm = TRUE)

# Convert grid cell center coordinates to WGS84 lon/lat for web maps
pts_3035 <- vect(df_3035[, c("x", "y")], geom = c("x", "y"), crs = "EPSG:3035")
pts_wgs84 <- project(pts_3035, "EPSG:4326")
lonlat <- crds(pts_wgs84)

json_data <- data.frame(
  lng = lonlat[, 1],
  lat = lonlat[, 2],
  GDD_now = df_3035$GDD_now,
  GDD_future = df_3035$GDD_future,
  GEN_now = df_3035$GEN_now,
  GEN_future = df_3035$GEN_future,
  GEN_diff = df_3035$GEN_diff
)

write_json(json_data, file.path(output_web_dir, "GDD_GEN_1km_points.json"), pretty = FALSE, auto_unbox = TRUE)

# ------------------------------
# 8) Web overlay export (aligned)
# ------------------------------
vals <- values(GEN_diff_wgs84, mat = FALSE)
nx <- ncol(GEN_diff_wgs84)
ny <- nrow(GEN_diff_wgs84)

m <- matrix(vals, nrow = ny, ncol = nx, byrow = TRUE)
rgba <- matrix("#00000000", nrow = ny, ncol = nx)
rgba[is.finite(m) & m <= -0.5] <- "#1565c0ff"
rgba[is.finite(m) & m > -0.5 & m <= 0.5] <- "#f0f4f8ff"
rgba[is.finite(m) & m > 0.5 & m <= 1.5] <- "#ff9800ff"
rgba[is.finite(m) & m > 1.5] <- "#c62828ff"

png(file.path(output_web_dir, "gdd_gen_diff_raster.png"), width = nx, height = ny, bg = "transparent")
par(mar = c(0, 0, 0, 0), xaxs = "i", yaxs = "i")
plot.new()
rasterImage(as.raster(rgba), 0, 0, 1, 1, interpolate = FALSE)
dev.off()

e_wgs84 <- ext(GEN_diff_wgs84)
meta <- list(
  image = "gdd_gen_diff_raster.png",
  bounds = list(
    south = unname(ymin(e_wgs84)),
    west = unname(xmin(e_wgs84)),
    north = unname(ymax(e_wgs84)),
    east = unname(xmax(e_wgs84))
  ),
  dimensions = list(
    width = nx,
    height = ny
  ),
  source = "gdd_gen_diff_1km_wgs84.tif",
  legend = list(
    labels = c("<= -0.5", "-0.5 to 0.5", "0.5 to 1.5", "> 1.5"),
    colors = c("#1565c0", "#f0f4f8", "#ff9800", "#c62828")
  )
)

write_json(meta, file.path(output_web_dir, "gdd_gen_diff_raster_meta.json"), pretty = TRUE, auto_unbox = TRUE)

cat("Done: 1 km aligned raster outputs, JSON export, and web overlay files created.\n")