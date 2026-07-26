#!/usr/bin/env Rscript
# Export the CHELSA/MPI-ESM1-2-HR bark-beetle risk rasters (built in
# FINAL_MODEL/CHELSA_DELIVERABLE/ by build_chelsa_predict.R + build_chelsa_
# apply_model.R) as a new "chelsa_mpi_ssp585" forecast-bundle model: same
# continuous-probability rendering as colorize_prob_continuous.R (shared
# RISK_RAMP), plus one classified risk_visible_t50 so the app's existing
# threshold-availability check (hydrateForecastControls) has something to
# match against.
#
# This is a single climatological-normal snapshot per period (not a 3-GCM,
# multi-year forecast series like the other models), applied with the
# already-trained v3 model -- no train/test evaluation exists for it, so
# metrics_confusion.json reports NA rather than invented skill numbers.
#
# Usage: Rscript webmap_data/export_chelsa_bb_forecast.R

suppressPackageStartupMessages({
  library(terra)
  library(jsonlite)
})

src_dir <- "/Users/jaroslavcepl/REENFOCE_LOCAL_MODEL_WIEN/FINAL_MODEL/CHELSA_DELIVERABLE"
out_root <- "webmap_data/forecast_bundle/chelsa_mpi_ssp585"
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

periods <- list(
  list(key = "1981_2010", src = file.path(src_dir, "bb_risk_reference.tif")),
  list(key = "2071_2100", src = file.path(src_dir, "bb_risk_mpi_ssp585.tif"))
)

# Unified Spectral ramp, shared app-wide (see index.html RISK_RAMP).
ramp_colors <- c("#5e4fa2", "#3288bd", "#66c2a5", "#abdda4", "#e6f598",
                 "#ffffbf", "#fee08b", "#fdae61", "#f46d43", "#d53e4f", "#9e0142")
color_ramp <- colorRamp(ramp_colors, space = "rgb")
fade_knee <- 0.3
fade_floor <- 90

risk_breaks <- c(0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9)
# 8 stops resampled off the 11-stop ramp (8 bins need 8 colors, not 11).
risk_colors <- c("#00000000", "#5e4fa2", "#48a1b3", "#a1d9a4", "#edf8a3",
                 "#fee99a", "#fca55d", "#e2524a", "#9e0142")

to_mercator_png <- function(r, out_png, colorfun, alpha_fun) {
  r <- project(r, "EPSG:3857", method = "bilinear")
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

for (p in periods) {
  r <- rast(p$src)
  out_dir <- file.path(out_root, p$key)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  # -- continuous probability (main view) --
  prob_tif <- file.path(out_dir, "probability_ensemble_mean.tif")
  writeRaster(r, prob_tif, overwrite = TRUE, datatype = "FLT4S", gdal = c("COMPRESS=LZW"))
  rp <- rast(prob_tif)
  merc_ext <- to_mercator_png(rp, file.path(out_dir, "probability_ensemble_mean.png"),
    colorfun = function(vv) color_ramp(pmin(pmax(vv, 0), 1)),
    alpha_fun = function(vv) round(fade_floor + (255 - fade_floor) * pmin(1, vv / fade_knee)))
  meta <- list(image = "probability_ensemble_mean.png", bounds = bounds_from_ext(merc_ext),
               dimensions = list(width = ncol(rp), height = nrow(rp)),
               source = "bb_risk (CHELSA-driven, v3 model, no retraining)",
               value_range = list(min = 0, max = 1))
  write_json(meta, file.path(out_dir, "probability_ensemble_mean_meta.json"), pretty = TRUE, auto_unbox = TRUE, digits = 8)

  # -- classified risk_visible_t50 (only threshold provided; satisfies the
  # app's threshold-availability check) --
  risk_tif <- file.path(out_dir, "risk_visible_t50.tif")
  writeRaster(r, risk_tif, overwrite = TRUE, datatype = "FLT4S", gdal = c("COMPRESS=LZW"))
  risk_colors_rgb <- substr(risk_colors, 1, 7)  # col2rgb can't parse 8-digit RGBA hex
  classify_color <- function(vv) {
    idx <- findInterval(vv, risk_breaks) + 1
    t(col2rgb(risk_colors_rgb[idx]))
  }
  classify_alpha <- function(vv) ifelse(vv < risk_breaks[1], 0, 255)
  rr <- rast(risk_tif)
  merc_ext2 <- to_mercator_png(rr, file.path(out_dir, "risk_visible_t50.png"),
    colorfun = classify_color, alpha_fun = classify_alpha)
  meta2 <- list(image = "risk_visible_t50.png", bounds = bounds_from_ext(merc_ext2),
                dimensions = list(width = ncol(rr), height = nrow(rr)),
                source = "bb_risk (CHELSA-driven, v3 model, no retraining)",
                value_range = list(min = 0.55, max = 1))
  write_json(meta2, file.path(out_dir, "risk_visible_t50_meta.json"), pretty = TRUE, auto_unbox = TRUE, digits = 8)

  cat("Exported", p$key, "\n")
}

metrics <- list(
  model = "chelsa_mpi_ssp585",
  notes = "Single climatological-normal snapshot per period (CHELSA 1981-2010 observed; MPI-ESM1-2-HR/SSP5-8.5 2071-2100), scored with the already-trained v3 geography-out model (no retraining). No train/test evaluation exists for this scenario -- metrics are NA, not invented skill numbers.",
  aggregate_by_cutoff = list(
    list(cutoff = "0.5", precision = "NA", recall = "NA", specificity = "NA", F1 = "NA", red_fraction = "NA")
  )
)
write_json(metrics, file.path(out_root, "metrics_confusion.json"), pretty = TRUE, auto_unbox = TRUE)
cat("Wrote metrics_confusion.json\n")
