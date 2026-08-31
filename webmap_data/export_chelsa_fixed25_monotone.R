#!/usr/bin/env Rscript
# Litmus-test export of the NEW, properly-validated CHELSA-native bark-beetle
# model (fixed25/CATALOGUE25 label, monotone-constrained XGBoost) as a
# "chelsa_fixed25_monotone" forecast-bundle model -- distinct from, and not
# replacing, the existing "chelsa_mpi_ssp585" model (the old CORDEX-trained
# v3 model just scored on CHELSA climate, no retraining, no held-out CV of
# its own). This one DOES have real held-out CV: see
# REENFOCE_LOCAL_MODEL_WIEN/CHELSA_NATIVE_TRAINING/
# PROJECT_REFERENCE_bark_beetle_model.md §5-9 for full provenance.
#
# Source rasters are CHELSA-native resolution (~1km, 5914x4449 domain-wide),
# much larger than the old coarse EUR-11/0.11 chelsa_mpi_ssp585 rasters, so
# this downsamples to max_dim before reprojecting (same pattern as
# colorize_prob_continuous.R / build_ensemble_rcp85.R), unlike
# export_chelsa_bb_forecast.R which projected at full resolution.
#
# Usage: Rscript webmap_data/export_chelsa_fixed25_monotone.R

suppressPackageStartupMessages({
  library(terra)
  library(jsonlite)
})
source("webmap_data/color_scales.R")

src_dir  <- "/Users/jaroslavcepl/REENFOCE_LOCAL_MODEL_WIEN/CHELSA_NATIVE_TRAINING/anna_deliverable"
out_root <- "webmap_data/forecast_bundle/chelsa_fixed25_monotone"
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
max_dim  <- 768

periods <- list(
  list(key = "1983_2010", src = file.path(src_dir, "barkbeetle_CZU_fixed25_monotone_reference_continuous_1983_2010.tif")),
  list(key = "2071_2100", src = file.path(src_dir, "barkbeetle_CZU_fixed25_monotone_future_continuous_2071_2100.tif"))
)

# Bark beetle's own scale -- see color_scales.R.
ramp_colors <- BARK_BEETLE_RAMP
color_ramp <- colorRamp(ramp_colors, space = "rgb")
fade_knee  <- 0.3
fade_floor <- 90

risk_breaks <- BARK_BEETLE_CLASSIFIED_BREAKS
risk_colors <- BARK_BEETLE_CLASSIFIED_COLORS

downsample_then_mercator <- function(r) {
  fact <- ceiling(max(ncol(r), nrow(r)) / max_dim)
  if (fact > 1) r <- aggregate(r, fact = fact, fun = "mean", na.rm = TRUE)
  project(r, "EPSG:3857", method = "bilinear")
}

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

for (p in periods) {
  cat("Reading", p$src, "...\n")
  r0 <- rast(p$src)
  out_dir <- file.path(out_root, p$key)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  r <- downsample_then_mercator(r0)

  # -- continuous probability (main view) --
  prob_tif <- file.path(out_dir, "probability_ensemble_mean.tif")
  writeRaster(r, prob_tif, overwrite = TRUE, datatype = "FLT4S", gdal = c("COMPRESS=LZW"))
  merc_ext <- to_png(r, file.path(out_dir, "probability_ensemble_mean.png"),
    colorfun = function(vv) color_ramp(pmin(pmax(vv, 0), 1)),
    alpha_fun = function(vv) round(fade_floor + (255 - fade_floor) * pmin(1, vv / fade_knee)))
  meta <- list(image = "probability_ensemble_mean.png", bounds = bounds_from_ext(merc_ext),
               dimensions = list(width = ncol(r), height = nrow(r)),
               source = "CHELSA-native fixed25/monotone-XGB (see PROJECT_REFERENCE_bark_beetle_model.md)",
               value_range = list(min = 0, max = 1))
  write_json(meta, file.path(out_dir, "probability_ensemble_mean_meta.json"), pretty = TRUE, auto_unbox = TRUE, digits = 8)

  # -- classified risk_visible_t50 (satisfies the app's threshold-availability
  # check; the model's own CV-calibrated operating threshold is 0.25, lower
  # than this display cutoff -- disclosed in metrics_confusion.json below,
  # not silently reconciled) --
  risk_colors_rgb <- substr(risk_colors, 1, 7)
  classify_color <- function(vv) {
    idx <- findInterval(vv, risk_breaks) + 1
    t(col2rgb(risk_colors_rgb[idx]))
  }
  classify_alpha <- function(vv) ifelse(vv < risk_breaks[1], 0, 255)
  risk_tif <- file.path(out_dir, "risk_visible_t50.tif")
  writeRaster(r, risk_tif, overwrite = TRUE, datatype = "FLT4S", gdal = c("COMPRESS=LZW"))
  merc_ext2 <- to_png(r, file.path(out_dir, "risk_visible_t50.png"),
    colorfun = classify_color, alpha_fun = classify_alpha)
  meta2 <- list(image = "risk_visible_t50.png", bounds = bounds_from_ext(merc_ext2),
                dimensions = list(width = ncol(r), height = nrow(r)),
                source = "CHELSA-native fixed25/monotone-XGB",
                value_range = list(min = 0.55, max = 1))
  write_json(meta2, file.path(out_dir, "risk_visible_t50_meta.json"), pretty = TRUE, auto_unbox = TRUE, digits = 8)

  cat("Exported", p$key, "\n")
}

# Real held-out CV metrics (PROJECT_REFERENCE_bark_beetle_model.md §6,
# fixed25/xgb_monotone, forward direction = the deployable direction, 5-seed
# mean, spatial-block tile CV). Unlike chelsa_mpi_ssp585's NA metrics, this
# model has genuine calibration -- reported at its own operating threshold
# (0.25), not the 0.50 the display raster above uses for UI consistency with
# other models.
metrics <- list(
  model = "chelsa_fixed25_monotone",
  # Threshold-independent, forward direction (= deployable; backward is a
  # confounded robustness check, see notes below and PROJECT_REFERENCE
  # bark_beetle_model.md). Displayed in the app's metrics grid alongside
  # the per-cutoff precision/recall/F1 below.
  auc = 0.663,
  notes = paste(
    "CHELSA-native training (observed climate, 1983-1997+2000-2020), fixed25",
    "(literature catalogue, 25km buffer) label, monotone-constrained XGBoost.",
    "Real held-out spatial-block CV (200km EPSG:3035 tiles, 80/20, 5 seeds,",
    "forward direction = deployable): AUC 0.663, F1@0.5 0.408. At the",
    "CV-calibrated operating threshold 0.25: precision 0.368, recall 0.580,",
    "F1 0.441 (forward); backward-direction AUC 0.782 (robustness check,",
    "confounded by unequal training-set size, see the reference doc).",
    "Future period is a single climatological-normal snapshot",
    "(SSP5-8.5/MPI-ESM1-2-HR 2071-2100), not a multi-year forecast series.",
    "Full provenance: REENFOCE_LOCAL_MODEL_WIEN/CHELSA_NATIVE_TRAINING/",
    "PROJECT_REFERENCE_bark_beetle_model.md"
  ),
  aggregate_by_cutoff = list(
    list(cutoff = "0.25", precision = 0.368, recall = 0.580, specificity = "NA", F1 = 0.441, red_fraction = "NA"),
    list(cutoff = "0.5", precision = "NA", recall = "NA", specificity = "NA", F1 = "NA", red_fraction = "NA")
  )
)
write_json(metrics, file.path(out_root, "metrics_confusion.json"), pretty = TRUE, auto_unbox = TRUE)
cat("Wrote metrics_confusion.json\n")
