#!/usr/bin/env Rscript
# Export the remaining CHELSA-native bark-beetle deployment models (all 4
# learners x both label sources, minus fixed25/xgb_monotone which
# export_chelsa_fixed25_monotone.R already covers) as forecast-bundle
# "model" entries, from the already-existing, already-exported GeoTIFFs in
# deployment_hazard_maps/geotiffs/ (built 2026-08-16; nothing recomputed
# here, just repackaged for the app in the same PNG+meta shape as
# export_chelsa_fixed25_monotone.R).
#
# CV metrics embedded per model/label are the real 5-seed spatial-block
# forward-direction numbers from CHELSA_NATIVE_TRAINING/
# PROJECT_REFERENCE_bark_beetle_model.md section 9's full 4-learner table
# (outputs/all_4learner_cv_results.csv) -- not invented, not NA.
#
# Usage: Rscript webmap_data/export_chelsa_all_learners.R

suppressPackageStartupMessages({
  library(terra)
  library(jsonlite)
})
source("webmap_data/color_scales.R")

src_root <- "/Users/jaroslavcepl/REENFOCE_LOCAL_MODEL_WIEN/CHELSA_NATIVE_TRAINING/deployment_hazard_maps/geotiffs"
max_dim  <- 768

# Bark beetle's own scale -- see color_scales.R.
ramp_colors <- BARK_BEETLE_RAMP
color_ramp <- colorRamp(ramp_colors, space = "rgb")
fade_knee  <- 0.3
fade_floor <- 90

# label_source, learner -> (skip fixed25/xgb_monotone, already exported)
# forward-direction, 5-seed-mean CV from PROJECT_REFERENCE section 9's table.
combos <- list(
  list(label_src = "fixed25", learner = "xgb",          auc = 0.703, f1_50 = 0.292, f1max = 0.478, recall = 0.664),
  list(label_src = "fixed25", learner = "rf",            auc = 0.780, f1_50 = 0.462, f1max = 0.549, recall = 0.797),
  list(label_src = "fixed25", learner = "gam",            auc = 0.601, f1_50 = 0.352, f1max = 0.424, recall = 0.497),
  list(label_src = "efda",    learner = "xgb",          auc = 0.599, f1_50 = 0.088, f1max = 0.436, recall = 0.820),
  list(label_src = "efda",    learner = "rf",            auc = 0.630, f1_50 = 0.137, f1max = 0.449, recall = 0.673),
  list(label_src = "efda",    learner = "gam",            auc = 0.527, f1_50 = 0.131, f1max = 0.415, recall = 0.948),
  list(label_src = "efda",    learner = "xgb_monotone", auc = 0.568, f1_50 = 0.223, f1max = 0.428, recall = 0.907)
)

learner_label <- c(xgb = "XGB", rf = "RF", gam = "GAM", xgb_monotone = "monotone XGB")
label_src_label <- c(fixed25 = "fixed25", efda = "EFDA")

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

risk_breaks <- BARK_BEETLE_CLASSIFIED_BREAKS
risk_colors <- BARK_BEETLE_CLASSIFIED_COLORS
risk_colors_rgb <- substr(risk_colors, 1, 7)
classify_color <- function(vv) { idx <- findInterval(vv, risk_breaks) + 1; t(col2rgb(risk_colors_rgb[idx])) }
classify_alpha <- function(vv) ifelse(vv < risk_breaks[1], 0, 255)

for (cmb in combos) {
  model_id <- paste0("chelsa_", cmb$label_src, "_", cmb$learner)
  out_root <- file.path("webmap_data/forecast_bundle", model_id)
  dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

  periods <- list(
    list(key = "1983_2010",
         src = file.path(src_root, paste0("reference_1983_2010_", cmb$label_src, "_spruce_only_", cmb$learner, "_gradient.tif"))),
    list(key = "2071_2100",
         src = file.path(src_root, paste0("future_ssp585_2071_2100_", cmb$label_src, "_spruce_only_", cmb$learner, "_gradient.tif")))
  )

  for (p in periods) {
    if (!file.exists(p$src)) { cat("MISSING", p$src, "-- skipping\n"); next }
    cat("Processing", model_id, p$key, "...\n")
    r0 <- rast(p$src)
    out_dir <- file.path(out_root, p$key)
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    r <- downsample_then_mercator(r0)

    prob_tif <- file.path(out_dir, "probability_ensemble_mean.tif")
    writeRaster(r, prob_tif, overwrite = TRUE, datatype = "FLT4S", gdal = c("COMPRESS=LZW"))
    merc_ext <- to_png(r, file.path(out_dir, "probability_ensemble_mean.png"),
      colorfun = function(vv) color_ramp(pmin(pmax(vv, 0), 1)),
      alpha_fun = function(vv) round(fade_floor + (255 - fade_floor) * pmin(1, vv / fade_knee)))
    meta <- list(image = "probability_ensemble_mean.png", bounds = bounds_from_ext(merc_ext),
                 dimensions = list(width = ncol(r), height = nrow(r)),
                 source = paste0("CHELSA-native ", label_src_label[[cmb$label_src]], "/", learner_label[[cmb$learner]]),
                 value_range = list(min = 0, max = 1))
    write_json(meta, file.path(out_dir, "probability_ensemble_mean_meta.json"), pretty = TRUE, auto_unbox = TRUE, digits = 8)

    risk_tif <- file.path(out_dir, "risk_visible_t50.tif")
    writeRaster(r, risk_tif, overwrite = TRUE, datatype = "FLT4S", gdal = c("COMPRESS=LZW"))
    merc_ext2 <- to_png(r, file.path(out_dir, "risk_visible_t50.png"),
      colorfun = classify_color, alpha_fun = classify_alpha)
    meta2 <- list(image = "risk_visible_t50.png", bounds = bounds_from_ext(merc_ext2),
                  dimensions = list(width = ncol(r), height = nrow(r)),
                  source = paste0("CHELSA-native ", label_src_label[[cmb$label_src]], "/", learner_label[[cmb$learner]]),
                  value_range = list(min = 0.55, max = 1))
    write_json(meta2, file.path(out_dir, "risk_visible_t50_meta.json"), pretty = TRUE, auto_unbox = TRUE, digits = 8)
  }

  metrics <- list(
    model = model_id,
    notes = paste0(
      "CHELSA-native training, ", label_src_label[[cmb$label_src]], " label, ",
      learner_label[[cmb$learner]], " learner. Real held-out spatial-block CV ",
      "(200km EPSG:3035 tiles, 80/20, 5 seeds, forward direction = deployable): ",
      "AUC ", cmb$auc, ", F1@0.5 ", cmb$f1_50, ", F1max ", cmb$f1max,
      " (recall@F1max ", cmb$recall, "). ",
      if (cmb$learner == "xgb_monotone")
        "Monotone-constrained (non-decreasing risk with warming, by construction)."
      else
        "NOT monotone-constrained -- per PROJECT_REFERENCE_bark_beetle_model.md section 4-5, this learner showed an out-of-training-envelope extrapolation artifact when projected to 2071-2100 (future risk understated relative to the monotone variant); compare against chelsa_fixed25_monotone before treating this future map as the primary read.",
      " Future period is a single climatological-normal snapshot (SSP5-8.5/MPI-ESM1-2-HR 2071-2100), not a multi-year forecast series. Full provenance: CHELSA_NATIVE_TRAINING/PROJECT_REFERENCE_bark_beetle_model.md section 9."
    ),
    aggregate_by_cutoff = list(
      list(cutoff = "0.5", precision = "NA", recall = "NA", specificity = "NA", F1 = cmb$f1_50, red_fraction = "NA"),
      list(cutoff = "F1max", precision = "NA", recall = cmb$recall, specificity = "NA", F1 = cmb$f1max, red_fraction = "NA")
    )
  )
  write_json(metrics, file.path(out_root, "metrics_confusion.json"), pretty = TRUE, auto_unbox = TRUE)
  cat("Wrote metrics for", model_id, "\n\n")
}
cat("DONE\n")
