#!/usr/bin/env Rscript
# Export the CHELSA-native PHENIPS GDD and generation-count rasters (built by
# scratchpad/build_phenips_chelsa_rasters.py -- scattered from the same
# future_predictor_stack.npz / assembled_fixed25_spruce_only.parquet the
# CHELSA-native bark-beetle model itself trains on) into the SAME
# forecast_bundle "model" shape used for the bark-beetle models
# (chelsa_fixed25_monotone etc), for both reference (CHELSA observed, mean
# of years <= 2010) and future (SSP5-8.5/MPI-ESM1-2-HR 2071-2100) periods.
#
# These are climate PREDICTORS, not a risk forecast -- reusing the forecast-
# bundle "model" shape (rather than inventing a new overlay/loader) is a
# deliberate low-risk choice: it needs zero new index.html JS (the model
# dropdown is already populated dynamically from manifest.json, proven this
# session), at the cost of appearing in the same "Predictions" model
# dropdown as the bark-beetle forecasts. Each variable gets its own real
# value_range (GDD degree-days, generations 0-3), not a 0-1 probability --
# same "layer has its own scale" precedent the wind layer already
# established (wind saturates at 0.5; see ARCHITECTURE.md). Labels/notes
# make the distinction from a risk forecast explicit so it's not
# mistaken for one.
#
# Usage: Rscript webmap_data/export_phenips_chelsa.R

suppressPackageStartupMessages({
  library(terra)
  library(jsonlite)
})

src_dir <- "webmap_data/phenips_chelsa"
max_dim <- 768

ramp_colors <- c("#5e4fa2", "#3288bd", "#66c2a5", "#abdda4", "#e6f598",
                 "#ffffbf", "#fee08b", "#fdae61", "#f46d43", "#d53e4f", "#9e0142")
color_ramp <- colorRamp(ramp_colors, space = "rgb")
fade_floor <- 90

vars <- list(
  list(id = "chelsa_phenips_gdd", label = "CHELSA PHENIPS GDD (bark-temperature degree-days)",
       vmin = 0, vmax = 6100, mid_cut = 2800,
       ref_src = file.path(src_dir, "phenips_gdd_reference_1984_2010.tif"),
       fut_src = file.path(src_dir, "phenips_gdd_future_2071_2100.tif")),
  list(id = "chelsa_phenips_gen", label = "CHELSA PHENIPS generations/year (capped at 3)",
       vmin = 0, vmax = 3, mid_cut = 1.5,
       ref_src = file.path(src_dir, "phenips_gen_reference_1984_2010.tif"),
       fut_src = file.path(src_dir, "phenips_gen_future_2071_2100.tif"))
)

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

for (v in vars) {
  out_root <- file.path("webmap_data/forecast_bundle", v$id)
  dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

  periods <- list(list(key = "reference_1984_2010", src = v$ref_src),
                   list(key = "future_2071_2100", src = v$fut_src))

  for (p in periods) {
    cat("Processing", v$id, p$key, "...\n")
    r0 <- rast(p$src)
    out_dir <- file.path(out_root, p$key)
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    r <- downsample_then_mercator(r0)
    norm_fun <- function(vv) (pmin(pmax(vv, v$vmin), v$vmax) - v$vmin) / (v$vmax - v$vmin)

    prob_tif <- file.path(out_dir, "probability_ensemble_mean.tif")
    writeRaster(r, prob_tif, overwrite = TRUE, datatype = "FLT4S", gdal = c("COMPRESS=LZW"))
    merc_ext <- to_png(r, file.path(out_dir, "probability_ensemble_mean.png"),
      colorfun = function(vv) color_ramp(norm_fun(vv)),
      alpha_fun = function(vv) round(fade_floor + (255 - fade_floor) * norm_fun(vv)))
    meta <- list(image = "probability_ensemble_mean.png", bounds = bounds_from_ext(merc_ext),
                 dimensions = list(width = ncol(r), height = nrow(r)),
                 source = v$label, value_range = list(min = v$vmin, max = v$vmax))
    write_json(meta, file.path(out_dir, "probability_ensemble_mean_meta.json"), pretty = TRUE, auto_unbox = TRUE, digits = 8)

    # Single "risk_visible_t50" stand-in threshold at each variable's own
    # midpoint (2800 GDD / 1.5 generations) -- satisfies the app's
    # threshold-availability schema check (same reason chelsa_fixed25_
    # monotone ships a t50 raster even though its real threshold is 0.25);
    # this is NOT a risk cutoff, just "above/below the midpoint".
    risk_tif <- file.path(out_dir, "risk_visible_t50.tif")
    rv <- r
    rv[rv < v$mid_cut] <- NA
    writeRaster(rv, risk_tif, overwrite = TRUE, datatype = "FLT4S", gdal = c("COMPRESS=LZW"))
    merc_ext2 <- to_png(rv, file.path(out_dir, "risk_visible_t50.png"),
      colorfun = function(vv) color_ramp(norm_fun(vv)),
      alpha_fun = function(vv) rep(255, length(vv)))
    meta2 <- list(image = "risk_visible_t50.png", bounds = bounds_from_ext(merc_ext2),
                  dimensions = list(width = ncol(rv), height = nrow(rv)),
                  source = paste0(v$label, " (above midpoint ", v$mid_cut, ")"),
                  value_range = list(min = v$mid_cut, max = v$vmax))
    write_json(meta2, file.path(out_dir, "risk_visible_t50_meta.json"), pretty = TRUE, auto_unbox = TRUE, digits = 8)

    cat("  wrote", p$key, "\n")
  }

  metrics <- list(
    model = v$id,
    notes = paste0(v$label, " -- a CHELSA-native climate PREDICTOR (input to the ",
                    "bark-beetle models), not a disturbance-risk forecast itself. No ",
                    "precision/recall/F1 applies. Same source pipeline as ",
                    "chelsa_fixed25_monotone / chelsa_*_efda; see ",
                    "CHELSA_NATIVE_TRAINING/PROJECT_REFERENCE_bark_beetle_model.md. ",
                    "Reference period is the CHELSA-observed mean of years 1984-2010 ",
                    "(not a single-year snapshot); future is the single SSP5-8.5/",
                    "MPI-ESM1-2-HR 2071-2100 climatological-normal snapshot."),
    aggregate_by_cutoff = list(
      list(cutoff = as.character(v$mid_cut), precision = "NA", recall = "NA",
           specificity = "NA", F1 = "NA", red_fraction = "NA")
    )
  )
  write_json(metrics, file.path(out_root, "metrics_confusion.json"), pretty = TRUE, auto_unbox = TRUE)
}
cat("\nDONE\n")
