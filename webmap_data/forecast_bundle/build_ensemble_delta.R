#!/usr/bin/env Rscript
# Bark-beetle ensemble forecast: "current vs future" delta layers, as
# RELATIVE percent change (not absolute probability points) -- more honest
# when baseline risk is low (a lot of this map is), where a small absolute
# gain is actually a large relative one.
#
# The ensemble forecast has no independent historical/observed baseline at
# this raster's resolution (bbtl/bbtpast are a different model+variable
# entirely, not comparable pixel-for-pixel) -- so, per the app's own existing
# convention for this ("near-future" already stands in for "now" throughout
# the RCP/period selectors), 2011-2040 is used as the "current" baseline.
# pct_change = (later-period - 2011_2040) / max(2011_2040, PCT_BASE_FLOOR) * 100.
#
# Output mirrors the existing ensemble/{rcp}/{period}/ layout:
#   delta_from_near.tif / .png / _meta.json
#
# IMPORTANT: this scenario's own probability_ensemble_mean.tif has embedded
# georeferencing that does NOT match its sibling probability_ensemble_mean_
# meta.json/.png (a pre-existing mismatch in this repo -- confirmed by
# comparing ext() of the tif to the meta.json's bounds: same west/north
# anchor, different south/east and a completely different pixel grid). The
# app has only ever displayed the .png at the meta.json's bounds, so THAT is
# the coordinate system every other layer is aligned to -- reusing the tif's
# own ext() here previously caused the delta layer to render visibly
# offset (color bleeding into the sea). Always copy bounds from the sibling
# meta.json, never recompute from the tif.
#
# Usage: Rscript webmap_data/forecast_bundle/build_ensemble_delta.R

suppressPackageStartupMessages(library(terra))

base_dir <- "/Users/jaroslavcepl/html_app/webmap_data/forecast_bundle/ensemble"
scenarios <- c("rcp26", "rcp45", "rcp85")
baseline_period <- "2011_2040"
future_periods <- c("2041_2070", "2071_2100")

# Diverging, colorblind-checked ColorBrewer RdBu (11-class), red = increase
# (more risk), blue = decrease -- deliberately a different family from the
# shared Spectral RISK_RAMP used for absolute risk everywhere else in the
# app, so a delta map is never mistaken for one. The two middle stops (near-
# zero change) additionally fade toward transparent -- see colorize_delta --
# so only areas of real change draw the eye; near-zero is barely there.
DELTA_RAMP <- c("#2166ac", "#4393c3", "#92c5de", "#d1e5f0", "#f7f7f7",
                "#f7f7f7", "#fddbc7", "#f4a582", "#d6604d", "#b2182b", "#67001f")
PCT_CAP <- 100        # symmetric clamp: |% change| >= 100 saturates the ramp
PCT_BASE_FLOOR <- 0.03 # floor on the baseline denominator (avoids divide-by-
                        # near-zero blowups where 2011-2040 risk is ~0)
FADE_BAND <- 0.12       # within +/-12% of no change, alpha fades toward FADE_FLOOR
FADE_FLOOR <- 40        # alpha (0-255) right at 0% change; full 255 outside FADE_BAND

log_step <- function(...) { cat(format(Sys.time(), "%H:%M:%S"), "|", ..., "\n"); flush.console() }

colorize_delta <- function(r, cap) {
  v <- values(r)[, 1]
  vc <- pmax(-cap, pmin(cap, v))
  frac <- (vc + cap) / (2 * cap)  # 0..1, 0.5 = no change
  idx <- round(frac * (length(DELTA_RAMP) - 1)) + 1
  idx[is.na(v)] <- NA
  cols <- rep(NA_character_, length(v))
  cols[!is.na(idx)] <- DELTA_RAMP[idx[!is.na(idx)]]
  hex2rgb <- function(h) list(
    r = strtoi(substr(h, 2, 3), 16), g = strtoi(substr(h, 4, 5), 16), b = strtoi(substr(h, 6, 7), 16)
  )
  rgb <- hex2rgb(cols)

  # Near-zero fade: distance from 0% as a fraction of the fade band, clamped
  # to [0,1], scaled between FADE_FLOOR and 255.
  dist_frac <- pmin(1, abs(vc) / (cap * FADE_BAND))
  alpha_val <- FADE_FLOOR + dist_frac * (255 - FADE_FLOOR)
  alpha <- ifelse(is.na(v), 0, alpha_val)

  nr <- nrow(r); nc <- ncol(r)
  img <- array(0, dim = c(nr, nc, 4))
  fill <- function(x, d = 0) matrix(replace(x, is.na(x), d), nr, nc, byrow = TRUE) / 255
  img[, , 1] <- fill(rgb$r); img[, , 2] <- fill(rgb$g)
  img[, , 3] <- fill(rgb$b); img[, , 4] <- fill(alpha)
  img
}

for (scen in scenarios) {
  base_file <- file.path(base_dir, scen, baseline_period, "probability_ensemble_mean.tif")
  base_meta_file <- file.path(base_dir, scen, baseline_period, "probability_ensemble_mean_meta.json")
  if (!file.exists(base_file) || !file.exists(base_meta_file)) { log_step("SKIP (no baseline):", scen); next }
  baseline <- rast(base_file)
  base_bounds <- jsonlite::read_json(base_meta_file)$bounds

  for (per in future_periods) {
    fut_file <- file.path(base_dir, scen, per, "probability_ensemble_mean.tif")
    if (!file.exists(fut_file)) { log_step("  SKIP (missing future):", scen, per); next }
    future_r <- rast(fut_file)

    base_floor_r <- clamp(baseline, lower = PCT_BASE_FLOOR, values = TRUE)
    pct_change <- (future_r - baseline) / base_floor_r * 100
    names(pct_change) <- "pct_change_from_near"

    out_dir <- file.path(base_dir, scen, per)
    writeRaster(pct_change, file.path(out_dir, "delta_from_near.tif"), overwrite = TRUE,
                wopt = list(gdal = c("COMPRESS=DEFLATE", "PREDICTOR=2", "ZLEVEL=6")))

    png::writePNG(colorize_delta(pct_change, PCT_CAP), file.path(out_dir, "delta_from_near.png"))

    v <- values(pct_change)[, 1]
    meta <- list(
      image = "delta_from_near.png",
      bounds = base_bounds,  # copied from the sibling meta.json -- see header note
      dimensions = list(width = ncol(pct_change), height = nrow(pct_change)),
      source = paste0("ensemble ", scen, ": % change ", per, " vs ", baseline_period, " (\"current\")"),
      value_range = list(min = -PCT_CAP, max = PCT_CAP, unit = "percent"),
      actual_range = list(min = round(min(v, na.rm = TRUE), 1), max = round(max(v, na.rm = TRUE), 1)),
      baseline_period = baseline_period,
      baseline_floor = PCT_BASE_FLOOR
    )
    jsonlite::write_json(meta, file.path(out_dir, "delta_from_near_meta.json"),
                          auto_unbox = TRUE, pretty = TRUE)

    log_step(scen, per, "-> % change range", round(min(v, na.rm=TRUE),1), "to", round(max(v, na.rm=TRUE),1))
  }
}

log_step("Done.")
