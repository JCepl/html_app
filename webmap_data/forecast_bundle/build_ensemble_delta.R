#!/usr/bin/env Rscript
# Bark-beetle ensemble forecast: "current vs future" delta layers.
#
# The ensemble forecast has no independent historical/observed baseline at
# this raster's resolution (bbtl/bbtpast are a different model+variable
# entirely, not comparable pixel-for-pixel) -- so, per the app's own existing
# convention for this ("near-future" already stands in for "now" throughout
# the RCP/period selectors), 2011-2040 is used as the "current" baseline.
# delta = later-period probability - 2011-2040 probability, per RCP scenario.
#
# Output mirrors the existing ensemble/{rcp}/{period}/ layout:
#   delta_from_near.tif / .png / _meta.json
#
# Usage: Rscript webmap_data/forecast_bundle/build_ensemble_delta.R

suppressPackageStartupMessages(library(terra))

base_dir <- "/Users/jaroslavcepl/html_app/webmap_data/forecast_bundle/ensemble"
scenarios <- c("rcp26", "rcp45", "rcp85")
baseline_period <- "2011_2040"
future_periods <- c("2041_2070", "2071_2100")

# Diverging, colorblind-checked ColorBrewer RdBu (11-class), red = increase
# (more risk), blue = decrease, pale grey-white = ~no change -- deliberately
# a different family from the shared Spectral RISK_RAMP used for absolute
# risk everywhere else in the app, so a delta map is never mistaken for one.
DELTA_RAMP <- c("#2166ac", "#4393c3", "#92c5de", "#d1e5f0", "#f7f7f7",
                "#f7f7f7", "#fddbc7", "#f4a582", "#d6604d", "#b2182b", "#67001f")
DELTA_CAP <- 0.30  # symmetric clamp: |delta| >= 0.30 saturates the ramp

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
  alpha <- ifelse(is.na(v), 0L, 255L)
  nr <- nrow(r); nc <- ncol(r)
  img <- array(0, dim = c(nr, nc, 4))
  fill <- function(x, d = 0L) matrix(replace(x, is.na(x), d), nr, nc, byrow = TRUE) / 255
  img[, , 1] <- fill(rgb$r); img[, , 2] <- fill(rgb$g)
  img[, , 3] <- fill(rgb$b); img[, , 4] <- fill(alpha)
  img
}

for (scen in scenarios) {
  base_file <- file.path(base_dir, scen, baseline_period, "probability_ensemble_mean.tif")
  if (!file.exists(base_file)) { log_step("SKIP (no baseline):", scen); next }
  baseline <- rast(base_file)

  for (per in future_periods) {
    fut_file <- file.path(base_dir, scen, per, "probability_ensemble_mean.tif")
    if (!file.exists(fut_file)) { log_step("  SKIP (missing future):", scen, per); next }
    future_r <- rast(fut_file)

    delta <- future_r - baseline
    names(delta) <- "delta_from_near"

    out_dir <- file.path(base_dir, scen, per)
    writeRaster(delta, file.path(out_dir, "delta_from_near.tif"), overwrite = TRUE,
                wopt = list(gdal = c("COMPRESS=DEFLATE", "PREDICTOR=2", "ZLEVEL=6")))

    b <- ext(delta)
    png::writePNG(colorize_delta(delta, DELTA_CAP), file.path(out_dir, "delta_from_near.png"))

    v <- values(delta)[, 1]
    meta <- list(
      image = "delta_from_near.png",
      bounds = list(south = ymin(b), west = xmin(b), north = ymax(b), east = xmax(b)),
      dimensions = list(width = ncol(delta), height = nrow(delta)),
      source = paste0("ensemble ", scen, ": ", per, " minus ", baseline_period, " (\"current\")"),
      value_range = list(min = -DELTA_CAP, max = DELTA_CAP),
      actual_range = list(min = round(min(v, na.rm = TRUE), 4), max = round(max(v, na.rm = TRUE), 4)),
      baseline_period = baseline_period
    )
    jsonlite::write_json(meta, file.path(out_dir, "delta_from_near_meta.json"),
                          auto_unbox = TRUE, pretty = TRUE)

    log_step(scen, per, "-> delta range", round(min(v, na.rm=TRUE),3), "to", round(max(v, na.rm=TRUE),3))
  }
}

log_step("Done.")
