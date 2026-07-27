#!/usr/bin/env Rscript
# Bark-beetle ensemble forecast: "current vs future" delta layers, as
# RELATIVE percent change (not absolute probability points) -- more honest
# when baseline risk is low (a lot of this map is), where a small absolute
# gain is actually a large relative one.
#
# Two comparisons per future period:
#   delta_from_near  = this period vs 2011-2040 ("current" -- the ensemble
#                       has no independent historical/observed baseline at
#                       this resolution, so near-future stands in for it,
#                       same convention the app already uses elsewhere)
#   delta_from_prior = this period vs the PREVIOUS period (sequential change:
#                       mid-century vs current for 2041-2070; end-century vs
#                       mid-century for 2071-2100) -- for 2041-2070 this is
#                       identical to delta_from_near since its "prior" IS
#                       2011-2040.
#
# IMPORTANT -- alignment: probability_ensemble_mean.tif's embedded
# georeferencing does NOT match its own sibling probability_ensemble_mean_
# meta.json/.png (pre-existing mismatch in this repo: same west/north
# anchor, different south/east, and a DIFFERENT PIXEL GRID ENTIRELY -- e.g.
# rcp85 is 383x331 in the tif vs 272x427 in the meta.json/png). The app has
# only ever displayed the .png at the meta.json's bounds+shape, so that's
# the coordinate system every layer must match. Just copying the meta.json's
# bounds onto the tif's native 383x331 array (first attempt) left the
# aspect ratio wrong -- Leaflet stretches whatever image into the given
# LatLngBounds rectangle regardless of the image's own pixel shape, so a
# 383x331 image forced into a 272x427-shaped box comes out visibly
# distorted/offset. Fixed by resampling onto a template that actually
# matches the sibling's grid (same ncol/nrow AND extent) before colorizing.
#
# Usage: Rscript webmap_data/forecast_bundle/build_ensemble_delta.R

suppressPackageStartupMessages(library(terra))

base_dir <- "/Users/jaroslavcepl/html_app/webmap_data/forecast_bundle/ensemble"
scenarios <- c("rcp26", "rcp45", "rcp85")
period_order <- c("2011_2040", "2041_2070", "2071_2100")

DELTA_RAMP <- c("#2166ac", "#4393c3", "#92c5de", "#d1e5f0", "#f7f7f7",
                "#f7f7f7", "#fddbc7", "#f4a582", "#d6604d", "#b2182b", "#67001f")
PCT_CAP <- 200
PCT_BASE_FLOOR <- 0.03
FADE_BAND <- 0.12
FADE_FLOOR <- 40

log_step <- function(...) { cat(format(Sys.time(), "%H:%M:%S"), "|", ..., "\n"); flush.console() }

colorize_delta <- function(r, cap) {
  v <- values(r)[, 1]
  vc <- pmax(-cap, pmin(cap, v))
  frac <- (vc + cap) / (2 * cap)
  idx <- round(frac * (length(DELTA_RAMP) - 1)) + 1
  idx[is.na(v)] <- NA
  cols <- rep(NA_character_, length(v))
  cols[!is.na(idx)] <- DELTA_RAMP[idx[!is.na(idx)]]
  hex2rgb <- function(h) list(
    r = strtoi(substr(h, 2, 3), 16), g = strtoi(substr(h, 4, 5), 16), b = strtoi(substr(h, 6, 7), 16)
  )
  rgb <- hex2rgb(cols)
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

# Read a period's probability raster AND resample it onto a template built
# from its OWN sibling meta.json (grid + bounds the app actually trusts).
read_aligned <- function(scen, per) {
  tif <- file.path(base_dir, scen, per, "probability_ensemble_mean.tif")
  meta_f <- file.path(base_dir, scen, per, "probability_ensemble_mean_meta.json")
  if (!file.exists(tif) || !file.exists(meta_f)) return(NULL)
  meta <- jsonlite::read_json(meta_f)
  b <- meta$bounds
  tmpl <- rast(xmin = b$west, xmax = b$east, ymin = b$south, ymax = b$north,
               ncols = meta$dimensions$width, nrows = meta$dimensions$height,
               crs = "EPSG:4326")
  r <- rast(tif)
  ext(r) <- ext(tmpl)  # same real-world footprint (confirmed via matching
                        # west/north anchors); only the pixel grid differs
  crs(r) <- crs(tmpl)
  resample(r, tmpl, method = "bilinear")
}

write_delta <- function(scen, later_per, earlier_per, out_name, label) {
  later_r <- read_aligned(scen, later_per)
  earlier_r <- read_aligned(scen, earlier_per)
  if (is.null(later_r) || is.null(earlier_r)) { log_step("  SKIP:", scen, out_name, "(missing input)"); return(invisible()) }

  base_floor_r <- clamp(earlier_r, lower = PCT_BASE_FLOOR, values = TRUE)
  pct_change <- (later_r - earlier_r) / base_floor_r * 100
  names(pct_change) <- out_name

  out_dir <- file.path(base_dir, scen, later_per)
  writeRaster(pct_change, file.path(out_dir, paste0(out_name, ".tif")), overwrite = TRUE,
              wopt = list(gdal = c("COMPRESS=DEFLATE", "PREDICTOR=2", "ZLEVEL=6")))
  png::writePNG(colorize_delta(pct_change, PCT_CAP), file.path(out_dir, paste0(out_name, ".png")))

  meta_f <- file.path(base_dir, scen, later_per, "probability_ensemble_mean_meta.json")
  bounds <- jsonlite::read_json(meta_f)$bounds
  v <- values(pct_change)[, 1]
  meta <- list(
    image = paste0(out_name, ".png"),
    bounds = bounds,
    dimensions = list(width = ncol(pct_change), height = nrow(pct_change)),
    source = paste0("ensemble ", scen, ": % change ", later_per, " vs ", earlier_per, " (", label, ")"),
    value_range = list(min = -PCT_CAP, max = PCT_CAP, unit = "percent"),
    actual_range = list(min = round(min(v, na.rm = TRUE), 1), max = round(max(v, na.rm = TRUE), 1)),
    baseline_period = earlier_per,
    baseline_floor = PCT_BASE_FLOOR
  )
  jsonlite::write_json(meta, file.path(out_dir, paste0(out_name, "_meta.json")),
                        auto_unbox = TRUE, pretty = TRUE)
  log_step(scen, out_name, "(", later_per, "vs", earlier_per, ") -> % change range",
            round(min(v, na.rm = TRUE), 1), "to", round(max(v, na.rm = TRUE), 1))
}

for (scen in scenarios) {
  # vs current (2011-2040), for both future periods
  write_delta(scen, "2041_2070", "2011_2040", "delta_from_near", "current")
  write_delta(scen, "2071_2100", "2011_2040", "delta_from_near", "current")
  # vs previous period (sequential change)
  write_delta(scen, "2041_2070", "2011_2040", "delta_from_prior", "previous period")
  write_delta(scen, "2071_2100", "2041_2070", "delta_from_prior", "previous period")
}

log_step("Done.")
