#!/usr/bin/env Rscript
# Bark-beetle ensemble forecast: "current vs future" delta layers, as
# RELATIVE percent change (not absolute probability points).
#
# delta_from_near = each future period (2041-2070, 2071-2100) vs 2011-2040
# ("current" -- no independent historical/observed baseline exists at this
# resolution, so near-future stands in for it, same convention used
# elsewhere in the app).
#
# CRITICAL: do NOT use probability_ensemble_mean.tif's own embedded
# georeferencing for anything, not even as a resampling source. It does not
# describe the same grid as its own sibling meta.json/.png -- confirmed
# repeatedly (matching a couple of corners was a coincidence, not proof of a
# shared frame; an earlier resample-based fix still came out shifted by
# ~200km). The .png is what the app actually displays, already exactly
# aligned with its meta.json bounds (that combination is what users have
# been looking at, correctly, all along). So: decode the EXISTING
# probability PNG's own pixels back into values via the same Spectral ramp
# the app colors them with, compute the delta directly on THAT array, and
# write the output on that IDENTICAL pixel grid. Zero resampling, zero
# reliance on the tif -- the output can't be misaligned because it's the
# same grid, pixel for pixel.
#
# The PNG's own colors were rendered from continuous float values via
# piecewise-LINEAR interpolation along the 11-stop ramp (same as the app's
# client-side rampColor()), not snapped to the 11 stops themselves -- an
# earlier version of this script decoded against only the 11 raw stops,
# which quantized everything to 11 discrete levels and produced blocky,
# uniform "plateau" regions in the delta (most visible over Spain/UK) that
# looked exactly like another alignment bug but wasn't one. Fixed by
# decoding against a 256-step LUT built the same way the PNG's colors were
# generated (interpolated between consecutive stops), recovering
# near-continuous precision.
#
# Usage: Rscript webmap_data/forecast_bundle/build_ensemble_delta.R

suppressPackageStartupMessages(library(png))

base_dir <- "/Users/jaroslavcepl/html_app/webmap_data/forecast_bundle/ensemble"
scenarios <- c("rcp26", "rcp45", "rcp85")

# Same Spectral ramp + domain (0..1 probability) the app colors
# probability_ensemble_mean.png with -- see RISK_RAMP in index.html.
PROB_RAMP <- c("#5e4fa2", "#3288bd", "#66c2a5", "#abdda4", "#e6f598",
               "#ffffbf", "#fee08b", "#fdae61", "#f46d43", "#d53e4f", "#9e0142")

DELTA_RAMP <- c("#2166ac", "#4393c3", "#92c5de", "#d1e5f0", "#f7f7f7",
                "#f7f7f7", "#fddbc7", "#f4a582", "#d6604d", "#b2182b", "#67001f")
PCT_CAP <- 200
PCT_BASE_FLOOR <- 0.03
FADE_BAND <- 0.12
FADE_FLOOR <- 40

log_step <- function(...) { cat(format(Sys.time(), "%H:%M:%S"), "|", ..., "\n"); flush.console() }

hex2rgb01 <- function(hexvec) {
  m <- t(grDevices::col2rgb(hexvec)) / 255
  m
}

# Piecewise-linear interpolation along a hex ramp at n samples -- mirrors
# the app's own rampColor()/buildRiskLut() (also n=256), so decoding against
# this recovers the same continuous value the ramp was built from.
build_fine_lut <- function(hex_ramp, n = 256) {
  stops <- hex2rgb01(hex_ramp)
  k <- nrow(stops)
  out <- matrix(0, nrow = n, ncol = 3)
  for (i in seq_len(n)) {
    t <- (i - 1) / (n - 1) * (k - 1)  # position in stop-space, 0..(k-1)
    lo <- floor(t); hi <- min(lo + 1, k - 1)
    frac <- t - lo
    out[i, ] <- stops[lo + 1, ] * (1 - frac) + stops[hi + 1, ] * frac
  }
  out
}
PROB_LUT <- build_fine_lut(PROB_RAMP, 256)

# Reverse-map an RGB pixel to its nearest ramp position (0..1), vectorized.
# alpha == 0 (fully transparent) means "no data" -> NA.
decode_prob_png <- function(path) {
  img <- readPNG(path)  # array [h, w, {3 or 4}]
  h <- dim(img)[1]; w <- dim(img)[2]
  has_alpha <- dim(img)[3] == 4
  r <- as.vector(img[, , 1]); g <- as.vector(img[, , 2]); b <- as.vector(img[, , 3])
  alpha <- if (has_alpha) as.vector(img[, , 4]) else rep(1, length(r))

  n <- length(r)
  best_t <- rep(NA_real_, n)
  best_d <- rep(Inf, n)
  for (i in seq_len(nrow(PROB_LUT))) {
    d <- (r - PROB_LUT[i, 1])^2 + (g - PROB_LUT[i, 2])^2 + (b - PROB_LUT[i, 3])^2
    better <- d < best_d
    best_d[better] <- d[better]
    best_t[better] <- (i - 1) / (nrow(PROB_LUT) - 1)
  }
  best_t[alpha < 0.01] <- NA
  list(values = matrix(best_t, nrow = h, ncol = w, byrow = FALSE), width = w, height = h)
}

colorize_delta_matrix <- function(vmat, cap) {
  h <- nrow(vmat); w <- ncol(vmat)
  v <- as.vector(vmat)
  vc <- pmax(-cap, pmin(cap, v))
  frac <- (vc + cap) / (2 * cap)
  idx <- round(frac * (length(DELTA_RAMP) - 1)) + 1
  idx[is.na(v)] <- NA
  rgb <- hex2rgb01(DELTA_RAMP)
  dist_frac <- pmin(1, abs(vc) / (cap * FADE_BAND))
  alpha_val <- (FADE_FLOOR + dist_frac * (255 - FADE_FLOOR)) / 255
  alpha <- ifelse(is.na(v), 0, alpha_val)

  img <- array(0, dim = c(h, w, 4))
  rch <- rep(0, length(v)); gch <- rep(0, length(v)); bch <- rep(0, length(v))
  ok <- !is.na(idx)
  rch[ok] <- rgb[idx[ok], 1]; gch[ok] <- rgb[idx[ok], 2]; bch[ok] <- rgb[idx[ok], 3]
  img[, , 1] <- matrix(rch, h, w); img[, , 2] <- matrix(gch, h, w)
  img[, , 3] <- matrix(bch, h, w); img[, , 4] <- matrix(alpha, h, w)
  img
}

write_delta <- function(scen, later_per, earlier_per, out_name, label) {
  later_png <- file.path(base_dir, scen, later_per, "probability_ensemble_mean.png")
  earlier_png <- file.path(base_dir, scen, earlier_per, "probability_ensemble_mean.png")
  if (!file.exists(later_png) || !file.exists(earlier_png)) {
    log_step("  SKIP:", scen, out_name, "(missing input)"); return(invisible())
  }

  later <- decode_prob_png(later_png)
  earlier <- decode_prob_png(earlier_png)
  if (later$width != earlier$width || later$height != earlier$height) {
    log_step("  SKIP:", scen, out_name, "(grid size mismatch between periods!)"); return(invisible())
  }

  base_floor <- pmax(earlier$values, PCT_BASE_FLOOR)
  pct_change <- (later$values - earlier$values) / base_floor * 100

  out_dir <- file.path(base_dir, scen, later_per)
  writePNG(colorize_delta_matrix(pct_change, PCT_CAP), file.path(out_dir, paste0(out_name, ".png")))

  # Bounds/dimensions copied VERBATIM from the later period's own meta.json
  # -- the same file the app already uses to place probability_ensemble_
  # mean.png, so this PNG (identical pixel grid) lands in exactly the same
  # place with no transform of any kind.
  meta_f <- file.path(base_dir, scen, later_per, "probability_ensemble_mean_meta.json")
  src_meta <- jsonlite::read_json(meta_f)
  v <- as.vector(pct_change)
  meta <- list(
    image = paste0(out_name, ".png"),
    bounds = src_meta$bounds,
    dimensions = src_meta$dimensions,
    source = paste0("ensemble ", scen, ": % change ", later_per, " vs ", earlier_per, " (", label, ")"),
    value_range = list(min = -PCT_CAP, max = PCT_CAP, unit = "percent"),
    actual_range = list(min = round(min(v, na.rm = TRUE), 1), max = round(max(v, na.rm = TRUE), 1)),
    baseline_period = earlier_per,
    baseline_floor = PCT_BASE_FLOOR,
    method = "decoded directly from sibling probability PNGs, no resampling"
  )
  jsonlite::write_json(meta, file.path(out_dir, paste0(out_name, "_meta.json")),
                        auto_unbox = TRUE, pretty = TRUE)
  log_step(scen, out_name, "(", later_per, "vs", earlier_per, ") -> % change range",
            round(min(v, na.rm = TRUE), 1), "to", round(max(v, na.rm = TRUE), 1))
}

for (scen in scenarios) {
  write_delta(scen, "2041_2070", "2011_2040", "delta_from_near", "current")
  write_delta(scen, "2071_2100", "2011_2040", "delta_from_near", "current")
}

log_step("Done.")
