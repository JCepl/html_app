#!/usr/bin/env Rscript

# Convert Tommaso's wind-damage-risk GeoTIFFs (EPSG:4326, prob 0-1, ~130 m/px,
# ~760 MB each) into lightweight web overlays for the DSS map:
#   - one RGBA PNG per period (smooth yellow->red gradient, saturated at 0.5)
#   - one shared wind_meta.json (Leaflet bounds + legend description)
#
# Visualization follows Tommaso's note: scale 0 (low risk) -> 0.5 (high risk);
# probabilities above 0.5 are clamped to the high-risk end. Low values fade in
# alpha so low-risk land does not blanket the basemap.
#
# Usage: Rscript webmap_data/export_wind_overlays.R

suppressPackageStartupMessages({
  library(terra)
  library(jsonlite)
})

# --- configuration ----------------------------------------------------------
src_dir   <- "/Users/jaroslavcepl/REENFOCE_LOCAL_MODEL_WIEN/Tommaso"
out_dir   <- "webmap_data/wind"
max_dim   <- 4000          # target longest side in pixels for the web PNG
risk_cap  <- 0.5           # saturate the color ramp at this probability
fade_knee <- 0.3           # low-risk end fades in across the bottom fraction
fade_floor <- 90           # minimum alpha at zero risk (stays faintly visible)

periods <- list(
  list(key = "historical", tif = "risk_of_wind_damage_historical_scenario_ssp585.tif"),
  list(key = "near",       tif = "risk_of_wind_damage_2041_2070_scenario_ssp585.tif"),
  list(key = "far",        tif = "risk_of_wind_damage_2071_2100_scenario_ssp585.tif")
)

# Intuitive sequential risk ramp (RdYlBu reversed): blue = low risk ->
# pale yellow = moderate -> red = high risk. Avoids green so it never clashes
# with the green land/forest on the basemap.
ramp_colors <- c("#14155f", "#2a45c2", "#3f7af0", "#6f74ee", "#9a5ee0",
                 "#c94fc0", "#ec4a86", "#ff7a1f", "#f5331a", "#8f0000")
# sRGB interpolation matches the CSS legend gradient and the JS ramp exactly.
color_ramp  <- colorRamp(ramp_colors, space = "rgb")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

bounds_out <- NULL

for (p in periods) {
  in_tif <- file.path(src_dir, p$tif)
  if (!file.exists(in_tif)) stop("Missing input raster: ", in_tif)

  cat("Processing", p$key, "->", p$tif, "\n")
  r <- rast(in_tif)

  # Downsample on the native grid first (cheap, windowed read) so we never
  # reproject the full 2.9-billion-cell raster. Use max so high-risk pixels
  # survive aggregation.
  fact <- ceiling(max(ncol(r), nrow(r)) / max_dim)
  if (fact > 1) {
    r <- aggregate(r, fact = fact, fun = "max", na.rm = TRUE)
  }

  # Reproject the small raster to Web Mercator so the imageOverlay aligns with
  # the EPSG:3857 basemap (same approach as export_web_overlay_from_tif.R).
  r <- project(r, "EPSG:3857", method = "bilinear")
  merc_ext <- ext(r)

  nx <- ncol(r); ny <- nrow(r)
  m  <- matrix(values(r, mat = FALSE), nrow = ny, ncol = nx, byrow = TRUE)
  finite_mask <- is.finite(m)

  # Normalize value to [0,1] over [0, risk_cap], clamp above.
  norm <- pmin(pmax(m, 0), risk_cap) / risk_cap

  rgba <- matrix("#00000000", nrow = ny, ncol = nx)
  if (any(finite_mask)) {
    vv  <- norm[finite_mask]
    rgb <- color_ramp(vv)                      # n x 3, 0-255
    # The lowest-risk (dark blue) end fades toward the basemap but never fully
    # disappears: alpha ramps from `fade_floor` at 0 up to 255 across the bottom
    # `fade_knee` of the scale, so even ~1% risk stays faintly visible.
    alpha <- round(fade_floor + (255 - fade_floor) * pmin(1, vv / fade_knee))
    hex <- sprintf("#%02X%02X%02X%02X",
                   round(rgb[, 1]), round(rgb[, 2]), round(rgb[, 3]), alpha)
    rgba[finite_mask] <- hex
  }

  out_png <- file.path(out_dir, paste0("wind_", p$key, ".png"))
  png(out_png, width = nx, height = ny, bg = "transparent")
  par(mar = c(0, 0, 0, 0), xaxs = "i", yaxs = "i")
  plot.new()
  rasterImage(as.raster(rgba), 0, 0, 1, 1, interpolate = FALSE)
  dev.off()

  if (is.null(bounds_out)) {
    corner_xy <- rbind(
      c(xmin(merc_ext), ymin(merc_ext)),
      c(xmin(merc_ext), ymax(merc_ext)),
      c(xmax(merc_ext), ymin(merc_ext)),
      c(xmax(merc_ext), ymax(merc_ext))
    )
    ll <- project(corner_xy, from = "EPSG:3857", to = "EPSG:4326")
    bounds_out <- list(
      south = min(ll[, 2]), west = min(ll[, 1]),
      north = max(ll[, 2]), east = max(ll[, 1])
    )
  }

  cat("  wrote", out_png, "(", nx, "x", ny, ")\n")
}

meta <- list(
  bounds = bounds_out,
  scenario = "RCP8.5 (SSP5-8.5)",
  risk_cap = risk_cap,
  images = list(
    historical = "wind/wind_historical.png",
    near = "wind/wind_near.png",
    far = "wind/wind_far.png"
  ),
  legend = list(
    min_label = "0 — low risk",
    max_label = "≥ 0.5 — high risk",
    colors = ramp_colors
  ),
  source = "Tommaso — risk_of_wind_damage_*_scenario_ssp585.tif",
  note = "Probability of wind damage 0-1 (threshold 0.5 = damage). Color saturates at 0.5."
)

write_json(meta, file.path(out_dir, "wind_meta.json"),
           pretty = TRUE, auto_unbox = TRUE, digits = 8)

cat("Done. Bounds:", bounds_out$west, bounds_out$south,
    bounds_out$east, bounds_out$north, "\n")
