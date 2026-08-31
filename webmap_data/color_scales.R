# One canonical color-scale definition per hazard driver.
#
# Source this file instead of retyping hex lists in each export script --
# change a driver's colors here and every script that builds its layers
# picks it up automatically. This replaces color definitions that used to
# be pasted independently into 8+ scripts (and had already drifted out of
# sync in at least one of them -- export_forecast_bundle_overlays.R's
# prob_colors no longer matched the "unified" ramp it claimed to mirror).
#
# Naming: <DRIVER>_RAMP = continuous probability ramp (sRGB-interpolated);
# <DRIVER>_CLASSIFIED_COLORS/_BREAKS = the discrete/thresholded variant.
# Drought (Stefan Ebner, BFW) and Past-disturbance already have their own
# independent classified palettes defined where they're used
# (export_drought_overlay.R, export_web_overlay_from_tif.R call sites) --
# not duplicated here, nothing to centralize for those.

# --- Wind: Tommaso's own QGIS rendering choice. This ramp originated as
# his wind color scheme and was, historically, reused for bark beetle too
# (see BARK_BEETLE_RAMP below) -- kept here as wind's own, canonical.
WIND_RAMP <- c(
  "#5e4fa2", "#3288bd", "#66c2a5", "#abdda4", "#e6f598",
  "#ffffbf", "#fee08b", "#fdae61", "#f46d43", "#d53e4f", "#9e0142"
)

# --- Bark beetle: continuous probability ramp for the forecast/ensemble
# layers. An independent copy of WIND_RAMP's values (visually identical for
# now) -- the point of splitting it out is that it can diverge from wind's
# scale without editing wind, and vice versa.
BARK_BEETLE_RAMP <- c(
  "#5e4fa2", "#3288bd", "#66c2a5", "#abdda4", "#e6f598",
  "#ffffbf", "#fee08b", "#fdae61", "#f46d43", "#d53e4f", "#9e0142"
)

# Bark beetle -- classified/thresholded view (8 discrete bins over the
# 0.55-1.0 "elevated risk" range), resampled from BARK_BEETLE_RAMP.
BARK_BEETLE_CLASSIFIED_BREAKS <- c(0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9)
BARK_BEETLE_CLASSIFIED_COLORS <- c(
  "#00000000", "#5e4fa2", "#48a1b3", "#a1d9a4",
  "#edf8a3", "#fee99a", "#fca55d", "#e2524a", "#9e0142"
)

# Helper: several scripts pass colors to another script/CLI as a comma
# string rather than an R vector.
ramp_to_csv <- function(x) paste(x, collapse = ",")
