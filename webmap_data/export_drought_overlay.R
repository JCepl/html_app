#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(terra)
  library(jsonlite)
})

# Convert Stefan Ebner's (BFW) two drought GeoTIFFs (CHELSA-native SPEI-12 /
# multi-year-drought indicators, ~1km, WGS84) into the app's standard
# web overlay assets (RGBA PNG + bounds/legend meta.json per band),
# same convention as export_web_overlay_from_tif.R, but classifying each
# band by EXACT categorical value (not continuous breaks) since these are
# already discrete class rasters.
#
# Definitions and layer selection per Stefan's 2026-08-25 notes
# (Stefan_DSS_data/Stefans_notes.docx): a multi-year-drought (MYD) event is
# SPEI-12 <= -1.1 for >=2 consecutive years, reference period 1981-2010.
# Stefan specified exactly 3 reference-period layers (not the 4 bands the
# source file has) -- "occurrence/presence" is deliberately dropped here,
# it's redundant with n_events_2yr >= 1 and not one of his 3 requested options.
#
# Source files (not served; sent by email, described in
# webmap_data/drought/drought_data_visualization.txt):
#   myd_SPEI-12_reference_DSS.tif        4 bands: presence_2yr (unused, see
#                                          above), n_events_2yr, max_duration,
#                                          total_duration
#   SPEI-equivalent-12_future_classes.tif 4 bands, all "class" (1-4), one per
#                                          future scenario/period. Band order
#                                          is NOT labelled in the source file;
#                                          inferred here from the matching
#                                          4-file drought bundle Stefan sent
#                                          in the same batch (band order =
#                                          2041-2070_SSP370, 2041-2070_SSP585,
#                                          2071-2100_SSP370, 2071-2100_SSP585).
#                                          Still flagged to Stefan to confirm.

SRC_DIR <- path.expand("~/REENFOCE_LOCAL_MODEL_WIEN/Stefan_DSS_data")
OUT_DIR <- "webmap_data/drought"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
MAX_DIM <- 900

to_rgba <- function(hex) {
  hex <- toupper(hex)
  if (grepl("^#[0-9A-F]{8}$", hex)) return(tolower(hex))
  if (grepl("^#[0-9A-F]{6}$", hex)) return(tolower(paste0(hex, "FF")))
  stop("bad color: ", hex)
}

# Render one band to RGBA PNG (exact value -> color lookup) + return bounds.
render_band <- function(r1, value_colors, out_png) {
  r1 <- project(r1, "EPSG:3857", method = "near")
  merc_ext <- ext(r1)
  nx <- ncol(r1); ny <- nrow(r1)
  scale_factor <- max(nx, ny) / MAX_DIM
  if (scale_factor > 1) {
    r1 <- aggregate(r1, fact = ceiling(scale_factor), fun = "modal", na.rm = TRUE)
    nx <- ncol(r1); ny <- nrow(r1)
  }
  vals <- values(r1, mat = FALSE)
  m <- matrix(vals, nrow = ny, ncol = nx, byrow = TRUE)

  rgba <- matrix("#00000000", nrow = ny, ncol = nx)
  for (nm in names(value_colors)) {
    v <- as.numeric(nm)
    class_mask <- !is.na(m) & abs(m - v) < 1e-6
    rgba[class_mask] <- to_rgba(value_colors[[nm]])
  }

  png(out_png, width = nx, height = ny, bg = "transparent")
  par(mar = c(0, 0, 0, 0), xaxs = "i", yaxs = "i")
  plot.new()
  rasterImage(as.raster(rgba), 0, 0, 1, 1, interpolate = FALSE)
  dev.off()

  corner_xy <- rbind(
    c(xmin(merc_ext), ymin(merc_ext)), c(xmin(merc_ext), ymax(merc_ext)),
    c(xmax(merc_ext), ymin(merc_ext)), c(xmax(merc_ext), ymax(merc_ext))
  )
  corner_ll <- project(corner_xy, from = "EPSG:3857", to = "EPSG:4326")
  list(
    south = unname(min(corner_ll[, 2], na.rm = TRUE)),
    west  = unname(min(corner_ll[, 1], na.rm = TRUE)),
    north = unname(max(corner_ll[, 2], na.rm = TRUE)),
    east  = unname(max(corner_ll[, 1], na.rm = TRUE))
  )
}

cat("== Reference (MYD) indicators ==\n")

ref_tif <- file.path(SRC_DIR, "myd_SPEI-12_reference_DSS.tif")
ref_images <- list()

ref_specs <- list(
  n_events_2yr = list(band = 2,
    colors = list(`0` = "#F7F7F7", `1` = "#fee391", `2` = "#fe9929", `3` = "#8c2d04"),
    labels = c("0 events", "1 event", "2 events", "3 events")),
  max_duration = list(band = 3,
    colors = list(`0` = "#F7F7F7", `2` = "#fee391", `3` = "#fe9929", `4` = "#cc4c02", `5` = "#8c2d04"),
    labels = c("No MYD", "2 years", "3 years", "4 years", "5 years")),
  total_duration = list(band = 4,
    colors = list(`0` = "#F7F7F7", `2` = "#fee391", `3` = "#fe9929", `4` = "#ec7014", `5` = "#cc4c02", `6` = "#8c2d04", `7` = "#8c2d04"),
    labels = c("No MYD", "2 years", "3 years", "4 years", "5 years", "6-7 years"))
)

for (key in names(ref_specs)) {
  spec <- ref_specs[[key]]
  cat(" -", key, "(band", spec$band, ")\n")
  r <- rast(ref_tif)[[spec$band]]
  out_png <- file.path(OUT_DIR, paste0("ref_", key, ".png"))
  bounds <- render_band(r, spec$colors, out_png)
  # de-dup color list for the legend (max/total_duration collapse 6&7 to one swatch)
  seen <- !duplicated(unlist(spec$colors, use.names = FALSE))
  ref_images[[key]] <- list(
    file = paste0("drought/ref_", key, ".png"),
    bounds = bounds,
    legend = list(
      labels = spec$labels,
      colors = unlist(spec$colors, use.names = FALSE)[seen]
    )
  )
}

cat("== Future SPEI-equivalent-12 classes ==\n")

fut_tif <- file.path(SRC_DIR, "SPEI-equivalent-12_future_classes.tif")
fut_keys <- c("2041_2070_ssp370", "2041_2070_ssp585", "2071_2100_ssp370", "2071_2100_ssp585")
fut_colors <- list(`1` = "#B2182B", `2` = "#EF8A62", `3` = "#FDDBC7", `4` = "#F7F7F7")
fut_labels <- c(
  "Extreme drought hazard (< -2.0)",
  "High drought hazard (-2.0 to -1.5)",
  "Moderate drought hazard (-1.5 to -1.0)",
  "Low drought hazard (>= -1.0)"
)
fut_images <- list()

for (i in seq_along(fut_keys)) {
  key <- fut_keys[[i]]
  cat(" -", key, "(band", i, ")\n")
  r <- rast(fut_tif)[[i]]
  out_png <- file.path(OUT_DIR, paste0("future_", key, ".png"))
  bounds <- render_band(r, fut_colors, out_png)
  fut_images[[key]] <- list(
    file = paste0("drought/future_", key, ".png"),
    bounds = bounds,
    legend = list(labels = fut_labels, colors = unlist(fut_colors, use.names = FALSE))
  )
}

write_json(list(
  source = "myd_SPEI-12_reference_DSS.tif (Stefan Ebner, BFW)",
  note = "Multi-year-drought (MYD) indicators, historical reference period 1981-2010, CHELSA-native ~1km. A MYD event is SPEI-12 <= -1.1 for >=2 consecutive years.",
  images = ref_images
), file.path(OUT_DIR, "drought_ref_meta.json"), pretty = TRUE, auto_unbox = TRUE)

write_json(list(
  source = "SPEI-equivalent-12_future_classes.tif (Stefan Ebner, BFW)",
  note = "Classified December SPEI-12-equivalent, future period, relative to the 1981-2010 reference climate. Band-to-scenario/period mapping is INFERRED (see script header) -- confirm with Stefan.",
  images = fut_images
), file.path(OUT_DIR, "drought_future_meta.json"), pretty = TRUE, auto_unbox = TRUE)

cat("\nDone. Wrote", length(ref_images) + length(fut_images), "PNG+meta layers to", OUT_DIR, "\n")
