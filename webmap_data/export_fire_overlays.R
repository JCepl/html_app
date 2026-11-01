#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(terra)
  library(jsonlite)
})

# Fire-weather layers for the DSS from the UNIPD-TESAF package (Flavio
# Taccaliti, 28 Oct 2026): Giannakopoulos & Karali 2020 (C3S), EURO-CORDEX
# EUR-11, RCA4 driven by MPI-ESM-LR, 0.15 degree grid, RCP2.6/4.5/8.5,
# historical 1981-2005, futures 2041-2070 and 2071-2098.
#
# Two views, both for summer (JJA):
#   change : relative change of the mean seasonal severity rating
#            (band delta_dsr_mean_JJA), classified in % bands. Cells where
#            |change| <= 1.96 * delta_dsr_se_rel_JJA (not distinguishable
#            from year-to-year noise, the partner's own rule) are drawn pale.
#   level  : mean Fire Weather Index (band fwi_mean_JJA), classified with
#            the EFFIS fire-danger classes.
# Bands are always selected BY NAME from the *_bands.csv (partner rule).
# Winter bands are not used (DJF year-count issue in the package).
#
# Usage: Rscript webmap_data/export_fire_overlays.R <folder with the unzipped geoTIFF.zip>
# Run from the repo root.

args <- commandArgs(trailingOnly = TRUE)
SRC_DIR <- if (length(args) >= 1) args[[1]] else
  "/Volumes/Extreme SSD/RE-ENFORCE_D211_DELIVERY_2026-10-31/2_Fire_Flavio/geoTIFF"
OUT_DIR <- "webmap_data/fire"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Europe window (same as the other driver layers; drops the Sahara corner).
EUROPE <- ext(-25, 45, 34, 72)

pick_band <- function(scen, name) {
  stem <- file.path(SRC_DIR, paste0("GiannakopoulosKarali2020_", scen))
  bands <- read.csv(paste0(stem, "_bands.csv"))
  idx <- bands$band[bands$name == name]
  if (length(idx) != 1) stop("band ", name, " not found once in ", scen)
  r <- crop(rast(paste0(stem, "_4326.tif"))[[idx]], EUROPE)
  names(r) <- name
  r
}

# Classes. Colours are listed low -> high; the app decodes clicks by exact colour.
LEVEL_BREAKS <- c(-Inf, 5.2, 11.2, 21.3, 38, 50, Inf)
LEVEL_COLORS <- c("#FFF7BC", "#FEE391", "#FEC44F", "#FB9A29", "#D95F0E", "#8C2D04")
LEVEL_LABELS <- c("Very low (< 5.2)", "Low (5.2-11.2)", "Moderate (11.2-21.3)",
                  "High (21.3-38)", "Very high (38-50)", "Extreme (> 50)")

CHANGE_BREAKS <- c(-Inf, -0.10, 0.10, 0.25, 0.50, 1.00, Inf)
CHANGE_COLORS <- c("#6BAED6", "#F0F0F0", "#FDD0A2", "#FDAE6B", "#E6550D", "#A63603")
CHANGE_LABELS <- c("Decrease (more than -10%)", "Little change (-10 to +10%)",
                   "+10 to +25%", "+25 to +50%", "+50 to +100%", "More than doubled (> +100%)")
MILD_COLOR <- "#D9D9D9"
MILD_LABEL <- "Fire weather too mild for a change to be computed"
PALE_ALPHA <- "66"   # alpha for changes within year-to-year noise

hex_alpha <- function(hex, a = "FF") paste0(toupper(hex), a)

to_merc <- function(r, method) {
  # 0.15 deg is ~10-17 km; project at ~5 km so cell edges stay crisp.
  project(r, "EPSG:3857", method = method, res = 5000)
}

write_png <- function(rgba, out_png) {
  png(out_png, width = ncol(rgba), height = nrow(rgba), bg = "transparent")
  par(mar = c(0, 0, 0, 0), xaxs = "i", yaxs = "i")
  plot.new()
  rasterImage(as.raster(rgba), 0, 0, 1, 1, interpolate = FALSE)
  dev.off()
}

merc_bounds <- function(r) {
  e <- ext(r)
  xy <- rbind(c(xmin(e), ymin(e)), c(xmax(e), ymax(e)))
  ll <- project(xy, from = "EPSG:3857", to = "EPSG:4326")
  list(south = ll[1, 2], west = ll[1, 1], north = ll[2, 2], east = ll[2, 1])
}

as_matrix <- function(r) matrix(values(r, mat = FALSE), nrow = nrow(r), ncol = ncol(r), byrow = TRUE)

render_level <- function(scen, out_png) {
  r <- to_merc(pick_band(scen, "fwi_mean_JJA"), "near")
  m <- as_matrix(r)
  cls <- cut(m, LEVEL_BREAKS, labels = FALSE, right = FALSE)
  rgba <- matrix("#00000000", nrow(m), ncol(m))
  ok <- !is.na(cls)
  rgba[ok] <- hex_alpha(LEVEL_COLORS[cls[ok]])
  write_png(rgba, out_png)
  list(bounds = merc_bounds(r), median = median(m, na.rm = TRUE))
}

render_change <- function(scen, out_png) {
  land <- pick_band("historical", "fwi_mean_JJA")
  d  <- pick_band(scen, "delta_dsr_mean_JJA")
  se <- pick_band(scen, "delta_dsr_se_rel_JJA")
  stack <- to_merc(c(land, d, se), "near")
  L <- as_matrix(stack[[1]]); D <- as_matrix(stack[[2]]); SE <- as_matrix(stack[[3]])
  rgba <- matrix("#00000000", nrow(D), ncol(D))
  mild <- !is.na(L) & is.na(D)
  rgba[mild] <- hex_alpha(MILD_COLOR)
  cls <- cut(D, CHANGE_BREAKS, labels = FALSE, right = FALSE)
  ok <- !is.na(cls)
  robust <- ok & !is.na(SE) & abs(D) > 1.96 * SE
  rgba[ok & robust]  <- hex_alpha(CHANGE_COLORS[cls[ok & robust]])
  rgba[ok & !robust] <- hex_alpha(CHANGE_COLORS[cls[ok & !robust]], PALE_ALPHA)
  write_png(rgba, out_png)
  list(bounds = merc_bounds(stack),
       robust_share = round(sum(robust) / sum(ok), 3),
       median_change = round(median(D, na.rm = TRUE), 3))
}

FUTURES <- c("rcp26_2041_2070", "rcp26_2071_2098", "rcp45_2041_2070",
             "rcp45_2071_2098", "rcp85_2041_2070", "rcp85_2071_2098")

images <- list()
cat("level historical\n")
info <- render_level("historical", file.path(OUT_DIR, "level_historical.png"))
images[["level::historical"]] <- list(file = "fire/level_historical.png", bounds = info$bounds)
for (s in FUTURES) {
  cat("level", s, "\n")
  info <- render_level(s, file.path(OUT_DIR, paste0("level_", s, ".png")))
  images[[paste0("level::", s)]] <- list(file = paste0("fire/level_", s, ".png"), bounds = info$bounds)
  cat("change", s, "\n")
  info <- render_change(s, file.path(OUT_DIR, paste0("change_", s, ".png")))
  images[[paste0("change::", s)]] <- list(file = paste0("fire/change_", s, ".png"), bounds = info$bounds,
                                          robust_share = info$robust_share, median_change = info$median_change)
}

write_json(list(
  source = "Giannakopoulos & Karali 2020 (C3S fire danger indicators), EURO-CORDEX EUR-11 RCA4 / MPI-ESM-LR; aggregated by F. Taccaliti (UNIPD-TESAF), package of 28 Oct 2026",
  grid = "0.15 degree (about 15 km)",
  reference_period = "1981-2005",
  legends = list(
    level  = list(labels = LEVEL_LABELS, colors = LEVEL_COLORS),
    change = list(labels = c(CHANGE_LABELS, MILD_LABEL), colors = c(CHANGE_COLORS, MILD_COLOR),
                  pale_note = "Pale colours: change smaller than the year-to-year noise (|change| <= 1.96 x standard error).")
  ),
  images = images
), file.path(OUT_DIR, "fire_meta.json"), pretty = TRUE, auto_unbox = TRUE, digits = 6)

cat("Done:", length(images), "layers in", OUT_DIR, "\n")
