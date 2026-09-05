#!/bin/bash
# Rebuilds fire/wind disturbance-pressure hex GeoJSONs the same way bbtl's
# were originally built: filter agent+metric='rate', drop null-mean
# (non-forest) cells, clamp tiny negative numerical-smoothing artifacts to
# 0, split by rcp x year, reproject EPSG:3035 (LAEA) -> EPSG:4326 (CRS84),
# keep only rcp/year/mean properties to match the existing bbtl_*.geojson
# schema. Source: Grunig et al. 2026, Science -- Dryad doi:10.5061/dryad.tb2rbp0dv,
# file 09_svd_simulations.tar.gz, processed_data/dist_rates_50km_shinyapp/
# (the paper's own pre-aggregated 50km hex product, not the raw 100m rasters
# in raw/ -- that's what bbtl was built from too; see git history commit
# e613dfe and REENFOCE_LOCAL_MODEL_WIEN memory for how this was traced).
# Requires that Dryad file to be downloaded and extracted locally first
# (Anubis bot-check blocks curl/wget -- use a real browser).
set -euo pipefail

SRC_BASE="/Volumes/Extreme SSD/Dryad_bark_beetle/09_svd_simulations/09_svd_simulations/processed_data/dist_rates_50km_shinyapp"
OUT_DIR="$(dirname "$0")/geojson"
YEARS=(2030 2040 2050 2060 2070 2080 2090 2100)
RCPS=(rcp26 rcp45 rcp85)

for agent in fire wind; do
  GPKG="$SRC_BASE/${agent}_dist_rate_50km.gpkg"
  LAYER="${agent}_dist_rate_50km"
  for rcp in "${RCPS[@]}"; do
    for year in "${YEARS[@]}"; do
      OUT="$OUT_DIR/${agent}_${rcp}_${year}.geojson"
      ogr2ogr -f GeoJSON -t_srs EPSG:4326 -dialect SQLite \
        -sql "SELECT rcp, year, MAX(mean,0) AS mean, geom FROM \"$LAYER\" WHERE rcp='${rcp}' AND year=${year} AND metric='rate' AND mean IS NOT NULL" \
        "$OUT" "$GPKG"
      echo "wrote $OUT"
    done
  done
done
