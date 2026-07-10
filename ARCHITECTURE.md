# RE-ENFORCE DSS — Architecture & Context

> Context document for developers and LLM agents working on this app.
> Last updated 2026-06-28. If you change structure, update this file.

## 1. What this is

A **single-file, static web map** — the **RE-ENFORCE Decision Support System** — that
visualises **forest-disturbance risk across Europe** under climate scenarios. It is
deployed to **GitHub Pages** (live at `jcepl.github.io/html_app`). No backend, no build
step: open `index.html` and it runs.

It is a **multi-stress** DSS. Stresses are grouped in the top ribbon:
- **Bark beetle** (the mature stress, several layers — see below)
- **Wind** (forest wind-damage probability)
- **Wildfire**, **Drought** — placeholders ("Coming soon")

## 2. Files

```
html_app/
  index.html              # THE APP — all HTML + CSS + JS in one file (~3.3k lines, ~116 KB)
  ARCHITECTURE.md         # this file
  README.md               # project blurb
  webmap_data/            # all map data + the R scripts that build it
    map_meta.js           # global: bbox, value range, bbtl (Science) hex geojson file list
    wind/                 # wind_<period>.png + wind_meta.json (RCP8.5, 3 periods)        ~11 MB
    phenips/              # phenips_*.png + phenips_data.json (GDD / generations)          ~4 MB
    forecast_bundle/      # bark-beetle forecast: manifest.json + per-model/scenario/period ~46 MB
    geojson/              # bbtl pressure hexes (Grünig et al. 2026 Science) per rcp+year
    bbt_past_disturbance.png(+meta)   # classified historical disturbance
    *.R                   # data-build scripts (see §6)
```

Source rasters (NOT served; inputs to the R scripts) live in the parent repo
`/Users/jaroslavcepl/REENFOCE_LOCAL_MODEL_WIEN/` (e.g. `Tommaso/` wind TIFFs,
`FINAL_MODEL/` forecast model outputs).

## 3. Tech stack (all via unpkg CDN, no bundler)

- **Leaflet 1.9.4** — the map.
- **leaflet.sync** — keeps the two split-screen maps' pan/zoom in lock-step.
- **leaflet-control-geocoder** — the place-search box (Nominatim).
- **@turf/turf** — point-in-polygon / bbox for click inspection & country filtering.
- `webmap_data/map_meta.js` — injects the global `MAP_META` constant.

(Earlier versions loaded `geotiff`/`georaster`/`leaflet-side-by-side`; these were removed
once everything moved to pre-rendered PNGs + the dual-map compare.)

## 4. The layers (`OVERLAY_CATALOG` in index.html)

| id | UI name | Type | Notes |
|----|---------|------|-------|
| `bbtforecast` | Predictions | **PNG** raster (per model·scenario·period·view) | bark-beetle forecast; the core layer. Defaults to the **ensemble** model. |
| `phenips` | PHENIPS climate | **PNG** raster | GDD / generations, 4 periods × scenarios × views |
| `bbtpast` | Past disturbance | **PNG** (classified) | historical 1985–2023 footprint |
| `bbtl` | Bark Beetle Pressure | **vector hexes (GeoJSON)** | Grünig et al. 2026, *Science*; per RCP × year |
| `wind` | Wind | **PNG** raster | wind-damage probability, RCP8.5, 3 periods |
| `wildfire`,`drought` | — | not available | placeholders |

All raster layers are **pre-rendered RGBA PNGs + a `*_meta.json`** (bounds, value range)
and drawn as `L.imageOverlay`. The only vector layer is `bbtl`.

In the ribbon, the bark-beetle layers are one group: only **Predictions** shows until the
user clicks the **"Bark beetle ▸"** label, which reveals the accessories (PHENIPS, Past,
Pressure). See `#bbToggle` / `#bbAccessories` and the `.bb-*` CSS.

## 5. JavaScript structure (inside the single `<script>`)

State lives in one object **`APP`**. Functions are grouped under banner comments
(`// ==== SECTION ====`). Sections, in order:

1. **STATE & CONSTANTS** — `APP`, `RISK_RAMP`/`RISK_GRADIENT` (the colour scale),
   `OVERLAY_CATALOG`, `DOM` (cached element refs), `ASSET_VERSION`.
2. **RISK COLOUR SCALE, THRESHOLD SLIDER & PNG MASKING** — `rampColor`, the
   `RISK_LUT` + `colorToRiskFraction` (recover a pixel's scale position from its colour),
   `loadThresholdedImage` (canvas masking), `riskFractionFor`/`effectiveRiskFraction`,
   `applyRiskThreshold`.
3. **OVERLAY REGISTRY HELPERS** — visibility (`visibleOverlays`/`showOverlay`/…), panes,
   `addImageOverlayLayer`, load tokens.
4. **BARK-BEETLE FORECAST BUNDLE** — manifest/entry/metrics loaders, `getForecastOverlayPaths`.
5. **PHENIPS & WIND RASTER LAYERS** — `updatePhenipsLayer`, `updateWindLayer`, etc.
6. **LEGEND** — `setLegendLabels` (per-layer legend + gradient + threshold marker).
7. **HEX / COUNTRY STYLING, CLICK MARKER & READOUT** — `getColor` (bbtl hex colour),
   `baseHexStyle`, country highlight, `fitDataBounds`.
8. **LAYER LOADERS & CONTROL-PANEL SYNC** — `OVERLAY_LOADERS`, `updateOverlayLayer`,
   `updateForecastOverlay`, `updateHexLayer`, `syncOverlayControls`.
9. **CLICK INSPECTION & POPUP** — `updatePanelFromClick`, `openClickPopup`.
10. **BASEMAPS** — `BASEMAPS` (osm/light/dark/plain), `setBasemap`, `addBasemapControl`.
11. **SPLIT SCREEN** — see §5.1.
12. **INIT, SHAREABLE LINK & CONTROL WIRING** — `initMap`, `buildShareUrl`/`applyStateFromUrl`,
    `initControls`, `loadCountryBorders`.

### 5.1 Split screen (the most important concept)

The "Split screen" button (top banner) shows **two real Leaflet maps** side by side.
The control panel always drives the **active** map; clicking a map makes it active
(amber border). Each map keeps its own independent settings.

Implementation = **state swapping**, not a second pipeline:
- All **per-map** fields are listed in `PANEL_FIELDS` (`map`, `layerInstances`, basemap,
  country layers, all `current*` settings, `riskThreshold`, …). They live directly on
  `APP` and describe the **active** map only.
- `savePanel(side)` parks the current `APP.*` fields into `APP.panels[side]`;
  `loadPanel(side)` copies the other side's parked state back onto `APP.*`.
- So when you switch sides, `APP.map`/`APP.layerInstances`/etc. point at the other map,
  and the **existing single-map rendering pipeline transparently drives whichever map is
  active**. No per-layer special-casing.
- `leaflet.sync` links the two physical maps' pan/zoom.

This is why "settings affect only the active map" works for every layer including the
forecast and the vector hexes.

## 6. Data pipeline (R scripts in `webmap_data/`)

All map rasters are **pre-rendered in R (terra)** to RGBA PNG + `*_meta.json`, then served
statically. The browser never reads TIFFs.

| script | builds |
|--------|--------|
| `export_wind_overlays.R` | wind PNGs from the big `Tommaso/*.tif` (downsample → Web-Mercator → RISK_RAMP, alpha fade) |
| `colorize_prob_continuous.R` | forecast **probability** PNGs (continuous RISK_RAMP) — re-run this for the probability view |
| `export_forecast_bundle_overlays.R` | forecast **thresholded risk** PNGs (classified). Its *probability* path is **superseded** by the script above |
| `export_web_overlay_from_tif.R` | generic classified TIF→PNG exporter (used by the others) |
| `build_ensemble_rcp85.R` | builds the RCP8.5 forecast **ensemble** (multi-model mean of xgb/rf/gam) — probability + threshold layers + entry.json |
| `build_gdd_*.R` | PHENIPS GDD/generation rasters |

### Forecast bundle layout
`forecast_bundle/manifest.json` defines `scenarios`, `periods`, `thresholds`, `models`
(each with `label`, `periods`, `metrics_file`), and `entries` keyed **`"<model>::<period>"`**.
Each entry points to `ensemble_probability_tif` and a `threshold_visible_tifs` map
(`"0.5"→…_t50.tif`). Folders: `forecast_bundle/<model_or_ensemble>/<scenario>/<period>/`.

**Models** (ensembles listed first; the app defaults to `ensemble_rcp45`):
`ensemble_rcp26/45/85`, `xgb_v2_rcp26/45`, `xgb_rcp85`, `rf_rcp{26,45,85}`, `gam_rcp{26,45,85}`.
The **ensemble** is a GCM-mean (rcp26/45) or a multi-model mean of available rcp85 models
(rcp85, built here). Evaluation metrics are historical/scenario-independent, so all
ensembles reuse `ensemble/metrics_confusion.json`.

### CRITICAL invariant: the colour ramp must stay in sync
`RISK_RAMP` (index.html JS), the `.legend-bar` CSS gradient (index.html), and `ramp_colors`
/ `prob_colors` in **every** R script must be the **same 10 hex stops**. The current ramp is
a "bold cool→hot": `#14155f #2a45c2 #3f7af0 #6f74ee #9a5ee0 #c94fc0 #ec4a86 #ff7a1f #f5331a #8f0000`.
Low-risk end fades but keeps a floor (`fade_floor`/`RISK_FADE_FLOOR`) so ~1% risk stays
faintly visible.

### Cache busting
Raster PNG URLs carry `?v=ASSET_VERSION`. **Bump `ASSET_VERSION` in index.html whenever you
regenerate PNGs**, or browsers serve stale images. (Current: `20260628-floor`.)

## 7. Notable UI features

- **Basemap switcher** (top-right, collapsed/subtle, expands on hover): Streets / Light /
  Dark / **Plain (borders only)**.
- **"Hide risk below X%"** slider (Legend card): masks low-risk areas. Threshold is an
  **absolute probability**; each layer's colour scale tops out at a different probability
  (wind saturates at 0.5 → its slider caps at 49%), so `riskFractionFor` converts the % to
  the right colour-scale fraction per layer. PNG layers are masked client-side via canvas
  (`loadThresholdedImage`); the hex layer via `fillOpacity`.
- **Click popup** — small bubble (`mini-popup`) above the borders; mirrors the side panel.
- **Share view** (banner) — encodes the full view into the URL hash (`buildShareUrl` /
  `applyStateFromUrl`): basemap, visible/focused layers, forecast & PHENIPS & wind settings,
  RCP/year, threshold.

## 8. Running & regenerating

- **Run locally:** `cd html_app && python3 -m http.server 8765` → open
  `http://localhost:8765/index.html`. (It's static; any web server works. GitHub Pages
  serves it as-is.)
- **Validate JS without a browser** (no Node here): extract the inline `<script>` and
  parse it, e.g. on macOS `osascript -l JavaScript -e 'Function(<code>)'` throws on syntax
  errors. Always do this after editing — it's one file and easy to break.
- **Regenerate data:** run the relevant R script in `webmap_data/` (terra + jsonlite
  required), then **bump `ASSET_VERSION`**.

## 9. Conventions / gotchas

- **One file.** No modules; functions are globals grouped by banner comments. Keep related
  code in its section.
- **Pre-rendered, not client-computed.** Don't reintroduce client-side raster reading;
  add layers as PNG + meta produced by an R script.
- **Ramp sync** (see §6) and **ASSET_VERSION bump** are the two easiest things to forget.
- The bark-beetle **pressure hexes (`bbtl`) can't be pinned into split-screen image
  panels** (vector); it works as a normal layer.
- Coordinates/extent: data is EPSG:4326; PNGs are reprojected to Web-Mercator at build time
  so `L.imageOverlay` aligns with the basemap.
