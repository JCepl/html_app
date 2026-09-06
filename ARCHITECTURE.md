# RE-ENFORCE DSS — Architecture & Context

> Context document for developers and LLM agents working on this app.
> Last updated 2026-06-28. If you change structure, update this file.

## 1. What this is

A **single-file, static web map** — the **RE-ENFORCE Decision Support System** — that
visualises **forest-disturbance indicators across Europe** under climate scenarios. It is
deployed to **GitHub Pages** (live at `jcepl.github.io/html_app`). No backend, no build
step: open `index.html` and it runs.

It is a **multi-stress** DSS. Stresses are grouped in the top ribbon:
- **Bark beetle** (the mature stress, several layers — see below)
- **Wind** (forest wind-damage susceptibility, + a Grunig et al. pressure-index accessory)
- **Drought** (SPEI-based, 2 layers — see below)
- **Wildfire** (Grunig et al. pressure index — this project has no in-house wildfire model)
- **Ash dieback** — placeholder ("Coming soon")

**Terminology (2026-09, consortium feedback):** avoid the word "risk"/"probability" as a
blanket label — most layers are a susceptibility index, a hazard-occurrence estimate, or
some blend of the two (varies per driver and per model; see each `OVERLAY_CATALOG[id]
.description`), not a single calibrated risk probability. Only the bark-beetle forecast
(a genuinely calibrated, spatial-block-CV-validated XGBoost classifier) earns the word
"probability" in its own copy. Everywhere else in the UI, prefer "pressure index",
"susceptibility", "model forecast", or "disturbance indicator". If you add copy, match
this — it's a live, disclosed decision (see §11), not an incidental style choice.

## 2. Files

```
html_app/
  index.html              # THE APP — all HTML + CSS + JS in one file (~3.3k lines, ~116 KB)
  ARCHITECTURE.md         # this file
  README.md               # project blurb
  webmap_data/            # all map data + the R scripts that build it
    map_meta.js           # MAP_META (bbtl) + MAP_META_FIRE + MAP_META_WIND: bbox, value
                           # range, years/rcps, hex geojson file list -- one object per
                           # hex layer, see HEX_OVERLAY_IDS in index.html
    wind/                 # wind_<period>.png + wind_meta.json (RCP8.5, 3 periods)        ~11 MB
    drought/              # drought_{future,ref}_meta.json + classified PNGs (Stefan Ebner, BFW) ~0.5 MB
    phenips/              # phenips_*.png + phenips_data.json (GDD / generations)          ~4 MB
    forecast_bundle/      # bark-beetle forecast: manifest.json + per-model/scenario/period ~46 MB
    geojson/              # pressure-index hexes (Grunig et al. 2026 Science), per rcp+year:
                           #   bbtl_*.geojson (bark beetle), fire_*.geojson (wildfire),
                           #   wind_*.geojson (wind) -- 24 files each (3 rcp x 8 years)
    bbt_past_disturbance.png(+meta)   # classified historical disturbance
    build_fire_wind_hex_geojsons.sh   # rebuilds fire_*/wind_*.geojson from the Dryad gpkg
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
| `windPressure` | Wind Pressure | **vector hexes (GeoJSON)** | Grünig et al. 2026, *Science*; per RCP × year; accessory of the Wind group |
| `droughtFuture` | Future SPEI classes | **PNG** raster (classified) | December SPEI-equivalent-12, 4 classes, 4 scenario/period combos |
| `droughtRef` | Reference indicators (MYD) | **PNG** raster (classified) | multi-year-drought indicators, 1981-2010 reference period, 3 selectable metrics |
| `wildfire` | Wildfire Pressure | **vector hexes (GeoJSON)** | Grünig et al. 2026, *Science*; per RCP × year; this project's only wildfire layer |
| `ashDieback` | — | not available | placeholder |

All raster layers are **pre-rendered RGBA PNGs + a `*_meta.json`** (bounds, value range)
and drawn as `L.imageOverlay`. The vector (hex) layers are `bbtl`, `wildfire`, and
`windPressure` — see `HEX_OVERLAY_IDS` in index.html; all three share one code path
(styling, click sampling, RCP/year control panel, citation box) parameterized by
overlay id, and more than one can be visible/cached at once (`APP.hexGeojson[id]`,
`APP.layerInstances[id]`). To add a fourth: add an `OVERLAY_CATALOG` entry, a
`MAP_META_*` object in map_meta.js, a ribbon card, and its id to `HEX_OVERLAY_IDS`
— no other JS changes needed.

In the ribbon, the bark-beetle layers are one group: only **Predictions** shows until the
user clicks the **"Bark beetle ▸"** label, which reveals the accessories (PHENIPS, Past,
Pressure). Wind works the same way: only **Wind** (the raster) shows until **"Wind ▸"**
(`#windToggle`/`#windAccessories`) reveals **Wind Pressure**. Wildfire has only the one
pressure-index card and no accessories to reveal, so its group has no toggle — same as
Wind's group had before Wind Pressure was added. See `.bb-*` CSS for the shared
toggle/accessories pattern.

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
| `export_drought_overlay.R` | Drought layers from Stefan Ebner's (BFW) 2 multi-band GeoTIFFs (`~/REENFOCE_LOCAL_MODEL_WIEN/Stefan_DSS_data/`) — 7 of the 8 bands used (occurrence/presence dropped per Stefan's spec), classified by exact categorical value (not continuous breaks) into `drought/`. **The future-classes band→scenario/period mapping is inferred, not embedded in the source file — confirm with Stefan before treating it as authoritative.** |

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

### Colour scales: one per driver, centralized (2026-08-28 refactor)
Each hazard driver owns its own scale, not a single scale shared app-wide. This
replaced an actual bug: the old single "shared" ramp was hand-copied into 9+ R
scripts, the JS, and the CSS separately, and two of those copies had already
silently drifted apart (the CSS `.legend-bar` gradient didn't match the JS/R
ramp it was rendering data for; `export_forecast_bundle_overlays.R`'s
`prob_colors` had drifted from `ramp_colors` elsewhere).

- **R side**: `webmap_data/color_scales.R` is the one canonical source —
  `WIND_RAMP` (Tommaso's own QGIS choice — the original ramp, historically
  reused everywhere else too), `BARK_BEETLE_RAMP` +
  `BARK_BEETLE_CLASSIFIED_COLORS`/`_BREAKS`. Every export/build script
  `source()`s this file instead of retyping hex lists. Drought (Stefan Ebner)
  and Past-disturbance already had their own independent palettes defined
  where they're used — nothing to centralize for those.
- **JS side**: `WIND_RAMP`/`BARK_BEETLE_RAMP` are separate arrays (same
  values as their R counterparts — keep them in sync by eye, there's no build
  step to do it automatically). `buildGradientCss(ramp)` generates the CSS
  gradient string from whichever ramp is active, so there is no second
  hand-typed gradient copy to drift. `WIND_LUT`/`colorToWindFraction` and
  `BARK_BEETLE_LUT`/`colorToRiskFraction` are separate reverse-color-match
  LUTs (see `sampleImageAt`/`sampleRasterOverlayFor`) — **use the one that
  matches the layer's actual palette**, or click-sampling (radar chart,
  threshold slider) silently nearest-matches against the wrong ramp.
- **`.legend-bar`'s CSS `background`** is a plain neutral placeholder,
  never actually shown — `setLegendLabels()` always sets an explicit inline
  background before the bar becomes visible.

Low-risk end fades but keeps a floor (`fade_floor`/`RISK_FADE_FLOOR`) so ~1% risk stays
faintly visible.

### Cache busting
Raster PNG URLs carry `?v=ASSET_VERSION`. **Bump `ASSET_VERSION` in index.html whenever you
regenerate PNGs**, or browsers serve stale images. (Current: `20260628-floor`.)

## 7. Notable UI features

- **Basemap switcher** (top-right, collapsed/subtle, expands on hover): Streets / Light /
  Dark / **Plain (borders only)**.
- **"Hide values below X%"** slider (Legend card): masks the low end of the displayed
  layer's own colour scale. Threshold is an **absolute value on that scale**; each layer's
  scale tops out at a different point (wind saturates at 0.5 → its slider caps at 49%), so
  `riskFractionFor` converts the % to the right colour-scale fraction per layer. PNG layers
  are masked client-side via canvas (`loadThresholdedImage`); hex layers via `fillOpacity`.
- **Click popup** — small bubble (`mini-popup`) above the borders; mirrors the side panel.
- **Share view** (banner) — encodes the full view into the URL hash (`buildShareUrl` /
  `applyStateFromUrl`): basemap, visible/focused layers, forecast & PHENIPS & wind settings,
  RCP/year, threshold, Simple/Advanced mode (`sv`).

## 8. Simple / Advanced view

**Default = Simple** (`APP.simpleView = true`), toggled by the banner's **"Show full
detail" / "Simple view"** button (`#detailModeBtn`, wired in `initControls`) and applied by
`applySimpleView()`. Simplifies the surface for a non-scientist audience, per 2026-09
consortium feedback (Anna Wöhlbrandt: *"think less of a scientist, more like a teacher or
journalist... I would totally omit picking climate models or ensemble models... include one
(or an average) and that's it"*). Nothing is deleted — Advanced restores the identical full
set instantly, no reload.

**Scope is deliberately narrow: model/GCM-ensemble complexity only, never time period.**
An earlier version of this also hid the far-future/end-of-century options (wind's `far`,
PHENIPS's `end`) and capped the pressure-index hex layers' year slider at 2070 — corrected
2026-09-06 per explicit user instruction ("DONT HIDE END OF CENTURY"): showing what changes
by end of century is the whole point of a climate-impact DSS, so every time period stays
available in **both** modes. Do not reintroduce a far-future/end-of-century cut here.

What Simple actually hides/forces:
- **Bark beetle forecast Model picker** (`#forecastModelField`) + **Historical Evaluation**
  metrics/confusion-matrix (`#forecastMetricsBox`) — hidden; `currentForecastModel` is
  pinned to `SIMPLE_VIEW_DEFAULT_FORECAST_MODEL` (`chelsa_fixed25_monotone`), the one
  variant with both a real held-out CV and a construction-guaranteed non-decreasing future
  response (see `CHELSA_NATIVE_TRAINING/PROJECT_REFERENCE_bark_beetle_model.md` §5–7). Its
  Future-period dropdown (reference / end-century) is untouched in both modes.
- **Wind Climate-model (GCM) picker** (`#windGcmField`) — hidden. Wind's Future-period
  dropdown (historical/near/far) is untouched in both modes.
- **RCP scenario, PHENIPS period, and the pressure-index hex layers' (`bbtl`/`wildfire`/
  `windPressure`) year slider are never restricted** — same reasoning as time period above:
  these are meaningful choices (policy scenario; what year to look at), not incidental
  model-ensemble clutter.

Toggling calls `applySimpleView()` again (idempotent, safe to call repeatedly).

## 9. Period-coverage strip

A 4-chip strip under every driver group's header (`.period-strip`/`.period-chip`),
added 2026-09 per user request for "clear structure -- stated empty if its empty":
Historic / Contemporary / Near future / Later future, one chip per slot, always
present so a gap reads as a gap rather than a silent absence.

- **Green chip** = a real layer covers that slot (hover for exactly which layer
  and years); **dashed grey `.is-empty`** = honestly no coverage.
- **Click-to-jump:** any non-empty chip is also a link. `data-jump-overlay` /
  `data-jump-period` attributes drive `jumpToOverlayPeriod()`, which shows that
  overlay, sets its period/year to match the chip, expands its accessory panel if
  collapsed (`expandAccessoryFor()` / `ACCESSORY_PARENT_GROUP`), and scrolls its
  card into view. `period` is a semantic slot (`historical`/`near`/`far`),
  translated per overlay since each one's own controls differ: drought's period
  selector also encodes a scenario (matched by year-range prefix, e.g.
  `2041_2070`), the pressure-index hex layers (`wildfire`/`bbtl`) use a literal
  year (2030/2100), and `bbtforecast`'s "later future" resolves to whichever
  period is chronologically last for the currently-selected model (fetched from
  the forecast manifest via `ensureForecastManifestLoaded()`) — the flagship
  model has no near-future run, so its strip only ever exposes Historic (via
  `bbtpast`, a different overlay) and Later future.
- Honest gaps this audit surfaced, worth not re-deriving: **no driver has a
  "Contemporary" (now) layer at all**; **Wildfire has no Historic layer**
  (Grünig's pre-aggregated 50km hex product ships future years only — a real one
  needs the raw 100m Dryad tarball, see the project memory's Dryad-archaeology
  note); **the bark-beetle forecast itself has no near-future run** (reference +
  end-century only — the Pressure hex stands in at the group level).

## 10. Running & regenerating

- **Run locally:** `cd html_app && python3 -m http.server 8765` → open
  `http://localhost:8765/index.html`. (It's static; any web server works. GitHub Pages
  serves it as-is.)
- **Validate JS without a browser** (no Node here): extract the inline `<script>` and
  parse it, e.g. on macOS `osascript -l JavaScript -e 'Function(<code>)'` throws on syntax
  errors. Always do this after editing — it's one file and easy to break.
- **Regenerate data:** run the relevant R script in `webmap_data/` (terra + jsonlite
  required), then **bump `ASSET_VERSION`**.

## 11. Conventions / gotchas

- **One file.** No modules; functions are globals grouped by banner comments. Keep related
  code in its section.
- **Pre-rendered, not client-computed.** Don't reintroduce client-side raster reading;
  add layers as PNG + meta produced by an R script.
- **Ramp sync** (see §6) and **ASSET_VERSION bump** are the two easiest things to forget.
- The pressure-index hex layers (`bbtl`, `wildfire`, `windPressure`) **can't be pinned
  into split-screen image panels** (vector); they work as normal layers.
- RCP/year are **global, shared controls**: changing either refreshes every currently
  *visible* hex layer (`HEX_OVERLAY_IDS.forEach` in the selector handlers), not just the
  focused one — so stacking e.g. `bbtl` + `wildfire` and dragging the year slider moves
  both.
- Coordinates/extent: data is EPSG:4326; PNGs are reprojected to Web-Mercator at build time
  so `L.imageOverlay` aligns with the basemap.
- **No "risk"/"probability" as a blanket label** in new copy — see §1's terminology note
  and §8 (Simple/Advanced). Say what the layer actually is: susceptibility, pressure
  index, model forecast, or (only where genuinely true) probability.
