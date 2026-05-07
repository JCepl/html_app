# RE-ENFORCE Decision Support Web Map

This repository contains a static Leaflet web app for exploring RE-ENFORCE forest disturbance layers on GitHub Pages.

Live app:

- https://jcepl.github.io/html_app/

The app is intentionally static: it has no build step, no server backend, and no package manager. `index.html` loads JavaScript libraries from CDNs and reads local data files from `webmap_data/`.

## Main Files

- `index.html`: the full app shell, styles, controls, and Leaflet logic.
- `webmap_data/map_meta.js`: metadata for the hex-pressure GeoJSON files, including years, RCP files, and value ranges.
- `webmap_data/geojson/`: projected bark beetle pressure by RCP scenario and year.
- `webmap_data/bbt_past_disturbance.png`: historical bark beetle disturbance image overlay.
- `webmap_data/bbt_past_disturbance_meta.json`: bounds and legend metadata for the historical disturbance overlay.
- `webmap_data/gdd_gen_diff_1km_wgs84_filled.tif`: GDD/GEN raster used by the app.
- `webmap_data/forecast_bundle/`: forecast model bundle for XGB, RF, GAM, and SVM.
- `webmap_data/forecast_bundle/README.md`: deeper provenance notes for forecast data generation.

## User-Facing Functionality

The map has four active layer cards in the ribbon below the header:

- `Hex pressure`: hexagon GeoJSON layer controlled by RCP scenario and year.
- `Predictions`: bark beetle forecast PNG overlays controlled by model, period, view, and threshold.
- `GDD raster`: 1 km GEN-difference GeoTIFF rendered directly in the browser.
- `Past disturbance`: classified historical bark beetle disturbance PNG overlay.

Each layer card has:

- A card body: clicking it selects that layer and opens its left-panel settings.
- A `Display` / `Hide` button: toggles whether that layer is drawn on the map.
- An opacity slider: controls only that layer.

Important behavior:

- Selecting a card and displaying a layer are separate actions.
- Multiple layers can be visible at the same time.
- The most recently displayed layer is drawn on top.
- The legend and click-readout follow the top visible layer.

## Code Structure In `index.html`

The JavaScript is organized around a few core ideas.

### `APP`

`APP` is the runtime state object. It stores the Leaflet map, loaded layer instances, selected features, current controls, cached manifests, and async load tokens.

Useful fields:

- `visibleOverlays`: ordered list of currently visible overlay ids. The last item is on top.
- `focusedOverlay`: the overlay whose settings are shown in the left panel.
- `layerInstances`: active Leaflet layers keyed by overlay id.
- `layerOpacity`: opacity values keyed by overlay id.
- `overlayLoadTokens`: guards against old async fetches replacing newer layer choices.

### `OVERLAY_CATALOG`

`OVERLAY_CATALOG` describes each possible layer:

- label
- description
- availability
- whether it uses RCP/year controls
- static metadata, where available

If a future layer is added, start by adding a catalog entry and a ribbon card in the HTML.

### Layer Loading

Layer loading is routed through `updateOverlayLayer(overlayId)`, which dispatches to `OVERLAY_LOADERS`.

Current loaders:

- `updateHexLayer()`: loads GeoJSON hexes from `MAP_META`.
- `updateForecastOverlay()`: loads PNG forecast overlays plus matching metadata.
- `updateGddRasterLayer()`: loads a GeoTIFF and renders it with `GeoRasterLayer`.
- `updateClassifiedImageOverlay()`: loads the historical bark beetle disturbance PNG.

The helper `beginOverlayLoad(overlayId)` increments a token and clears the old layer. After each async fetch, the loader checks `isCurrentOverlayLoad(...)`. This prevents race bugs when users click controls quickly.

### Layer Visibility And Order

Visibility is controlled by:

- `showOverlay(overlayId)`
- `hideOverlay(overlayId)`
- `setLayerVisibilityFromRibbon(overlayId, shouldShow)`

Layer order is controlled by `visibleOverlays`. When a layer is displayed, it is moved to the end of that array. `updateOverlayPaneOrder()` then updates Leaflet pane z-indexes.

### Left Panel Sync

`syncOverlayControls(refreshTarget)` updates the left panel after card selection or display changes.

It keeps these concepts separate:

- Focus: which layer is being edited.
- Visibility: which layers are currently drawn.
- Refresh target: which layer should be reloaded after a setting change.

This separation is why the user can edit one layer while several layers remain visible.

### Click Readout

All map clicks go through `updatePanelFromClick(latlng)`.

The app intentionally uses one map-level click handler instead of click handlers on every layer. This avoids duplicate events when overlays are stacked.

Click behavior:

- Always updates latitude, longitude, marker, and country.
- If GDD raster is top visible layer, samples the raster value.
- If forecast or past disturbance is top visible layer, reports that the PNG is visualization-only.
- Otherwise, tries to find and highlight a hex feature.

## Data Contracts

### Hex Pressure

Hex files are discovered through `MAP_META` in `webmap_data/map_meta.js`.

Expected structure:

- `MAP_META.years`
- `MAP_META.rcps`
- `MAP_META.files[rcp][year]`
- `MAP_META.valueRange`
- `MAP_META.bbox`

GeoJSON features are expected to contain:

- `properties.mean`

### Forecast Bundle

Forecast options are discovered from:

- `webmap_data/forecast_bundle/manifest.json`

The app expects:

- `models.{model}.periods`
- `models.{model}.metrics_file`
- `entries.{model::period}`
- `entries.{model::period}.ensemble_probability_tif`
- `entries.{model::period}.threshold_visible_tifs`

The app uses PNG sidecars for display:

- `probability_ensemble_mean.png`
- `risk_visible_tXX.png`

The metrics panel uses:

- `metrics_confusion.json`
- `aggregate_by_cutoff`

The threshold selector only shows thresholds that exist in all three places:

- manifest thresholds
- forecast entry threshold rasters
- model metrics rows

This avoids UI options that cannot render or cannot show performance metrics.

### Raster/Image Overlays

PNG overlays need companion metadata with:

- `bounds.south`
- `bounds.west`
- `bounds.north`
- `bounds.east`
- optional `legend`

GeoTIFF overlays are parsed in the browser, so very large rasters can be slow. For public GitHub Pages deployment, prefer PNG overlays when the raster is only used visually.

## Development Workflow

Because the app is static, local testing can be done with a simple HTTP server from the repo root:

```bash
python3 -m http.server 8000
```

Then open:

```text
http://localhost:8000/
```

Avoid opening `index.html` directly with `file://`, because browser security rules can block local fetches.

Before pushing, useful checks are:

```bash
git diff --check -- index.html README.md
git status --short
```

If Node.js is available, a JavaScript syntax check can be added, but this repository currently does not require Node.

## How To Add A New Layer

1. Add a card in the layer ribbon HTML.
2. Add an entry to `OVERLAY_CATALOG`.
3. Add a key to `APP.layerOpacity`, `APP.overlayLoadTokens`, and `APP.layerInstances`.
4. If it is a raster/image layer, add its id to `RASTER_OVERLAYS`.
5. Write a loader function that creates the Leaflet layer.
6. Add the loader to `OVERLAY_LOADERS`.
7. Update `setLegendLabels()` if the layer needs a custom legend.
8. Update `updatePanelFromClick()` if the layer needs custom click readout.

Keep layer ids short and stable. They are used in DOM `data-*` attributes, state keys, Leaflet pane names, and loader dispatch.

## Review Notes

The code favors plain JavaScript and explicit functions over clever abstractions. That is deliberate: the app is maintained as a single deployable HTML file, and most future editing will likely happen during data/model iteration rather than front-end framework work.

When changing behavior, protect these invariants:

- Card selection must not automatically hide other visible layers.
- Displaying a layer should move it to the top.
- Opacity sliders should affect only their own layer.
- Async layer loads must ignore stale responses.
- GitHub Pages paths should stay relative to the repository root.
