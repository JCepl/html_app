    // Static Leaflet DSS prototype for GitHub Pages deployment.
    //
    // ARCHITECTURE
    //  - One big state object `APP`. Per-map fields (map, layerInstances, all
    //    settings) describe the ACTIVE map only.
    //  - Overlays: bark-beetle pressure (vector hexes) + image overlays (PNG)
    //    for the forecast, PHENIPS, past disturbance and wind layers.
    //  - SPLIT SCREEN keeps two real Leaflet maps; switching the active side
    //    parks the current per-map state in `APP.panels` and loads the other
    //    side's, so the existing single-map pipeline only ever drives the
    //    active map (see the SPLIT SCREEN section).
    //  - Data is pre-rendered by the R scripts in webmap_data/ into PNGs +
    //    *_meta.json; bump ASSET_VERSION when those are regenerated.
    // Keep most runtime state in one object so the app remains easy to inspect in DevTools.
    //
    // SECTIONS IN THIS FILE, IN ORDER (search for the "====" banners):
    //  1. RISK COLOUR SCALE, THRESHOLD SLIDER & PNG MASKING
    //  2. OVERLAY REGISTRY HELPERS (visibility, panes, image layers)
    //  3. BARK-BEETLE FORECAST BUNDLE (manifest, entries, metrics)
    //  4. PHENIPS & WIND RASTER LAYERS
    //  5. LEGEND
    //  6. HEX / COUNTRY STYLING, CLICK MARKER & READOUT
    //  7. LAYER LOADERS & CONTROL-PANEL SYNC
    //  8. CLICK INSPECTION & POPUP
    //  9. RISK RADAR CHART -- 5-hazard profile at a clicked point
    // 10. BASEMAPS
    // 11. SPLIT SCREEN -- two synced maps; the panel drives the ACTIVE one
    // 12. INIT, SHAREABLE LINK & CONTROL WIRING
    // --- Optional, privacy-friendly usage statistics + feedback link --------
    // Both stay OFF until a value is filled in:
    //   goatcounterCode : the code of a GoatCounter site (https://CODE.goatcounter.com).
    //                     Cookie-free, no personal data; skipped when the browser sends
    //                     Do Not Track / Global Privacy Control.
    //   feedbackUrl     : link of a short feedback form (a quiet link in the menu, never a pop-up).
    const SITE_SETTINGS = { goatcounterCode: "", feedbackUrl: "" };
    let analyticsOn = false;

    function track(name) {
      if (!analyticsOn) return;
      try { window.goatcounter.count({ path: "event/" + name, title: name, event: true }); } catch (e) { /* analytics must never break the app */ }
    }

    function initSiteSettings() {
      const feedback = document.getElementById("feedbackLink");
      if (feedback && SITE_SETTINGS.feedbackUrl) { feedback.href = SITE_SETTINGS.feedbackUrl; feedback.hidden = false; }

      const code = SITE_SETTINGS.goatcounterCode;
      if (!code || navigator.doNotTrack === "1" || navigator.globalPrivacyControl) return;
      window.goatcounter = { no_onload: false };
      const tag = document.createElement("script");
      tag.async = true;
      tag.src = "https://gc.zgo.at/count.js";
      tag.setAttribute("data-goatcounter", "https://" + code + ".goatcounter.com/count");
      document.head.appendChild(tag);
      analyticsOn = true;
      const note = document.getElementById("privacyNote");
      if (note) note.hidden = false;

      // One delegated listener: what people use, never who they are.
      document.addEventListener("click", function(e) {
        const t = e.target.closest && e.target.closest("[data-toggle-overlay], .period-chip[data-jump-overlay], .dss-pillar, #shareBtn, #detailModeBtn, #feedbackLink, #mobLayersBtn, #mobLegendBtn");
        if (!t) return;
        if (t.dataset.toggleOverlay) {
          const id = t.dataset.toggleOverlay;
          setTimeout(function() { if (isOverlayVisible(id)) track("layer-" + id); }, 600);
        } else if (t.dataset.jumpOverlay) track("chip-" + t.dataset.jumpOverlay + "-" + (t.dataset.jumpPeriod || "historic"));
        else if (t.classList.contains("dss-pillar")) track("pillar-" + (t.id || t.dataset.pillar || "x"));
        else track("click-" + t.id);
      }, true);  // capture: the cards stop clicks from bubbling
      tag.addEventListener("load", function() { track(Math.min(screen.width, screen.height) < 600 ? "visit-phone" : "visit-desktop"); });
    }

    // Final bark-beetle model (P1C4, SSP5-8.5): its own sequential ramp (light to dark) over 0..P1C4_CAP, flag at 0.5.
    // The ramp equals SEQ in webmap_data/build_p1c4_bundle.py; keep them in step.
    const P1C4_MODEL = "p1c4_ssp585", P1C4_CAP = 0.8, P1C4_FLAG = 0.5;
    const P1C4_RAMP = ["#ffffcc", "#ffeda0", "#fed976", "#feb24c", "#fd8d3c", "#fc4e2a", "#e31a1c", "#bd0026", "#800026"];
    const EXTRAP_RAMP = ["#f0f0f0", "#c7c7c7", "#9f9f9f", "#767676", "#4e4e4e", "#252525"];

    const APP = {
      map: null,
      // Stable reference to the original (left) Leaflet map instance, set
      // once in initMap() and never reassigned -- unlike `map`, which swaps
      // to whichever side is active. Pilot-action markers key off this so
      // split screen never cross-adds a side's pins to the other map.
      leftMap: null,
      showPilotActions: false,
      showSeed4Forest: false,
      activeTab: "drivers",
      lastClick: null,
      pilotLayerLeft: null,
      pilotLayerRight: null,
      // Natura 2000: fetched once (11MB, lazy) and cached here; each side
      // gets its own Leaflet layer built from the same cached GeoJSON.
      showNatura2000: false,
      natura2000Opacity: 1,
      natura2000Data: null,
      natura2000LayerLeft: null,
      natura2000LayerRight: null,
      // Derived from bbtpast (Past disturbance) visibility, not user-toggled
      // directly -- see syncCatalogueEventsFromPastDisturbance().
      showCatalogueEvents: false,
      catalogueLayerLeft: null,
      catalogueLayerRight: null,
      radarPanelEl: null,
      radarRequestToken: 0,
      basemapLayer: null,
      currentBasemap: "osm",
      clickMarker: null,
      countryHaloLayer: null,
      selectedHexLayer: null,
      countryLayer: null,
      selectedCountryLayer: null,
      currentYear: null,
      currentRcp: "rcp85",
      // Hide map values whose risk-scale fraction is below this (0..1). Applies
      // to all risk layers (bbtl hexes, wind + forecast probability PNGs).
      riskThreshold: 0,
      // Default = simplified surface (one model, no GCM picker, near-future
      // only). "Show full detail" flips this off. Not split-screen-per-side
      // state -- deliberately global, same for both maps.
      simpleView: true,
      // Split screen: two real maps. The active one is what `map`/`layerInstances`
      // (and all per-map fields) point to, so the whole pipeline drives only it.
      // The inactive side's state is parked in `panels`.
      split: false,
      activeSide: "left",
      panels: { left: null, right: null },
      // `visibleOverlays` controls what is drawn. The last item is drawn on top.
      visibleOverlays: [],   // nothing on at start; the card of focusedOverlay is just the one shown in the panel
      // `focusedOverlay` controls which settings are shown in the left panel.
      focusedOverlay: "droughtFuture",
      layerOpacity: {
        bbtl: 0.85,
        bbtforecast: 0.85,
        phenips: 0.85,
        bbtpast: 0.85,
        wind: 0.85,
        droughtFuture: 0.85,
        droughtRef: 0.85,
        fireWeather: 0.85,
        wildfire: 0.85,
        windPressure: 0.85
      },
      currentPhenipsScen: "obs",
      currentPhenisPeriod: "m2",
      currentPhenipsVar: "gdd",
      currentPhenipsView: "abs",
      currentPhenipsExp: "shade",
      phenipsV5: null,
      phenipsData: null,
      currentForecastModel: P1C4_MODEL,
      currentForecastPeriod: "2041_2070",
      currentForecastView: "probability",
      currentForecastThreshold: "0.5",
      currentWindPeriod: "near",
      currentWindGcm: "mpi",
      currentDroughtFuturePeriod: "2041_2070_ssp370",
      currentDroughtRefMetric: "n_events_2yr",
      currentFireView: "change",
      currentFireScen: "rcp85_2041_2070",
      currentForecastRasterRange: null,
      forecastManifest: null,
      forecastMetricsCache: {},
      // Incremented before each async fetch. Old responses are ignored if a newer load started.
      overlayLoadTokens: {
        bbtl: 0,
        bbtforecast: 0,
        phenips: 0,
        bbtpast: 0,
        wind: 0,
        droughtFuture: 0,
        droughtRef: 0,
        fireWeather: 0,
        wildfire: 0,
        windPressure: 0
      },
      overlayVisible: false,
      // One cached GeoJSON per hex overlay (see HEX_OVERLAY_IDS) -- keyed so
      // two hex layers can be visible/cached at once without fighting over
      // a single slot.
      hexGeojson: { bbtl: null, wildfire: null, windPressure: null },
      countryGeojson: null,
      layerInstances: {
        bbtl: null,
        bbtforecast: null,
        phenips: null,
        bbtpast: null,
        wind: null,
        droughtFuture: null,
        droughtRef: null,
        fireWeather: null,
        wildfire: null,
        windPressure: null
      }
    };

    const RASTER_OVERLAYS = ["phenips", "bbtforecast", "bbtpast", "wind", "droughtFuture", "droughtRef", "fireWeather"];
    const WIND_DATA_URL = "webmap_data/wind/wind_meta.json";
    const DROUGHT_FUTURE_DATA_URL = "webmap_data/drought/drought_future_meta.json";
    const DROUGHT_REF_DATA_URL = "webmap_data/drought/drought_ref_meta.json";
    // Wind colour scale saturates at this probability (risk_cap in the R export),
    // so the top of the wind legend (t=1) corresponds to p=0.5.
    const WIND_RISK_CAP = 0.5;
    // Bump to bust browser cache when overlay PNGs are regenerated.
    const ASSET_VERSION = "20261103-p1c4";

    // Unified risk palette, shared by every hazard/probability layer (bark
    // beetle history/forecast, wind, PHENIPS, CHELSA): Spectral, sampled
    // directly from Tommaso's own QGIS wind-risk rendering (2026-07 RE-
    // ENFORCE meeting screenshot). Purple -> blue -> teal -> green -> pale
    // yellow -> orange -> red -> dark red. Deliberately chosen over a
    // perceptually-monotonic warm ramp for full color variety and to match
    // the consortium's shared wind-risk convention everywhere, not just wind.
    // Builds a CSS linear-gradient string from a ramp array so there is only
    // ONE place gradient strings get generated -- previously RISK_GRADIENT
    // was a second, hand-typed copy of RISK_RAMP that had already drifted
    // out of sync with it (and with the .legend-bar CSS default below).
    function buildGradientCss(ramp) {
      const n = ramp.length - 1;
      return "linear-gradient(to right, " + ramp.map(function(hex, i) {
        return hex + " " + Math.round((i / n) * 100) + "%";
      }).join(", ") + ")";
    }

    // Bark beetle's own scale -- see webmap_data/color_scales.R
    // (BARK_BEETLE_RAMP); PNGs are pre-rendered from that same file so these
    // must stay in sync with it, not with WIND_RAMP below.
    const BARK_BEETLE_RAMP = ["#5e4fa2", "#3288bd", "#66c2a5", "#abdda4", "#e6f598",
                       "#ffffbf", "#fee08b", "#fdae61", "#f46d43", "#d53e4f", "#9e0142"];
    const BARK_BEETLE_GRADIENT = buildGradientCss(BARK_BEETLE_RAMP);
    const P1C4_GRADIENT = buildGradientCss(P1C4_RAMP);
    const EXTRAP_GRADIENT = buildGradientCss(EXTRAP_RAMP);
    const isP1c4 = function(model) { return (model || APP.currentForecastModel) === P1C4_MODEL; };
    const P1C4_CLASSES = [[0.25, "low"], [P1C4_FLAG, "moderate"], [Infinity, "flagged, 0.5 or more"]];

    // Every vector hex ("pressure index") layer -- bark beetle, wildfire,
    // wind -- shares one code path (styling, click sampling, RCP/year
    // control panel, citation). Add an id here (plus an OVERLAY_CATALOG
    // entry, a MAP_META_* object, and a ribbon card) to add another.
    const HEX_OVERLAY_IDS = ["bbtl", "wildfire", "windPressure"];

    // Simple/Advanced view (2026-09 consortium feedback: default to the
    // least-confusing surface for a non-scientist audience -- one model, no
    // GCM picker -- and let "Show full detail" restore everything). Every
    // time period, including end-of-century, stays available in BOTH modes
    // -- explicitly not something to hide (2026-09-06 correction: an
    // earlier version of this also capped the far-future/end-of-century
    // options, which is wrong -- the whole point of a climate-impact DSS is
    // showing what changes by end of century). See applySimpleView().
    const SIMPLE_VIEW_DEFAULT_FORECAST_MODEL = P1C4_MODEL;

    // Wind's own scale (Tommaso's QGIS rendering) -- see WIND_RAMP in
    // webmap_data/color_scales.R. An independent copy: historically this and
    // BARK_BEETLE_RAMP were the same array (wind's ramp reused everywhere),
    // now split so either can change without touching the other.
    const WIND_RAMP = ["#5e4fa2", "#3288bd", "#66c2a5", "#abdda4", "#e6f598",
                       "#ffffbf", "#fee08b", "#fdae61", "#f46d43", "#d53e4f", "#9e0142"];
    const WIND_GRADIENT = buildGradientCss(WIND_RAMP);
    // Diverging (ColorBrewer RdBu) ramp for the bark-beetle ensemble forecast's
    // "delta vs. current" view -- deliberately a different family from the
    // Spectral BARK_BEETLE_GRADIENT so a change map is never mistaken for an
    // absolute-risk one. Matches DELTA_RAMP/PCT_CAP in build_ensemble_delta.R.
    // Stop positions are NOT evenly spaced -- they mirror the signed-sqrt
    // transform the R script applies before coloring (0/+-18/+-32/+-42/+-48%
    // land at 0/18/32/42/48% etc.), so the legend bar's own color spacing
    // shows the same "fast near zero, compressed at the extremes"
    // sensitivity as the map -- a straight 0/10/20.../100% spacing would
    // make the legend lie about how visible a small change actually is.
    const DELTA_GRADIENT = "linear-gradient(to right, #053061 0%, #2166ac 18%, #4393c3 32%, #92c5de 42%, #d1e5f0 48%, #f7f7f7 50%, #fddbc7 52%, #f4a582 58%, #d6604d 68%, #b2182b 82%, #67001f 100%)";
    const DELTA_CAP = 200; // percent change, saturates the ramp at +/-200%
    // Low-risk end fades across the bottom RISK_FADE_KNEE of the scale, but
    // never below RISK_FADE_FLOOR opacity, so even ~1% risk stays faintly
    // visible (kept in sync with fade_knee / fade_floor in the R exports).
    const RISK_FADE_KNEE = 0.3;
    const RISK_FADE_FLOOR = 90 / 255;


    // ====================================================================
    // RISK COLOUR SCALE, THRESHOLD SLIDER & PNG MASKING
    // ====================================================================
    function lerpHexColor(a, b, t) {
      const ah = parseInt(a.slice(1), 16), bh = parseInt(b.slice(1), 16);
      const r = Math.round(((ah >> 16) & 255) + (((bh >> 16) & 255) - ((ah >> 16) & 255)) * t);
      const g = Math.round(((ah >> 8) & 255) + (((bh >> 8) & 255) - ((ah >> 8) & 255)) * t);
      const bl = Math.round((ah & 255) + ((bh & 255) - (ah & 255)) * t);
      return "#" + ((1 << 24) + (r << 16) + (g << 8) + bl).toString(16).slice(1);
    }

    function rampColor(ramp, t) {
      t = Math.max(0, Math.min(1, t));
      const seg = t * (ramp.length - 1);
      const i = Math.min(ramp.length - 2, Math.floor(seg));
      return lerpHexColor(ramp[i], ramp[i + 1], seg - i);
    }

    // 256-step sRGB lookup of a risk ramp, used to recover the scale fraction
    // (0..1) of an already-colored overlay pixel so the threshold slider can
    // hide low-risk pixels in the wind / forecast probability PNGs. Built per
    // ramp since the wind PNGs and forecast-probability PNGs use different
    // palettes and must each be matched against their own ramp.
    function buildRiskLut(ramp) {
      const n = 256, lut = new Array(n);
      for (let i = 0; i < n; i++) {
        const c = rampColor(ramp, i / (n - 1));
        lut[i] = {
          r: parseInt(c.slice(1, 3), 16),
          g: parseInt(c.slice(3, 5), 16),
          b: parseInt(c.slice(5, 7), 16),
          t: i / (n - 1)
        };
      }
      return lut;
    }
    const BARK_BEETLE_LUT = buildRiskLut(BARK_BEETLE_RAMP);

    function makeColorToRiskFraction(lut) {
      const cache = new Map();
      return function colorToRiskFraction(r, g, b) {
        const key = (r << 16) | (g << 8) | b;
        const cached = cache.get(key);
        if (cached !== undefined) return cached;
        let bt = 0, bd = Infinity;
        for (let i = 0; i < lut.length; i++) {
          const e = lut[i];
          const dr = e.r - r, dg = e.g - g, db = e.b - b;
          const d = dr * dr + dg * dg + db * db;
          if (d < bd) { bd = d; bt = e.t; }
        }
        cache.set(key, bt);
        return bt;
      };
    }
    // Default/bark-beetle sampler (bbtl, bbtforecast, bbtpast). Wind gets its
    // own colorToWindFraction below since it now has its own palette -- do
    // not use this one for wind, it would nearest-match against the wrong LUT.
    const colorToRiskFraction = makeColorToRiskFraction(BARK_BEETLE_LUT);
    const WIND_LUT = buildRiskLut(WIND_RAMP);
    const colorToWindFraction = makeColorToRiskFraction(WIND_LUT);
    const colorToP1c4Fraction = makeColorToRiskFraction(buildRiskLut(P1C4_RAMP));
    const colorToExtrapFraction = makeColorToRiskFraction(buildRiskLut(EXTRAP_RAMP));

    // The drought future-classes layer uses its own 4-colour classified
    // palette (see export_drought_overlay.R fut_colors), unrelated to either
    // ramp above -- matching its pixels against colorToRiskFraction's LUT
    // would produce meaningless nearest-neighbour noise. This is the same
    // "match a sampled pixel back to a scale position" trick, but built
    // against drought's own 4 known colours mapped onto severity (Extreme=1
    // .. Low=0), so the radar chart can plot a real severity reading.
    const DROUGHT_SEVERITY_LUT = [
      { r: 0xB2, g: 0x18, b: 0x2B, t: 1.00 }, // Extreme drought hazard
      { r: 0xEF, g: 0x8A, b: 0x62, t: 0.67 }, // High
      { r: 0xFD, g: 0xDB, b: 0xC7, t: 0.33 }, // Moderate
      { r: 0xF7, g: 0xF7, b: 0xF7, t: 0.00 }  // Low
    ];
    const colorToDroughtSeverity = makeColorToRiskFraction(DROUGHT_SEVERITY_LUT);

    // Load a same-origin overlay PNG, make pixels below the risk threshold
    // transparent, and return a data URL. Falls back to the original URL on
    // any failure (e.g. a tainted canvas).
    function loadThresholdedImage(url, threshold, onReady, onError, colorFractionFn) {
      const fractionOf = colorFractionFn || colorToRiskFraction;
      const img = new Image();
      img.onload = function() {
        try {
          const c = document.createElement("canvas");
          c.width = img.naturalWidth;
          c.height = img.naturalHeight;
          const ctx = c.getContext("2d");
          ctx.drawImage(img, 0, 0);
          const id = ctx.getImageData(0, 0, c.width, c.height);
          const d = id.data;
          for (let i = 0; i < d.length; i += 4) {
            if (d[i + 3] === 0) continue;
            if (fractionOf(d[i], d[i + 1], d[i + 2]) < threshold) {
              d[i + 3] = 0;
            }
          }
          ctx.putImageData(id, 0, 0);
          onReady(c.toDataURL());
        } catch (e) {
          if (onError) onError(e); else onReady(url);
        }
      };
      img.onerror = function(e) { if (onError) onError(e); else onReady(url); };
      img.src = url;
    }

    // Absolute probability at the top (t=1) of each layer's colour scale.
    function overlayScaleCap(overlayId) {
      if (overlayId === "wind") return WIND_RISK_CAP;
      if (HEX_OVERLAY_IDS.includes(overlayId)) {
        const m = OVERLAY_CATALOG[overlayId].meta;
        return (m && m.valueRange) ? Number(m.valueRange.zmax) : 1;
      }
      if (overlayId === "bbtforecast" && isP1c4()) return P1C4_CAP;
      return 1; // earlier bbtforecast models: probability-like 0..1
    }

    // Convert the slider's absolute probability threshold into a colour-scale
    // fraction (0..1) for the given layer, so "50%" always means p = 0.5.
    function riskFractionFor(overlayId) {
      if (APP.riskThreshold <= 0) return 0;
      return Math.max(0, Math.min(1, APP.riskThreshold / overlayScaleCap(overlayId)));
    }

    // The slider only filters the layer that is actually on top (the one the
    // legend and slider range describe); stacked layers underneath render full.
    // Views drawn on their own colour scale (change, share of years, extrapolation): the hide-below slider does not apply.
    const FORECAST_UNFILTERED_VIEWS = ["delta", "share05", "extrapolation"];
    function effectiveRiskFraction(overlayId) {
      if (overlayId !== getTopVisibleOverlayId()) return 0;
      // The delta view is colored with a diverging (blue/white/red) ramp, not
      // the shared Spectral risk ramp the slider's color reversal assumes --
      // thresholding it would reverse-match delta colors against the wrong
      // LUT. Never mask it; it's not a "risk" surface to begin with.
      if (overlayId === "bbtforecast" && FORECAST_UNFILTERED_VIEWS.includes(APP.currentForecastView)) return 0;
      return riskFractionFor(overlayId);
    }

    // Add an image overlay, applying the risk threshold for PNG layers. `token`
    // guards against a stale async result overwriting a newer load.
    function addRiskImageOverlay(overlayId, url, bounds, token) {
      const place = function(src) {
        if (token !== undefined && !isCurrentOverlayLoad(overlayId, token)) return;
        addImageOverlayLayer(overlayId, src, bounds);
        updateOverlayPaneOrder();
      };
      const fraction = effectiveRiskFraction(overlayId);
      if (fraction > 0) {
        // Wind now has its own palette (WIND_RAMP) distinct from bark
        // beetle's -- loadThresholdedImage's default colorFractionFn is the
        // bark-beetle LUT, so wind must explicitly pass its own or the
        // threshold mask reverse-matches wind's pixels against the wrong ramp.
        const colorFn = overlayId === "wind" ? colorToWindFraction : (overlayId === "bbtforecast" && isP1c4()) ? colorToP1c4Fraction : undefined;
        loadThresholdedImage(url, fraction, place, function() { place(url); }, colorFn);
      } else {
        place(url);
      }
    }

    // Clear the risk slider and re-render any visible risk layer unfiltered.
    // No-op when already at zero so it is cheap to call on every toggle.
    function resetRiskThreshold() {
      if (APP.riskThreshold <= 0) return;
      APP.riskThreshold = 0;
      if (DOM.riskSlider) DOM.riskSlider.value = 0;
      if (DOM.riskValue) DOM.riskValue.textContent = "0%";
      applyRiskThreshold();
    }

    // Re-apply the current risk threshold to the top (displayed) risk layer
    // only — that is the surface the slider and legend describe.
    function applyRiskThreshold() {
      const top = getTopVisibleOverlayId();
      if (HEX_OVERLAY_IDS.includes(top) && APP.layerInstances[top]) {
        resetSelectedHex();
        APP.layerInstances[top].setStyle(baseHexStyle(top));
      } else if (top === "wind") {
        updateWindLayer();
      } else if (top === "bbtforecast") {
        updateForecastOverlay();
      }
    }

    // Wind probability saturates at 0.5, so its slider caps at 50%; the beetle
    // layers are full 0..1 probability and run to 100%.
    function updateRiskSliderRange(topOverlayId) {
      const max = topOverlayId === "wind" ? 49 : (topOverlayId === "bbtforecast" && isP1c4()) ? 79 : 100;
      if (Number(DOM.riskSlider.max) !== max) {
        DOM.riskSlider.max = max;
        if (Number(DOM.riskSlider.value) > max) {
          DOM.riskSlider.value = max;
          APP.riskThreshold = max / 100;
          DOM.riskValue.textContent = max + "%";
        }
      }
    }

    // The "hide values below" slider only does something for continuous-scale
    // layers (wind, bark-beetle susceptibility, the pressure indices). For the
    // classed layers (drought, fire weather, past disturbance, PHENIPS) it did
    // nothing, so it is removed there instead of sitting on the panel idle.
    function riskFilterApplies(topOverlayId) {
      if (topOverlayId === "wind") return true;
      if (topOverlayId === "bbtforecast") return !FORECAST_UNFILTERED_VIEWS.includes(APP.currentForecastView);
      return HEX_OVERLAY_IDS.includes(topOverlayId);
    }
    function syncRiskFilterBox(topOverlayId) {
      const box = document.getElementById("riskFilterBox");
      if (!box) return;
      const applies = APP.overlayVisible && riskFilterApplies(topOverlayId);
      box.style.display = applies ? "" : "none";
      if (!applies && APP.riskThreshold > 0) {   // never leave a hidden cut-off behind
        APP.riskThreshold = 0;
        if (DOM.riskSlider) DOM.riskSlider.value = 0;
        if (DOM.riskValue) DOM.riskValue.textContent = "0%";
      }
    }

    // Shade the hidden (left) portion of the legend bar up to the threshold.
    function updateLegendThreshMarker(topOverlayId) {
      const marker = DOM.legendThreshMarker;
      if (!marker) return;
      const isRisk = topOverlayId === "wind" || topOverlayId === "bbtforecast" || HEX_OVERLAY_IDS.includes(topOverlayId);
      if (isRisk && APP.riskThreshold > 0) {
        marker.style.left = (riskFractionFor(topOverlayId) * 100) + "%";
        marker.style.display = "block";
      } else {
        marker.style.display = "none";
      }
    }
    const PHENIPS_BOUNDS = { south: 34.5, west: -10.7, north: 71.57, east: 38.58 };
    const PHENIPS_DATA_URL = "webmap_data/phenips/phenips_data.json";

    const OVERLAY_CATALOG = {
      bbtl: {
        label: "Bark beetle pressure index (Grunig et al.)",
        description: "Modelled bark beetle pressure index by scenario and year (Grunig et al. 2026, Science). Select RCP scenario and year using the controls below.",
        available: true,
        usesRcp: true,
        usesYear: true,
        meta: MAP_META
      },
      phenips: {
        label: "PHENIPS degree-days & generations",
        description: "Bark beetle (Ips typographus) development from the PHENIPS phenology model: effective degree-days (bark-temperature corrected) and completed generations per year. Observed climate 1981-2024 (CHELSA, per year; shaded stand or sun-exposed) and projections to 2100 (EURO-CORDEX, 3 climate models, photoperiod-gated, max 3 generations). Development potential, not an infestation forecast.",
        available: true,
        usesRcp: false,
        usesYear: false,
        meta: null
      },
      bbtforecast: {
        label: "Bark beetle hazard",
        description: "Hazard for Norway spruce forest from bark-beetle and wind disturbance, learnt from documented outbreak events (223 events: the EFI DFDE database plus a literature search) overlaid with the European Forest Disturbance Atlas, 1985–2024. Boosted trees on 16 predictors: terrain, climate normals and 4-year climate memory (autumn precipitation, beetle generations from PHENIPS, evaporative demand against precipitation, multi-year drought run length). Future climate: CHELSA-CMIP6 MPI-ESM1-2-HR, SSP5-8.5 only so far (SSP1-2.6 and SSP3-7.0 to follow); each future year replays an observed 1981–2010 year shifted to the future monthly mean, so year-to-year variability is today's. Validation: AUC 0.87 with a region and a period both held out, 0.82 trained to 2014 and tested on 2017–2024. Shown only where the Norway spruce probability is at least 0.2, the model's domain. Share of all land cells (not only spruce areas) with hazard of 0.5 or more (the flag): about 6.1% in 1981–2010, 6.7% in 2011–2040, 7.9% in 2041–2070 and 9.4% in 2071–2100. The value is a susceptibility index, a ranking calibrated at the sample prevalence: not an annual attack probability and not a forecast of single years. Future maps are a climate shift on today's fixed spruce area, not a forecast of where forest will be.",
        available: true,
        usesRcp: false,
        usesYear: false,
        meta: null
      },
      bbtpast: {
        label: "Past wind / bark-beetle disturbance (EFDA v3)",
        description: "Number of years between 1985 and 2024 in which satellite data (European Forest Disturbance Atlas v3.0, Landsat, 1 km) detected forest disturbance attributed to the wind / bark-beetle complex (cells with more than 1 % of their area disturbed). Wind and bark beetle cannot be separated, so large storm events such as the 2009 Landes storm also show; fire and harvest are excluded. Dots are individually documented bark-beetle outbreaks (event catalogue). Zero is transparent.",
        available: true,
        usesRcp: false,
        usesYear: false,
        meta: {
          image: "webmap_data/efda_v3/efda_v3_wbb_hit_years.webp",
          metaUrl: "webmap_data/efda_v3/efda_v3_meta.json",
          bounds: null,
          legend: null
        }
      },
      fireWeather: {
        label: "Fire-weather severity",
        description: "Summer (June-August) fire-weather severity from one climate-model chain (RCA4 / MPI-ESM-LR, 0.15 degree grid): how hard a fire would be to contain, not how likely one is to start. Show the change relative to 1981-2005 or the absolute Fire Weather Index. Data: F. Taccaliti (UNIPD-TESAF).",
        available: true,
        usesRcp: false,
        usesYear: false,
        meta: { metaUrl: "webmap_data/fire/fire_meta.json", bounds: null, images: null }
      },
      wildfire: {
        label: "Fire pressure index (Grunig et al.)",
        description: "Simulated fire pressure index of an independent study (Grunig et al. 2026, Science), by scenario and year. Context only, not a RE-ENFORCE model. Select RCP scenario and year using the controls below.",
        available: true,
        usesRcp: true,
        usesYear: true,
        meta: MAP_META_FIRE
      },
      ashDieback: {
        label: "Ash dieback",
        description: "Planned overlay for ash dieback disturbance. Currently only referenced as an event-marker tag (Pilot actions), no model/layer exists yet.",
        available: false,
        usesRcp: true,
        usesYear: true,
        meta: null
      },
      droughtFuture: {
        label: "Future SPEI drought classes",
        description: "December SPEI-12-equivalent values, indicating projected climatic water balance relative to the historical reference climate (1981-2010), across 4 future scenario/period combinations. Negative = drier, positive = wetter than reference. Data: Stefan Ebner (BFW).",
        available: true,
        usesRcp: false,
        usesYear: false,
        meta: {
          metaUrl: "webmap_data/drought/drought_future_meta.json",
          bounds: null,
          images: null
        }
      },
      droughtRef: {
        label: "Reference drought indicators (MYD)",
        description: "Multi-year-drought (MYD) indicators, historical reference period 1981-2010: number of events, max event duration, cumulative duration. A MYD event is SPEI-12 <= -1.1 for >=2 consecutive years. CHELSA-native ~1km. Data: Stefan Ebner (BFW).",
        available: true,
        usesRcp: false,
        usesYear: false,
        meta: {
          metaUrl: "webmap_data/drought/drought_ref_meta.json",
          bounds: null,
          images: null
        }
      },
      wind: {
        label: "Wind damage susceptibility",
        description: "Probability of stand damage under a 100-year wind speed (0-1; 0.5 and above = damage expected). Combines the wind hazard with the vulnerability of the 2020 forest stands, so it is the one layer that includes vulnerability. Future scenario SSP5-8.5. The change between periods is small because the stand map is fixed at 2020. Data: T. Baggio (UNIPD).",
        available: true,
        usesRcp: false,
        usesYear: false,
        meta: {
          metaUrl: "webmap_data/wind/wind_meta.json",
          bounds: null,
          images: null
        }
      },
      windPressure: {
        label: "Wind pressure index",
        description: "Modelled wind disturbance pressure index by scenario and year (Grunig et al. 2026, Science). Select RCP scenario and year using the controls below.",
        available: true,
        usesRcp: true,
        usesYear: true,
        meta: MAP_META_WIND
      }
    };

    // Local copy of GISCO CNTR_RG_01M_2024 clipped to Europe (see webmap_data/build_country_borders.py).
    const COUNTRY_GEOJSON_URL = "webmap_data/borders/countries_europe.geojson";
    const FORECAST_MANIFEST_URL = "webmap_data/forecast_bundle/manifest.json";

    // RE-ENFORCE's 6 official numbered Pilot Actions (PA1-PA6), one pin per
    // action -- sourced from the project's own authoritative deliverable AND
    // its own overview map, not news articles or guesses:
    //   Report on Pilot Actions (Deliverable 2.3.3, Output 2.1, v1, 05 2026)
    //   https://www.interreg-central.eu/wp-content/uploads/2026/05/Report_Pilot-actions_D.2.3.3.pdf
    //   Figure 1 (table) -> which driver/PA-code per location.
    //   Figure 2 (map, page 7) -> exactly which icon per driver, AND that
    //   Stefan's own map shows 10 pins, not 9: Duna-Dráva (Hungary) and
    //   Agordino (Italy) each get TWO adjacent pins (one per action, not one
    //   pin listing two codes), while Austria's two PA5 sites (Japons AND
    //   Horn, both Lower Austria) are shown as a SINGLE bark-beetle pin on
    //   his map -- matched exactly here rather than our own 9-pin grouping.
    // `hazard` selects the icon (webmap_data/pilot_icons/*.png) -- these are
    // Stefan's own icon graphics, cropped from the embedded map image in the
    // PDF above (page 7), not a redrawn approximation.
    // This superseded an earlier set built from general news articles, which
    // had real errors: Pilsen is PA-2 ONLY (no PA-5 there -- PA-5's second
    // site is Horn, Austria, missing entirely before); PA-6's Polish site is
    // Toruń, not Rytel (Rytel isn't in the official deliverable at all);
    // Goriška Karst, Slovenia was mislabeled "PA-3" -- Slovenia hosts no
    // numbered Pilot Action in this deliverable, so it's dropped rather than
    // shown with a fabricated code.
    // Only Pilsen has a published exact coordinate; the rest are approximate
    // (town/region centroid, or a small offset off a shared site so two pins
    // there don't sit exactly on top of each other) -- marked `approx: true`
    // so the popup is honest about precision. `link` points to each action's
    // own dedicated news/event page where one exists (verified live
    // 2026-09-04); otherwise to the official PDF above, deep-linked to that
    // action's page.
    const PILOT_ACTIONS = [
      { id: "hu-dunadrava-pa1", pa: "PA-1", hazard: "ash_dieback", label: "Duna-Dráva National Park", country: "Hungary",
        type: "Ash dieback (hardwood gallery forests)", lat: 45.95, lng: 18.08, approx: true,
        link: "https://www.interreg-central.eu/wp-content/uploads/2026/05/Report_Pilot-actions_D.2.3.3.pdf#page=8" },
      { id: "hu-dunadrava-pa4", pa: "PA-4", hazard: "riparian", label: "Duna-Dráva National Park", country: "Hungary",
        type: "Improper management of riparian forests", lat: 45.95, lng: 18.14, approx: true,
        link: "https://www.interreg-central.eu/wp-content/uploads/2026/05/Report_Pilot-actions_D.2.3.3.pdf#page=25" },
      { id: "de-pomerania", pa: "PA-1", hazard: "ash_dieback", label: "Western Pomerania", country: "Germany",
        type: "Ash dieback (since 2008; storm-thrown 2022)", lat: 53.70, lng: 12.64, approx: true,
        link: "https://www.interreg-central.eu/news/ash-dieback-and-multiple-reforestation-approaches-in-germany/" },
      { id: "hr-jastrebarsko", pa: "PA-2", hazard: "drought", label: "Jastrebarsko", country: "Croatia",
        type: "Drought (European beech)", lat: 45.6667, lng: 15.6500, approx: true,
        link: "https://www.interreg-central.eu/wp-content/uploads/2026/05/Report_Pilot-actions_D.2.3.3.pdf#page=12" },
      { id: "cz-pilsen", pa: "PA-2", hazard: "drought", label: "Near Pilsen (Plzeň)", country: "Czech Republic",
        type: "Drought (Scots pine, since 2018)", lat: 49.8516714, lng: 13.4085067, approx: false,
        link: "https://www.youtube.com/watch?v=D6szlpYw4Vc" },
      { id: "pl-myszyniec", pa: "PA-3", hazard: "fire", label: "Myszyniec Forest District", country: "Poland",
        type: "Fire (stand-replacing, 2014)", lat: 53.2833, lng: 21.35, approx: true,
        link: "https://www.interreg-central.eu/news/uav-poland/" },
      { id: "it-agordino-pa3", pa: "PA-3", hazard: "fire", label: "Agordino (Dolomites)", country: "Italy",
        type: "Fire (2018, 2011)", lat: 46.3167, lng: 12.03, approx: true,
        link: "https://www.interreg-central.eu/wp-content/uploads/2026/05/Report_Pilot-actions_D.2.3.3.pdf#page=22" },
      { id: "it-agordino-pa6", pa: "PA-6", hazard: "wind", label: "Agordino (Dolomites)", country: "Italy",
        type: "Storm (Vaia storm, 2018)", lat: 46.3167, lng: 12.09, approx: true,
        link: "https://www.interreg-central.eu/wp-content/uploads/2026/05/Report_Pilot-actions_D.2.3.3.pdf#page=34" },
      { id: "at-japons-horn", pa: "PA-5", hazard: "bark_beetle", label: "Japons & Horn", country: "Austria",
        type: "Spruce bark beetle & drought (Lower Austria)", lat: 48.643, lng: 15.428, approx: true,
        link: "https://www.interreg-central.eu/event/open-day-in-austria/" },
      { id: "pl-torun", pa: "PA-6", hazard: "wind", label: "Toruń", country: "Poland",
        type: "Storm (2017 hurricane, 18,425 ha lost)", lat: 53.0138, lng: 18.5981, approx: true,
        link: "https://www.interreg-central.eu/wp-content/uploads/2026/05/Report_Pilot-actions_D.2.3.3.pdf#page=32" }
    ];
    // European Bark-Beetle Event Catalogue (~/barkbeetle-catalogue, own repo).
    // 220 events, 378 source links (enhanced catalogue incl. DFDE additions; webmap_data/build_catalogue_events.py) Full multi-source citations per event
    // (not just one link) -- see the sources array on each event. Independent
    // point-reference layer, same simple non-split-screen-registered pattern as
    // PILOT_ACTIONS above -- always present on both map instances, toggled by a
    // single boolean, not wired into OVERLAY_CATALOG/PANEL_FIELDS. Regenerate
    // by re-running the join in barkbeetle-catalogue whenever events.csv/
    // event_sources.csv/sources.csv change.
    const CATALOGUE_EVENTS = [{"id":"BB0001","name":"Southern Norway spruce bark-beetle outbreak (1971-1981)","country":"Norway","y0":1971,"y1":1981,"lat":60.5,"lon":11.0,"sev":"high","conf":"high","note":"Five million m3 killed in the documented 1971-1981 southern Norway outbreak.","sources":[{"title":"Host tree and bark beetle interaction during a mass outbreak of Ips typographus in Norway","url":"https://onlinelibrary.wiley.com/doi/10.1111/j.1439-0418.1983.tb03651.x","role":"primary","citation":"Bakke, A. (1983). Host tree and bark beetle interaction during a mass outbreak of Ips typographus in Norway. Zeitschrift fuer Angewandte Entomologie 96:118-125. https://doi.org/10.1111/j.1439-0418.1983.tb03651.x"},{"title":"Natural disturbances in the European forests in the 19th and 20th centuries","url":"https://onlinelibrary.wiley.com/doi/abs/10.1046/j.1365-2486.2003.00684.x","role":"secondary","citation":"Schelhaas, M.-J., Nabuurs, G.-J. & Schuck, A. (2003). Natural disturbances in the European forests in the 19th and 20th centuries. Global Change Biology 9(11):1620-1633. https://doi.org/10.1046/j.1365-2486.2003.00684.x"}]},{"id":"BB0003","name":"Kronoberg and Jonkoping counties bark-beetle outbreak (2005–2010)","country":"Sweden","y0":2005,"y1":2010,"lat":57.0,"lon":14.5,"sev":"high","conf":"high","note":"Storm Gudrun 2005; large windthrow volume.","sources":[{"title":"The Storm 'Gudrun' and the Spruce Bark Beetle in Sweden","url":"https://publications.slu.se/?file=publ/show&id=20775","role":"primary","citation":"Lindelöw, Å. & Schroeder, L.M. (2008). The Storm 'Gudrun' and the Spruce Bark Beetle in Sweden. Forstschutz Aktuell 44:5-7."},{"title":"Rise and decline of the southern Swedish storm-triggered bark-beetle outbreak","url":"","role":"secondary","citation":"Schroeder, M. et al. Recent synthesis of the rise and decline of the southern Swedish storm-triggered outbreak."}]},{"id":"BB0004","name":"Southern Swedish storm Per zone bark-beetle outbreak (2007–2011)","country":"Sweden","y0":2007,"y1":2011,"lat":56.8,"lon":14.0,"sev":"moderate","conf":"medium","note":"Storm Per 2007 added breeding substrate.","sources":[{"title":"The Storm 'Gudrun' and the Spruce Bark Beetle in Sweden","url":"https://publications.slu.se/?file=publ/show&id=20775","role":"primary","citation":"Lindelöw, Å. & Schroeder, L.M. (2008). The Storm 'Gudrun' and the Spruce Bark Beetle in Sweden. Forstschutz Aktuell 44:5-7."},{"title":"Åtgärder vid granbarkborreangrepp efter stormfällning inom skyddade områden i Jönköpings län 2005-2015 (Measures against spruce bark beetle attacks after storm-felling in protected areas in Jönköping county 2005-2015)","url":"https://www.lansstyrelsen.se/download/18.1b1d393819324610c374acad/1732520357149/%C3%85tg%C3%A4rder%20vid%20granbarkborreangrepp%20efter%20stormf%C3%A4llning%20inom%20skyddade%20omr%C3%A5den%20i%20J%C3%B6nk%C3%B6pings%20l%C3%A4n%202005-2015.pdf","role":"primary","citation":"Länsstyrelsen Jönköping (2015). Åtgärder vid granbarkborreangrepp efter stormfällning inom skyddade områden i Jönköpings län 2005-2015."},{"title":"Aerial inventory of beetle-attacked spruce, Skaraborg, April 2007","url":"https://www.sverigesradio.se/artikel/1317790","role":"secondary","citation":"Sveriges Radio P4 Skaraborg (2007, April)."}]},{"id":"BB0005","name":"Southern Finland storm-damaged spruce forests bark-beetle outbreak (2012–2013)","country":"Finland","y0":2012,"y1":2013,"lat":61.0,"lon":25.0,"sev":"minor","conf":"medium","note":"Warm summers and storm damage.","sources":[{"title":"Changing climatic drivers of European spruce bark beetle outbreaks: a comparison of locations around the Northern Baltic Sea","url":"https://www.silvafennica.fi/article/23003","role":"primary","citation":"Tikkanen, O.-P. & Lehtonen, I. (2023). Changing climatic drivers of European spruce bark beetle outbreaks: a comparison of locations around the Northern Baltic Sea. Silva Fennica 57(3):23003. https://doi.org/10.14214/sf.23003"},{"title":"Analyzing the environmental risk factors of European spruce bark beetle damage at the local scale","url":"https://link.springer.com/article/10.1007/s10342-024-01662-4","role":"secondary","citation":"Huo, L., Persson, H.J. & Lindberg, E. (2024). Analyzing the environmental risk factors of European spruce bark beetle damage at the local scale. European Journal of Forest Research 143:985-1000. https://doi.org/10.1007/s10342-024-01662-4"}]},{"id":"BB0006","name":"East-central Sweden bark-beetle outbreak (2017–2024)","country":"Sweden","y0":2017,"y1":2024,"lat":59.5,"lon":16.5,"sev":"extreme","conf":"high","note":"2018 drought and heat; susceptible spruce.","sources":[{"title":"Changing climatic drivers of European spruce bark beetle outbreaks: a comparison of locations around the Northern Baltic Sea","url":"https://www.silvafennica.fi/article/23003","role":"primary","citation":"Tikkanen, O.-P. & Lehtonen, I. (2023). Changing climatic drivers of European spruce bark beetle outbreaks: a comparison of locations around the Northern Baltic Sea. Silva Fennica 57(3):23003. https://doi.org/10.14214/sf.23003"},{"title":"Granbarkborre (spruce bark beetle) risk maps and monitoring","url":"https://www.skogsstyrelsen.se/granbarkborre","role":"secondary","citation":"Skogsstyrelsen. Granbarkborre risk maps and monitoring. https://www.skogsstyrelsen.se/granbarkborre and https://www.skogsstyrelsen.se/sjalvservice/karttjanster/kartor-over-skador-pa-skog/ (CC0 open geodata)."}]},{"id":"BB0007","name":"Gotaland interior bark-beetle outbreak (2018–2023)","country":"Sweden","y0":2018,"y1":2023,"lat":57.5,"lon":14.0,"sev":"high","conf":"high","note":"Drought and high beetle pressure.","sources":[{"title":"Changing climatic drivers of European spruce bark beetle outbreaks: a comparison of locations around the Northern Baltic Sea","url":"https://www.silvafennica.fi/article/23003","role":"primary","citation":"Tikkanen, O.-P. & Lehtonen, I. (2023). Changing climatic drivers of European spruce bark beetle outbreaks: a comparison of locations around the Northern Baltic Sea. Silva Fennica 57(3):23003. https://doi.org/10.14214/sf.23003"},{"title":"Granbarkborre (spruce bark beetle) risk maps and monitoring","url":"https://www.skogsstyrelsen.se/granbarkborre","role":"secondary","citation":"Skogsstyrelsen. Granbarkborre risk maps and monitoring. https://www.skogsstyrelsen.se/granbarkborre and https://www.skogsstyrelsen.se/sjalvservice/karttjanster/kartor-over-skador-pa-skog/ (CC0 open geodata)."}]},{"id":"BB0008","name":"Svealand spruce belt bark-beetle outbreak (2018–2023)","country":"Sweden","y0":2018,"y1":2023,"lat":59.5,"lon":15.5,"sev":"high","conf":"medium","note":"Drought and warm summers.","sources":[{"title":"Changing climatic drivers of European spruce bark beetle outbreaks: a comparison of locations around the Northern Baltic Sea","url":"https://www.silvafennica.fi/article/23003","role":"primary","citation":"Tikkanen, O.-P. & Lehtonen, I. (2023). Changing climatic drivers of European spruce bark beetle outbreaks: a comparison of locations around the Northern Baltic Sea. Silva Fennica 57(3):23003. https://doi.org/10.14214/sf.23003"},{"title":"Skogsvård och insats mot granbarkborren (forest care and action against spruce bark beetle)","url":"https://morakommun.se/arkiv/nyhetsarkiv/nyheter/2022-12-15-skogsvard-och-insats-mot-granbarkborren.html","role":"primary","citation":"Mora kommun (2022, 15 December). Skogsvård och insats mot granbarkborren."},{"title":"Granbarkborre (spruce bark beetle) risk maps and monitoring","url":"https://www.skogsstyrelsen.se/granbarkborre","role":"secondary","citation":"Skogsstyrelsen. Granbarkborre risk maps and monitoring. https://www.skogsstyrelsen.se/granbarkborre and https://www.skogsstyrelsen.se/sjalvservice/karttjanster/kartor-over-skador-pa-skog/ (CC0 open geodata)."},{"title":"Skogsskador i Sverige 2022 / 2023 (Forest damage in Sweden 2022/2023)","url":"https://www.skogsstyrelsen.se/globalassets/om-oss/rapporter/rapporter-2024/rapport-2024-06-skogsskador-i-sverige-2023.pdf","role":"secondary","citation":"Skogsstyrelsen (2023, 2024). Skogsskador i Sverige 2022; Skogsskador i Sverige 2023."}]},{"id":"BB0009","name":"Karula and Otepaa post-windthrow bark-beetle spread (2016-2020)","country":"Estonia","y0":2016,"y1":2020,"lat":58.6,"lon":25.0,"sev":"moderate","conf":"high","note":"Spread from 2016 windthrows assessed by forest inventory in 2020.","sources":[{"title":"Drivers behind the spatial dispersion of European spruce bark beetle (Ips typographus) infestation in protected areas in Estonia, four years after a major storm","url":"https://doi.org/10.1016/j.foreco.2024.122469","role":"primary","citation":"Vodde, F., Ait, K., Orumaa, A., Jõgiste, K., Kaart, T., Kiviste, A., Õunap, H., Stanturf, J.A. & Metslaid, M. (2025). Drivers behind the spatial dispersion of European spruce bark beetle (Ips typographus) infestation in protected areas in Estonia, four years after a major storm. Forest Ecology and Management 578:122469. https://doi.org/10.1016/j.foreco.2024.122469"},{"title":"Kuuse-kooreürask vallutab järk-järgult Eesti kuusikuid (Spruce bark beetle gradually conquers Estonian spruce stands)","url":"https://novaator.err.ee/1089240/kuuse-kooreurask-vallutab-jark-jargult-eesti-kuusikuid","role":"secondary","citation":"ERR Novaator (2020). Kuuse-kooreürask vallutab järk-järgult Eesti kuusikuid."},{"title":"Estonian Environment Agency 2023 monitoring result; spruce bark beetle damage increasingly visible in Virumaa","url":"https://keskkonnaagentuur.ee/node/1647","role":"secondary","citation":"Estonian Environment Agency (2023). 2023 monitoring result. ERR News (2023). Spruce bark beetle damage increasingly visible in Virumaa. https://news.err.ee/1609093430/"}]},{"id":"BB0010","name":"Latvia spruce bark-beetle emergency (2023)","country":"Latvia","y0":2023,"y1":2023,"lat":56.8,"lon":25.0,"sev":"high","conf":"high","note":"Official emergency declared in 32 municipalities for April-June 2023.","sources":[{"title":"Emergency situation declared in 32 municipalities to control Ips typographus","url":"https://www.zm.gov.lv/lv/jaunums/32-novados-latvija-no-1-aprila-lidz-30-junijam-izsludinata-arkarteja-situacija-eglu-mezaudzem-loti-bistama-eglu-astonzobu-mizgrauza-ierobezosanai","role":"primary","citation":"Latvian Ministry of Agriculture (2023). Emergency situation declared in 32 municipalities to control the spruce bark beetle."}]},{"id":"BB0011","name":"Lithuania spruce bark-beetle outbreak (2021-2023)","country":"Lithuania","y0":2021,"y1":2023,"lat":55.3,"lon":23.5,"sev":"high","conf":"high","note":"Official monitoring records a steep rise from 2021 and more than 1 million m3 of dead spruce in 2023.","sources":[{"title":"Spruce bark beetle damage and population forecast for 2024","url":"https://am.lrv.lt/lt/naujienos/zievegrauzio-tipografo-daromos-zalos-ir-populiacijos-kitimo-prognoze-2024-metams","role":"primary","citation":"Ministry of Environment of Lithuania (2024). Spruce bark beetle damage and population forecast for 2024."}]},{"id":"BB0012","name":"Sudetes bark-beetle outbreak (1981–1987)","country":"Poland","y0":1981,"y1":1987,"lat":50.8,"lon":15.7,"sev":"high","conf":"medium","note":"Air pollution, drought and weakened spruce.","sources":[{"title":"Ips typographus (IPSXTY) datasheet","url":"https://gd.eppo.int/taxon/IPSXTY/datasheet","role":"primary","citation":"EPPO Global Database. Ips typographus (IPSXTY) datasheet. https://gd.eppo.int/taxon/IPSXTY/datasheet (last updated 2020-09-01)."},{"title":"Natural disturbances in the European forests in the 19th and 20th centuries","url":"https://onlinelibrary.wiley.com/doi/abs/10.1046/j.1365-2486.2003.00684.x","role":"secondary","citation":"Schelhaas, M.-J., Nabuurs, G.-J. & Schuck, A. (2003). Natural disturbances in the European forests in the 19th and 20th centuries. Global Change Biology 9(11):1620-1633. https://doi.org/10.1046/j.1365-2486.2003.00684.x"},{"title":"Polish forest-health literature and European historical reviews of Sudetes/Carpathian spruce decline","url":"","role":"secondary","citation":"Polish forest-health literature and European historical reviews of Sudetes and Carpathian spruce decline."}]},{"id":"BB0014","name":"Polish Tatras bark-beetle outbreak (1990–2000)","country":"Poland","y0":1990,"y1":2000,"lat":49.25,"lon":19.95,"sev":"moderate","conf":"high","note":"Wind disturbance and susceptible spruce.","sources":[{"title":"Factors influencing the wind-bark beetles' disturbance system in the course of an Ips typographus outbreak in the Tatra Mountains","url":"https://doi.org/10.1016/j.foreco.2013.10.020","role":"primary","citation":"Mezei, P., Grodzki, W., Blaženec, M. & Jakuš, R. (2014). Factors influencing the wind-bark beetles' disturbance system in the course of an Ips typographus outbreak in the Tatra Mountains. Forest Ecology and Management 312:67-77. https://doi.org/10.1016/j.foreco.2013.10.020"},{"title":"Storms, temperature maxima and the Eurasian spruce bark beetle Ips typographus - an infernal trio in Norway spruce forests of the Central European High Tatra Mountains","url":"https://doi.org/10.1016/j.agrformet.2017.04.004","role":"secondary","citation":"Mezei, P. et al. (2017). Storms, temperature maxima and the Eurasian spruce bark beetle Ips typographus - an infernal trio in Norway spruce forests of the Central European High Tatra Mountains. Agricultural and Forest Meteorology 242:85-95. https://doi.org/10.1016/j.agrformet.2017.04.004"}]},{"id":"BB0015","name":"Belarus spruce bark-beetle mass mortality episode (1993-1997)","country":"Belarus","y0":1993,"y1":1997,"lat":53.5,"lon":27.5,"sev":"high","conf":"low","note":"Country-level episode reported in Belarusian forestry literature for 1993-1997.","sources":[{"title":"Specific features of the mass dying of spruce trees in forests of Belarus","url":"https://gd.eppo.int/taxon/IPSXTY/distribution/BY","role":"primary","citation":"Fedorov, N. I., Sarnatskii, V. V. & Kharitonov, V. A. (1998). Specific features of the mass dying of spruce trees in forests of Belarus. Lesovedenie 6:12-23."},{"title":"Natural disturbances in the European forests in the 19th and 20th centuries","url":"https://onlinelibrary.wiley.com/doi/abs/10.1046/j.1365-2486.2003.00684.x","role":"secondary","citation":"Schelhaas, M.-J., Nabuurs, G.-J. & Schuck, A. (2003). Natural disturbances in the European forests in the 19th and 20th centuries. Global Change Biology 9(11):1620-1633. https://doi.org/10.1046/j.1365-2486.2003.00684.x"}]},{"id":"BB0017","name":"Bialowieza Forest bark-beetle outbreak (2012–2022)","country":"Poland","y0":2012,"y1":2022,"lat":52.7,"lon":23.85,"sev":"high","conf":"high","note":"Drought, mature spruce and population buildup.","sources":[{"title":"Changes of Norway Spruce Health in the Białowieża Forest (CE Europe) in 2013-2019 during a Bark Beetle Infestation, Studied with Landsat Imagery","url":"https://doi.org/10.3390/f12010034","role":"primary","citation":"Bryk, M., Kołodziej, B. & Pliszka, R. (2021). Changes of Norway Spruce Health in the Białowieża Forest (CE Europe) in 2013-2019 during a Bark Beetle Infestation, Studied with Landsat Imagery. Forests 12:34. https://doi.org/10.3390/f12010034"},{"title":"Mass outbreaks and factors related to the spatial dynamics of spruce bark beetle (Ips typographus) dieback considering diverse management regimes in the Białowieża forest","url":"https://doi.org/10.1016/j.foreco.2021.119530","role":"secondary","citation":"Kamińska, A. et al. (2021). Mass outbreaks and factors related to the spatial dynamics of spruce bark beetle (Ips typographus) dieback considering diverse management regimes in the Białowieża forest. Forest Ecology and Management 498:119530. https://doi.org/10.1016/j.foreco.2021.119530"},{"title":"Impact of forest stands dieback as a result of bark beetle (Ips typographus) outbreak on the fire risk in the Bialowieza Forest","url":"https://doi.org/10.48538/lpb-2023-0001","role":"secondary","citation":"Szczygiel, R., Kwiatkowski, M. & Tyburski, L. (2023). Impact of forest stands dieback as a result of bark beetle (Ips typographus) outbreak on the fire risk in the Bialowieza Forest. Forest Research Papers 83(1):1-15. https://doi.org/10.48538/lpb-2023-0001"}]},{"id":"BB0019","name":"Lower Silesian and Opole spruce forests bark-beetle outbreak (2015–2022)","country":"Poland","y0":2015,"y1":2022,"lat":50.9,"lon":17.5,"sev":"high","conf":"medium","note":"Drought and spillover from Czechia/Germany.","sources":[{"title":"Living with bark beetles: impacts, outlook and management options","url":"https://efi.int/sites/default/files/files/publication-bank/2019/efi_fstp_8_2019.pdf","role":"primary","citation":"Hlásny, T., Krokene, P., Liebhold, A. et al. (2019). Living with bark beetles: impacts, outlook and management options. European Forest Institute, From Science to Policy 8."},{"title":"Forest Europe and national forest-health reports on the post-2015 Central European epidemic","url":"","role":"secondary","citation":"Forest Europe and national forest-health reports describing the post-2015 Central European epidemic."}]},{"id":"BB0020","name":"Malopolskie and Podkarpackie Carpathians bark-beetle outbreak (2015–2022)","country":"Poland","y0":2015,"y1":2022,"lat":49.6,"lon":21.5,"sev":"high","conf":"medium","note":"Drought and regional population pressure.","sources":[{"title":"Living with bark beetles: impacts, outlook and management options","url":"https://efi.int/sites/default/files/files/publication-bank/2019/efi_fstp_8_2019.pdf","role":"primary","citation":"Hlásny, T., Krokene, P., Liebhold, A. et al. (2019). Living with bark beetles: impacts, outlook and management options. European Forest Institute, From Science to Policy 8."},{"title":"Forest Europe and national forest-health reports on the post-2015 Central European epidemic","url":"","role":"secondary","citation":"Forest Europe and national forest-health reports describing the post-2015 Central European epidemic."}]},{"id":"BB0022","name":"Czech Republic bark-beetle outbreak (1992-1995)","country":"Czechia","y0":1992,"y1":1995,"lat":49.05,"lon":13.65,"sev":"high","conf":"high","note":"Documented 1992-1995 Czech episode with 5.59 million m3 of infested wood.","sources":[{"title":"Meteorological and climatological triggers of notable past and present bark beetle outbreaks in the Czech Republic","url":"https://doi.org/10.5194/cp-18-2155-2022","role":"primary","citation":"Brazdil, R. et al. (2022). Meteorological and climatological triggers of notable past and present bark beetle outbreaks in the Czech Republic. Climate of the Past 18:2155-2180. https://doi.org/10.5194/cp-18-2155-2022"}]},{"id":"BB0025","name":"Jeseniky and northern Moravia bark-beetle outbreak (1991–1999)","country":"Czechia","y0":1991,"y1":1999,"lat":50.1,"lon":17.2,"sev":"moderate","conf":"medium","note":"Wind damage and warming.","sources":[{"title":"Ips typographus (IPSXTY) datasheet","url":"https://gd.eppo.int/taxon/IPSXTY/datasheet","role":"primary","citation":"EPPO Global Database. Ips typographus (IPSXTY) datasheet. https://gd.eppo.int/taxon/IPSXTY/datasheet (last updated 2020-09-01)."},{"title":"Devastating outbreak of bark beetles in the Czech Republic: drivers, impacts, and management implications","url":"https://doi.org/10.1016/j.foreco.2021.119075","role":"secondary","citation":"Hlásny, T. et al. (2021). Devastating outbreak of bark beetles in the Czech Republic: drivers, impacts, and management implications. Forest Ecology and Management 490:119075. https://doi.org/10.1016/j.foreco.2021.119075"}]},{"id":"BB0026","name":"Sumava-Bavarian Forest border zone bark-beetle outbreak (2003–2016)","country":"Czechia","y0":2003,"y1":2016,"lat":49.05,"lon":13.35,"sev":"moderate","conf":"high","note":"Windthrow, mature spruce and non-intervention zones.","sources":[{"title":"Spatio-temporal infestation patterns of Ips typographus (L.) in the Bavarian Forest National Park, Germany","url":"https://doi.org/10.1016/j.ecolind.2012.07.026","role":"primary","citation":"Lausch, A., Heurich, M. & Fahse, L. (2013). Spatio-temporal infestation patterns of Ips typographus (L.) in the Bavarian Forest National Park, Germany. Ecological Indicators 31:73-81. https://doi.org/10.1016/j.ecolind.2012.07.026"},{"title":"Do bark beetle outbreaks amplify or dampen future bark beetle disturbances in Central Europe?","url":"https://doi.org/10.1111/1365-2745.13502","role":"secondary","citation":"Sommerfeld, A. et al. (2021). Do bark beetle outbreaks amplify or dampen future bark beetle disturbances in Central Europe? Journal of Ecology 109(6). https://doi.org/10.1111/1365-2745.13502"}]},{"id":"BB0027","name":"Northern Moravia and Silesia bark-beetle outbreak (2003–2016)","country":"Czechia","y0":2003,"y1":2016,"lat":49.8,"lon":18.2,"sev":"moderate","conf":"medium","note":"Recurrent drought and host stress.","sources":[{"title":"Devastating outbreak of bark beetles in the Czech Republic: drivers, impacts, and management implications","url":"https://doi.org/10.1016/j.foreco.2021.119075","role":"primary","citation":"Hlásny, T. et al. (2021). Devastating outbreak of bark beetles in the Czech Republic: drivers, impacts, and management implications. Forest Ecology and Management 490:119075. https://doi.org/10.1016/j.foreco.2021.119075"},{"title":"Patterns and impacts of an unprecedented outbreak of bark beetles in Central Europe: a glimpse into the future?","url":"https://doi.org/10.1016/j.fecs.2024.100243","role":"secondary","citation":"Washaya, P., Modlinger, R., Tyšer, D. et al. (2024). Patterns and impacts of an unprecedented outbreak of bark beetles in Central Europe: a glimpse into the future? Forest Ecosystems 11(6):100243. https://doi.org/10.1016/j.fecs.2024.100243"}]},{"id":"BB0028","name":"Northern Moravia-Silesia epicentre bark-beetle outbreak (2015–2017)","country":"Czechia","y0":2015,"y1":2017,"lat":49.75,"lon":18.1,"sev":"moderate","conf":"high","note":"Drought and accumulated beetle populations.","sources":[{"title":"Devastating outbreak of bark beetles in the Czech Republic: drivers, impacts, and management implications","url":"https://doi.org/10.1016/j.foreco.2021.119075","role":"primary","citation":"Hlásny, T. et al. (2021). Devastating outbreak of bark beetles in the Czech Republic: drivers, impacts, and management implications. Forest Ecology and Management 490:119075. https://doi.org/10.1016/j.foreco.2021.119075"},{"title":"Patterns and impacts of an unprecedented outbreak of bark beetles in Central Europe: a glimpse into the future?","url":"https://doi.org/10.1016/j.fecs.2024.100243","role":"secondary","citation":"Washaya, P., Modlinger, R., Tyšer, D. et al. (2024). Patterns and impacts of an unprecedented outbreak of bark beetles in Central Europe: a glimpse into the future? Forest Ecosystems 11(6):100243. https://doi.org/10.1016/j.fecs.2024.100243"}]},{"id":"BB0029","name":"Vysocina expansion zone bark-beetle outbreak (2017–2021)","country":"Czechia","y0":2017,"y1":2021,"lat":49.4,"lon":15.6,"sev":"extreme","conf":"high","note":"Drought, heat and rapid population spread.","sources":[{"title":"Devastating outbreak of bark beetles in the Czech Republic: drivers, impacts, and management implications","url":"https://doi.org/10.1016/j.foreco.2021.119075","role":"primary","citation":"Hlásny, T. et al. (2021). Devastating outbreak of bark beetles in the Czech Republic: drivers, impacts, and management implications. Forest Ecology and Management 490:119075. https://doi.org/10.1016/j.foreco.2021.119075"},{"title":"Patterns and impacts of an unprecedented outbreak of bark beetles in Central Europe: a glimpse into the future?","url":"https://doi.org/10.1016/j.fecs.2024.100243","role":"secondary","citation":"Washaya, P., Modlinger, R., Tyšer, D. et al. (2024). Patterns and impacts of an unprecedented outbreak of bark beetles in Central Europe: a glimpse into the future? Forest Ecosystems 11(6):100243. https://doi.org/10.1016/j.fecs.2024.100243"},{"title":"Kůrovcová mapa (Czech bark-beetle risk/damage map)","url":"https://www.kurovcovamapa.cz/","role":"secondary","citation":"ÚHÚL. Kůrovcová mapa. https://www.kurovcovamapa.cz/"},{"title":"DEFID2: Database of European Forest Insect and Disease Disturbances","url":"https://onlinelibrary.wiley.com/doi/full/10.1111/gcb.16912","role":"secondary","citation":"Forzieri, G. et al. (2023). The Database of European Forest Insect and Disease Disturbances: DEFID2. Global Change Biology 29. https://doi.org/10.1111/gcb.16912 (data hosted via JRC Forest Information System for Europe)."}]},{"id":"BB0030","name":"Olomouc and Zlin regions bark-beetle outbreak (2017–2021)","country":"Czechia","y0":2017,"y1":2021,"lat":49.5,"lon":17.5,"sev":"extreme","conf":"high","note":"Drought and east-to-west spread.","sources":[{"title":"Devastating outbreak of bark beetles in the Czech Republic: drivers, impacts, and management implications","url":"https://doi.org/10.1016/j.foreco.2021.119075","role":"primary","citation":"Hlásny, T. et al. (2021). Devastating outbreak of bark beetles in the Czech Republic: drivers, impacts, and management implications. Forest Ecology and Management 490:119075. https://doi.org/10.1016/j.foreco.2021.119075"},{"title":"Patterns and impacts of an unprecedented outbreak of bark beetles in Central Europe: a glimpse into the future?","url":"https://doi.org/10.1016/j.fecs.2024.100243","role":"secondary","citation":"Washaya, P., Modlinger, R., Tyšer, D. et al. (2024). Patterns and impacts of an unprecedented outbreak of bark beetles in Central Europe: a glimpse into the future? Forest Ecosystems 11(6):100243. https://doi.org/10.1016/j.fecs.2024.100243"},{"title":"DEFID2: Database of European Forest Insect and Disease Disturbances","url":"https://onlinelibrary.wiley.com/doi/full/10.1111/gcb.16912","role":"secondary","citation":"Forzieri, G. et al. (2023). The Database of European Forest Insect and Disease Disturbances: DEFID2. Global Change Biology 29. https://doi.org/10.1111/gcb.16912 (data hosted via JRC Forest Information System for Europe)."}]},{"id":"BB0031","name":"Bohemian-Moravian uplands and central Czechia bark-beetle outbreak (2018–2023)","country":"Czechia","y0":2018,"y1":2023,"lat":49.6,"lon":15.9,"sev":"high","conf":"high","note":"Drought, heat and depleted host resistance.","sources":[{"title":"Patterns and impacts of an unprecedented outbreak of bark beetles in Central Europe: a glimpse into the future?","url":"https://doi.org/10.1016/j.fecs.2024.100243","role":"primary","citation":"Washaya, P., Modlinger, R., Tyšer, D. et al. (2024). Patterns and impacts of an unprecedented outbreak of bark beetles in Central Europe: a glimpse into the future? Forest Ecosystems 11(6):100243. https://doi.org/10.1016/j.fecs.2024.100243"},{"title":"Kůrovcová mapa (Czech bark-beetle risk/damage map)","url":"https://www.kurovcovamapa.cz/","role":"secondary","citation":"ÚHÚL. Kůrovcová mapa. https://www.kurovcovamapa.cz/"},{"title":"DEFID2: Database of European Forest Insect and Disease Disturbances","url":"https://onlinelibrary.wiley.com/doi/full/10.1111/gcb.16912","role":"secondary","citation":"Forzieri, G. et al. (2023). The Database of European Forest Insect and Disease Disturbances: DEFID2. Global Change Biology 29. https://doi.org/10.1111/gcb.16912 (data hosted via JRC Forest Information System for Europe)."}]},{"id":"BB0032","name":"Residual Czech hotspots bark-beetle outbreak (2021–2025)","country":"Czechia","y0":2021,"y1":2025,"lat":49.6,"lon":15.9,"sev":"minor","conf":"medium","note":"Residual populations and local drought.","sources":[{"title":"Patterns and impacts of an unprecedented outbreak of bark beetles in Central Europe: a glimpse into the future?","url":"https://doi.org/10.1016/j.fecs.2024.100243","role":"primary","citation":"Washaya, P., Modlinger, R., Tyšer, D. et al. (2024). Patterns and impacts of an unprecedented outbreak of bark beetles in Central Europe: a glimpse into the future? Forest Ecosystems 11(6):100243. https://doi.org/10.1016/j.fecs.2024.100243"},{"title":"Kůrovcová mapa (Czech bark-beetle risk/damage map)","url":"https://www.kurovcovamapa.cz/","role":"secondary","citation":"ÚHÚL. Kůrovcová mapa. https://www.kurovcovamapa.cz/"}]},{"id":"BB0033","name":"High Tatras bark-beetle outbreak (1990–2000)","country":"Slovakia","y0":1990,"y1":2000,"lat":49.15,"lon":20.05,"sev":"moderate","conf":"high","note":"Wind disturbance and susceptible spruce.","sources":[{"title":"Factors influencing the wind-bark beetles' disturbance system in the course of an Ips typographus outbreak in the Tatra Mountains","url":"https://doi.org/10.1016/j.foreco.2013.10.020","role":"primary","citation":"Mezei, P., Grodzki, W., Blaženec, M. & Jakuš, R. (2014). Factors influencing the wind-bark beetles' disturbance system in the course of an Ips typographus outbreak in the Tatra Mountains. Forest Ecology and Management 312:67-77. https://doi.org/10.1016/j.foreco.2013.10.020"},{"title":"Storms, temperature maxima and the Eurasian spruce bark beetle Ips typographus - an infernal trio in Norway spruce forests of the Central European High Tatra Mountains","url":"https://doi.org/10.1016/j.agrformet.2017.04.004","role":"secondary","citation":"Mezei, P. et al. (2017). Storms, temperature maxima and the Eurasian spruce bark beetle Ips typographus - an infernal trio in Norway spruce forests of the Central European High Tatra Mountains. Agricultural and Forest Meteorology 242:85-95. https://doi.org/10.1016/j.agrformet.2017.04.004"}]},{"id":"BB0034","name":"High Tatras after windstorm Alzbeta bark-beetle outbreak (2004–2014)","country":"Slovakia","y0":2004,"y1":2014,"lat":49.15,"lon":20.05,"sev":"high","conf":"high","note":"Large 2004 windthrow and warm years.","sources":[{"title":"Post-disaster Forest Management and Bark Beetle Outbreak in Tatra National Park, Slovakia","url":"https://doi.org/10.1659/MRD-JOURNAL-D-13-00017.1","role":"primary","citation":"Nikolov, C. et al. (2014). Post-disaster Forest Management and Bark Beetle Outbreak in Tatra National Park, Slovakia. Mountain Research and Development 34(4):326-335. https://doi.org/10.1659/MRD-JOURNAL-D-13-00017.1"},{"title":"Slovak and Tatra disturbance literature on the 2004 windstorm and subsequent Ips typographus outbreak","url":"https://web.nlcsk.org/","role":"secondary","citation":"Slovak and Tatra disturbance literature on the 2004 windstorm and subsequent Ips typographus outbreak; institutional home: Národné lesnícke centrum, https://web.nlcsk.org/"}]},{"id":"BB0035","name":"High and Low Tatras bark-beetle outbreak (2018–2023)","country":"Slovakia","y0":2018,"y1":2023,"lat":49.0,"lon":19.8,"sev":"high","conf":"medium","note":"Drought and disturbance legacies.","sources":[{"title":"Living with bark beetles: impacts, outlook and management options","url":"https://efi.int/sites/default/files/files/publication-bank/2019/efi_fstp_8_2019.pdf","role":"primary","citation":"Hlásny, T., Krokene, P., Liebhold, A. et al. (2019). Living with bark beetles: impacts, outlook and management options. European Forest Institute, From Science to Policy 8."},{"title":"Slovak and Tatra disturbance literature on the 2004 windstorm and subsequent Ips typographus outbreak","url":"https://web.nlcsk.org/","role":"secondary","citation":"Slovak and Tatra disturbance literature on the 2004 windstorm and subsequent Ips typographus outbreak; institutional home: Národné lesnícke centrum, https://web.nlcsk.org/"}]},{"id":"BB0036","name":"Orava and Kysuce bark-beetle outbreak (2018–2023)","country":"Slovakia","y0":2018,"y1":2023,"lat":49.35,"lon":19.3,"sev":"moderate","conf":"medium","note":"Drought and regional spread.","sources":[{"title":"Living with bark beetles: impacts, outlook and management options","url":"https://efi.int/sites/default/files/files/publication-bank/2019/efi_fstp_8_2019.pdf","role":"primary","citation":"Hlásny, T., Krokene, P., Liebhold, A. et al. (2019). Living with bark beetles: impacts, outlook and management options. European Forest Institute, From Science to Policy 8."}]},{"id":"BB0037","name":"Slovenske rudohorie bark-beetle outbreak (2018–2023)","country":"Slovakia","y0":2018,"y1":2023,"lat":48.7,"lon":20.2,"sev":"moderate","conf":"medium","note":"Drought and warm summers.","sources":[{"title":"Living with bark beetles: impacts, outlook and management options","url":"https://efi.int/sites/default/files/files/publication-bank/2019/efi_fstp_8_2019.pdf","role":"primary","citation":"Hlásny, T., Krokene, P., Liebhold, A. et al. (2019). Living with bark beetles: impacts, outlook and management options. European Forest Institute, From Science to Policy 8."}]},{"id":"BB0038","name":"Swiss post-Vivian bark-beetle epidemic (1990-1997)","country":"Switzerland","y0":1990,"y1":1997,"lat":46.9,"lon":8.3,"sev":"high","conf":"high","note":"Direct Swiss chronology: Vivian-triggered epidemic peaked in 1992/93 and ended in 1997.","sources":[{"title":"Ursachen und Verlauf der Buchdrucker-Epidemien (Ips typographus L.) in der Schweiz von 1984 bis 1999 [Causes and progress of the eight-toothed spruce bark beetle epidemics (Ips typographus) in Switzerland from 1984 to 1999]","url":"https://doi.org/10.3188/szf.2003.0437","role":"primary","citation":"Meier, F., Gall, R. & Forster, B. (2003). Ursachen und Verlauf der Buchdrucker-Epidemien (Ips typographus L.) in der Schweiz von 1984 bis 1999. Schweizerische Zeitschrift fur Forstwesen 154(11):437-441. https://doi.org/10.3188/szf.2003.0437"},{"title":"Ecology and management of the spruce bark beetle Ips typographus - a review of recent research","url":"https://doi.org/10.1016/j.foreco.2004.07.018","role":"secondary","citation":"Wermelinger, B. (2004). Ecology and management of the spruce bark beetle Ips typographus - a review of recent research. Forest Ecology and Management 202(1-3):67-82. https://doi.org/10.1016/j.foreco.2004.07.018"}]},{"id":"BB0039","name":"Gandberg bark-beetle outbreak (1992–1997)","country":"Switzerland","y0":1992,"y1":1997,"lat":47.05,"lon":9.15,"sev":"minor","conf":"high","note":"Uncleared Vivian windthrow.","sources":[{"title":"Tree regeneration in a Norway spruce snag stand after tree die-back caused by Ips typographus","url":"https://www.researchgate.net/publication/262974333_Tree_regeneration_in_a_Norway_spruce_snag_after_tree_die-back_caused_by_Ips_typographus","role":"primary","citation":"Kupferschmid, A.D., Schönenberger, W. & Wasem, U. (2002). Tree regeneration in a Norway spruce snag stand after tree die-back caused by Ips typographus. Forest Snow and Landscape Research 77(1/2):149-160."}]},{"id":"BB0040","name":"Bavaria and Baden-Wurttemberg after Vivian bark-beetle outbreak (1990–1994)","country":"Germany","y0":1990,"y1":1994,"lat":48.5,"lon":10.0,"sev":"moderate","conf":"medium","note":"Storm Vivian and warm summers.","sources":[{"title":"Ips typographus (IPSXTY) datasheet","url":"https://gd.eppo.int/taxon/IPSXTY/datasheet","role":"primary","citation":"EPPO Global Database. Ips typographus (IPSXTY) datasheet. https://gd.eppo.int/taxon/IPSXTY/datasheet (last updated 2020-09-01)."},{"title":"Ecology and management of the spruce bark beetle Ips typographus - a review of recent research","url":"https://doi.org/10.1016/j.foreco.2004.07.018","role":"secondary","citation":"Wermelinger, B. (2004). Ecology and management of the spruce bark beetle Ips typographus - a review of recent research. Forest Ecology and Management 202(1-3):67-82. https://doi.org/10.1016/j.foreco.2004.07.018"}]},{"id":"BB0041","name":"Switzerland after storm Lothar bark-beetle outbreak (1999–2007)","country":"Switzerland","y0":1999,"y1":2007,"lat":47.0,"lon":8.0,"sev":"high","conf":"high","note":"Storm Lothar, then 2003 heatwave.","sources":[{"title":"Ecology and management of the spruce bark beetle Ips typographus - a review of recent research","url":"https://doi.org/10.1016/j.foreco.2004.07.018","role":"primary","citation":"Wermelinger, B. (2004). Ecology and management of the spruce bark beetle Ips typographus - a review of recent research. Forest Ecology and Management 202(1-3):67-82. https://doi.org/10.1016/j.foreco.2004.07.018"},{"title":"Swiss, French and German post-Lothar bark-beetle studies","url":"","role":"secondary","citation":"Swiss, French and German post-Lothar studies summarized in Wermelinger (2004) and later reviews."}]},{"id":"BB0042","name":"Baden-Wurttemberg after Lothar bark-beetle outbreak (1999–2007)","country":"Germany","y0":1999,"y1":2007,"lat":48.5,"lon":9.0,"sev":"high","conf":"high","note":"Storm Lothar and 2003 drought.","sources":[{"title":"Ecology and management of the spruce bark beetle Ips typographus - a review of recent research","url":"https://doi.org/10.1016/j.foreco.2004.07.018","role":"primary","citation":"Wermelinger, B. (2004). Ecology and management of the spruce bark beetle Ips typographus - a review of recent research. Forest Ecology and Management 202(1-3):67-82. https://doi.org/10.1016/j.foreco.2004.07.018"},{"title":"Swiss, French and German post-Lothar bark-beetle studies","url":"","role":"secondary","citation":"Swiss, French and German post-Lothar studies summarized in Wermelinger (2004) and later reviews."}]},{"id":"BB0043","name":"Bavaria after Lothar bark-beetle outbreak (1999–2007)","country":"Germany","y0":1999,"y1":2007,"lat":48.8,"lon":11.5,"sev":"high","conf":"medium","note":"Storm Lothar and hot summers.","sources":[{"title":"Ecology and management of the spruce bark beetle Ips typographus - a review of recent research","url":"https://doi.org/10.1016/j.foreco.2004.07.018","role":"primary","citation":"Wermelinger, B. (2004). Ecology and management of the spruce bark beetle Ips typographus - a review of recent research. Forest Ecology and Management 202(1-3):67-82. https://doi.org/10.1016/j.foreco.2004.07.018"},{"title":"Swiss, French and German post-Lothar bark-beetle studies","url":"","role":"secondary","citation":"Swiss, French and German post-Lothar studies summarized in Wermelinger (2004) and later reviews."}]},{"id":"BB0044","name":"Alsace bark-beetle outbreak (1999–2007)","country":"France","y0":1999,"y1":2007,"lat":48.5,"lon":7.0,"sev":"high","conf":"high","note":"Storm Lothar and 2003 heatwave.","sources":[{"title":"Ecology and management of the spruce bark beetle Ips typographus - a review of recent research","url":"https://doi.org/10.1016/j.foreco.2004.07.018","role":"primary","citation":"Wermelinger, B. (2004). Ecology and management of the spruce bark beetle Ips typographus - a review of recent research. Forest Ecology and Management 202(1-3):67-82. https://doi.org/10.1016/j.foreco.2004.07.018"},{"title":"Swiss, French and German post-Lothar bark-beetle studies","url":"","role":"secondary","citation":"Swiss, French and German post-Lothar studies summarized in Wermelinger (2004) and later reviews."}]},{"id":"BB0046","name":"Austria bark-beetle disturbance regime (1992-2012)","country":"Austria","y0":1992,"y1":2012,"lat":47.3,"lon":13.5,"sev":"moderate","conf":"high","note":"Direct Austrian analysis of the 1992-2012 bark-beetle disturbance regime.","sources":[{"title":"Drivers of the bark beetle disturbance regime in Alpine forests in Austria","url":"https://doi.org/10.1016/j.foreco.2014.01.044","role":"primary","citation":"Pasztor, F., Matulla, C., Rammer, W. & Lexer, M.J. (2014). Drivers of the bark beetle disturbance regime in Alpine forests in Austria. Forest Ecology and Management 318:349-358. https://doi.org/10.1016/j.foreco.2014.01.044"}]},{"id":"BB0050","name":"Waldviertel bark-beetle outbreak (2015–2020)","country":"Austria","y0":2015,"y1":2020,"lat":48.7,"lon":15.3,"sev":"high","conf":"high","note":"Hot dry summers; spruce outside optimum.","sources":[{"title":"Two unprecedented outbreaks of the European spruce bark beetle, Ips typographus L. (Col., Scolytinae) in Austria since 2015: different causes and different impacts on forests","url":"https://doi.org/10.2478/forj-2024-0014","role":"primary","citation":"Hallas, T., Steyrer, G., Laaha, G. & Hoch, G. (2024). Two unprecedented outbreaks of the European spruce bark beetle, Ips typographus L. (Col., Scolytinae) in Austria since 2015: different causes and different impacts on forests. Central European Forestry Journal 70:263-274. https://doi.org/10.2478/forj-2024-0014"}]},{"id":"BB0051","name":"Muhlviertel bark-beetle outbreak (2015–2020)","country":"Austria","y0":2015,"y1":2020,"lat":48.5,"lon":14.4,"sev":"high","conf":"high","note":"Hot dry summers; spruce outside optimum.","sources":[{"title":"Two unprecedented outbreaks of the European spruce bark beetle, Ips typographus L. (Col., Scolytinae) in Austria since 2015: different causes and different impacts on forests","url":"https://doi.org/10.2478/forj-2024-0014","role":"primary","citation":"Hallas, T., Steyrer, G., Laaha, G. & Hoch, G. (2024). Two unprecedented outbreaks of the European spruce bark beetle, Ips typographus L. (Col., Scolytinae) in Austria since 2015: different causes and different impacts on forests. Central European Forestry Journal 70:263-274. https://doi.org/10.2478/forj-2024-0014"}]},{"id":"BB0052","name":"Northern and eastern Austria bark-beetle outbreak (2018–2023)","country":"Austria","y0":2018,"y1":2023,"lat":48.3,"lon":15.5,"sev":"high","conf":"high","note":"2018 drought and heat.","sources":[{"title":"Living with bark beetles: impacts, outlook and management options","url":"https://efi.int/sites/default/files/files/publication-bank/2019/efi_fstp_8_2019.pdf","role":"primary","citation":"Hlásny, T., Krokene, P., Liebhold, A. et al. (2019). Living with bark beetles: impacts, outlook and management options. European Forest Institute, From Science to Policy 8."},{"title":"Two unprecedented outbreaks of the European spruce bark beetle, Ips typographus L. (Col., Scolytinae) in Austria since 2015: different causes and different impacts on forests","url":"https://doi.org/10.2478/forj-2024-0014","role":"secondary","citation":"Hallas, T., Steyrer, G., Laaha, G. & Hoch, G. (2024). Two unprecedented outbreaks of the European spruce bark beetle, Ips typographus L. (Col., Scolytinae) in Austria since 2015: different causes and different impacts on forests. Central European Forestry Journal 70:263-274. https://doi.org/10.2478/forj-2024-0014"}]},{"id":"BB0053","name":"East Tyrol bark-beetle outbreak (2021–2024)","country":"Austria","y0":2021,"y1":2024,"lat":46.85,"lon":12.6,"sev":"high","conf":"high","note":"Storm and snow-break damage, then warm years.","sources":[{"title":"Two unprecedented outbreaks of the European spruce bark beetle, Ips typographus L. (Col., Scolytinae) in Austria since 2015: different causes and different impacts on forests","url":"https://doi.org/10.2478/forj-2024-0014","role":"primary","citation":"Hallas, T., Steyrer, G., Laaha, G. & Hoch, G. (2024). Two unprecedented outbreaks of the European spruce bark beetle, Ips typographus L. (Col., Scolytinae) in Austria since 2015: different causes and different impacts on forests. Central European Forestry Journal 70:263-274. https://doi.org/10.2478/forj-2024-0014"},{"title":"Waldschäden (Waldzustand Tirol)","url":"https://tirol.gv.at/umwelt/wald/waldzustand/waldschaeden","role":"primary","citation":"Land Tirol (n.d., accessed 2026-10-27). Waldschäden. tirol.gv.at."},{"title":"Waldschutz: Im Angesicht des Buchdruckers [Forest protection: facing the spruce bark beetle]","url":"https://www.bfw.gv.at/pressemeldungen/waldschutz-im-angesicht-des-buchdruckers/","role":"secondary","citation":"Hoch, G., Hallas, T. & Schreck, M. (2024, 26 June). Waldschutz: Im Angesicht des Buchdruckers. Austrian Research Centre for Forests (BFW). https://www.bfw.gv.at/pressemeldungen/waldschutz-im-angesicht-des-buchdruckers/"}]},{"id":"BB0054","name":"Upper Carinthia bark-beetle outbreak (2021–2023)","country":"Austria","y0":2021,"y1":2023,"lat":46.85,"lon":13.2,"sev":"high","conf":"high","note":"Storm and snow-break damage.","sources":[{"title":"Two unprecedented outbreaks of the European spruce bark beetle, Ips typographus L. (Col., Scolytinae) in Austria since 2015: different causes and different impacts on forests","url":"https://doi.org/10.2478/forj-2024-0014","role":"primary","citation":"Hallas, T., Steyrer, G., Laaha, G. & Hoch, G. (2024). Two unprecedented outbreaks of the European spruce bark beetle, Ips typographus L. (Col., Scolytinae) in Austria since 2015: different causes and different impacts on forests. Central European Forestry Journal 70:263-274. https://doi.org/10.2478/forj-2024-0014"},{"title":"Waldschutz: Im Angesicht des Buchdruckers [Forest protection: facing the spruce bark beetle]","url":"https://www.bfw.gv.at/pressemeldungen/waldschutz-im-angesicht-des-buchdruckers/","role":"secondary","citation":"Hoch, G., Hallas, T. & Schreck, M. (2024, 26 June). Waldschutz: Im Angesicht des Buchdruckers. Austrian Research Centre for Forests (BFW). https://www.bfw.gv.at/pressemeldungen/waldschutz-im-angesicht-des-buchdruckers/"}]},{"id":"BB0055","name":"Belgian Ardennes spruce bark-beetle outbreak (2017-2022)","country":"Belgium","y0":2017,"y1":2022,"lat":50.2,"lon":5.8,"sev":"high","conf":"high","note":"Direct 2025 academic analysis of the 2017-2022 Belgian outbreak.","sources":[{"title":"Mitigating bark beetle damage in Norway spruce stands: Insights from Belgium's latest outbreak","url":"https://doi.org/10.1016/j.foreco.2025.123136","role":"primary","citation":"Gilles, A., Lisein, J., Latte, N. & Claessens, H. (2025). Mitigating bark beetle damage in Norway spruce stands: Insights from Belgium's latest outbreak. Forest Ecology and Management 597:123136. https://doi.org/10.1016/j.foreco.2025.123136"}]},{"id":"BB0056","name":"Grand Est and Bourgogne-Franche-Comte spruce bark-beetle crisis (2018-2021)","country":"France","y0":2018,"y1":2021,"lat":48.7,"lon":6.2,"sev":"extreme","conf":"high","note":"Official French assessment reports 19 million m3 and about 55,000 ha across the two regions by end-2021.","sources":[{"title":"Crise scolytes sur epiceas: quel est le bilan fin 2021? [Spruce bark-beetle crisis: assessment at the end of 2021]","url":"https://agriculture.gouv.fr/crise-scolytes-sur-epiceas-quel-est-le-bilan-fin-2021","role":"primary","citation":"Departement de la sante des forets (2022, 21 January). Crise scolytes sur epiceas: quel est le bilan fin 2021? French Ministry of Agriculture. https://agriculture.gouv.fr/crise-scolytes-sur-epiceas-quel-est-le-bilan-fin-2021"}]},{"id":"BB0058","name":"Bavaria bark-beetle outbreak (2018–2024)","country":"Germany","y0":2018,"y1":2024,"lat":48.9,"lon":11.4,"sev":"extreme","conf":"high","note":"Drought, heat and windthrow.","sources":[{"title":"Living with bark beetles: impacts, outlook and management options","url":"https://efi.int/sites/default/files/files/publication-bank/2019/efi_fstp_8_2019.pdf","role":"primary","citation":"Hlásny, T., Krokene, P., Liebhold, A. et al. (2019). Living with bark beetles: impacts, outlook and management options. European Forest Institute, From Science to Policy 8."},{"title":"Schriftliche Anfrage, Drucksache 19/8060 (Antwort des StMELF)","url":"https://www.bayern.landtag.de/www/ElanTextAblage_WP19/Drucksachen/Schriftliche%20Anfragen/19_0008060.pdf","role":"primary","citation":"Bayerischer Landtag (2025, 3 October). Drucksache 19/8060, Schriftliche Anfrage mit Antwort des StMELF."},{"title":"German national and regional forest-health reports and Ips typographus climate studies","url":"","role":"secondary","citation":"German national and regional forest-health reports and recent Ips typographus climate studies."}]},{"id":"BB0059","name":"Saxony and Thuringia bark-beetle outbreak (2018-2023)","country":"Germany","y0":2018,"y1":2023,"lat":50.9,"lon":12.5,"sev":"extreme","conf":"high","note":"Recurrent drought since 2018 made bark beetle the largest damage factor; over 19 million solid m3 of beetle-damaged wood since 2018, 5 million m3 in 2023 through September (official report).","sources":[{"title":"Living with bark beetles: impacts, outlook and management options","url":"https://efi.int/sites/default/files/files/publication-bank/2019/efi_fstp_8_2019.pdf","role":"primary","citation":"Hlásny, T., Krokene, P., Liebhold, A. et al. (2019). Living with bark beetles: impacts, outlook and management options. European Forest Institute, From Science to Policy 8."},{"title":"Experimental evidence of reduced Ips typographus damage in mixed spruce plantations","url":"https://link.springer.com/article/10.1007/s10340-025-01995-y","role":"primary","citation":"Scherer-Lorenzen, M. et al. (2025). Experimental evidence of reduced Ips typographus damage in mixed spruce plantations. Journal of Pest Science, 99, 32."},{"title":"Waldzustandsbericht Thüringen 2023 (Thuringia Forest Condition Report 2023)","url":"https://digitales-infrastruktur.thueringen.de/fileadmin/Forst_und_Jagd_Fischerei/Forstwirtschaft/2023_Waldzustandsbericht.pdf","role":"primary","citation":"Thüringer Ministerium für Infrastruktur und Landwirtschaft (2023). Waldzustandsbericht Thüringen 2023."},{"title":"German national and regional forest-health reports and Ips typographus climate studies","url":"","role":"secondary","citation":"German national and regional forest-health reports and recent Ips typographus climate studies."}]},{"id":"BB0060","name":"Harz Mountains bark-beetle outbreak (2018–2024)","country":"Germany","y0":2018,"y1":2024,"lat":51.75,"lon":10.7,"sev":"extreme","conf":"high","note":"Drought-stressed spruce and population buildup.","sources":[{"title":"Es entsteht ein klimaresistenter Wald: der niedersaechsische Umweltminister Meyer besuchte den Nationalpark Harz [A climate-resilient forest is emerging]","url":"https://www.nationalpark-harz.de/de/aktuelles/2024/2024_07_12_Umweltminister_Meyer_besucht_NLP.php","role":"primary","citation":"Nationalpark Harz (2024, 12 July). Es entsteht ein klimaresistenter Wald: der niedersaechsische Umweltminister Meyer besuchte den Nationalpark Harz. https://www.nationalpark-harz.de/de/aktuelles/2024/2024_07_12_Umweltminister_Meyer_besucht_NLP.php"}]},{"id":"BB0061","name":"Sauerland bark-beetle outbreak (2018-2023)","country":"Germany","y0":2018,"y1":2023,"lat":51.1,"lon":8.5,"sev":"high","conf":"high","note":"Reframed to Sauerland/NRW because the former source did not document Hesse.","sources":[{"title":"Forest Condition Report North Rhine-Westphalia 2024","url":"https://www.wald-und-holz.nrw.de/fileadmin/Wald_in_NRW/250102_waldzustandsbericht_nrw_2024_lang.pdf","role":"primary","citation":"Wald und Holz Nordrhein-Westfalen (2024). Forest Condition Report North Rhine-Westphalia 2024."},{"title":"Waldzustandsbericht 2022 (Vorlage 18/538)","url":"https://www.landtag.nrw.de/portal/WWW/dokumentenarchiv/Dokument/MMV18-538.pdf","role":"primary","citation":"Ministerium für Landwirtschaft und Verbraucherschutz NRW (2022, December). Waldzustandsbericht 2022. Landtag NRW, Vorlage 18/538."}]},{"id":"BB0062","name":"Black Forest and Swabian uplands bark-beetle outbreak (2018–2024)","country":"Germany","y0":2018,"y1":2024,"lat":48.2,"lon":8.2,"sev":"moderate","conf":"medium","note":"Drought and heat.","sources":[{"title":"Forest-condition monitoring results for Baden-Wuerttemberg","url":"https://www.fva-bw.de/daten-tools/monitoring/die-waldzustandserhebung/ergebnisse","role":"primary","citation":"Forest Research Institute Baden-Wuerttemberg. Forest-condition monitoring results."},{"title":"Waldschutzsituation 2024/2025 in Baden-Württemberg","url":"https://www.waldwissen.net/de/waldwirtschaft/schadensmanagement/waldschutzsituation-bw-24-25","role":"primary","citation":"Delb H et al. (2025). Waldschutzsituation 2024/2025 in Baden-Württemberg. AFZ-DerWald 09/2025: 16-20."}]},{"id":"BB0063","name":"Slovenian spruce regions bark-beetle outbreak (2003–2009)","country":"Slovenia","y0":2003,"y1":2009,"lat":46.2,"lon":14.6,"sev":"high","conf":"medium","note":"Drought and recurrent sanitary-felling peaks.","sources":[{"title":"Sanitary felling of Norway spruce due to spruce bark beetles in Slovenia: a model and projections for various climate change scenarios","url":"https://doi.org/10.1016/j.ecolmodel.2009.05.015","role":"primary","citation":"Ogris, N. & Jurc, M. (2010). Sanitary felling of Norway spruce due to spruce bark beetles in Slovenia: a model and projections for various climate change scenarios. Ecological Modelling 221:290-302. https://doi.org/10.1016/j.ecolmodel.2009.05.015"},{"title":"Slovenia Forest Service (ZGS) bark-beetle monitoring reports","url":"https://www.zdravgozd.si/","role":"secondary","citation":"Zavod za gozdove Slovenije. Bark-beetle (podlubniki) monitoring and prognosis reports. https://www.zdravgozd.si/"}]},{"id":"BB0064","name":"Western and central Slovenia ice-storm zone bark-beetle outbreak (2014–2018)","country":"Slovenia","y0":2014,"y1":2018,"lat":46.0,"lon":14.3,"sev":"high","conf":"high","note":"2014 ice storm and damaged spruce.","sources":[{"title":"Slovenia Forest Service (ZGS) bark-beetle monitoring reports","url":"https://www.zdravgozd.si/","role":"primary","citation":"Zavod za gozdove Slovenije. Bark-beetle (podlubniki) monitoring and prognosis reports. https://www.zdravgozd.si/"},{"title":"Studies of the 2014 Slovenian ice storm and subsequent bark-beetle outbreak","url":"","role":"secondary","citation":"Studies of the 2014 Slovenian ice storm and subsequent bark-beetle outbreak."}]},{"id":"BB0065","name":"Eastern Slovenia bark-beetle outbreak (2020–2023)","country":"Slovenia","y0":2020,"y1":2023,"lat":46.4,"lon":15.5,"sev":"minor","conf":"medium","note":"High trap catches and localized sanitary felling.","sources":[{"title":"Slovenia Forest Service (ZGS) bark-beetle monitoring reports","url":"https://www.zdravgozd.si/","role":"primary","citation":"Zavod za gozdove Slovenije. Bark-beetle (podlubniki) monitoring and prognosis reports. https://www.zdravgozd.si/"}]},{"id":"BB0066","name":"Northeastern Romania Ips duplicatus outbreak (2005-2014)","country":"Romania","y0":2005,"y1":2014,"lat":47.65,"lon":26.25,"sev":"high","conf":"high","note":"Direct study reports a 2005-2014 outbreak with more than 1.3 million trees killed.","sources":[{"title":"Geographical Distribution of Three Forest Invasive Beetle Species in Romania","url":"https://www.mdpi.com/2075-4450/13/7/621","role":"primary","citation":"Olenici, N. et al. (2022). Geographical Distribution of Three Forest Invasive Beetle Species in Romania. Insects 13:621. https://doi.org/10.3390/insects13070621"}]},{"id":"BB0071","name":"South Tyrol bark-beetle outbreak (2020–2025)","country":"Italy","y0":2020,"y1":2025,"lat":46.65,"lon":11.35,"sev":"high","conf":"high","note":"Storm Vaia, drought and heat.","sources":[{"title":"Spatial and temporal dynamics of a bark beetle outbreak in the Eastern European Alps","url":"https://doi.org/10.1002/ecs2.70603","role":"primary","citation":"Candotti, A., Carvalhais, N., Wang, S. & Tomelleri, E. (2026). Spatial and temporal dynamics of a bark beetle outbreak in the Eastern European Alps. Ecosphere 17(4):e70603. https://doi.org/10.1002/ecs2.70603"}]},{"id":"BB0072","name":"Trentino bark-beetle outbreak (2020–2025)","country":"Italy","y0":2020,"y1":2025,"lat":46.15,"lon":11.0,"sev":"high","conf":"high","note":"Storm Vaia and warm dry summers.","sources":[{"title":"Bostrico: forests and monitoring in Trentino","url":"https://forestefauna.provincia.tn.it/Foreste/Foreste-in-Trentino/Bostrico","role":"primary","citation":"Autonomous Province of Trento. Bostrico: forests and monitoring in Trentino."}]},{"id":"BB0073","name":"Northern Veneto Alpine zone bark-beetle outbreak (2020–2025)","country":"Italy","y0":2020,"y1":2025,"lat":46.4,"lon":12.1,"sev":"high","conf":"medium","note":"Storm Vaia and drought.","sources":[{"title":"Bostrico dell'abete rosso","url":"https://www.regione.veneto.it/web/agricoltura-e-foreste/bostrico","role":"primary","citation":"Regione del Veneto. Bostrico dell'abete rosso."}]},{"id":"BB0074","name":"Friuli-Venezia Giulia Alpine zone bark-beetle outbreak (2020–2025)","country":"Italy","y0":2020,"y1":2025,"lat":46.4,"lon":13.0,"sev":"moderate","conf":"medium","note":"Vaia-related damage and drought.","sources":[{"title":"Forest resources: compensation criteria established for bark-beetle damage","url":"https://www.regione.fvg.it/rafvg/giunta/dettaglio.act%3Bjsessionid%3DCFDE65601E6FB2BB275D0D36314F55D2?WT.ti=Ricerca+comunicati+stampa&ass=C10&dir=%2Frafvg%2Fcms%2FRAFVG%2FGiunta%2FZannier%2Fcomunicati%2F&id=128618","role":"primary","citation":"Regione Autonoma Friuli Venezia Giulia (2023, 29 December). Forest resources: compensation criteria established for bark-beetle damage."}]},{"id":"BB0075","name":"Moscow Oblast bark-beetle outbreak (2010–2014)","country":"Russia","y0":2010,"y1":2014,"lat":55.5,"lon":37.8,"sev":"extreme","conf":"high","note":"2010 heatwave and drought.","sources":[{"title":"Ips typographus (IPSXTY) datasheet","url":"https://gd.eppo.int/taxon/IPSXTY/datasheet","role":"primary","citation":"EPPO Global Database. Ips typographus (IPSXTY) datasheet. https://gd.eppo.int/taxon/IPSXTY/datasheet (last updated 2020-09-01)."},{"title":"Russian forest-health records and European syntheses of the Moscow-region outbreak after the 2010 heatwave","url":"","role":"secondary","citation":"Russian forest-health records and European syntheses of the Moscow-region outbreak after the 2010 heatwave."}]},{"id":"BB0077","name":"Dvinsko-Pinegskiy / Arkhangelsk bark-beetle outbreak (2001-2014)","country":"Russia","y0":2001,"y1":2014,"lat":64.5,"lon":40.5,"sev":"moderate","conf":"high","note":"High-latitude outbreak documented for 2001-2014 on Picea obovata.","sources":[{"title":"Northernmost European spruce bark beetle Ips typographus outbreak: modelling tree mortality using remote sensing and climate data","url":"https://doi.org/10.1016/j.foreco.2021.119829","role":"primary","citation":"Trubin, A., Mezei, P., Zabihi, K., Surový, P. & Jakuš, R. (2022). Northernmost European spruce bark beetle Ips typographus outbreak: modelling tree mortality using remote sensing and climate data. Forest Ecology and Management 505:119829. https://doi.org/10.1016/j.foreco.2021.119829"}]},{"id":"BB0078","name":"Leningrad Oblast bark-beetle outbreak (2010-2016)","country":"Russia","y0":2010,"y1":2016,"lat":59.5,"lon":33.0,"sev":"moderate","conf":"high","note":"Dedicated regional chronology identifies a mass outbreak during 2010-2016.","sources":[{"title":"Historical outbreaks of Ips typographus in Leningrad Oblast","url":"https://sciencejournals.ru/view-article/?a=LesVed2303009Selikhovkin&j=lesved&n=3&v=0&y=2023","role":"primary","citation":"Selikhovkin, A. V. et al. (2023). Historical outbreaks of Ips typographus in Leningrad Oblast. Lesovedenie 3."}]},{"id":"BB0079","name":"Kent initial breeding site bark-beetle outbreak (2018–2018)","country":"England, UK","y0":2018,"y1":2018,"lat":51.28,"lon":1.08,"sev":"minor","conf":"high","note":"Likely wind-assisted or imported arrival into stressed spruce.","sources":[{"title":"Larger eight-toothed European spruce bark beetle (Ips typographus) pest resource","url":"https://www.forestresearch.gov.uk/tools-and-resources/fthr/pest-and-disease-resources/larger-eight-toothed-european-spruce-bark-beetle-ips-typographus/","role":"primary","citation":"Forest Research (UK). Larger eight-toothed European spruce bark beetle (Ips typographus). https://www.forestresearch.gov.uk/tools-and-resources/fthr/pest-and-disease-resources/larger-eight-toothed-european-spruce-bark-beetle-ips-typographus/"},{"title":"Eight-toothed spruce bark beetle (Ips typographus) guidance","url":"https://www.gov.uk/guidance/eight-toothed-european-spruce-bark-beetle-ips-typographus","role":"secondary","citation":"UK Government / Forestry Commission. Eight-toothed spruce bark beetle (Ips typographus) guidance. https://www.gov.uk/guidance/eight-toothed-european-spruce-bark-beetle-ips-typographus"}]},{"id":"BB0080","name":"Kent outbreak sites bark-beetle outbreak (2021–2021)","country":"England, UK","y0":2021,"y1":2021,"lat":51.2,"lon":0.9,"sev":"minor","conf":"high","note":"Cross-Channel dispersal during large flight event.","sources":[{"title":"Larger eight-toothed European spruce bark beetle (Ips typographus) pest resource","url":"https://www.forestresearch.gov.uk/tools-and-resources/fthr/pest-and-disease-resources/larger-eight-toothed-european-spruce-bark-beetle-ips-typographus/","role":"primary","citation":"Forest Research (UK). Larger eight-toothed European spruce bark beetle (Ips typographus). https://www.forestresearch.gov.uk/tools-and-resources/fthr/pest-and-disease-resources/larger-eight-toothed-european-spruce-bark-beetle-ips-typographus/"},{"title":"Eight-toothed spruce bark beetle (Ips typographus) guidance","url":"https://www.gov.uk/guidance/eight-toothed-european-spruce-bark-beetle-ips-typographus","role":"secondary","citation":"UK Government / Forestry Commission. Eight-toothed spruce bark beetle (Ips typographus) guidance. https://www.gov.uk/guidance/eight-toothed-european-spruce-bark-beetle-ips-typographus"},{"title":"New action to protect against impacts of bark beetle tree pest / Ips typographus grid summaries","url":"https://www.gov.uk/government/news/new-action-to-protect-against-impacts-of-bark-beetle-tree-pest","role":"secondary","citation":"Forestry Commission (2024). New action to protect against impacts of bark beetle tree pest; periodic Ips typographus grid summaries, e.g. https://assets.publishing.service.gov.uk/media/66e9490ac069f68b7681bc08/20240925_Ips_typographus_25km_grid_summary.pdf"}]},{"id":"BB0081","name":"East Sussex outbreak sites bark-beetle outbreak (2021–2021)","country":"England, UK","y0":2021,"y1":2021,"lat":50.9,"lon":0.3,"sev":"minor","conf":"high","note":"Cross-Channel dispersal.","sources":[{"title":"Larger eight-toothed European spruce bark beetle (Ips typographus) pest resource","url":"https://www.forestresearch.gov.uk/tools-and-resources/fthr/pest-and-disease-resources/larger-eight-toothed-european-spruce-bark-beetle-ips-typographus/","role":"primary","citation":"Forest Research (UK). Larger eight-toothed European spruce bark beetle (Ips typographus). https://www.forestresearch.gov.uk/tools-and-resources/fthr/pest-and-disease-resources/larger-eight-toothed-european-spruce-bark-beetle-ips-typographus/"},{"title":"Eight-toothed spruce bark beetle (Ips typographus) guidance","url":"https://www.gov.uk/guidance/eight-toothed-european-spruce-bark-beetle-ips-typographus","role":"secondary","citation":"UK Government / Forestry Commission. Eight-toothed spruce bark beetle (Ips typographus) guidance. https://www.gov.uk/guidance/eight-toothed-european-spruce-bark-beetle-ips-typographus"},{"title":"New action to protect against impacts of bark beetle tree pest / Ips typographus grid summaries","url":"https://www.gov.uk/government/news/new-action-to-protect-against-impacts-of-bark-beetle-tree-pest","role":"secondary","citation":"Forestry Commission (2024). New action to protect against impacts of bark beetle tree pest; periodic Ips typographus grid summaries, e.g. https://assets.publishing.service.gov.uk/media/66e9490ac069f68b7681bc08/20240925_Ips_typographus_25km_grid_summary.pdf"}]},{"id":"BB0082","name":"Southeast England confirmed sites bark-beetle outbreak (2022–2022)","country":"England, UK","y0":2022,"y1":2022,"lat":51.1,"lon":0.5,"sev":"minor","conf":"high","note":"Repeated continental arrivals.","sources":[{"title":"Larger eight-toothed European spruce bark beetle (Ips typographus) pest resource","url":"https://www.forestresearch.gov.uk/tools-and-resources/fthr/pest-and-disease-resources/larger-eight-toothed-european-spruce-bark-beetle-ips-typographus/","role":"primary","citation":"Forest Research (UK). Larger eight-toothed European spruce bark beetle (Ips typographus). https://www.forestresearch.gov.uk/tools-and-resources/fthr/pest-and-disease-resources/larger-eight-toothed-european-spruce-bark-beetle-ips-typographus/"},{"title":"Eight-toothed spruce bark beetle (Ips typographus) guidance","url":"https://www.gov.uk/guidance/eight-toothed-european-spruce-bark-beetle-ips-typographus","role":"secondary","citation":"UK Government / Forestry Commission. Eight-toothed spruce bark beetle (Ips typographus) guidance. https://www.gov.uk/guidance/eight-toothed-european-spruce-bark-beetle-ips-typographus"}]},{"id":"BB0083","name":"Southeast England confirmed sites bark-beetle outbreak (2023–2023)","country":"England, UK","y0":2023,"y1":2023,"lat":51.1,"lon":0.5,"sev":"minor","conf":"high","note":"Repeated continental arrivals.","sources":[{"title":"Larger eight-toothed European spruce bark beetle (Ips typographus) pest resource","url":"https://www.forestresearch.gov.uk/tools-and-resources/fthr/pest-and-disease-resources/larger-eight-toothed-european-spruce-bark-beetle-ips-typographus/","role":"primary","citation":"Forest Research (UK). Larger eight-toothed European spruce bark beetle (Ips typographus). https://www.forestresearch.gov.uk/tools-and-resources/fthr/pest-and-disease-resources/larger-eight-toothed-european-spruce-bark-beetle-ips-typographus/"},{"title":"Eight-toothed spruce bark beetle (Ips typographus) guidance","url":"https://www.gov.uk/guidance/eight-toothed-european-spruce-bark-beetle-ips-typographus","role":"secondary","citation":"UK Government / Forestry Commission. Eight-toothed spruce bark beetle (Ips typographus) guidance. https://www.gov.uk/guidance/eight-toothed-european-spruce-bark-beetle-ips-typographus"}]},{"id":"BB0084","name":"Southeast England confirmed sites bark-beetle outbreak (2024–2024)","country":"England, UK","y0":2024,"y1":2024,"lat":51.1,"lon":0.5,"sev":"minor","conf":"high","note":"Repeated continental arrivals.","sources":[{"title":"Eight-toothed spruce bark beetle (Ips typographus) guidance","url":"https://www.gov.uk/guidance/eight-toothed-european-spruce-bark-beetle-ips-typographus","role":"primary","citation":"UK Government / Forestry Commission. Eight-toothed spruce bark beetle (Ips typographus) guidance. https://www.gov.uk/guidance/eight-toothed-european-spruce-bark-beetle-ips-typographus"},{"title":"New action to protect against impacts of bark beetle tree pest / Ips typographus grid summaries","url":"https://www.gov.uk/government/news/new-action-to-protect-against-impacts-of-bark-beetle-tree-pest","role":"secondary","citation":"Forestry Commission (2024). New action to protect against impacts of bark beetle tree pest; periodic Ips typographus grid summaries, e.g. https://assets.publishing.service.gov.uk/media/66e9490ac069f68b7681bc08/20240925_Ips_typographus_25km_grid_summary.pdf"}]},{"id":"BB0085","name":"East Anglia localized establishments bark-beetle outbreak (2024–2024)","country":"England, UK","y0":2024,"y1":2024,"lat":52.3,"lon":1.0,"sev":"minor","conf":"high","note":"Wind-assisted dispersal from continental Europe.","sources":[{"title":"New action to protect against impacts of bark beetle tree pest / Ips typographus grid summaries","url":"https://www.gov.uk/government/news/new-action-to-protect-against-impacts-of-bark-beetle-tree-pest","role":"primary","citation":"Forestry Commission (2024). New action to protect against impacts of bark beetle tree pest; periodic Ips typographus grid summaries, e.g. https://assets.publishing.service.gov.uk/media/66e9490ac069f68b7681bc08/20240925_Ips_typographus_25km_grid_summary.pdf"}]},{"id":"BB0086","name":"West Sussex Sitka spruce finding bark-beetle outbreak (2024–2024)","country":"England, UK","y0":2024,"y1":2024,"lat":50.95,"lon":-0.5,"sev":"minor","conf":"high","note":"Breeding galleries on cut/fallen Sitka spruce near infested Norway spruce.","sources":[{"title":"New action to protect against impacts of bark beetle tree pest / Ips typographus grid summaries","url":"https://www.gov.uk/government/news/new-action-to-protect-against-impacts-of-bark-beetle-tree-pest","role":"primary","citation":"Forestry Commission (2024). New action to protect against impacts of bark beetle tree pest; periodic Ips typographus grid summaries, e.g. https://assets.publishing.service.gov.uk/media/66e9490ac069f68b7681bc08/20240925_Ips_typographus_25km_grid_summary.pdf"}]},{"id":"BB0087","name":"Southeast England confirmed outbreak grid cells bark-beetle outbreak (2025–2025)","country":"England, UK","y0":2025,"y1":2025,"lat":51.1,"lon":0.5,"sev":"minor","conf":"high","note":"Continued incursions and surveillance detections.","sources":[{"title":"Ips typographus confirmed outbreak-site grid summary (10km grid, SE England)","url":"https://assets.publishing.service.gov.uk/media/6903b862fabc9f10a832a8b8/20251022_Ips_typographus_10km_grid_summary_SE.pdf","role":"primary","citation":"UK Government. Ips typographus confirmed outbreak-site grid summaries, data through November 2025. https://assets.publishing.service.gov.uk/media/6903b862fabc9f10a832a8b8/20251022_Ips_typographus_10km_grid_summary_SE.pdf"}]},{"id":"BB0088","name":"Karelian Isthmus bark-beetle outbreak (2021–2023)","country":"Russia","y0":2021,"y1":2023,"lat":60.5,"lon":29.7,"sev":"moderate","conf":"medium","note":"Elevated June-July 2021 temperatures enabled two beetle generations from a large pre-existing pest stock.","sources":[{"title":"A New Outbreak of the European Spruce Bark Beetle, Ips typographus (L.) (Coleoptera, Curculionidae) in Leningrad Province","url":"https://doi.org/10.1134/S0013873822030034","role":"primary","citation":"Selikhovkin, A.V. et al. (2022). A New Outbreak of the European Spruce Bark Beetle, Ips typographus (L.) (Coleoptera, Curculionidae) in Leningrad Province. Entomological Review 102(3):1-11. https://doi.org/10.1134/S0013873822030034"}]},{"id":"BB0089","name":"Rold Skov bark-beetle outbreak (1982–1985)","country":"Denmark","y0":1982,"y1":1985,"lat":56.83,"lon":9.98,"sev":"high","conf":"high","note":"24 November 1981 severe windthrow gale, one of the most severe in Danish forestry history.","sources":[{"title":"The spread of Ips typographus (L.) (Coleoptera, Scolytidae) attacks following heavy windthrow in Denmark, analysed using GIS","url":"https://doi.org/10.1016/S0378-1127(00)00477-1","role":"primary","citation":"Wichmann, L. & Ravn, H.P. (2001). The spread of Ips typographus (L.) (Coleoptera, Scolytidae) attacks following heavy windthrow in Denmark, analysed using GIS. Forest Ecology and Management 148:31-39. https://doi.org/10.1016/S0378-1127(00)00477-1"}]},{"id":"BB0090","name":"Bistrishko Branishte reserve bark-beetle outbreak (2001–2006)","country":"Bulgaria","y0":2001,"y1":2006,"lat":42.56,"lon":23.28,"sev":"moderate","conf":"medium","note":"Post-windthrow bark-beetle population buildup in a biosphere reserve.","sources":[{"title":"Outbreak of Ips typographus (L.) (Coleoptera: Scolytidae) and drying of Norway spruce (Picea abies L. Karst.) on Vitosha Mountain","url":"https://www.researchgate.net/publication/259645173_Outbreak_of_Ips_typographus_L_Coleoptera_Scolitidae_and_drying_of_norway_spruse_on_Vitosha_mountain","role":"primary","citation":"Georgiev, G., Mirchev, P., Tsankov, G., Rosnev, B. & Petkov, P. (2006). Outbreak of Ips typographus (L.) (Coleoptera: Scolytidae) and drying of Norway spruce (Picea abies L. Karst.) on Vitosha Mountain. Proceedings of FORMEC 2006, 24-28 September 2006, Sofia, Bulgaria, Expressprint Ltd., 218-220."}]},{"id":"BB0091","name":"Chuprene Reserve bark-beetle outbreak (2017–2021)","country":"Bulgaria","y0":2017,"y1":2021,"lat":43.63,"lon":22.8,"sev":"minor","conf":"medium","note":"Remote-sensing-documented expansion of attacks from a mid-1980s establishment.","sources":[{"title":"Remote Sensing Assessment of the Expansion of Ips typographus Attacks in the Chuprene Reserve, Western Balkan Range","url":"https://doi.org/10.3390/f13010039","role":"primary","citation":"Georgiev, G. et al. (2022). Remote Sensing Assessment of the Expansion of Ips typographus Attacks in the Chuprene Reserve, Western Balkan Range. Forests 13(1):39. https://doi.org/10.3390/f13010039"}]},{"id":"BB0093","name":"Vidzeme and Selija bark-beetle outbreak (2022–2023)","country":"Latvia","y0":2022,"y1":2023,"lat":56.8,"lon":25.7,"sev":"high","conf":"high","note":"Damage rose from a historical 550-850 ha/yr to ~2,300 ha in 2022; national state of emergency declared 1 Apr-30 Jun 2023.","sources":[{"title":"Eight-toothed spruce bark beetle emergency in Latvia","url":"https://eng.lsm.lv/article/society/environment/28.03.2023-eight-toothed-spruce-bark-beetle-emergency-to-be-announced-in-latvia.a502715/","role":"primary","citation":"LSM.lv (2023). Eight-toothed spruce bark beetle emergency in Latvia. https://eng.lsm.lv/article/society/environment/28.03.2023-eight-toothed-spruce-bark-beetle-emergency-to-be-announced-in-latvia.a502715/"},{"title":"2023 forest-pest and disease monitoring report (MBRF monitoring results)","url":"https://silava.lv/images/Petijumi/Nacionalais-meza-monitorings/MBRF-monitorings-rezultati/2023-MBRF-monitorings-Kaitekli-slimibas-Parskats.pdf","role":"primary","citation":"Silava (2023). 2023 forest-pest and disease monitoring report."}]},{"id":"BB0094","name":"Kopaonik National Park spruce decline and bark-beetle episode (2014-2017)","country":"Serbia","y0":2014,"y1":2017,"lat":43.28,"lon":20.82,"sev":"moderate","conf":"high","note":"Pronounced decline in 2014-2017 following drought, heat and earlier disturbance.","sources":[{"title":"Unfavorable climatic factors and their impact on the decline of spruce at the Kopaonik National Park (central Serbia)","url":"https://www.researchgate.net/publication/360963654_Unfavorable_climatic_factors_and_their_impact_on_the_decline_of_spruce_at_the_Kopaonik_national_park_Central_Serbia","role":"primary","citation":"Cesljar, G. et al. (2022). Unfavorable climatic factors and their impact on the decline of spruce at the Kopaonik National Park (central Serbia). Fresenius Environmental Bulletin 31:5238-5249."}]},{"id":"BB0095","name":"Golija Nature Park spruce decline and bark-beetle episode (2012-2016)","country":"Serbia","y0":2012,"y1":2016,"lat":43.45,"lon":20.28,"sev":"moderate","conf":"high","note":"Mortality peaked in 2015 and bark-beetle abundance in 2016.","sources":[{"title":"The Correlation between the Changes in Climate, the Intensity of Spruce Decline and the Abundance of Spruce Bark Beetles in 'Golija' Nature Park","url":"https://doi.org/10.5937/sustfor1673059t","role":"primary","citation":"Tabakovic-Tosic, M. & Milosavljevic, M. (2016). The Correlation between the Changes in Climate, the Intensity of Spruce Decline and the Abundance of Spruce Bark Beetles in 'Golija' Nature Park. Sustainable Forestry Collection 73-74:59-68. https://doi.org/10.5937/sustfor1673059t"}]},{"id":"BB0096","name":"Tara National Park bark-beetle phenology monitoring (2016-2018)","country":"Serbia","y0":2016,"y1":2018,"lat":43.9,"lon":19.45,"sev":"minor","conf":"high","note":"Three-year phenology and trap-monitoring study; not a mortality-footprint survey.","sources":[{"title":"Phenological Study of the Spruce Bark Beetles Ips typographus (L.) and Pityogenes chalcographus (L.) (Coleoptera: Curculionidae) on Tara Mountain, Serbia","url":"https://doi.org/10.15177/seefor.25-16","role":"primary","citation":"Milosavljevic, M. et al. (2025). Phenological Study of the Spruce Bark Beetles Ips typographus (L.) and Pityogenes chalcographus (L.) (Coleoptera: Curculionidae) on Tara Mountain, Serbia. South-east European Forestry 16(2):123-130. https://doi.org/10.15177/seefor.25-16"}]},{"id":"BB0097","name":"Teberda Nature Reserve Picea orientalis mass drying-out (2012-2016)","country":"Russia","y0":2012,"y1":2016,"lat":43.45,"lon":41.75,"sev":"high","conf":"high","note":"Mass Ips typographus reproduction in the hot seasons of 2012 and 2015.","sources":[{"title":"Reconstruction of the Dynamics of Dark Coniferous Forests of Teberdinskii Nature Reserve and Prospects for Their Natural Recovery after Mass Drying-Out","url":"https://doi.org/10.1134/S1062359023100126","role":"primary","citation":"Pukinskaya, M.Y. (2023). Reconstruction of the Dynamics of Dark Coniferous Forests of Teberdinskii Nature Reserve and Prospects for Their Natural Recovery after Mass Drying-Out. Biology Bulletin 50:2647-2658. https://doi.org/10.1134/S1062359023100126"}]},{"id":"BB0098","name":"Ips typographus detection in Valle de Roncal (2022-2023)","country":"Spain","y0":2022,"y1":2023,"lat":42.95,"lon":-1.1,"sev":"minor","conf":"high","note":"Range-edge detection in a planted Norway-spruce stand; not a landscape outbreak.","sources":[{"title":"Deteccion de Ips typographus (Linne, 1758) (Coleoptera: Curculionidae: Scolytinae) en una plantacion forestal de Picea abies en el Valle de Roncal (Navarra, Espana)","url":"https://www.researchgate.net/publication/371133437_Deteccion_de_Ips_typographus_Linne_1758_Coleoptera_Curculionidae_Scolytinae_en_una_plantacion_forestal_de_Picea_abies_en_el_Valle_de_Roncal_Navarra_Espana_Detection_of_Ips_typographus_Linne_1758_Coleo","role":"primary","citation":"Riba-Flinch, J.M. (2023). Deteccion de Ips typographus (Linne, 1758) en una plantacion forestal de Picea abies en el Valle de Roncal (Navarra, Espana). Revista Gaditana de Entomologia XIV:43-53."}]},{"id":"BB0099","name":"Ips typographus confirmed presence, Valle de Aran plantations (2020)","country":"Spain","y0":2020,"y1":2020,"lat":42.7,"lon":0.85,"sev":"minor","conf":"high","note":"Confirmed presence in Picea abies plantations, specimens collected July 2020.","sources":[{"title":"Deteccion de Dendroctonus micans (Kugelann, 1794) (Coleoptera: Curculionidae: Scolytinae) y confirmacion de la presencia de Ips typographus (Linne, 1758) en plantaciones forestales de Picea abies en el Valle de Aran (Lleida, Espana)","url":"https://dialnet.unirioja.es/servlet/articulo?codigo=7997324","role":"primary","citation":"Riba Flinch, J.M. (2021). Deteccion de Dendroctonus micans (Kugelann, 1794) y confirmacion de la presencia de Ips typographus (Linne, 1758) en plantaciones forestales de Picea abies en el Valle de Aran (Lleida, Espana). Boletin de la SEA 68."}]},{"id":"BB0101","name":"Smolyan region windthrow-associated bark-beetle infestation, Rhodope Mountains (2018-2020)","country":"Bulgaria","y0":2018,"y1":2020,"lat":41.58,"lon":24.7,"sev":"moderate","conf":"high","note":"January 2018 windthrow followed by bark-beetle attacks on healthy spruce by August 2020.","sources":[{"title":"Early detection of Ips typographus infestations by using Sentinel-2 satellite images in windthrow affected Norway spruce forests in Smolyan region, Bulgaria","url":"https://doi.org/10.3897/silvabalcanica.22.e98314","role":"primary","citation":"Georgiev, G. et al. (2022). Early detection of Ips typographus infestations by using Sentinel-2 satellite images in windthrow affected Norway spruce forests in Smolyan region, Bulgaria. Silva Balcanica 23(2):27-34. https://doi.org/10.3897/silvabalcanica.22.e98314"}]},{"id":"BB0102","name":"Gorski Kotar bark-beetle outbreak and sanitary felling (2016-2021)","country":"Croatia","y0":2016,"y1":2021,"lat":45.35,"lon":14.8,"sev":"extreme","conf":"high","note":"Ice-storm-driven population increase; natural disaster declared 2017.","sources":[{"title":"A Decade of Sanitary Fellings Followed by Climate Extremes in Croatian Managed Forests","url":"https://doi.org/10.3390/land14040766","role":"primary","citation":"Djuka, A. et al. (2025). A Decade of Sanitary Fellings Followed by Climate Extremes in Croatian Managed Forests. Land 14(4):766. https://doi.org/10.3390/land14040766"}]},{"id":"BB0103","name":"Sarajevo Canton protected-area Ips typographus monitoring (2018-2021)","country":"Bosnia and Herzegovina","y0":2018,"y1":2021,"lat":43.78,"lon":18.35,"sev":"minor","conf":"high","note":"Pheromone-trap monitoring record; not a mortality footprint.","sources":[{"title":"Influence of Forest Type and Climate Factors on the Number of Caught Ips typographus (Coleoptera, Curculionidae) Bark Beetles in Pheromone Traps in Protected Areas of Bosnia and Herzegovina","url":"https://doi.org/10.15177/seefor.23-16","role":"primary","citation":"Mujezinovic, O. et al. (2023). Influence of Forest Type and Climate Factors on the Number of Caught Ips typographus (Coleoptera, Curculionidae) Bark Beetles in Pheromone Traps in Protected Areas of Bosnia and Herzegovina. South-east European Forestry 14(2). https://doi.org/10.15177/seefor.23-16"}]},{"id":"BB0104","name":"Sumava (Bohemian Forest) pre-park storm-triggered outbreak (1982-1985)","country":"Czechia","y0":1982,"y1":1985,"lat":49.05,"lon":13.65,"sev":"high","conf":"high","note":"1982-83 predisposing drought followed by a major 1984 windstorm.","sources":[{"title":"Meteorological and climatological triggers of notable past and present bark beetle outbreaks in the Czech Republic","url":"https://doi.org/10.5194/cp-18-2155-2022","role":"primary","citation":"Brazdil, R. et al. (2022). Meteorological and climatological triggers of notable past and present bark beetle outbreaks in the Czech Republic. Climate of the Past 18:2155-2180. https://doi.org/10.5194/cp-18-2155-2022"}]},{"id":"BB0105","name":"Switzerland pre-Vivian windthrow/drought-triggered epidemic (1982-1989)","country":"Switzerland","y0":1982,"y1":1989,"lat":46.8,"lon":8.2,"sev":"high","conf":"high","note":"1982-83 windthrow and drought; epidemic peaked 1984-1985.","sources":[{"title":"Ursachen und Verlauf der Buchdrucker-Epidemien (Ips typographus L.) in der Schweiz von 1984 bis 1999 [Causes and progress of the eight-toothed spruce bark beetle epidemics (Ips typographus) in Switzerland from 1984 to 1999]","url":"https://doi.org/10.3188/szf.2003.0437","role":"primary","citation":"Meier, F., Gall, R. & Forster, B. (2003). Ursachen und Verlauf der Buchdrucker-Epidemien (Ips typographus L.) in der Schweiz von 1984 bis 1999. Schweizerische Zeitschrift fur Forstwesen 154(11):437-441. https://doi.org/10.3188/szf.2003.0437"}]},{"id":"BB0106","name":"Brasov and Sibiu counties real-data disturbance cluster (2014-2017)","country":"Romania","y0":2014,"y1":2017,"lat":45.68,"lon":24.75,"sev":"minor","conf":"medium","note":"Satellite-detected (DEFID2) disturbance polygons, concentrated in 2017.","sources":[{"title":"DEFID2: Database of European Forest Insect and Disease Disturbances","url":"https://onlinelibrary.wiley.com/doi/full/10.1111/gcb.16912","role":"primary","citation":"Forzieri, G. et al. (2023). The Database of European Forest Insect and Disease Disturbances: DEFID2. Global Change Biology 29. https://doi.org/10.1111/gcb.16912 (data hosted via JRC Forest Information System for Europe)."}]},{"id":"BB0107","name":"South Tyrol / Trentino pre-Vaia minor activity (2010-2015)","country":"Italy","y0":2010,"y1":2015,"lat":46.8,"lon":11.6,"sev":"minor","conf":"high","note":"Satellite-detected minor disturbance activity, predating storm Vaia.","sources":[{"title":"DEFID2: Database of European Forest Insect and Disease Disturbances","url":"https://onlinelibrary.wiley.com/doi/full/10.1111/gcb.16912","role":"primary","citation":"Forzieri, G. et al. (2023). The Database of European Forest Insect and Disease Disturbances: DEFID2. Global Change Biology 29. https://doi.org/10.1111/gcb.16912 (data hosted via JRC Forest Information System for Europe)."}]},{"id":"BB0108","name":"Dinaric mountain forests (Karst plateau) bark-beetle monitoring (1986-2000)","country":"Slovenia","y0":1986,"y1":2000,"lat":45.6,"lon":14.3,"sev":"moderate","conf":"high","note":"Long-term phenological monitoring and modeling in Dinaric Fir-Beech forest.","sources":[{"title":"Spruce bark beetles (Ips typographus, Pityogenes chalcographus, Col.: Scolytidae) in the Dinaric mountain forests of Slovenia: monitoring and modeling","url":"https://doi.org/10.1016/j.ecolmodel.2005.10.014","role":"primary","citation":"Jurc, M. et al. (2006). Spruce bark beetles (Ips typographus, Pityogenes chalcographus, Col.: Scolytidae) in the Dinaric mountain forests of Slovenia: monitoring and modeling. Ecological Modelling 194(1-3):219-226. https://doi.org/10.1016/j.ecolmodel.2005.10.014"}]},{"id":"BB0109","name":"Borjomi Gorge Ips typographus sampling record (2007-2008)","country":"Georgia","y0":2007,"y1":2008,"lat":41.84,"lon":43.39,"sev":"minor","conf":"high","note":"Pathogen/nematode study sampled Ips typographus at three Borjomi sites in 2007-2008.","sources":[{"title":"Occurrence of pathogens and nematodes in the spruce bark beetle, Ips typographus (Col., Scolytidae), in Borjomi Gorge","url":"https://citeseerx.ist.psu.edu/document?doi=44800e9a18213bcc952dce565760aa258f69e067&repid=rep1&type=pdf","role":"primary","citation":"Burjanadze, M. & Goginashvili, N. (2009). Occurrence of pathogens and nematodes in the spruce bark beetle, Ips typographus (Col., Scolytidae), in Borjomi Gorge. Bulletin of the Georgian National Academy of Sciences 3(1):145-150."}]},{"id":"BB0110","name":"Rheinland-Pfalz historically unprecedented bark-beetle losses (2018-2019)","country":"Germany","y0":2018,"y1":2019,"lat":49.9,"lon":7.4,"sev":"extreme","conf":"high","note":"Official state release reports about 4 million m3 in 2018-2019.","sources":[{"title":"About four million cubic metres of bark-beetle damage in Rheinland-Pfalz in 2018-2019","url":"https://mlwuf.rlp.de/service/pressemitteilungen/detail/hoefken-unser-deutscher-wald-ist-klimakrank-die-treibhausgase-muessen-dringend-reduziert-werden","role":"primary","citation":"Ministry for Environment, Rhineland-Palatinate (2020). Official forest-health press release."}]},{"id":"BB0112","name":"Southeastern Norway post-storm bark-beetle crisis (2021-2024)","country":"Norway","y0":2021,"y1":2024,"lat":60.3,"lon":9.5,"sev":"high","conf":"medium","note":"2021-2022 storms created ideal breeding conditions; industry describes a crisis by 2024.","sources":[{"title":"Skognaeringen slar alarm om invasjon av granbarkbiller etter kraftig uvaer tilbake i 2021 [The forestry industry sounds the alarm over an invasion of spruce bark beetles following severe weather back in 2021]","url":"https://www.nrk.no/innlandet/skognaeringen-slar-alarm-om-invasjon-av-granbarkbiller-etter-kraftig-uvaer-tilbake-i-2021.-1.16985875","role":"primary","citation":"NRK Innlandet (2024, 12 August). Skognaeringen slar alarm om invasjon av granbarkbiller etter kraftig uvaer tilbake i 2021. https://www.nrk.no/innlandet/skognaeringen-slar-alarm-om-invasjon-av-granbarkbiller-etter-kraftig-uvaer-tilbake-i-2021.-1.16985875"}]},{"id":"BB0113","name":"Graubunden, Glarus and Wallis snow-damage-triggered increase (2025)","country":"Switzerland","y0":2025,"y1":2025,"lat":46.7,"lon":9.3,"sev":"moderate","conf":"high","note":"Heavy wet-snow damage created breeding material, against a national downward trend.","sources":[{"title":"Buchdruckersituation entspannt sich weiter [Spruce bark beetle situation continues to ease]","url":"https://www.wsl.ch/de/news/buchdruckersituation-entspannt-sich-weiter/","role":"primary","citation":"WSL (2026, 31 March). Buchdruckersituation entspannt sich weiter. https://www.wsl.ch/de/news/buchdruckersituation-entspannt-sich-weiter/"}]},{"id":"BB0114","name":"Lviv mountain forests bark-beetle outbreak (2003–2017)","country":"Ukraine","y0":2003,"y1":2017,"lat":49.05,"lon":23.5,"sev":"high","conf":"high","note":"A documented regional outbreak expanded after 2003 and reached about 1,900 ha of active foci in 2009.","sources":[{"title":"The role of biotic factors in spruce decline in the Ukrainian Carpathians","url":"https://doi.org/10.15421/411827","role":"primary","citation":"Kramarets, V. & Matsiakh, I. (2018). The role of biotic factors in spruce decline in the Ukrainian Carpathians. Proceedings of the Forestry Academy of Sciences of Ukraine 17:121-132. https://doi.org/10.15421/411827"},{"title":"State of the bark beetle (Ips typographus L.) outbreak foci in different regions of Ukraine","url":"https://doi.org/10.33220/1026-3365.137.2020.120","role":"secondary","citation":"Kukina, O. M. & Zinchenko, O. V. (2020). State of the bark beetle (Ips typographus L.) outbreak foci in different regions of Ukraine. Forestry and Forest Melioration 137:120-126. https://doi.org/10.33220/1026-3365.137.2020.120"}]},{"id":"BB0116","name":"Estonia bark-beetle outbreak (2006-2007)","country":"Estonia","y0":2006,"y1":2007,"lat":58.7,"lon":25.5,"sev":"high","conf":"medium","note":"Listed as a major national outbreak (2006-2007) in the Estonian Environment Agency's official bark-beetle chronology.","sources":[{"title":"Kuuse-kooreüraskist Eestis (About the spruce bark beetle in Estonia)","url":"https://www.keskkonnaagentuur.ee/sites/default/files/kuuse_kooreyraskist_eestis.pdf","role":"primary","citation":"Õunap, H. Kuuse-kooreüraskist Eestis. Keskkonnaagentuur."}]},{"id":"BB0117","name":"Goetaland-Svealand spruce bark-beetle episode (1992-2000)","country":"Sweden","y0":1992,"y1":2000,"lat":58.5,"lon":14.0,"sev":"high","conf":"high","note":"Damage increased from 1992, culminated around 1997 (approx. 500,000-600,000 m3 killed that year), then generally declined by 2000.","sources":[{"title":"Skador på skog (Damage to forest), Rapport 8O/2001","url":"https://shop.skogsstyrelsen.se/shop/9098/art61/4646061-6b588c-1709.pdf","role":"primary","citation":"Samuelsson, H. & Örlander, G. (2001). Skador på skog. Skogsstyrelsen Rapport 8O/2001."}]},{"id":"BB0118","name":"Bohemian-Saxon Switzerland bark-beetle outbreak (2017-2022)","country":"Czechia","y0":2017,"y1":2022,"lat":50.9,"lon":14.25,"sev":"high","conf":"high","note":"Peer-reviewed study: bark-beetle-infested forest covered 17.6% of a 700 km2 transboundary study region.","sources":[{"title":"Effects of bark beetle disturbance and fuel types on fire radiative power and burn severity in the Bohemian-Saxon Switzerland","url":"https://academic.oup.com/forestry/article/98/1/59/7688820","role":"primary","citation":"Beetz, K. et al. (2025). Effects of bark beetle disturbance and fuel types on fire radiative power and burn severity in the Bohemian-Saxon Switzerland. Forestry 98(1):59-70. https://doi.org/10.1093/forestry/cpae024"},{"title":"Bohemian Switzerland National Park bark-beetle outbreak study","url":"https://www.frontiersin.org/journals/forests-and-global-change/articles/10.3389/ffgc.2025.1492622/full","role":"secondary","citation":"Frontiers in Forests and Global Change (2025). Bohemian Switzerland National Park bark-beetle outbreak study."},{"title":"Rückblick Borkenkäfer 2020 (Bark beetle 2020 review)","url":"https://www.wald.sachsen.de/rueckblick-borkenkaefer-2020.html","role":"secondary","citation":"Sachsenforst (2020). Rückblick Borkenkäfer 2020."}]},{"id":"BB0119","name":"Budafapuszta Arboretum bark-beetle outbreak (2003-2004)","country":"Hungary","y0":2003,"y1":2004,"lat":46.510503,"lon":16.699948,"sev":"moderate","conf":"high","note":"Direct sources report 7-8 ha of destroyed conifer stands by 2003; field and laboratory work identified Ips typographus as the principal agent during the 2003-2004 spread.","sources":[{"title":"Bark beetle outbreak in the arboretum of Budafapuszta","url":"https://real-j.mtak.hu/23761/1/novenyvedelem_2006_1.pdf","role":"primary","citation":"Lakatos, F. & Kovács, K. (2006). Bark beetle outbreak in the arboretum of Budafapuszta. Növényvédelem 42(1):25-30."},{"title":"Pusztuló arborétum","url":"https://erdeszetilapok.oszk.hu/00294/pdf/00294.pdf","role":"secondary","citation":"Nemecz, F. (2003). Pusztuló arborétum. Erdészeti Lapok 138(12):372."}]},{"id":"BB0120","name":"Hungary reported spruce bark-beetle damage with Central Bükk hotspot (2022)","country":"Hungary","y0":2022,"y1":2022,"lat":48.04,"lon":20.53,"sev":"high","conf":"high","note":"OENyR recorded 258.63 ha of Ips typographus damage nationally in 2022, including 87.86 ha in Central Bükk.","sources":[{"title":"Magyarországi erdőkárok 2022-ben, a kárjelentések alapján","url":"https://real-j.mtak.hu/28908/1/2024_01_Novenyvedelem_boritoval.pdf","role":"primary","citation":"Hirka, A., Koltay, A. & Csóka, G. (2024). Magyarországi erdőkárok 2022-ben, a kárjelentések alapján. Növényvédelem 85 [N.S. 60](1):1-7."},{"title":"A 2022. évi biotikus és abiotikus erdőgazdasági károk, valamint a 2023-ban várható károsítások","url":"https://publicatio.uni-sopron.hu/2827/1/PROGNOZIS_2022_2023.pdf","role":"primary","citation":"Hirka, A. (2023). A 2022. évi biotikus és abiotikus erdőgazdasági károk, valamint a 2023-ban várható károsítások. SOE ERTI; NFK EFO. 212 pp. ISSN 1418-4427."}]},{"id":"BB0121","name":"Sopron Mountains spruce bark-beetle outbreak (2017-2020)","country":"Hungary","y0":2017,"y1":2020,"lat":47.67,"lon":16.5,"sev":"high","conf":"high","note":"Peer-reviewed Sentinel-2 and OENyR study reports 450 ha of spruce-dominated forest affected during the 2017-2020 outbreak wave.","sources":[{"title":"Sentinel-2 satellite-based analysis of bark beetle damage in Sopron Mountains Hungary","url":"https://ojs.lib.unideb.hu/landsenv/article/download/12920/11883/29799","role":"primary","citation":"Molnár, T. & Király, G. (2023). Sentinel-2 satellite-based analysis of bark beetle damage in Sopron Mountains, Hungary. Landscape & Environment 17(1):33-40. https://doi.org/10.21120/LE/17/1/3"}]},{"id":"BB0122","name":"Estonia drought-triggered spruce bark-beetle outbreak (1992-1995)","country":"Estonia","y0":1992,"y1":1995,"lat":58.6,"lon":25.01,"sev":"high","conf":"high","note":"The Estonian Environment Agency's historical chronology directly lists a major 1992-1995 outbreak triggered by the 1992 drought.","sources":[{"title":"Kuuse-kooreüraskist Eestis (About the spruce bark beetle in Estonia)","url":"https://www.keskkonnaagentuur.ee/sites/default/files/kuuse_kooreyraskist_eestis.pdf","role":"primary","citation":"Õunap, H. Kuuse-kooreüraskist Eestis. Keskkonnaagentuur."}]},{"id":"BB0123","name":"Sopron Mountains historical spruce bark-beetle catastrophe (1993-2003)","country":"Hungary","y0":1993,"y1":2003,"lat":47.67,"lon":16.5,"sev":"high","conf":"high","note":"A contemporaneous forestry-enterprise account reports about 130,000 m3 of beetle-damaged spruce removed from 1993 to 1997; a later operational review records about 180,000 m3 removed through 2003.","sources":[{"title":"A „Soproni fenyveserdő” végnapjai ...?","url":"https://erdeszetilapok.oszk.hu/01394/pdf/EL_132_1997_10_297-320.pdf","role":"primary","citation":"Jurásek, F. & Román, J. (1997). A „Soproni fenyveserdő” végnapjai ...? Erdészeti Lapok 132(10):313-315."},{"title":"Erdőfelújítás szúkatasztrófát követően","url":"https://erdeszetilapok.oszk.hu/01794/pdf/EPA01192_erdeszeti_lapok_2014_07-08_225-248.pdf","role":"secondary","citation":"Csapó, J. & Lakatos, F. (2014). Erdőfelújítás szúkatasztrófát követően. Erdészeti Lapok 149(7-8):228."}]},{"id":"BB0125","name":"Belgian Ardennes post-windthrow spruce bark-beetle outbreak (1990-1993)","country":"Belgium","y0":1990,"y1":1993,"lat":50.2,"lon":5.8,"sev":"high","conf":"high","note":"A Belgian outbreak is explicitly dated to 1990-1993 after the severe 1990 windthrows; contemporary Walloon ministry publications document the explosive attack and forest-protection response.","sources":[{"title":"Protection des forêts contre l'Ips typographus","url":"https://bibliotheques.wallonie.be/index.php?id=1785&l_typdoc=a%2Ce%2Cw%2Cx&lvl=publisher_see&nb_per_page_custom=200&nbr_lignes=230&page=1","role":"primary","citation":"Grégoire, J.-C. & De Proft, M. (1996). Protection des forêts contre l'Ips typographus. Fiche technique no. 6. Ministère de la Région Wallonne, DGRNE. 16 pp."},{"title":"Ips typographus contre épicéa: hygiène forestière et arbre-piège","url":"https://bibliotheques.wallonie.be/index.php?id=150&l_typdoc=&lvl=publisher_see&nb_per_page_custom=50&nbr_lignes=1303&page=27","role":"primary","citation":"Ministère de la Région wallonne, Division de la Nature et des Forêts (1992). Ips typographus contre épicéa: hygiène forestière et arbre-piège."},{"title":"Mitigating bark beetle damage in Norway spruce stands: Insights from Belgium's latest outbreak","url":"https://doi.org/10.1016/j.foreco.2025.123136","role":"secondary","citation":"Gilles, A., Lisein, J., Latte, N. & Claessens, H. (2025). Mitigating bark beetle damage in Norway spruce stands: Insights from Belgium's latest outbreak. Forest Ecology and Management 597:123136. https://doi.org/10.1016/j.foreco.2025.123136"}]},{"id":"BB0126","name":"Lithuania nationwide spruce bark-beetle outbreak (1994-1997)","country":"Lithuania","y0":1994,"y1":1997,"lat":55.3,"lon":23.5,"sev":"extreme","conf":"high","note":"Lithuania's Environment Ministry identifies 1994-1997 as the country's largest recorded Ips typographus mass outbreak, with more than 8 million m3 of spruce killed.","sources":[{"title":"Žievėgraužio tipografo daromos žalos ir populiacijos kitimo prognozė 2024 metams","url":"https://am.lrv.lt/lt/naujienos/zievegrauzio-tipografo-daromos-zalos-ir-populiacijos-kitimo-prognoze-2024-metams/","role":"primary","citation":"Lietuvos Respublikos aplinkos ministerija (2024). Žievėgraužio tipografo daromos žalos ir populiacijos kitimo prognozė 2024 metams."},{"title":"Outbreak of bark beetles in Lithuania","url":"https://www.euforgen.org/uploads/tx_news/493_Picea_abies_network__report_of_the_2nd_meeting__5-7_September_1996__Hyytia%CC%88la%CC%88__Finland.PA02.pdf","role":"primary","citation":"Cirtautas, V. G. (1997). Outbreak of bark beetles in Lithuania. In: Turok, J. & Koski, V. (compilers), Picea abies Network, Report of the Second Meeting, 5-7 September 1996, Hyytiälä, Finland, pp. 23-24. International Plant Genetic Resources Institute, Rome."}]},{"id":"BB0127","name":"Friuli-Venezia Giulia first monitored spruce bark-beetle outbreak (1995-2001)","country":"Italy","y0":1995,"y1":2001,"lat":46.4,"lon":13.0,"sev":"moderate","conf":"high","note":"The regional forest-health inventory identifies Friuli-Venezia Giulia's first Ips typographus outbreak in 1995-1996; six years of field monitoring show damage declining through 2001.","sources":[{"title":"Stato fitosanitario delle foreste del Friuli Venezia Giulia nel 2005","url":"https://www.ersa.fvg.it/export/sites/ersa/aziende/monitoraggi/allegati/BAUSINVE-2005.pdf","role":"primary","citation":"Stergulc, F., Frigimelica, G., Zandigiacomo, P., Osler, R. & Carpanelli, A. (2007). Stato fitosanitario delle foreste del Friuli Venezia Giulia nel 2005. Regione Autonoma Friuli Venezia Giulia / Università degli Studi di Udine. 56 pp."},{"title":"Ips typographus (L.) (Coleoptera: Scolytidae) in Southeastern Alps: Results of a Six-Year-Long Monitoring Program","url":"https://research.fs.usda.gov/treesearch/download/6071.pdf","role":"primary","citation":"Stergulc, F. & Faccoli, M. (2003). Ips typographus (L.) (Coleoptera: Scolytidae) in Southeastern Alps: Results of a Six-Year-Long Monitoring Program. In: Proceedings: Ecology, Survey and Management of Forest Insects. USDA Forest Service General Technical Report NE-311, pp. 168-169."}]},{"id":"BB0128","name":"Nazarje forest-management region spruce bark-beetle gradation (1992-1994)","country":"Slovenia","y0":1992,"y1":1994,"lat":46.32,"lon":14.95,"sev":"high","conf":"medium","note":"Slovenia's government-adopted Nazarje forest plan records a three-year spruce bark-beetle gradation in 1992-1994, with 73,578 m3 affected and intensity class 2-3.","sources":[{"title":"Gozdnogospodarski načrt gozdnogospodarskega območja Nazarje 2011-2020","url":"https://www.gov.si/assets/Ministrstva/MKGP/DOKUMENTI/GOZDARSTVO/Gozdnogospodarski-nacrti/Nazarje/aa9b57a4fe/10_NAZARJE_2011-2020.pdf","role":"primary","citation":"Zavod za gozdove Slovenije, Območna enota Nazarje (2012). Gozdnogospodarski načrt gozdnogospodarskega območja Nazarje 2011-2020, št. 10/11. Adopted by the Government of the Republic of Slovenia, 8 November 2012."},{"title":"Spruce bark beetles (Ips typographus, Pityogenes chalcographus, Col.: Scolytidae) in the Dinaric mountain forests of Slovenia: monitoring and modeling","url":"https://doi.org/10.1016/j.ecolmodel.2005.10.014","role":"secondary","citation":"Jurc, M. et al. (2006). Spruce bark beetles (Ips typographus, Pityogenes chalcographus, Col.: Scolytidae) in the Dinaric mountain forests of Slovenia: monitoring and modeling. Ecological Modelling 194(1-3):219-226. https://doi.org/10.1016/j.ecolmodel.2005.10.014"}]},{"id":"BB0129","name":"Great Britain great spruce bark-beetle establishment and outbreak response (1982-1989)","country":"England and Wales, UK","y0":1982,"y1":1989,"lat":52.5,"lon":-3.0,"sev":"moderate","conf":"high","note":"The Forestry Commission identified Dendroctonus micans as the cause of spruce deaths in 1982 and documented a well-established regional infestation through 1989.","sources":[{"title":"Dendroctonus micans in Britain - its Biology and Control","url":"https://cdn.forestresearch.gov.uk/1989/03/fcbu085.pdf","role":"primary","citation":"King, C. J. & Fielding, N. J. (1989). Dendroctonus micans in Britain - its Biology and Control. Forestry Commission Bulletin 85. HMSO, London."},{"title":"Distribution and Spread of the Great European Spruce Bark Beetle, Dendroctonus micans, in Britain - 1982 to 1989","url":"https://academic.oup.com/forestry/article/64/4/345/565308","role":"primary","citation":"Fielding, N. J., Evans, H. F., Williams, J. M. & Evans, B. (1991). Distribution and Spread of the Great European Spruce Bark Beetle, Dendroctonus micans, in Britain - 1982 to 1989. Forestry 64(4):345-358. https://doi.org/10.1093/forestry/64.4.345"},{"title":"Pest categorisation of Dendroctonus micans","url":"https://efsa.onlinelibrary.wiley.com/doi/10.2903/j.efsa.2017.4880","role":"secondary","citation":"EFSA Panel on Plant Health (2017). Pest categorisation of Dendroctonus micans. EFSA Journal 15(7):e04880. https://doi.org/10.2903/j.efsa.2017.4880"}]},{"id":"BB0130","name":"Southern Massif Central great spruce bark-beetle outbreak (1979-1985 documented phase)","country":"France","y0":1979,"y1":1985,"lat":44.9,"lon":3.9,"sev":"high","conf":"high","note":"Contemporaneous research describes outbreak-level Dendroctonus micans damage in many Massif Central spruce stands; dendrochronology dates the Ardèche outbreak onset to 1979.","sources":[{"title":"Effects of climate on radial growth of Norway spruce and interactions with attacks by the bark beetle Dendroctonus micans: a dendroecological study in the French Massif Central","url":"https://www.sciencedirect.com/science/article/abs/pii/S037811270400516X","role":"primary","citation":"Rolland, C. & Lempérière, G. (2004). Effects of climate on radial growth of Norway spruce and interactions with attacks by the bark beetle Dendroctonus micans: a dendroecological study in the French Massif Central. Forest Ecology and Management 201:89-104. https://doi.org/10.1016/j.foreco.2004.05.059"},{"title":"Biocontrol of Dendroctonus micans by Rhizophagus grandis in Massif Central (France): a first appraisal of the mass-rearing and release methods","url":"https://doi.org/10.1111/j.1439-0418.1985.tb01977.x","role":"primary","citation":"Grégoire, J.-C., Merlin, J., Pasteels, J. M., Jaffuel, R., Vouland, G. & Schvester, D. (1985). Biocontrol of Dendroctonus micans by Rhizophagus grandis in Massif Central (France): a first appraisal of the mass-rearing and release methods. Zeitschrift für Angewandte Entomologie 99:182-190. https://doi.org/10.1111/j.1439-0418.1985.tb01977.x"},{"title":"Pest categorisation of Dendroctonus micans","url":"https://efsa.onlinelibrary.wiley.com/doi/10.2903/j.efsa.2017.4880","role":"secondary","citation":"EFSA Panel on Plant Health (2017). Pest categorisation of Dendroctonus micans. EFSA Journal 15(7):e04880. https://doi.org/10.2903/j.efsa.2017.4880"}]},{"id":"BB0131","name":"Dalarna-Varmland post-1969-storm spruce bark-beetle outbreak (1976-1979 conservative window)","country":"Sweden","y0":1976,"y1":1979,"lat":60.5,"lon":14.0,"sev":"high","conf":"medium","note":"EPPO reports 2.0 million m3 of wood loss in Sweden during 1976-1979; an SLU historical synthesis places most damage from the broader post-1969 episode in Dalarna and Varmland.","sources":[{"title":"The spruce bark beetle Ips typographus in a changing climate - Effects of weather conditions on the biology of Ips typographus","url":"https://pub.epsilon.slu.se/id/document/12933","role":"secondary","citation":"Öhrn, P. (2012). The spruce bark beetle Ips typographus in a changing climate - Effects of weather conditions on the biology of Ips typographus. Introductory Research Essay No. 18. Department of Ecology, Swedish University of Agricultural Sciences, Uppsala. 27 pp."},{"title":"Ips typographus (IPSXTY) datasheet","url":"https://gd.eppo.int/taxon/IPSXTY/datasheet","role":"secondary","citation":"EPPO Global Database. Ips typographus (IPSXTY) datasheet. https://gd.eppo.int/taxon/IPSXTY/datasheet (last updated 2020-09-01)."},{"title":"Destructive Storms in European Forests: Past and Forthcoming Impacts","url":"https://op.europa.eu/en/publication-detail/-/publication/a2488068-31a8-4c90-a724-e563472c053c","role":"secondary","citation":"Gardiner, B., Blennow, K., Carnus, J.-M. et al. (2010). Destructive Storms in European Forests: Past and Forthcoming Impacts. Final report to the European Commission, DG Environment. European Forest Institute. 138 pp."}]},{"id":"BB0132","name":"Germany post-1972-storm spruce bark-beetle outbreak (1972-1978)","country":"Germany","y0":1972,"y1":1978,"lat":51.2,"lon":10.4,"sev":"high","conf":"medium","note":"A review chapter reports more than 710,370 m3 of Ips typographus damage in Germany during 1972-1978 after 9.2 million m3 of spruce windthrow.","sources":[{"title":"Economics and Politics of Bark Beetles","url":"https://doi.org/10.1016/B978-0-12-417156-5.00015-0","role":"secondary","citation":"Grégoire, J.-C., Raffa, K. F. & Lindgren, B. S. (2015). Economics and Politics of Bark Beetles. In: Vega, F. E. & Hofstetter, R. W. (eds), Bark Beetles: Biology and Ecology of Native and Invasive Species, pp. 585-613. Academic Press. https://doi.org/10.1016/B978-0-12-417156-5.00015-0"},{"title":"Destructive Storms in European Forests: Past and Forthcoming Impacts","url":"https://op.europa.eu/en/publication-detail/-/publication/a2488068-31a8-4c90-a724-e563472c053c","role":"secondary","citation":"Gardiner, B., Blennow, K., Carnus, J.-M. et al. (2010). Destructive Storms in European Forests: Past and Forthcoming Impacts. Final report to the European Commission, DG Environment. European Forest Institute. 138 pp."}]},{"id":"BB0133","name":"Bialowieza Forest spruce bark-beetle outbreak (1983-1988)","country":"Poland","y0":1983,"y1":1988,"lat":52.7,"lon":23.85,"sev":"high","conf":"medium","note":"A retrospective forest synthesis identifies 1983-1988 as a separate large-scale Ips typographus mass outbreak in Bialowieza Forest.","sources":[{"title":"Zmiany udziału świerka pospolitego w drzewostanach Puszczy Białowieskiej w kontekście dynamiki liczebności kornika drukarza Ips typographus (L.)","url":"https://www.lasy.gov.pl/pl/publikacje/prace-naukowe/stan-ekosystemow-lesnych-puszczy-bialowieskiej/@@download/file/Stan%20ekosystem%C3%B3w%20le%C5%9Bnych%20PB.pdf","role":"secondary","citation":"Gutowski, J. M. & Jaroszewicz, B. (2016). Zmiany udziału świerka pospolitego w drzewostanach Puszczy Białowieskiej w kontekście dynamiki liczebności kornika drukarza Ips typographus (L.). In A. Wikło (ed.), Stan ekosystemów leśnych Puszczy Białowieskiej, pp. 87-108. Centrum Informacyjne Lasów Państwowych, Warszawa."},{"title":"Ten years (2012-2021) of spruce bark beetle Ips typographus (L.) activity in the Białowieża Forest District of the Białowieża Primeval Forest","url":"https://doi.org/10.26202/sylwan.2022005","role":"secondary","citation":"Ginszt, T. & Laskowska-Ginszt, A. (2022). Ten years (2012-2021) of spruce bark beetle Ips typographus (L.) activity in the Białowieża Forest District of the Białowieża Primeval Forest. Sylwan, 166(3), 183-193. https://doi.org/10.26202/sylwan.2022005"}]},{"id":"BB0134","name":"Bialowieza Forest spruce bark-beetle outbreak (1994-1997)","country":"Poland","y0":1994,"y1":1997,"lat":52.7,"lon":23.85,"sev":"moderate","conf":"medium","note":"Separate secondary syntheses identify a 1994-1997 outbreak that culminated in 1995 at about 40,000 m3 of infested trees.","sources":[{"title":"Zmiany udziału świerka pospolitego w drzewostanach Puszczy Białowieskiej w kontekście dynamiki liczebności kornika drukarza Ips typographus (L.)","url":"https://www.lasy.gov.pl/pl/publikacje/prace-naukowe/stan-ekosystemow-lesnych-puszczy-bialowieskiej/@@download/file/Stan%20ekosystem%C3%B3w%20le%C5%9Bnych%20PB.pdf","role":"secondary","citation":"Gutowski, J. M. & Jaroszewicz, B. (2016). Zmiany udziału świerka pospolitego w drzewostanach Puszczy Białowieskiej w kontekście dynamiki liczebności kornika drukarza Ips typographus (L.). In A. Wikło (ed.), Stan ekosystemów leśnych Puszczy Białowieskiej, pp. 87-108. Centrum Informacyjne Lasów Państwowych, Warszawa."},{"title":"Mass outbreaks of the spruce bark beetle Ips typographus in the context of the controversies around the Białowieża Primeval Forest","url":"https://www.ibles.pl/wp-content/uploads/2023/05/DOI-10.1515-frp-2016-0033-en.pdf","role":"secondary","citation":"Grodzki, W. (2016). Mass outbreaks of the spruce bark beetle Ips typographus in the context of the controversies around the Białowieża Primeval Forest. Forest Research Papers, 77(4), 324-331. https://doi.org/10.1515/frp-2016-0033"},{"title":"Ten years (2012-2021) of spruce bark beetle Ips typographus (L.) activity in the Białowieża Forest District of the Białowieża Primeval Forest","url":"https://doi.org/10.26202/sylwan.2022005","role":"secondary","citation":"Ginszt, T. & Laskowska-Ginszt, A. (2022). Ten years (2012-2021) of spruce bark beetle Ips typographus (L.) activity in the Białowieża Forest District of the Białowieża Primeval Forest. Sylwan, 166(3), 183-193. https://doi.org/10.26202/sylwan.2022005"}]},{"id":"BB0135","name":"Bialowieza Forest spruce bark-beetle outbreak (2001-2004)","country":"Poland","y0":2001,"y1":2004,"lat":52.7,"lon":23.85,"sev":"high","conf":"medium","note":"Secondary syntheses distinguish a large-scale 2001-2003/2004 outbreak that culminated in 2003 at about 75,000 m3 of infested trees.","sources":[{"title":"Zmiany udziału świerka pospolitego w drzewostanach Puszczy Białowieskiej w kontekście dynamiki liczebności kornika drukarza Ips typographus (L.)","url":"https://www.lasy.gov.pl/pl/publikacje/prace-naukowe/stan-ekosystemow-lesnych-puszczy-bialowieskiej/@@download/file/Stan%20ekosystem%C3%B3w%20le%C5%9Bnych%20PB.pdf","role":"secondary","citation":"Gutowski, J. M. & Jaroszewicz, B. (2016). Zmiany udziału świerka pospolitego w drzewostanach Puszczy Białowieskiej w kontekście dynamiki liczebności kornika drukarza Ips typographus (L.). In A. Wikło (ed.), Stan ekosystemów leśnych Puszczy Białowieskiej, pp. 87-108. Centrum Informacyjne Lasów Państwowych, Warszawa."},{"title":"Mass outbreaks of the spruce bark beetle Ips typographus in the context of the controversies around the Białowieża Primeval Forest","url":"https://www.ibles.pl/wp-content/uploads/2023/05/DOI-10.1515-frp-2016-0033-en.pdf","role":"secondary","citation":"Grodzki, W. (2016). Mass outbreaks of the spruce bark beetle Ips typographus in the context of the controversies around the Białowieża Primeval Forest. Forest Research Papers, 77(4), 324-331. https://doi.org/10.1515/frp-2016-0033"},{"title":"Ten years (2012-2021) of spruce bark beetle Ips typographus (L.) activity in the Białowieża Forest District of the Białowieża Primeval Forest","url":"https://doi.org/10.26202/sylwan.2022005","role":"secondary","citation":"Ginszt, T. & Laskowska-Ginszt, A. (2022). Ten years (2012-2021) of spruce bark beetle Ips typographus (L.) activity in the Białowieża Forest District of the Białowieża Primeval Forest. Sylwan, 166(3), 183-193. https://doi.org/10.26202/sylwan.2022005"}]},{"id":"BB0136","name":"Bialowieza Forest spruce bark-beetle outbreak (2006-2009)","country":"Poland","y0":2006,"y1":2009,"lat":52.7,"lon":23.85,"sev":"moderate","conf":"medium","note":"A recent historical chronology identifies a 2006-2009 outbreak; a separate review records its distinct 2008 culmination at about 40,000 m3 of infested trees.","sources":[{"title":"Mass outbreaks of the spruce bark beetle Ips typographus in the context of the controversies around the Białowieża Primeval Forest","url":"https://www.ibles.pl/wp-content/uploads/2023/05/DOI-10.1515-frp-2016-0033-en.pdf","role":"secondary","citation":"Grodzki, W. (2016). Mass outbreaks of the spruce bark beetle Ips typographus in the context of the controversies around the Białowieża Primeval Forest. Forest Research Papers, 77(4), 324-331. https://doi.org/10.1515/frp-2016-0033"},{"title":"Ten years (2012-2021) of spruce bark beetle Ips typographus (L.) activity in the Białowieża Forest District of the Białowieża Primeval Forest","url":"https://doi.org/10.26202/sylwan.2022005","role":"secondary","citation":"Ginszt, T. & Laskowska-Ginszt, A. (2022). Ten years (2012-2021) of spruce bark beetle Ips typographus (L.) activity in the Białowieża Forest District of the Białowieża Primeval Forest. Sylwan, 166(3), 183-193. https://doi.org/10.26202/sylwan.2022005"}]},{"id":"BB0137","name":"Hatila Valley National Park post-windthrow spruce bark-beetle outbreak (2003-2009)","country":"Türkiye","y0":2003,"y1":2009,"lat":41.11,"lon":41.66,"sev":"high","conf":"medium","note":"A retrospective peer-reviewed account reports a severe 2003-2009 outbreak over 7,476 ha, with 260,012 m3 of Oriental spruce cut.","sources":[{"title":"Analyzing the Local Spread of Ips typographus (L.) (Coleoptera: Curculionidae Scolytinae) by Pheromone Catches in Turkey's Hatila Valley National Park","url":"https://www.pjoes.com/pdf-76363-24320?filename=Analyzing-the-Local-Sprea.pdf","role":"secondary","citation":"Alkan Akinci, H. & Aksu, Y. (2018). Analyzing the Local Spread of Ips typographus (L.) (Coleoptera: Curculionidae, Scolytinae) by Pheromone Catches in Turkey's Hatila Valley National Park. Polish Journal of Environmental Studies, 27(3), 979-985. https://doi.org/10.15244/pjoes/76363"}]},{"id":"BB0138","name":"Eastern Carpathians post-1995-windfall bark-beetle outbreak (1996-2001 documented phase)","country":"Romania","y0":1996,"y1":2001,"lat":46.6,"lon":25.5,"sev":"high","conf":"medium","note":"Field monitoring found increasing bark-beetle infestation after the 1995 windfall, reaching high and very high attack levels across all damaged areas in 2000.","sources":[{"title":"The Frequency and Intensity of Bark Beetle Infestations Recorded in Romanian Forests Affected by Windfall in 1995","url":"https://research.fs.usda.gov/treesearch/19355","role":"primary","citation":"Mihalciuc, V., Danci, A., Bujila, M. & Chira, D. (2003). The Frequency and Intensity of Bark Beetle Infestations Recorded in Romanian Forests Affected by Windfall in 1995. In M. L. McManus & A. M. Liebhold (eds), Proceedings: Ecology, Survey and Management of Forest Insects. USDA Forest Service General Technical Report NE-311, pp. 162-164."},{"title":"Report on Current State of Forest Resources in the Carpathians","url":"https://doczz.net/doc/519452/activity-2.7--forestry-and-timber-industry-current-state-...","role":"secondary","citation":"Anfodillo, T., Carrer, M., Dalla Valle, E., Giacoma, E., Lamedica, S. & Pettenella, D. (2008). Report on Current State of Forest Resources in the Carpathians. Activity 2.7: Forestry and Timber Industry. INTERREG III B CADSES Programme Carpathian Project, University of Padova."}]},{"id":"BB0139","name":"Ormoz great spruce bark-beetle gradation (2000-2001)","country":"Slovenia","y0":2000,"y1":2001,"lat":46.41,"lon":16.15,"sev":"minor","conf":"medium","note":"A Slovenian bark-beetle overview records a Dendroctonus micans gradation near Ormoz in 2000 and about 80 m3 of dying spruce in 24 small groups in 2001.","sources":[{"title":"Bark Beetles (Scolytidae Coleoptera) in Slovenia with Special Regard to Species in Burnt Pine Forests","url":"https://research.fs.usda.gov/treesearch/19342","role":"secondary","citation":"Jurc, M. (2003). Bark Beetles (Scolytidae, Coleoptera) in Slovenia with Special Regard to Species in Burnt Pine Forests. In M. L. McManus & A. M. Liebhold (eds), Proceedings: Ecology, Survey and Management of Forest Insects. USDA Forest Service General Technical Report NE-311, pp. 157-159."}]},{"id":"BB0140","name":"Denmark post-1947-drought great spruce bark-beetle outbreak (1947-1949)","country":"Denmark","y0":1947,"y1":1949,"lat":56.1,"lon":9.4,"sev":"moderate","conf":"medium","note":"A contemporaneous Forestry Commission synthesis records frequent Dendroctonus micans reports across Denmark from 1947 and identifies a 1947-1949 outbreak period, especially in Sitka-spruce plantations.","sources":[{"title":"The Great Spruce Bark-Beetle Dendroctonus micans in North West Europe","url":"https://cdn.forestresearch.gov.uk/1966/03/fcbu038.pdf","role":"secondary","citation":"Brown, J. M. B. & Bevan, D. (1966). The Great Spruce Bark-Beetle Dendroctonus micans in North West Europe. Forestry Commission Bulletin 38. London: Her Majesty's Stationery Office. 41 pp."},{"title":"The Greater European Spruce Beetle","url":"https://www.researchgate.net/publication/312987008_The_greater_European_spruce_beetle_In_AA_Berryman_ed","role":"secondary","citation":"Gregoire, J.-C. (1988). The Greater European Spruce Beetle. In A. A. Berryman (ed.), Dynamics of Forest Insect Populations: Patterns, Causes, Implications, pp. 455-478. New York: Plenum Press."}]},{"id":"BB0141","name":"Suderlugum Sitka-spruce great bark-beetle outbreak (1947-early 1950s)","country":"Germany","y0":1947,"y1":1952,"lat":54.87,"lon":8.91,"sev":"high","conf":"medium","note":"A historical forestry survey records severe Dendroctonus micans invasion after the 1947 drought at Suderlugum, followed by disappearance of the affected Sitka-spruce stand within a few years.","sources":[{"title":"The Great Spruce Bark-Beetle Dendroctonus micans in North West Europe","url":"https://cdn.forestresearch.gov.uk/1966/03/fcbu038.pdf","role":"secondary","citation":"Brown, J. M. B. & Bevan, D. (1966). The Great Spruce Bark-Beetle Dendroctonus micans in North West Europe. Forestry Commission Bulletin 38. London: Her Majesty's Stationery Office. 41 pp."},{"title":"The Greater European Spruce Beetle","url":"https://www.researchgate.net/publication/312987008_The_greater_European_spruce_beetle_In_AA_Berryman_ed","role":"secondary","citation":"Gregoire, J.-C. (1988). The Greater European Spruce Beetle. In A. A. Berryman (ed.), Dynamics of Forest Insect Populations: Patterns, Causes, Implications, pp. 455-478. New York: Plenum Press."}]},{"id":"BB0142","name":"Saar great spruce bark-beetle outbreak (1952-1954)","country":"Germany","y0":1952,"y1":1954,"lat":49.4,"lon":6.9,"sev":"minor","conf":"medium","note":"A major historical synthesis identifies a localized Dendroctonus micans outbreak in the Saar during 1952-1954 and states that it rapidly subsided.","sources":[{"title":"The Greater European Spruce Beetle","url":"https://www.researchgate.net/publication/312987008_The_greater_European_spruce_beetle_In_AA_Berryman_ed","role":"secondary","citation":"Gregoire, J.-C. (1988). The Greater European Spruce Beetle. In A. A. Berryman (ed.), Dynamics of Forest Insect Populations: Patterns, Causes, Implications, pp. 455-478. New York: Plenum Press."}]},{"id":"BB0143","name":"Morvan great spruce bark-beetle outbreak (1959-1960)","country":"France","y0":1959,"y1":1960,"lat":47.08,"lon":4.08,"sev":"moderate","conf":"medium","note":"Historical syntheses identify a Dendroctonus micans outbreak in Morvan Norway-spruce forests during 1959-1960 that subsided soon afterward.","sources":[{"title":"The Greater European Spruce Beetle","url":"https://www.researchgate.net/publication/312987008_The_greater_European_spruce_beetle_In_AA_Berryman_ed","role":"secondary","citation":"Gregoire, J.-C. (1988). The Greater European Spruce Beetle. In A. A. Berryman (ed.), Dynamics of Forest Insect Populations: Patterns, Causes, Implications, pp. 455-478. New York: Plenum Press."},{"title":"Deteccion de Dendroctonus micans (Kugelann, 1794) (Coleoptera: Curculionidae: Scolytinae) y confirmacion de la presencia de Ips typographus (Linne, 1758) en plantaciones forestales de Picea abies en el Valle de Aran (Lleida, Espana)","url":"https://dialnet.unirioja.es/servlet/articulo?codigo=7997324","role":"secondary","citation":"Riba Flinch, J.M. (2021). Deteccion de Dendroctonus micans (Kugelann, 1794) y confirmacion de la presencia de Ips typographus (Linne, 1758) en plantaciones forestales de Picea abies en el Valle de Aran (Lleida, Espana). Boletin de la SEA 68."}]},{"id":"BB0144","name":"Georgian Caucasus great spruce bark-beetle invasion outbreak (1956-1979)","country":"Georgia","y0":1956,"y1":1979,"lat":42.2,"lon":43.3,"sev":"extreme","conf":"medium","note":"A historical synthesis reports rapid invasion after 1956, with about 135,000 ha of Georgian Oriental spruce invaded and seriously damaged by 1976 before decline late in the decade.","sources":[{"title":"The Greater European Spruce Beetle","url":"https://www.researchgate.net/publication/312987008_The_greater_European_spruce_beetle_In_AA_Berryman_ed","role":"secondary","citation":"Gregoire, J.-C. (1988). The Greater European Spruce Beetle. In A. A. Berryman (ed.), Dynamics of Forest Insect Populations: Patterns, Causes, Implications, pp. 455-478. New York: Plenum Press."},{"title":"Pest categorisation of Dendroctonus micans","url":"https://efsa.onlinelibrary.wiley.com/doi/10.2903/j.efsa.2017.4880","role":"secondary","citation":"EFSA Panel on Plant Health (2017). Pest categorisation of Dendroctonus micans. EFSA Journal 15(7):e04880. https://doi.org/10.2903/j.efsa.2017.4880"},{"title":"Late Glacial and Holocene Records of Tree-Killing Conifer Bark Beetles in Europe and North America: Implications for Forest Disturbance Dynamics","url":"https://research.fs.usda.gov/download/treesearch/59981.pdf","role":"secondary","citation":"Schafstall, N., Kuosmanen, N., Fettig, C. J., Knizek, M. & Clear, J. L. (2020). Late Glacial and Holocene records of tree-killing conifer bark beetles in Europe and North America: Implications for forest disturbance dynamics. The Holocene 30(6):847-857. https://doi.org/10.1177/0959683620902214"},{"title":"Deteccion de Dendroctonus micans (Kugelann, 1794) (Coleoptera: Curculionidae: Scolytinae) y confirmacion de la presencia de Ips typographus (Linne, 1758) en plantaciones forestales de Picea abies en el Valle de Aran (Lleida, Espana)","url":"https://dialnet.unirioja.es/servlet/articulo?codigo=7997324","role":"secondary","citation":"Riba Flinch, J.M. (2021). Deteccion de Dendroctonus micans (Kugelann, 1794) y confirmacion de la presencia de Ips typographus (Linne, 1758) en plantaciones forestales de Picea abies en el Valle de Aran (Lleida, Espana). Boletin de la SEA 68."}]},{"id":"BB0145","name":"Northeastern Turkiye great spruce bark-beetle expansion outbreak (1966-1984 documented phase)","country":"Türkiye","y0":1966,"y1":1984,"lat":40.95,"lon":41.3,"sev":"extreme","conf":"medium","note":"Historical sources document invasion from 1966, 326,244 felled or burned trees in 1967-1972, and an infested area reaching about 55,000 ha in the early 1980s.","sources":[{"title":"The Greater European Spruce Beetle","url":"https://www.researchgate.net/publication/312987008_The_greater_European_spruce_beetle_In_AA_Berryman_ed","role":"secondary","citation":"Gregoire, J.-C. (1988). The Greater European Spruce Beetle. In A. A. Berryman (ed.), Dynamics of Forest Insect Populations: Patterns, Causes, Implications, pp. 455-478. New York: Plenum Press."},{"title":"Attack Strategy and Development of Dendroctonus micans on Oriental Spruce in Turkey","url":"https://dergipark.org.tr/en/download/article-file/65254","role":"secondary","citation":"Alkan-Akinci, H., Eroglu, M. & Ozcan, G. E. (2014). Attack Strategy and Development of Dendroctonus micans (Kug.) (Coleoptera: Curculionidae) on Oriental Spruce in Turkey. Turkish Journal of Entomology 38(1):31-41. https://doi.org/10.16970/ted.68484"},{"title":"Pest categorisation of Dendroctonus micans","url":"https://efsa.onlinelibrary.wiley.com/doi/10.2903/j.efsa.2017.4880","role":"secondary","citation":"EFSA Panel on Plant Health (2017). Pest categorisation of Dendroctonus micans. EFSA Journal 15(7):e04880. https://doi.org/10.2903/j.efsa.2017.4880"},{"title":"Late Glacial and Holocene Records of Tree-Killing Conifer Bark Beetles in Europe and North America: Implications for Forest Disturbance Dynamics","url":"https://research.fs.usda.gov/download/treesearch/59981.pdf","role":"secondary","citation":"Schafstall, N., Kuosmanen, N., Fettig, C. J., Knizek, M. & Clear, J. L. (2020). Late Glacial and Holocene records of tree-killing conifer bark beetles in Europe and North America: Implications for forest disturbance dynamics. The Holocene 30(6):847-857. https://doi.org/10.1177/0959683620902214"},{"title":"Deteccion de Dendroctonus micans (Kugelann, 1794) (Coleoptera: Curculionidae: Scolytinae) y confirmacion de la presencia de Ips typographus (Linne, 1758) en plantaciones forestales de Picea abies en el Valle de Aran (Lleida, Espana)","url":"https://dialnet.unirioja.es/servlet/articulo?codigo=7997324","role":"secondary","citation":"Riba Flinch, J.M. (2021). Deteccion de Dendroctonus micans (Kugelann, 1794) y confirmacion de la presencia de Ips typographus (Linne, 1758) en plantaciones forestales de Picea abies en el Valle de Aran (Lleida, Espana). Boletin de la SEA 68."}]},{"id":"BB0146","name":"Western Estonia great spruce bark-beetle pine outbreak (1971-1979 documented phase)","country":"Estonia","y0":1971,"y1":1979,"lat":58.6,"lon":23.7,"sev":"high","conf":"medium","note":"Primary Estonian research documents a Dendroctonus micans damage outbreak in young western-Estonian Scots-pine stands during the second half of the 1970s; a historical review dates severe attacks from 1971.","sources":[{"title":"Hibernation Peculiarities and Cold-Hardiness of the Great Spruce Bark Beetle Dendroctonus micans Kug.","url":"https://kirj.ee/wp-content/plugins/kirj/pub/biol-3-1990-214-218_20230304162114.pdf","role":"primary","citation":"Luik, A. & Voolma, K. (1990). Hibernation Peculiarities and Cold-Hardiness of the Great Spruce Bark Beetle, Dendroctonus micans Kug. Proceedings of the Estonian Academy of Sciences, Biology 39(3):214-218. https://doi.org/10.3176/biol.1990.3.07"},{"title":"The Greater European Spruce Beetle","url":"https://www.researchgate.net/publication/312987008_The_greater_European_spruce_beetle_In_AA_Berryman_ed","role":"secondary","citation":"Gregoire, J.-C. (1988). The Greater European Spruce Beetle. In A. A. Berryman (ed.), Dynamics of Forest Insect Populations: Patterns, Causes, Implications, pp. 455-478. New York: Plenum Press."},{"title":"Deteccion de Dendroctonus micans (Kugelann, 1794) (Coleoptera: Curculionidae: Scolytinae) y confirmacion de la presencia de Ips typographus (Linne, 1758) en plantaciones forestales de Picea abies en el Valle de Aran (Lleida, Espana)","url":"https://dialnet.unirioja.es/servlet/articulo?codigo=7997324","role":"secondary","citation":"Riba Flinch, J.M. (2021). Deteccion de Dendroctonus micans (Kugelann, 1794) y confirmacion de la presencia de Ips typographus (Linne, 1758) en plantaciones forestales de Picea abies en el Valle de Aran (Lleida, Espana). Boletin de la SEA 68."}]},{"id":"BB0147","name":"Northeastern France post-1990-storm spruce bark-beetle outbreak (1990-1997)","country":"France","y0":1990,"y1":1997,"lat":48.25,"lon":6.75,"sev":"high","conf":"medium","note":"French and EU sources document a 1990-1997 post-storm bark-beetle outbreak, including 425,500 m3 of Ips typographus damage in 1991-1993 and more than 500,000 m3 killed by Ips typographus and Pityogenes chalcographus in 1992-1994.","sources":[{"title":"Ophiostomatoid Fungi Associated with the Spruce Bark Beetle, Ips typographus, in Three Areas in France","url":"https://www.afs-journal.org/articles/forest/pdf/2004/03/F4302.pdf","role":"secondary","citation":"Viiri, H. & Lieutier, F. (2004). Ophiostomatoid fungi associated with the spruce bark beetle, Ips typographus, in three areas in France. Annals of Forest Science 61:215-219. https://doi.org/10.1051/forest:2004013"},{"title":"Bark and Wood Boring Insects in Living Trees","url":"https://cordis.europa.eu/project/id/E16","role":"secondary","citation":"European Commission (2003). Bark and Wood Boring Insects in Living Trees. CORDIS fact sheet for COST Action E16, project period 1998-2002, last updated 5 March 2003."},{"title":"La pullulation de scolytes dans le nord-est de la France de 1990 a 1997 suite aux tempetes de fevrier 1990","url":"https://agriculture.gouv.fr/sites/default/files/documents/pdf/Publications_cle01b31d.pdf","role":"secondary","citation":"Nageleisen, L.-M. (2000). La pullulation de scolytes dans le nord-est de la France de 1990 a 1997 suite aux tempetes de fevrier 1990. Les Cahiers du DSF 1-2000 (La sante des forets en 1999), pp. 45-49."}]},{"id":"BB0148","name":"Eastern Belgium great spruce bark-beetle invasion outbreak (1897-1900)","country":"Belgium","y0":1897,"y1":1900,"lat":50.58,"lon":6.06,"sev":"high","conf":"medium","note":"A historical Forestry Commission synthesis records heavy cutting of Dendroctonus-damaged spruce in eastern Belgium throughout 1897-1900.","sources":[{"title":"The Great Spruce Bark-Beetle Dendroctonus micans in North West Europe","url":"https://cdn.forestresearch.gov.uk/1966/03/fcbu038.pdf","role":"primary","citation":"Brown, J. M. B. & Bevan, D. (1966). The Great Spruce Bark-Beetle Dendroctonus micans in North West Europe. Forestry Commission Bulletin 38. London: Her Majesty's Stationery Office. 41 pp."},{"title":"The Greater European Spruce Beetle","url":"https://www.researchgate.net/publication/312987008_The_greater_European_spruce_beetle_In_AA_Berryman_ed","role":"secondary","citation":"Gregoire, J.-C. (1988). The Greater European Spruce Beetle. In A. A. Berryman (ed.), Dynamics of Forest Insect Populations: Patterns, Causes, Implications, pp. 455-478. New York: Plenum Press."}]},{"id":"BB0149","name":"Northern Uppland post-1931-storm spruce bark-beetle outbreak (1932-1935)","country":"Sweden","y0":1932,"y1":1935,"lat":60.25,"lon":17.7,"sev":"moderate","conf":"medium","note":"A Swedish research synthesis records an Ips typographus mass propagation in northern Uppland during 1932-1935 after the December 1931 storm.","sources":[{"title":"The spruce bark beetle Ips typographus in a changing climate - Effects of weather conditions on the biology of Ips typographus","url":"https://pub.epsilon.slu.se/id/document/12933","role":"secondary","citation":"Öhrn, P. (2012). The spruce bark beetle Ips typographus in a changing climate - Effects of weather conditions on the biology of Ips typographus. Introductory Research Essay No. 18. Department of Ecology, Swedish University of Agricultural Sciences, Uppsala. 27 pp."}]},{"id":"BB0150","name":"Nas Dalarna spruce bark-beetle outbreak (1935-1938)","country":"Sweden","y0":1935,"y1":1938,"lat":60.45,"lon":14.49,"sev":"moderate","conf":"medium","note":"A Swedish research synthesis identifies a distinct Ips typographus outbreak at Nas in Dalarna from 1935 through 1938.","sources":[{"title":"The spruce bark beetle Ips typographus in a changing climate - Effects of weather conditions on the biology of Ips typographus","url":"https://pub.epsilon.slu.se/id/document/12933","role":"secondary","citation":"Öhrn, P. (2012). The spruce bark beetle Ips typographus in a changing climate - Effects of weather conditions on the biology of Ips typographus. Introductory Research Essay No. 18. Department of Ecology, Swedish University of Agricultural Sciences, Uppsala. 27 pp."}]},{"id":"BB0152","name":"Gaestrikland and southeastern Dalarna post-1954-storm outbreak (1956-1957)","country":"Sweden","y0":1956,"y1":1957,"lat":60.7,"lon":16.0,"sev":"high","conf":"medium","note":"Localized Ips typographus outbreaks in 1956 prompted a 1957 campaign using 4,514 treated standing trap trees.","sources":[{"title":"The spruce bark beetle Ips typographus in a changing climate - Effects of weather conditions on the biology of Ips typographus","url":"https://pub.epsilon.slu.se/id/document/12933","role":"secondary","citation":"Öhrn, P. (2012). The spruce bark beetle Ips typographus in a changing climate - Effects of weather conditions on the biology of Ips typographus. Introductory Research Essay No. 18. Department of Ecology, Swedish University of Agricultural Sciences, Uppsala. 27 pp."}]},{"id":"BB0153","name":"Southern Sweden six-toothed spruce bark-beetle outbreak (1970-1971)","country":"Sweden","y0":1970,"y1":1971,"lat":57.2,"lon":14.5,"sev":"moderate","conf":"medium","note":"A peer-reviewed pest review reports that many young spruces in southern Sweden were attacked and killed by Pityogenes chalcographus in 1970-1971.","sources":[{"title":"Pityogenes chalcographus (Coleoptera: Curculionidae): biology pest status and current management options","url":"https://academic.oup.com/jipm/article/16/1/11/8112036","role":"secondary","citation":"Holusa, J. & Fiala, T. (2025). Pityogenes chalcographus (Coleoptera: Curculionidae): biology, pest status, and current management options. Journal of Integrated Pest Management 16(1):11. https://doi.org/10.1093/jipm/pmaf007"}]},{"id":"BB0154","name":"Czech border-mountain post-windstorm bark-beetle calamity (1868-1877)","country":"Czechia","y0":1868,"y1":1877,"lat":49.1,"lon":13.5,"sev":"extreme","conf":"medium","note":"A documentary reconstruction identifies the most disastrous Czech bark-beetle outbreak of the nineteenth century, continuing from the 1868 and 1870 storms until 1875-1877.","sources":[{"title":"Meteorological and climatological triggers of notable past and present bark beetle outbreaks in the Czech Republic","url":"https://doi.org/10.5194/cp-18-2155-2022","role":"primary","citation":"Brazdil, R. et al. (2022). Meteorological and climatological triggers of notable past and present bark beetle outbreaks in the Czech Republic. Climate of the Past 18:2155-2180. https://doi.org/10.5194/cp-18-2155-2022"}]},{"id":"BB0156","name":"Southwestern Germany postwar spruce bark-beetle calamity (1944-1951)","country":"Germany","y0":1944,"y1":1951,"lat":48.3,"lon":8.3,"sev":"extreme","conf":"medium","note":"Historical syntheses identify a 1944-1951 southwest-German calamity; a 1947 report records 2.4 million festmeters of new Bavarian damage in 1946 and major damage in Baden and Wuerttemberg.","sources":[{"title":"Insect Pests in North-West Germany During 1946-1947","url":"https://cdn.forestresearch.gov.uk/1948/04/fcjo019.pdf","role":"primary","citation":"Hanson, H. S. (1948). Insect Pests in North-West Germany During 1946-1947. Journal of the Forestry Commission 19:104-107."},{"title":"Meteorological and climatological triggers of notable past and present bark beetle outbreaks in the Czech Republic","url":"https://doi.org/10.5194/cp-18-2155-2022","role":"secondary","citation":"Brazdil, R. et al. (2022). Meteorological and climatological triggers of notable past and present bark beetle outbreaks in the Czech Republic. Climate of the Past 18:2155-2180. https://doi.org/10.5194/cp-18-2155-2022"}]},{"id":"BB0157","name":"Swiss spruce bark-beetle mass outbreak (1946-1949)","country":"Switzerland","y0":1946,"y1":1949,"lat":47.1,"lon":7.5,"sev":"high","conf":"medium","note":"A contemporaneous Swiss monograph is dedicated to the 1946-1949 Ips typographus mass outbreak; later WSL synthesis documents severe lowland and western-Swiss bark-beetle damage.","sources":[{"title":"Das Massenauftreten des achtzaehnigen Fichtenborkenkaefers Ips typographus L. nach Untersuchungen in schweizerischen Waldungen 1946 bis 1949","url":"https://www.dora.lib4ri.ch/wsl/islandora/object/wsl%3A24080","role":"primary","citation":"Kuhn, W. (1949). Das Massenauftreten des achtzaehnigen Fichtenborkenkaefers Ips typographus L. nach Untersuchungen in schweizerischen Waldungen 1946 bis 1949. Mitteilungen der Schweizerischen Anstalt fuer das Forstliche Versuchswesen 26(1):245-330. https://doi.org/10.3929/ethz-a-000187269"},{"title":"Effects of the 1947 and 2018 droughts on Swiss forests","url":"https://www.dora.lib4ri.ch/wsl/islandora/object/wsl%3A32101/datastream/PDF/Wohlgemuth-2022-Effects_of_the_1947_and-%28published_version%29.pdf","role":"secondary","citation":"Wohlgemuth, T. & Buergi, M. (2022). Effects of the 1947 and 2018 droughts on Swiss forests. In Hot and Dry Summers in Switzerland: Causes and Impacts of the Record Summers 1947, 2003, and 2018, pp. 18-19."}]},{"id":"BB0158","name":"Swiss crooked-toothed fir bark-beetle outbreak (1947-1949)","country":"Switzerland","y0":1947,"y1":1949,"lat":47.05,"lon":7.5,"sev":"high","conf":"high","note":"WSL reports more than 200,000 m3 of fir attacked nationwide in the two years following the hot and dry summer of 1947, especially in the Swiss Plateau and Jura.","sources":[{"title":"Der Krummzaehnige Weisstannenborkenkaefer im Klimawandel","url":"https://www.waldwissen.net/assets/waldwirtschaft/schaden/insekten/wsl_weisstannenborkenkaefer_klimawandel/Blaser-2022-Der_Krummzaehnige_Weisstannenborkenkaefer_im_Klimawandel-_published_version_.pdf","role":"primary","citation":"Blaser, S. (2022). Der Krummzaehnige Weisstannenborkenkaefer im Klimawandel. Buendner Wald 75(6):32-35."},{"title":"Effects of the 1947 and 2018 droughts on Swiss forests","url":"https://www.dora.lib4ri.ch/wsl/islandora/object/wsl%3A32101/datastream/PDF/Wohlgemuth-2022-Effects_of_the_1947_and-%28published_version%29.pdf","role":"secondary","citation":"Wohlgemuth, T. & Buergi, M. (2022). Effects of the 1947 and 2018 droughts on Swiss forests. In Hot and Dry Summers in Switzerland: Causes and Impacts of the Record Summers 1947, 2003, and 2018, pp. 18-19."}]},{"id":"BB0159","name":"Bialowieza Forest post-windstorm spruce bark-beetle outbreak (1882-1883)","country":"Poland","y0":1882,"y1":1883,"lat":52.7,"lon":23.85,"sev":"moderate","conf":"medium","note":"Peer-reviewed historical chronology identifies a documented Ips typographus mass outbreak in Bialowieza Forest during 1882-1883.","sources":[{"title":"Ten years (2012-2021) of spruce bark beetle Ips typographus (L.) activity in the Białowieża Forest District of the Białowieża Primeval Forest","url":"https://doi.org/10.26202/sylwan.2022005","role":"secondary","citation":"Ginszt, T. & Laskowska-Ginszt, A. (2022). Ten years (2012-2021) of spruce bark beetle Ips typographus (L.) activity in the Białowieża Forest District of the Białowieża Primeval Forest. Sylwan, 166(3), 183-193. https://doi.org/10.26202/sylwan.2022005"}]},{"id":"BB0160","name":"Bialowieza Forest interwar spruce bark-beetle outbreak (1919-1922)","country":"Poland","y0":1919,"y1":1922,"lat":52.7,"lon":23.85,"sev":"high","conf":"medium","note":"Peer-reviewed historical chronology identifies an Ips typographus mass outbreak in Bialowieza Forest from 1919 through 1922.","sources":[{"title":"Ten years (2012-2021) of spruce bark beetle Ips typographus (L.) activity in the Białowieża Forest District of the Białowieża Primeval Forest","url":"https://doi.org/10.26202/sylwan.2022005","role":"secondary","citation":"Ginszt, T. & Laskowska-Ginszt, A. (2022). Ten years (2012-2021) of spruce bark beetle Ips typographus (L.) activity in the Białowieża Forest District of the Białowieża Primeval Forest. Sylwan, 166(3), 183-193. https://doi.org/10.26202/sylwan.2022005"}]},{"id":"BB0161","name":"Bialowieza Forest spruce bark-beetle outbreak (1950-1957)","country":"Poland","y0":1950,"y1":1957,"lat":52.7,"lon":23.85,"sev":"moderate","conf":"medium","note":"Peer-reviewed historical chronology identifies a documented Ips typographus mass outbreak in Bialowieza Forest during 1950-1957.","sources":[{"title":"Ten years (2012-2021) of spruce bark beetle Ips typographus (L.) activity in the Białowieża Forest District of the Białowieża Primeval Forest","url":"https://doi.org/10.26202/sylwan.2022005","role":"secondary","citation":"Ginszt, T. & Laskowska-Ginszt, A. (2022). Ten years (2012-2021) of spruce bark beetle Ips typographus (L.) activity in the Białowieża Forest District of the Białowieża Primeval Forest. Sylwan, 166(3), 183-193. https://doi.org/10.26202/sylwan.2022005"}]},{"id":"BB0162","name":"Bialowieza Forest spruce bark-beetle outbreak (1960-1966)","country":"Poland","y0":1960,"y1":1966,"lat":52.7,"lon":23.85,"sev":"moderate","conf":"medium","note":"Peer-reviewed historical chronology identifies a documented Ips typographus mass outbreak in Bialowieza Forest during 1960-1966.","sources":[{"title":"Ten years (2012-2021) of spruce bark beetle Ips typographus (L.) activity in the Białowieża Forest District of the Białowieża Primeval Forest","url":"https://doi.org/10.26202/sylwan.2022005","role":"secondary","citation":"Ginszt, T. & Laskowska-Ginszt, A. (2022). Ten years (2012-2021) of spruce bark beetle Ips typographus (L.) activity in the Białowieża Forest District of the Białowieża Primeval Forest. Sylwan, 166(3), 183-193. https://doi.org/10.26202/sylwan.2022005"}]},{"id":"BB0163","name":"Central Slovakia postwar spruce bark-beetle outbreak (1947-1950)","country":"Slovakia","y0":1947,"y1":1950,"lat":48.75,"lon":19.5,"sev":"high","conf":"medium","note":"A Slovak National Forest Centre synthesis identifies the first serious mass spruce-mortality episode in Central Slovakia at the turn of the 1940s and 1950s as an Ips typographus outbreak.","sources":[{"title":"Bark Beetle - the Main Driving Force Behind Spruce Forests Decline","url":"https://www.researchgate.net/publication/293156483_Bark_beetle_-_The_main_driving_force_behind_spruce_forests_decline","role":"secondary","citation":"Hlasny, T., Grodzki, W., Turcani, M., Sitkova, Z. & Konopka, B. (2010). Bark Beetle - the Main Driving Force Behind Spruce Forests Decline. In T. Hlasny & Z. Sitkova (eds.), Spruce Forests Decline in the Beskids, pp. 69-91. Zvolen: National Forest Centre."},{"title":"Meteorological and climatological triggers of notable past and present bark beetle outbreaks in the Czech Republic","url":"https://doi.org/10.5194/cp-18-2155-2022","role":"secondary","citation":"Brazdil, R. et al. (2022). Meteorological and climatological triggers of notable past and present bark beetle outbreaks in the Czech Republic. Climate of the Past 18:2155-2180. https://doi.org/10.5194/cp-18-2155-2022"}]},{"id":"BB0164","name":"Krkonose National Park (Czech side) bark-beetle outbreak (2017-2019 documented phase)","country":"Czechia","y0":2017,"y1":2019,"lat":50.73,"lon":15.6,"sev":"high","conf":"high","note":"Bark-beetle-infested timber harvest rose from ~5,000 m3 (2017) to ~33,000 m3 (2018) to a ~40,000 m3 2019 target; about one-third of the Czech-side park's forest area affected.","sources":[{"title":"Kurovec se siri vys a vys. V ohrozeni jsou uz i Krkonose","url":"https://www.seznamzpravy.cz/clanek/kurovec-se-siri-vys-a-vys-krkonose-letos-vytezi-rekordni-mnozstvi-napadeneho-drivi-80894","role":"primary","citation":"Seznam Zpravy (2019, 16 October). Kurovec se siri vys a vys. V ohrozeni jsou uz i Krkonose."},{"title":"Kurovec dobyva Krkonose. A vitezi","url":"https://novinky.cz/domaci/clanek/kurovec-dobyva-krkonose-a-vitezi-40293080","role":"primary","citation":"Novinky.cz (2019, 14 August). Kurovec dobyva Krkonose. A vitezi."}]},{"id":"BB0165","name":"Krkonose National Park (Czech side) emission-compounded bark-beetle dieback (1984-1994 documented phase)","country":"Czechia","y0":1984,"y1":1994,"lat":50.79,"lon":15.54,"sev":"extreme","conf":"medium","note":"KRNAP's own 50th-anniversary retrospective states spruce stands weakened by long-range sulfur emissions began additionally declining from bark-beetle gradation around 1984; the park was listed among IUCN's 12 most endangered national parks worldwide that same year.","sources":[{"title":"50 let Krkonosskeho narodniho parku","url":"https://www.krnap.cz/media/ag5hgnpy/krnap_kniha_50_let_vnitrek_final_verze-pro-web.pdf","role":"primary","citation":"Sprava Krkonosskeho narodniho parku (2013). 50 let Krkonosskeho narodniho parku."}]},{"id":"BB0166","name":"Luxembourg nationwide spruce bark-beetle outbreak (2018-2019 documented phase)","country":"Luxembourg","y0":2018,"y1":2019,"lat":49.85,"lon":6.05,"sev":"extreme","conf":"high","note":"About 50% of Luxembourg's spruce stands affected; 55,000 m3 already felled in public forests alone; wood price collapsed from 80 to 20 EUR/m3.","sources":[{"title":"Epiceas attaques par le scolyte, une catastrophe economique","url":"https://lequotidien.lu/a-la-une/epiceas-attaques-par-le-scolyte-une-catastrophe-economique/","role":"primary","citation":"Le Quotidien (2019, 11 September). Epiceas attaques par le scolyte, une catastrophe economique."},{"title":"Les forets luxembourgeoises, martyres de la crise climatique","url":"https://paperjam.lu/article/forets-luxembourgeoises-martyr","role":"secondary","citation":"Paperjam (2020). Les forets luxembourgeoises, martyres de la crise climatique."}]},{"id":"BB0167","name":"Styria bark-beetle outbreak (2023)","country":"Austria","y0":2023,"y1":2023,"lat":47.35,"lon":14.75,"sev":"extreme","conf":"high","note":"964,000 Vorratsfestmeter of bark-beetle damage in 2023, +48% versus 2022 and the third-highest value since detailed records began in the late 1980s.","sources":[{"title":"Waldschutzsituation 2023","url":"https://www.forstzeitung.at/markt/2024/06/waldschutzsituation-2023.html","role":"primary","citation":"Forstzeitung (2024, June). Waldschutzsituation 2023."}]},{"id":"BB0168","name":"Ticino canton-wide bark-beetle response program (2017-2020)","country":"Switzerland","y0":2017,"y1":2020,"lat":46.35,"lon":8.85,"sev":"high","conf":"high","note":"Outbreak foci tripled from 2015 to 2016; canton launched a 4-year, CHF 12 million monitoring/felling program (2017-2020, 80% cantonal/federal funded).","sources":[{"title":"Il Ticino dichiara guerra al bostrico","url":"https://www.tio.ch/ticino/attualita/1160790/il-ticino-dichiara-guerra-al-bostrico","role":"primary","citation":"Ticino Online (2017, 24 August). Il Ticino dichiara guerra al bostrico."}]},{"id":"BB0169","name":"St. Gallen bark-beetle outbreak (2018)","country":"Switzerland","y0":2018,"y1":2018,"lat":47.2,"lon":9.3,"sev":"high","conf":"high","note":"About 20,000 m3 of beetle-affected wood in 2018, double the 2017 figure (~10,000 m3).","sources":[{"title":"Der Borkenkaefer ist dieses Jahr unersaettlich: in der Ostschweiz ist er besonders gefraessig","url":"https://www.tagblatt.ch/ostschweiz/der-borkenkaefer-ist-dieses-jahr-unersaettlich-in-der-ostschweiz-ist-er-besonders-gefraessig-ld.1116247","role":"primary","citation":"Tagblatt (2019, 6 May). Der Borkenkaefer ist dieses Jahr unersaettlich: in der Ostschweiz ist er besonders gefraessig."}]},{"id":"BB0170","name":"Thurgau bark-beetle outbreak (2018)","country":"Switzerland","y0":2018,"y1":2018,"lat":47.55,"lon":9.1,"sev":"high","conf":"high","note":"30,000-35,000 m3 of beetle-affected wood in 2018, the highest level since 2006.","sources":[{"title":"Der Borkenkaefer ist dieses Jahr unersaettlich: in der Ostschweiz ist er besonders gefraessig","url":"https://www.tagblatt.ch/ostschweiz/der-borkenkaefer-ist-dieses-jahr-unersaettlich-in-der-ostschweiz-ist-er-besonders-gefraessig-ld.1116247","role":"primary","citation":"Tagblatt (2019, 6 May). Der Borkenkaefer ist dieses Jahr unersaettlich: in der Ostschweiz ist er besonders gefraessig."}]},{"id":"BB0171","name":"Prokletije National Park (Plav, Montenegro) bark-beetle outbreak (2020-2024 documented phase)","country":"Montenegro","y0":2020,"y1":2024,"lat":42.58,"lon":19.9,"sev":"extreme","conf":"high","note":"Damaged timber volume rose from ~6,000 m3 (April 2022) to ~60,000 m3 (Nov 2023) to ~100,000 m3 (by publication, Sept 2024); ~5,000 ha affected; EUR 20 million estimated economic loss.","sources":[{"title":"Potkornjak osusio sume, nazire kraj drame drvecem u nacionalnom parku Prokletije koje napadaju stetocine","url":"https://www.novosti.rs/crna-gora/vesti/1408313/potkornjak-osusio-sume-nazire-kraj-drame-drvecem-nacionalnom-parku-prokletije-koje-napadaju-stetocine","role":"primary","citation":"Novosti (2024, 16 September). Potkornjak osusio sume, nazire kraj drame drvecem u nacionalnom parku Prokletije koje napadaju stetocine."}]},{"id":"BB0172","name":"Zakarpattia (Transcarpathian) Oblast spruce dieback and bark-beetle outbreak (2016 documented phase)","country":"Ukraine","y0":2016,"y1":2016,"lat":48.05,"lon":24.2,"sev":"high","conf":"high","note":"Zakarpattia held about one-third of the 426,000 ha of dying Ukrainian Carpathian spruce reported in 2016 (roughly 142,000 ha); the Rahiv forest economy example shows 2,000 of its 5,776 ha affected.","sources":[{"title":"Smerekova katastrofa. Chomu v Karpatakh masovo vsykhaiut yalynnyky?","url":"https://www.ukrinform.ua/rubric-economy/2027668-smerekova-katastrofa-comu-v-karpatah-masovo-vsihaut-alinniki.html","role":"primary","citation":"Ukrinform (2016, 4 June). Smerekova katastrofa. Chomu v Karpatakh masovo vsykhaiut yalynnyky?"}]},{"id":"BB0173","name":"Alpes du Nord and Ain spruce bark-beetle escalation (2019-2021 documented phase)","country":"France","y0":2019,"y1":2021,"lat":45.6,"lon":6.3,"sev":"high","conf":"high","note":"Official French forest-health service reports damage rising sharply from 2019, stabilizing at a high level in 2020, then rising a further ~50% on average across departments in 2021; damage moving upslope past 1,000 m.","sources":[{"title":"Dommages lies aux scolytes des epiceas sur les Alpes du Nord et l'Ain -- Point en janvier 2022","url":"https://draaf.auvergne-rhone-alpes.agriculture.gouv.fr/IMG/pdf/Bilan_scolytes_alpes_du_nord_janv_2022_-Vdef_cle8a6feb.pdf","role":"primary","citation":"Departement de la Sante des Forets (2022, January). Dommages lies aux scolytes des epiceas sur les Alpes du Nord et l'Ain -- Point en janvier 2022."}]},{"id":"BB0174","name":"Swiss Jura bark-beetle increase (2022–2024)","country":"Switzerland","y0":2022,"y1":2024,"lat":46.8,"lon":6.5,"sev":"moderate","conf":"medium","note":"Swiss Ips typographus infested spruce 713,210 m3 in 2023 (+10%) and 622,850 m3 in 2024; Neuchâtel +138% in 2023; Vaud Jura reported a very large amount of beetle wood in 2024, worsening over the last three years.","sources":[{"title":"Buchdrucker: gebietsweise Zunahme (Waldschutz aktuell 1/2024)","url":"https://waldschutz.wsl.ch/fileadmin/user_upload/News/2024/04_Borkenkaefer/WSaktuell_D_1_2024_opt.pdf","role":"primary","citation":"Stroheker S, Blaser S, Queloz V (2024). Buchdrucker: gebietsweise Zunahme. Waldschutz aktuell 1/2024. WSL, Birmensdorf."},{"title":"Bark beetle situation has eased somewhat in Switzerland","url":"https://www.bluewin.ch/en/news/international/bark-beetle-situation-has-eased-somewhat-in-switzerland-2614953.html","role":"secondary","citation":"Keystone-SDA (2025, 20 March). Bark beetle situation has eased somewhat in Switzerland. blue News."}]},{"id":"BB0175","name":"Saarland western Hunsrück bark-beetle increase (2023–2024)","country":"Germany","y0":2023,"y1":2024,"lat":49.55,"lon":6.9,"sev":"moderate","conf":"high","note":"Saarland bark-beetle wood (spruce and fir) rose to about 0.18 million Fm in 2024 from 0.13 million Fm in 2023, against the general decline.","sources":[{"title":"Waldschutzsituation 2024/2025 in Rheinland-Pfalz und Saarland","url":"https://www.waldwissen.net/de/waldwirtschaft/schadensmanagement/waldschutzsituation-rp-sl-24-25","role":"primary","citation":"Delb H et al. (2025). Waldschutzsituation 2024/2025 in Rheinland-Pfalz und Saarland. AFZ-DerWald 09/2025: 26-29."}]},{"id":"BB0176","name":"Thuringia (Thüringer Wald and East Thuringia) bark-beetle calamity (2018–2024)","country":"Germany","y0":2018,"y1":2024,"lat":50.55,"lon":11.15,"sev":"extreme","conf":"high","note":"ThüringenForst: about 6 million Fm of spruce damage timber expected for 2023 and about 3.5 million Fm for 2024; focus in the Thüringer Wald and East Thuringia.","sources":[{"title":"Borkenkäfer aktuell (Alle Borkenkäfer-News für Thüringen)","url":"https://www.thueringenforst.de/startseite/borkenkaefer-aktuell","role":"primary","citation":"ThüringenForst (n.d., accessed 2026-10-28). Borkenkäfer aktuell. thueringenforst.de."},{"title":"Waldzustandsbericht Thüringen 2023 (Thuringia Forest Condition Report 2023)","url":"https://digitales-infrastruktur.thueringen.de/fileadmin/Forst_und_Jagd_Fischerei/Forstwirtschaft/2023_Waldzustandsbericht.pdf","role":"primary","citation":"Thüringer Ministerium für Infrastruktur und Landwirtschaft (2023). Waldzustandsbericht Thüringen 2023."}]},{"id":"BB0177","name":"South-eastern Finland (South Karelia and Kymenlaakso) bark-beetle increase (2022–2023)","country":"Finland","y0":2022,"y1":2023,"lat":61.0,"lon":27.4,"sev":"moderate","conf":"medium","note":"Luke (National Forest Inventory): national bark-beetle damage area 64,300 ha in 2023 against typically about 20,000 ha; the MetuKKa project (2022–2024) documents rising damage in South Karelia and Kymenlaakso (for example Lappeenranta and Savitaipale).","sources":[{"title":"Kirjanpainajan aiheuttama kuusikuolleisuus kolminkertaistunut – Metsätuhot vuonna 2023 -raportti on julkaistu","url":"https://luke.fi/fi/seurannat/metsatuhojen-seuranta/kirjanpainajan-aiheuttama-kuusikuolleisuus-kolminkertaistunut-metsatuhot-vuonna-2023-raportti-on-julkaistu","role":"primary","citation":"Luonnonvarakeskus (2024, 22 April). Kirjanpainajan aiheuttama kuusikuolleisuus kolminkertaistunut – Metsätuhot vuonna 2023 -raportti on julkaistu. luke.fi."},{"title":"Metsätuhot kuriin Kaakkois-Suomessa – loppuraportti","url":"https://www.metsakeskus.fi/sites/default/files/document/metukka-hanke-2022-2024-loppuraportti.pdf","role":"primary","citation":"Kokkonen J, Ylioja T (n.d., written at the end of 2024). Metsätuhot kuriin Kaakkois-Suomessa – loppuraportti. Suomen metsäkeskus and Luonnonvarakeskus."}]},{"id":"BB0178","name":"Mecklenburg-Vorpommern spruce bark-beetle outbreak (1992–1997), state reports via DFDE","country":"Germany","y0":1992,"y1":1997,"lat":53.797,"lon":12.518,"sev":"high","conf":"medium","note":"Infested wood peaked at 222,360 m3 in June 1994 - May 1995.","sources":[{"title":"Waldschutzinformation 3/14: Stehendbefallserfassung der Rindenbrüter Januar/Februar 2014","url":"https://www.wald-mv.de/serviceassistent/download?id=1647582","role":"secondary","citation":"Landesforst Mecklenburg-Vorpommern (2014). Waldschutzinformation 3/14: Stehendbefallserfassung der Rindenbrüter Januar/Februar 2014. Waldschutzinformation 3/14."},{"title":"Bemmann, M. and e. al. (2000). “Waldschutzsituation 1999/2000 in Mecklenburg-Vorpommern.” Allgemeine Forstzeitschrift 55(7): 350-351.","url":"https://dfde.efi.int/db/dfde_app.php","role":"secondary","citation":"Bemmann, M. and e. al. (2000). “Waldschutzsituation 1999/2000 in Mecklenburg-Vorpommern.” Allgemeine Forstzeitschrift 55(7): 350-351. (accessed through the DFDE, EFI; original not seen)"},{"title":"Majunke, C., C. Walter, et al. (1996). “Waldschutzsituation 1995/96 in Brandenburg, Mecklenburg-Vorpommern und Berlin.” Allgemeine Forstzeitschrift 51(7): 378-382.","url":"https://dfde.efi.int/db/dfde_app.php","role":"secondary","citation":"Majunke, C., C. Walter, et al. (1996). “Waldschutzsituation 1995/96 in Brandenburg, Mecklenburg-Vorpommern und Berlin.” Allgemeine Forstzeitschrift 51(7): 378-382. (accessed through the DFDE, EFI; original not seen)"}]},{"id":"BB0179","name":"Hessen spruce bark-beetle outbreak (1991–1995), Waldschutz reports via DFDE","country":"Germany","y0":1991,"y1":1995,"lat":50.526,"lon":9.019,"sev":"high","conf":"medium","note":"Sanitation fellings 1994 sum to 737,551 m3 in the three regional districts (state and private forest) and about 1,300-1,700 ha of 'threatening' damage area 1991-1993.","sources":[{"title":"Dimitri, L. (1996). “Waldschutzsituation 1995/96 in Hessen.” Allgemeine Forstzeitschrift 51(7): 354-358.","url":"https://dfde.efi.int/db/dfde_app.php","role":"secondary","citation":"Dimitri, L. (1996). “Waldschutzsituation 1995/96 in Hessen.” Allgemeine Forstzeitschrift 51(7): 354-358. (accessed through the DFDE, EFI; original not seen)"}]},{"id":"BB0180","name":"Thuringia (Mühlhausen) bark-beetle damage (1994–1995), via DFDE","country":"Germany","y0":1994,"y1":1995,"lat":51.21,"lon":10.46,"sev":"moderate","conf":"medium","note":"60,033 m3 in 1994 and 55,719 m3 in 1995 attributed to the first beetle generation.","sources":[{"title":"Baier, U. (1996). “Waldschutzsituation 1995/96 in Thüringen.” Allgemeine Forstzeitschrift 51(7): 365-368.","url":"https://dfde.efi.int/db/dfde_app.php","role":"secondary","citation":"Baier, U. (1996). “Waldschutzsituation 1995/96 in Thüringen.” Allgemeine Forstzeitschrift 51(7): 365-368. (accessed through the DFDE, EFI; original not seen)"}]},{"id":"BB0181","name":"Western Germany (Federal Republic before 1990 boundaries) spruce bark-beetle epidemic (1990–1998), via DFDE","country":"Germany","y0":1990,"y1":1998,"lat":51.2,"lon":10.4,"sev":"high","conf":"medium","note":"Attacked volume 5,901,250 m3 on 192,760 ha, threatened area 323,300 ha, 1990-1998 (Grégoire and Evans 2004).","sources":[{"title":"J-C.Gregoire & H.F. Evans, 2004, 'Damage and control of BAWBILT organisms, an overview', F.lieutier et al., (eds), Bark and Wood Boring Insects in Living Trees in Europe, A synthesis, 2004, 19-37","url":"https://dfde.efi.int/db/dfde_app.php","role":"secondary","citation":"J-C.Gregoire & H.F. Evans, 2004, 'Damage and control of BAWBILT organisms, an overview', F.lieutier et al., (eds), Bark and Wood Boring Insects in Living Trees in Europe, A synthesis, 2004, 19-37 (accessed through the DFDE, EFI; original not seen)"}]},{"id":"BB0182","name":"France national Ips typographus damage (1990–2001), via DFDE","country":"France","y0":1990,"y1":2001,"lat":46.6,"lon":2.5,"sev":"moderate","conf":"medium","note":"1,500,000 m3 attacked 1990-2001 (Grégoire and Evans 2004); INRA 2003 notes insect damage increased by 15-35% after the storms.","sources":[{"title":"J-C.Gregoire & H.F. Evans, 2004, 'Damage and control of BAWBILT organisms, an overview', F.lieutier et al., (eds), Bark and Wood Boring Insects in Living Trees in Europe, A synthesis, 2004, 19-37","url":"https://dfde.efi.int/db/dfde_app.php","role":"secondary","citation":"J-C.Gregoire & H.F. Evans, 2004, 'Damage and control of BAWBILT organisms, an overview', F.lieutier et al., (eds), Bark and Wood Boring Insects in Living Trees in Europe, A synthesis, 2004, 19-37 (accessed through the DFDE, EFI; original not seen)"},{"title":"Le dossiers de l'environnement de l'INRA 2003, Forets et tempete, Mission environnement-societe, Institut National de la Recherce Agronomique","url":"https://dfde.efi.int/db/dfde_app.php","role":"secondary","citation":"Le dossiers de l'environnement de l'INRA 2003, Forets et tempete, Mission environnement-societe, Institut National de la Recherce Agronomique (accessed through the DFDE, EFI; original not seen)"}]},{"id":"BB0183","name":"Rhône-Alpes (Isère valley, Savoie, Haute-Savoie) bark-beetle damage (1984–1985), via DFDE","country":"France","y0":1984,"y1":1985,"lat":45.6,"lon":6.2,"sev":"high","conf":"medium","note":"235,000 m3 damaged and salvaged; 90% of the damage in the whole Rhône-Alpes region was in these districts.","sources":[{"title":"Le dossiers de l'environnement de l'INRA 2003, Forets et tempete, Mission environnement-societe, Institut National de la Recherce Agronomique","url":"https://dfde.efi.int/db/dfde_app.php","role":"secondary","citation":"Le dossiers de l'environnement de l'INRA 2003, Forets et tempete, Mission environnement-societe, Institut National de la Recherce Agronomique (accessed through the DFDE, EFI; original not seen)"}]},{"id":"BB0184","name":"Hungary national Ips typographus damage (1990–1997), via DFDE","country":"Hungary","y0":1990,"y1":1997,"lat":47.2,"lon":19.5,"sev":"moderate","conf":"medium","note":"420,000 m3 attacked 1990-1997, 24,000 ha threatened (Grégoire and Evans 2004).","sources":[{"title":"J-C.Gregoire & H.F. Evans, 2004, 'Damage and control of BAWBILT organisms, an overview', F.lieutier et al., (eds), Bark and Wood Boring Insects in Living Trees in Europe, A synthesis, 2004, 19-37","url":"https://dfde.efi.int/db/dfde_app.php","role":"secondary","citation":"J-C.Gregoire & H.F. Evans, 2004, 'Damage and control of BAWBILT organisms, an overview', F.lieutier et al., (eds), Bark and Wood Boring Insects in Living Trees in Europe, A synthesis, 2004, 19-37 (accessed through the DFDE, EFI; original not seen)"}]},{"id":"BB0185","name":"Valle d'Aosta spruce bark-beetle damage after Lymantria monacha and storm Vivian (1991–1998), via DFDE","country":"Italy","y0":1991,"y1":1998,"lat":45.74,"lon":7.4,"sev":"moderate","conf":"medium","note":"Groups of 10-20 spruce first, later groups of hundreds of trees killed (Negro et al. 1999).","sources":[{"title":"Valle d'Aosta regional chronicle (French/Italian), entry for 1990","url":"https://gestionewww.regione.vda.it/allegato.aspx?pk=77520","role":"secondary","citation":"Regione Autonoma Valle d'Aosta (2020). Valle d'Aosta regional chronicle (French/Italian), entry for 1990. regional government document."},{"title":"Negro, M., L. Dovigo, et al. (1999). Survey methods on insect forest pests in the Aosta Valley. Methodology of Forest Insect and Disease Survey in Central Europe. Second Workshop of the IUFRO WP 7.03.10, Sion-Châteauneuf, Switzerland, Swiss Federal Institu","url":"https://dfde.efi.int/db/dfde_app.php","role":"secondary","citation":"Negro, M., L. Dovigo, et al. (1999). Survey methods on insect forest pests in the Aosta Valley. Methodology of Forest Insect and Disease Survey in Central Europe. Second Workshop of the IUFRO WP 7.03.10, Sion-Châteauneuf, Switzerland, Swiss Federal Institu (accessed through the DFDE, EFI; original not seen)"}]},{"id":"BB0186","name":"Lithuania Ips typographus attack (1998–2001), via DFDE","country":"Lithuania","y0":1998,"y1":2001,"lat":55.3,"lon":23.9,"sev":"moderate","conf":"medium","note":"27,462 ha attacked and 221,155 ha threatened, 1998-2001 (Grégoire and Evans 2004); no volume.","sources":[{"title":"J-C.Gregoire & H.F. Evans, 2004, 'Damage and control of BAWBILT organisms, an overview', F.lieutier et al., (eds), Bark and Wood Boring Insects in Living Trees in Europe, A synthesis, 2004, 19-37","url":"https://dfde.efi.int/db/dfde_app.php","role":"secondary","citation":"J-C.Gregoire & H.F. Evans, 2004, 'Damage and control of BAWBILT organisms, an overview', F.lieutier et al., (eds), Bark and Wood Boring Insects in Living Trees in Europe, A synthesis, 2004, 19-37 (accessed through the DFDE, EFI; original not seen)"}]},{"id":"BB0187","name":"Lithuania state forests spruce bark-beetle damage (2015–2019), via DFDE","country":"Lithuania","y0":2015,"y1":2019,"lat":55.3,"lon":23.9,"sev":"moderate","conf":"medium","note":"State forest service yearbook: 34,224 m3 (2015), 27,321 (2016), 19,976 (2017), 24,772 (2018), 55,367 m3 on 1,315 ha (2019) in state forests only.","sources":[{"title":"Vienas svarbiausių miškininkų rūpesčių – eglynų sveikata","url":"https://alkas.lt/2016/04/02/vienas-svarbiausiu-miskininku-rupesciu-eglynu-sveikata","role":"secondary","citation":"alkas.lt (quoting the State Forest Service 2015 sanitary review) (2016). Vienas svarbiausių miškininkų rūpesčių – eglynų sveikata. news article, 2 April 2016."},{"title":"Lithuanian statistical yearbook of forestry. ","url":"https://dfde.efi.int/db/dfde_app.php","role":"secondary","citation":"Lithuanian statistical yearbook of forestry.  (accessed through the DFDE, EFI; original not seen)"}]},{"id":"BB0188","name":"Poland national Ips typographus attack (1990–1999), via DFDE","country":"Poland","y0":1990,"y1":1999,"lat":52.0,"lon":19.4,"sev":"high","conf":"medium","note":"6,270,000 m3 attacked 1990-1999 and 230,000 ha threatened (Grégoire and Evans 2004).","sources":[{"title":"J-C.Gregoire & H.F. Evans, 2004, 'Damage and control of BAWBILT organisms, an overview', F.lieutier et al., (eds), Bark and Wood Boring Insects in Living Trees in Europe, A synthesis, 2004, 19-37","url":"https://dfde.efi.int/db/dfde_app.php","role":"secondary","citation":"J-C.Gregoire & H.F. Evans, 2004, 'Damage and control of BAWBILT organisms, an overview', F.lieutier et al., (eds), Bark and Wood Boring Insects in Living Trees in Europe, A synthesis, 2004, 19-37 (accessed through the DFDE, EFI; original not seen)"}]},{"id":"BB0189","name":"Romania national Ips typographus attack (1990–2000), via DFDE","country":"Romania","y0":1990,"y1":2000,"lat":45.9,"lon":24.9,"sev":"moderate","conf":"medium","note":"811,195 m3 on 259,883 ha attacked 1990-2000; 600,000 ha threatened (Grégoire and Evans 2004).","sources":[{"title":"J-C.Gregoire & H.F. Evans, 2004, 'Damage and control of BAWBILT organisms, an overview', F.lieutier et al., (eds), Bark and Wood Boring Insects in Living Trees in Europe, A synthesis, 2004, 19-37","url":"https://dfde.efi.int/db/dfde_app.php","role":"secondary","citation":"J-C.Gregoire & H.F. Evans, 2004, 'Damage and control of BAWBILT organisms, an overview', F.lieutier et al., (eds), Bark and Wood Boring Insects in Living Trees in Europe, A synthesis, 2004, 19-37 (accessed through the DFDE, EFI; original not seen)"}]},{"id":"BB0190","name":"Slovakia national Ips typographus attack (1990–1999), via DFDE","country":"Slovakia","y0":1990,"y1":1999,"lat":48.7,"lon":19.5,"sev":"high","conf":"medium","note":"4,490,090 m3 attacked 1990-1999 (Grégoire and Evans 2004); 618,700 m3 processed and 700,000 m3 infested in 1996; 537,613 m3 infested and 477,678 m3 processed in 1998.","sources":[{"title":"J-C.Gregoire & H.F. Evans, 2004, 'Damage and control of BAWBILT organisms, an overview', F.lieutier et al., (eds), Bark and Wood Boring Insects in Living Trees in Europe, A synthesis, 2004, 19-37","url":"https://dfde.efi.int/db/dfde_app.php","role":"secondary","citation":"J-C.Gregoire & H.F. Evans, 2004, 'Damage and control of BAWBILT organisms, an overview', F.lieutier et al., (eds), Bark and Wood Boring Insects in Living Trees in Europe, A synthesis, 2004, 19-37 (accessed through the DFDE, EFI; original not seen)"},{"title":"Novotný, J. and J. Varínsky (1997). “The most significant injurious factors in the forests of Slovakia - the present situation.” Slovak Forest 53(9): 14-17.","url":"https://dfde.efi.int/db/dfde_app.php","role":"secondary","citation":"Novotný, J. and J. Varínsky (1997). “The most significant injurious factors in the forests of Slovakia - the present situation.” Slovak Forest 53(9): 14-17. (accessed through the DFDE, EFI; original not seen)"},{"title":"Zúbrik, M., J. Novotný, et al. (1999). The dominant pest and injurious agents in Slovakia. Methodology of Forest Insect and Disease Survey in Central Europe. Second Workshop of the IUFRO WP 7.03.10, Sion-Châteauneuf, Switzerland, Swiss Federal Institute fo","url":"https://dfde.efi.int/db/dfde_app.php","role":"secondary","citation":"Zúbrik, M., J. Novotný, et al. (1999). The dominant pest and injurious agents in Slovakia. Methodology of Forest Insect and Disease Survey in Central Europe. Second Workshop of the IUFRO WP 7.03.10, Sion-Châteauneuf, Switzerland, Swiss Federal Institute fo (accessed through the DFDE, EFI; original not seen)"}]},{"id":"BB0191","name":"Kysuce (north-western Slovakia) spruce bark-beetle damage (1999–2004), via DFDE","country":"Slovakia","y0":1999,"y1":2004,"lat":49.44,"lon":18.79,"sev":"moderate","conf":"medium","note":"427,000 m3 damaged over the preceding five years (Novotny 2004).","sources":[{"title":"Influence of selected factors on bark beetle outbreak dynamics in the Western Carpathians","url":"https://reference-global.com/article/10.1515/forj-2015-0023","role":"secondary","citation":"Vakula J., Zúbrik M., Galko J., Gubka A., Kunca A., Nikolov C., Bošeľa M. (2015). Influence of selected factors on bark beetle outbreak dynamics in the Western Carpathians. Central European Forestry Journal 61(3): 149-156."},{"title":"J. Novotny, 2004, 'Community involvement in forest protection in Slovakia', Unasylva 217, vol. 55, 2004","url":"https://dfde.efi.int/db/dfde_app.php","role":"secondary","citation":"J. Novotny, 2004, 'Community involvement in forest protection in Slovakia', Unasylva 217, vol. 55, 2004 (accessed through the DFDE, EFI; original not seen)"}]},{"id":"BB0192","name":"Šumava National Park spruce bark-beetle outbreak (1993–1998)","country":"Czechia","y0":1993,"y1":1998,"lat":49.05,"lon":13.5,"sev":"high","conf":"high","note":"Park administration: sudden onset of outbreaks in 1994, 1996 and 1997, 187,000 m3 processed in 1996 and about 165,000 m3 of new dead trees in 1997; DFDE yearly series (Zatloukal 1998): 50,996 m3 (1993), 45,484 (1994), 57,052 (1995), 181,683","sources":[{"title":"Vývoj situace v lesních porostech s ohledem na napadení lýkožroutem smrkovým v lesích ve správě Správy NP a CHKO Šumava","url":"https://mzp.gov.cz/sites/default/files/tiskove-zpravy/tz_nps090825kurovec_sumava.pdf","role":"primary","citation":"Správa NP a CHKO Šumava (2009). Vývoj situace v lesních porostech s ohledem na napadení lýkožroutem smrkovým v lesích ve správě Správy NP a CHKO Šumava. press release of 25 August 2009."},{"title":"Liška, J., P. Kapitola, et al. (1999). Occurrence of Forest Insect Pests and Fungal Diseases in Czechia in 1998. Methodology of Forest Insect and Disease Survey in Central Europe. Second Workshop of the IUFRO WP 7.03.10, Sion-Châteauneuf, Switzerland, Swis","url":"https://dfde.efi.int/db/dfde_app.php","role":"secondary","citation":"Liška, J., P. Kapitola, et al. (1999). Occurrence of Forest Insect Pests and Fungal Diseases in Czechia in 1998. Methodology of Forest Insect and Disease Survey in Central Europe. Second Workshop of the IUFRO WP 7.03.10, Sion-Châteauneuf, Switzerland, Swis (accessed through the DFDE, EFI; original not seen)"},{"title":"Zatloukal, V. (1998). “Historické a soucasné príciny kurovcové kalamity v Národním parku Šumava. [Historical and current factors of the bark beetle calamity in the Šumava National Park].” Silva Gabreta 2: 327-357.","url":"https://dfde.efi.int/db/dfde_app.php","role":"secondary","citation":"Zatloukal, V. (1998). “Historické a soucasné príciny kurovcové kalamity v Národním parku Šumava. [Historical and current factors of the bark beetle calamity in the Šumava National Park].” Silva Gabreta 2: 327-357. (accessed through the DFDE, EFI; original not seen)"}]},{"id":"BB0193","name":"Denmark Ips typographus attacks after the storms of 1981, via DFDE","country":"Denmark","y0":1981,"y1":1981,"lat":56.0,"lon":9.5,"sev":"moderate","conf":"medium","note":"Severe attacks following storms in 1981 (Franklin et al. 1999).","sources":[{"title":"Franklin, A., L. Raty, et al. (1999). Monitoring Ips typographus populations using logging and trapping data: two complementary approaches. Methodology of Forest Insect and Disease Survey in Central Europe. Second Workshop of the IUFRO WP 7.03.10, Sion-Châ","url":"https://dfde.efi.int/db/dfde_app.php","role":"secondary","citation":"Franklin, A., L. Raty, et al. (1999). Monitoring Ips typographus populations using logging and trapping data: two complementary approaches. Methodology of Forest Insect and Disease Survey in Central Europe. Second Workshop of the IUFRO WP 7.03.10, Sion-Châ (accessed through the DFDE, EFI; original not seen)"}]},{"id":"BB0194","name":"Austria national bark-beetle damage (1980–1991), Donaubauer et al. 1996 via DFDE","country":"Austria","y0":1980,"y1":1991,"lat":47.3,"lon":13.5,"sev":"high","conf":"medium","note":"Yearly volume 270,000-460,000 m3 (peak 460,000 m3 in 1984) before the jump to 1.5-1.9 million m3 in 1992-1995.","sources":[{"title":"Donaubauer, E.; Krehan, H.; Tomiczek, Ch. 1996. Forstschadenssituation in Österreich. AFZ/Der Wald, München, 51(7): 386-388.","url":"https://dfde.efi.int/db/dfde_app.php","role":"secondary","citation":"Donaubauer, E.; Krehan, H.; Tomiczek, Ch. 1996. Forstschadenssituation in Österreich. AFZ/Der Wald, München, 51(7): 386-388. (accessed through the DFDE, EFI; original not seen)"}]},{"id":"BB0195","name":"Czech Republic national Ips typographus damage (1997–1998), via DFDE","country":"Czechia","y0":1997,"y1":1998,"lat":49.8,"lon":15.5,"sev":"high","conf":"medium","note":"390,000 m3 in 1997 and 330,000 m3 in 1998, mainly Ips typographus (Liška et al. 1999; the 1998 figure covers about 75% of the forest area).","sources":[{"title":"Liška, J., P. Kapitola, et al. (1999). Occurrence of Forest Insect Pests and Fungal Diseases in Czechia in 1998. Methodology of Forest Insect and Disease Survey in Central Europe. Second Workshop of the IUFRO WP 7.03.10, Sion-Châteauneuf, Switzerland, Swis","url":"https://dfde.efi.int/db/dfde_app.php","role":"secondary","citation":"Liška, J., P. Kapitola, et al. (1999). Occurrence of Forest Insect Pests and Fungal Diseases in Czechia in 1998. Methodology of Forest Insect and Disease Survey in Central Europe. Second Workshop of the IUFRO WP 7.03.10, Sion-Châteauneuf, Switzerland, Swis (accessed through the DFDE, EFI; original not seen)"}]},{"id":"BB0196","name":"Spiš region (Poprad, Spišská Nová Ves, Gelnica, Kežmarok; Slovakia) Ips typographus damage (1996), via DFDE","country":"Slovakia","y0":1996,"y1":1996,"lat":49.06,"lon":20.3,"sev":"high","conf":"medium","note":"Summed damage of 276,100 m3 in four districts in 1996 (94,200; 75,200; 60,000; 46,700 m3), probably I. typographus and I. amitinus (Novotný and Varínsky 1997).","sources":[{"title":"Novotný, J. and J. Varínsky (1997). “The most significant injurious factors in the forests of Slovakia - the present situation.” Slovak Forest 53(9): 14-17.","url":"https://dfde.efi.int/db/dfde_app.php","role":"secondary","citation":"Novotný, J. and J. Varínsky (1997). “The most significant injurious factors in the forests of Slovakia - the present situation.” Slovak Forest 53(9): 14-17. (accessed through the DFDE, EFI; original not seen)"}]},{"id":"BB0197","name":"Northwestern Slovakia (Žilina, Púchov, Považská Bystrica, Čadca, Bytča) Ips typographus damage (1996), via DFDE","country":"Slovakia","y0":1996,"y1":1996,"lat":49.22,"lon":18.74,"sev":"moderate","conf":"medium","note":"Summed damage of 155,200 m3 in five districts in 1996, probably I. typographus and I. amitinus (Novotný and Varínsky 1997).","sources":[{"title":"Novotný, J. and J. Varínsky (1997). “The most significant injurious factors in the forests of Slovakia - the present situation.” Slovak Forest 53(9): 14-17.","url":"https://dfde.efi.int/db/dfde_app.php","role":"secondary","citation":"Novotný, J. and J. Varínsky (1997). “The most significant injurious factors in the forests of Slovakia - the present situation.” Slovak Forest 53(9): 14-17. (accessed through the DFDE, EFI; original not seen)"}]},{"id":"BB0198","name":"Kočevje (Slovenia) bark-beetle outbreak after the 2003 drought (2003–2007), via DFDE","country":"Slovenia","y0":2003,"y1":2007,"lat":45.64,"lon":14.86,"sev":"high","conf":"medium","note":"820,000 m3 over 2003-2007, 227,000 m3 in 2005 alone (Poljanec et al. 2014).","sources":[{"title":"Poljanec A., Ščap Š., Bončina A. 2014. Količina, struktura in razporeditev sanitarnega poseka v Sloveniji v obdobju 1995–2012 - Volume, Structure and Distribution of Salvage Logging in Slovenia in the Period\n1995-2012","url":"https://dfde.efi.int/db/dfde_app.php","role":"secondary","citation":"Poljanec A., Ščap Š., Bončina A. 2014. Količina, struktura in razporeditev sanitarnega poseka v Sloveniji v obdobju 1995–2012 - Volume, Structure and Distribution of Salvage Logging in Slovenia in the Period\n1995-2012 (accessed through the DFDE, EFI; original not seen)"}]},{"id":"BB0199","name":"Kranj (Slovenia) bark-beetle outbreak after snow and ice breaks (1998), via DFDE","country":"Slovenia","y0":1998,"y1":1998,"lat":46.24,"lon":14.36,"sev":"moderate","conf":"medium","note":"52,004 m3 in 1998 after the 1995-1997 snow and ice breaks (Poljanec et al. 2014).","sources":[{"title":"Poljanec A., Ščap Š., Bončina A. 2014. Količina, struktura in razporeditev sanitarnega poseka v Sloveniji v obdobju 1995–2012 - Volume, Structure and Distribution of Salvage Logging in Slovenia in the Period\n1995-2012","url":"https://dfde.efi.int/db/dfde_app.php","role":"secondary","citation":"Poljanec A., Ščap Š., Bončina A. 2014. Količina, struktura in razporeditev sanitarnega poseka v Sloveniji v obdobju 1995–2012 - Volume, Structure and Distribution of Salvage Logging in Slovenia in the Period\n1995-2012 (accessed through the DFDE, EFI; original not seen)"}]},{"id":"BB0200","name":"Burgenland-Süd (Güssing) bark-beetle damage (2002–2008), regional series via DFDE","country":"Austria","y0":2002,"y1":2008,"lat":47.07,"lon":16.32,"sev":"moderate","conf":"medium","note":"Reported damage in Burgenland-Süd (Güssing) of at least 50,000 m3 per year in 2002-2008, peak 167,530 m3 in 2005 (Thom, Seidl et al. 2013).","sources":[{"title":"Thom, D, Seidl, R, Steyrer, G, Krehan, H, Formayer, H (2013): Slow and fast drivers of the natural disturbance regime in Central European forest ecosystems. Forest Ecology and Management, 307: 293-302.","url":"https://dfde.efi.int/db/dfde_app.php","role":"secondary","citation":"Thom, D, Seidl, R, Steyrer, G, Krehan, H, Formayer, H (2013): Slow and fast drivers of the natural disturbance regime in Central European forest ecosystems. Forest Ecology and Management, 307: 293-302. (accessed through the DFDE, EFI; original not seen)"}]},{"id":"BB0201","name":"Flachgau (Salzburg) bark-beetle damage (2004), regional series via DFDE","country":"Austria","y0":2004,"y1":2004,"lat":47.8,"lon":13.04,"sev":"moderate","conf":"medium","note":"Reported damage in Flachgau (Salzburg) of at least 50,000 m3 per year in 2004-2004, peak 66,075 m3 in 2004 (Thom, Seidl et al. 2013).","sources":[{"title":"Thom, D, Seidl, R, Steyrer, G, Krehan, H, Formayer, H (2013): Slow and fast drivers of the natural disturbance regime in Central European forest ecosystems. Forest Ecology and Management, 307: 293-302.","url":"https://dfde.efi.int/db/dfde_app.php","role":"secondary","citation":"Thom, D, Seidl, R, Steyrer, G, Krehan, H, Formayer, H (2013): Slow and fast drivers of the natural disturbance regime in Central European forest ecosystems. Forest Ecology and Management, 307: 293-302. (accessed through the DFDE, EFI; original not seen)"}]},{"id":"BB0202","name":"Flachgau (Salzburg) bark-beetle damage (2007–2010), regional series via DFDE","country":"Austria","y0":2007,"y1":2010,"lat":47.8,"lon":13.04,"sev":"moderate","conf":"medium","note":"Reported damage in Flachgau (Salzburg) of at least 50,000 m3 per year in 2007-2010, peak 164,900 m3 in 2010 (Thom, Seidl et al. 2013).","sources":[{"title":"Thom, D, Seidl, R, Steyrer, G, Krehan, H, Formayer, H (2013): Slow and fast drivers of the natural disturbance regime in Central European forest ecosystems. Forest Ecology and Management, 307: 293-302.","url":"https://dfde.efi.int/db/dfde_app.php","role":"secondary","citation":"Thom, D, Seidl, R, Steyrer, G, Krehan, H, Formayer, H (2013): Slow and fast drivers of the natural disturbance regime in Central European forest ecosystems. Forest Ecology and Management, 307: 293-302. (accessed through the DFDE, EFI; original not seen)"}]},{"id":"BB0203","name":"Lungau (Tamsweg) bark-beetle damage (2004–2005), regional series via DFDE","country":"Austria","y0":2004,"y1":2005,"lat":47.12,"lon":13.81,"sev":"moderate","conf":"medium","note":"Reported damage in Lungau (Tamsweg) of at least 50,000 m3 per year in 2004-2005, peak 104,590 m3 in 2005 (Thom, Seidl et al. 2013).","sources":[{"title":"Thom, D, Seidl, R, Steyrer, G, Krehan, H, Formayer, H (2013): Slow and fast drivers of the natural disturbance regime in Central European forest ecosystems. Forest Ecology and Management, 307: 293-302.","url":"https://dfde.efi.int/db/dfde_app.php","role":"secondary","citation":"Thom, D, Seidl, R, Steyrer, G, Krehan, H, Formayer, H (2013): Slow and fast drivers of the natural disturbance regime in Central European forest ecosystems. Forest Ecology and Management, 307: 293-302. (accessed through the DFDE, EFI; original not seen)"}]},{"id":"BB0204","name":"Pinzgau (Zell am See) bark-beetle damage (2004–2007), regional series via DFDE","country":"Austria","y0":2004,"y1":2007,"lat":47.32,"lon":12.8,"sev":"high","conf":"medium","note":"Reported damage in Pinzgau (Zell am See) of at least 50,000 m3 per year in 2004-2007, peak 324,500 m3 in 2005 (Thom, Seidl et al. 2013).","sources":[{"title":"Thom, D, Seidl, R, Steyrer, G, Krehan, H, Formayer, H (2013): Slow and fast drivers of the natural disturbance regime in Central European forest ecosystems. Forest Ecology and Management, 307: 293-302.","url":"https://dfde.efi.int/db/dfde_app.php","role":"secondary","citation":"Thom, D, Seidl, R, Steyrer, G, Krehan, H, Formayer, H (2013): Slow and fast drivers of the natural disturbance regime in Central European forest ecosystems. Forest Ecology and Management, 307: 293-302. (accessed through the DFDE, EFI; original not seen)"}]},{"id":"BB0205","name":"Baden-Württemberg spruce bark-beetle outbreak (2013–2014), state reports via DFDE","country":"Germany","y0":2013,"y1":2014,"lat":48.647,"lon":9.058,"sev":"high","conf":"medium","note":"Reported Ips typographus damage in Baden-Württemberg of at least 50,000 m3 per year in 2013-2014, peak 242,000 m3 in 2013 (Waldschutz reports, AFZ Der Wald).","sources":[{"title":"AFZ Der Wald (2014). \"Waldschutzsituation: Tabellarischer Überblick\". AFZ Der Wald - Allgemeine Forst Zeitschrift für Waldwirtschaft und Umweltvorsorge, 69(7). BLV Verlagsgesellschaft, Berlin.","url":"https://dfde.efi.int/db/dfde_app.php","role":"secondary","citation":"AFZ Der Wald (2014). \"Waldschutzsituation: Tabellarischer Überblick\". AFZ Der Wald - Allgemeine Forst Zeitschrift für Waldwirtschaft und Umweltvorsorge, 69(7). BLV Verlagsgesellschaft, Berlin. (accessed through the DFDE, EFI; original not seen)"},{"title":"AFZ Der Wald (2015). \"Übersicht Waldschutz 2014/2015\". AFZ Der Wald - Allgemeine Forst Zeitschrift für Waldwirtschaft und Umweltvorsorge, 70(7). BLV Verlagsgesellschaft, Berlin.","url":"https://dfde.efi.int/db/dfde_app.php","role":"secondary","citation":"AFZ Der Wald (2015). \"Übersicht Waldschutz 2014/2015\". AFZ Der Wald - Allgemeine Forst Zeitschrift für Waldwirtschaft und Umweltvorsorge, 70(7). BLV Verlagsgesellschaft, Berlin. (accessed through the DFDE, EFI; original not seen)"}]},{"id":"BB0206","name":"Brandenburg spruce bark-beetle outbreak (2019–2021), state reports via DFDE","country":"Germany","y0":2019,"y1":2021,"lat":52.459,"lon":14.159,"sev":"high","conf":"medium","note":"Reported Ips typographus damage in Brandenburg of at least 50,000 m3 per year in 2019-2021, peak 297,598 m3 in 2020 (Waldschutz reports, AFZ Der Wald).","sources":[{"title":"AFZ Der Wald (2020). \"Übersicht Waldschutz 2019/2020\". AFZ Der Wald - Allgemeine Forst Zeitschrift für Waldwirtschaft und Umweltvorsorge, 75(11). BLV Verlagsgesellschaft, Berlin.","url":"https://dfde.efi.int/db/dfde_app.php","role":"secondary","citation":"AFZ Der Wald (2020). \"Übersicht Waldschutz 2019/2020\". AFZ Der Wald - Allgemeine Forst Zeitschrift für Waldwirtschaft und Umweltvorsorge, 75(11). BLV Verlagsgesellschaft, Berlin. (accessed through the DFDE, EFI; original not seen)"},{"title":"AFZ Der Wald (2021). \"Übersicht Waldschutz 2020/2021\". AFZ Der Wald - Allgemeine Forst Zeitschrift für Waldwirtschaft und Umweltvorsorge, 76(9). BLV Verlagsgesellschaft, Berlin.","url":"https://dfde.efi.int/db/dfde_app.php","role":"secondary","citation":"AFZ Der Wald (2021). \"Übersicht Waldschutz 2020/2021\". AFZ Der Wald - Allgemeine Forst Zeitschrift für Waldwirtschaft und Umweltvorsorge, 76(9). BLV Verlagsgesellschaft, Berlin. (accessed through the DFDE, EFI; original not seen)"},{"title":"AFZ Der Wald (2022)","url":"https://dfde.efi.int/db/dfde_app.php","role":"secondary","citation":"AFZ Der Wald (2022) (accessed through the DFDE, EFI; original not seen)"}]},{"id":"BB0207","name":"Mecklenburg-Vorpommern spruce bark-beetle outbreak (2018–2021), state reports via DFDE","country":"Germany","y0":2018,"y1":2021,"lat":53.797,"lon":12.518,"sev":"moderate","conf":"medium","note":"Reported Ips typographus damage in Mecklenburg-Vorpommern of at least 50,000 m3 per year in 2018-2021, peak 167,000 m3 in 2019 (Waldschutz reports, AFZ Der Wald).","sources":[{"title":"Waldschutzinformation 3/14: Stehendbefallserfassung der Rindenbrüter Januar/Februar 2014","url":"https://www.wald-mv.de/serviceassistent/download?id=1647582","role":"secondary","citation":"Landesforst Mecklenburg-Vorpommern (2014). Waldschutzinformation 3/14: Stehendbefallserfassung der Rindenbrüter Januar/Februar 2014. Waldschutzinformation 3/14."},{"title":"AFZ Der Wald (2019). \"Übersicht Waldschutz 2018/2019\". AFZ Der Wald - Allgemeine Forst Zeitschrift für Waldwirtschaft und Umweltvorsorge, 74(7). BLV Verlagsgesellschaft, Berlin.","url":"https://dfde.efi.int/db/dfde_app.php","role":"secondary","citation":"AFZ Der Wald (2019). \"Übersicht Waldschutz 2018/2019\". AFZ Der Wald - Allgemeine Forst Zeitschrift für Waldwirtschaft und Umweltvorsorge, 74(7). BLV Verlagsgesellschaft, Berlin. (accessed through the DFDE, EFI; original not seen)"},{"title":"AFZ Der Wald (2020). \"Übersicht Waldschutz 2019/2020\". AFZ Der Wald - Allgemeine Forst Zeitschrift für Waldwirtschaft und Umweltvorsorge, 75(11). BLV Verlagsgesellschaft, Berlin.","url":"https://dfde.efi.int/db/dfde_app.php","role":"secondary","citation":"AFZ Der Wald (2020). \"Übersicht Waldschutz 2019/2020\". AFZ Der Wald - Allgemeine Forst Zeitschrift für Waldwirtschaft und Umweltvorsorge, 75(11). BLV Verlagsgesellschaft, Berlin. (accessed through the DFDE, EFI; original not seen)"},{"title":"AFZ Der Wald (2021). \"Übersicht Waldschutz 2020/2021\". AFZ Der Wald - Allgemeine Forst Zeitschrift für Waldwirtschaft und Umweltvorsorge, 76(9). BLV Verlagsgesellschaft, Berlin.","url":"https://dfde.efi.int/db/dfde_app.php","role":"secondary","citation":"AFZ Der Wald (2021). \"Übersicht Waldschutz 2020/2021\". AFZ Der Wald - Allgemeine Forst Zeitschrift für Waldwirtschaft und Umweltvorsorge, 76(9). BLV Verlagsgesellschaft, Berlin. (accessed through the DFDE, EFI; original not seen)"},{"title":"AFZ Der Wald (2022)","url":"https://dfde.efi.int/db/dfde_app.php","role":"secondary","citation":"AFZ Der Wald (2022) (accessed through the DFDE, EFI; original not seen)"}]},{"id":"BB0208","name":"Nordrhein-Westfalen spruce bark-beetle outbreak (2003–2007), state reports via DFDE","country":"Germany","y0":2003,"y1":2007,"lat":51.426,"lon":7.58,"sev":"moderate","conf":"medium","note":"Reported Ips typographus damage in Nordrhein-Westfalen of at least 50,000 m3 per year in 2003-2007, peak 151,957 m3 in 2006 (Waldschutz reports, AFZ Der Wald).","sources":[{"title":"AFZ Der Wald (2004). \"Waldschutzsituation 2003/2004: Tabellarischer Überblick\". AFZ Der Wald - Allgemeine Forst Zeitschrift für Waldwirtschaft und Umweltvorsorge, 59(7). BLV Verlagsgesellschaft, Berlin.","url":"https://dfde.efi.int/db/dfde_app.php","role":"secondary","citation":"AFZ Der Wald (2004). \"Waldschutzsituation 2003/2004: Tabellarischer Überblick\". AFZ Der Wald - Allgemeine Forst Zeitschrift für Waldwirtschaft und Umweltvorsorge, 59(7). BLV Verlagsgesellschaft, Berlin. (accessed through the DFDE, EFI; original not seen)"},{"title":"AFZ Der Wald (2005). \"Tabellarischer Überblick: Waldschutz 2004/2005\". AFZ Der Wald - Allgemeine Forst Zeitschrift für Waldwirtschaft und Umweltvorsorge, 60(7). BLV Verlagsgesellschaft, Berlin.","url":"https://dfde.efi.int/db/dfde_app.php","role":"secondary","citation":"AFZ Der Wald (2005). \"Tabellarischer Überblick: Waldschutz 2004/2005\". AFZ Der Wald - Allgemeine Forst Zeitschrift für Waldwirtschaft und Umweltvorsorge, 60(7). BLV Verlagsgesellschaft, Berlin. (accessed through the DFDE, EFI; original not seen)"},{"title":"AFZ Der Wald (2006). \"Tabellarischer Überblick: Waldschutz 2005/2006\". AFZ Der Wald - Allgemeine Forst Zeitschrift für Waldwirtschaft und Umweltvorsorge, 61(7). BLV Verlagsgesellschaft, Berlin.","url":"https://dfde.efi.int/db/dfde_app.php","role":"secondary","citation":"AFZ Der Wald (2006). \"Tabellarischer Überblick: Waldschutz 2005/2006\". AFZ Der Wald - Allgemeine Forst Zeitschrift für Waldwirtschaft und Umweltvorsorge, 61(7). BLV Verlagsgesellschaft, Berlin. (accessed through the DFDE, EFI; original not seen)"},{"title":"AFZ Der Wald (2007). \"Tabellarischer Überblick: Waldschutz 2006/2007\". AFZ Der Wald - Allgemeine Forst Zeitschrift für Waldwirtschaft und Umweltvorsorge, 62(7). BLV Verlagsgesellschaft, Berlin.","url":"https://dfde.efi.int/db/dfde_app.php","role":"secondary","citation":"AFZ Der Wald (2007). \"Tabellarischer Überblick: Waldschutz 2006/2007\". AFZ Der Wald - Allgemeine Forst Zeitschrift für Waldwirtschaft und Umweltvorsorge, 62(7). BLV Verlagsgesellschaft, Berlin. (accessed through the DFDE, EFI; original not seen)"},{"title":"AFZ Der Wald (2008). \"Überblick: Waldschutz 2007/2008\". AFZ Der Wald - Allgemeine Forst Zeitschrift für Waldwirtschaft und Umweltvorsorge, 63(7). BLV Verlagsgesellschaft, Berlin.","url":"https://dfde.efi.int/db/dfde_app.php","role":"secondary","citation":"AFZ Der Wald (2008). \"Überblick: Waldschutz 2007/2008\". AFZ Der Wald - Allgemeine Forst Zeitschrift für Waldwirtschaft und Umweltvorsorge, 63(7). BLV Verlagsgesellschaft, Berlin. (accessed through the DFDE, EFI; original not seen)"}]},{"id":"BB0209","name":"Nordrhein-Westfalen spruce bark-beetle outbreak (2012), state reports via DFDE","country":"Germany","y0":2012,"y1":2012,"lat":51.426,"lon":7.58,"sev":"moderate","conf":"medium","note":"Reported Ips typographus damage in Nordrhein-Westfalen of at least 50,000 m3 per year in 2012-2012, peak 62,983 m3 in 2012 (Waldschutz reports, AFZ Der Wald).","sources":[{"title":"AFZ Der Wald (2013). \"Tabellarischer Überblick: Waldschutzsituation 2012/2013\". AFZ Der Wald - Allgemeine Forst Zeitschrift für Waldwirtschaft und Umweltvorsorge, 68(7). BLV Verlagsgesellschaft, Berlin.","url":"https://dfde.efi.int/db/dfde_app.php","role":"secondary","citation":"AFZ Der Wald (2013). \"Tabellarischer Überblick: Waldschutzsituation 2012/2013\". AFZ Der Wald - Allgemeine Forst Zeitschrift für Waldwirtschaft und Umweltvorsorge, 68(7). BLV Verlagsgesellschaft, Berlin. (accessed through the DFDE, EFI; original not seen)"}]},{"id":"BB0210","name":"Rheinland-Pfalz spruce bark-beetle outbreak (2003–2009), state reports via DFDE","country":"Germany","y0":2003,"y1":2009,"lat":49.954,"lon":7.274,"sev":"high","conf":"medium","note":"Reported Ips typographus damage in Rheinland-Pfalz of at least 50,000 m3 per year in 2003-2009, peak 315,000 m3 in 2006 (Waldschutz reports, AFZ Der Wald).","sources":[{"title":"AFZ Der Wald (2004). \"Waldschutzsituation 2003/2004: Tabellarischer Überblick\". AFZ Der Wald - Allgemeine Forst Zeitschrift für Waldwirtschaft und Umweltvorsorge, 59(7). BLV Verlagsgesellschaft, Berlin.","url":"https://dfde.efi.int/db/dfde_app.php","role":"secondary","citation":"AFZ Der Wald (2004). \"Waldschutzsituation 2003/2004: Tabellarischer Überblick\". AFZ Der Wald - Allgemeine Forst Zeitschrift für Waldwirtschaft und Umweltvorsorge, 59(7). BLV Verlagsgesellschaft, Berlin. (accessed through the DFDE, EFI; original not seen)"},{"title":"AFZ Der Wald (2005). \"Tabellarischer Überblick: Waldschutz 2004/2005\". AFZ Der Wald - Allgemeine Forst Zeitschrift für Waldwirtschaft und Umweltvorsorge, 60(7). BLV Verlagsgesellschaft, Berlin.","url":"https://dfde.efi.int/db/dfde_app.php","role":"secondary","citation":"AFZ Der Wald (2005). \"Tabellarischer Überblick: Waldschutz 2004/2005\". AFZ Der Wald - Allgemeine Forst Zeitschrift für Waldwirtschaft und Umweltvorsorge, 60(7). BLV Verlagsgesellschaft, Berlin. (accessed through the DFDE, EFI; original not seen)"},{"title":"AFZ Der Wald (2006). \"Tabellarischer Überblick: Waldschutz 2005/2006\". AFZ Der Wald - Allgemeine Forst Zeitschrift für Waldwirtschaft und Umweltvorsorge, 61(7). BLV Verlagsgesellschaft, Berlin.","url":"https://dfde.efi.int/db/dfde_app.php","role":"secondary","citation":"AFZ Der Wald (2006). \"Tabellarischer Überblick: Waldschutz 2005/2006\". AFZ Der Wald - Allgemeine Forst Zeitschrift für Waldwirtschaft und Umweltvorsorge, 61(7). BLV Verlagsgesellschaft, Berlin. (accessed through the DFDE, EFI; original not seen)"},{"title":"AFZ Der Wald (2007). \"Tabellarischer Überblick: Waldschutz 2006/2007\". AFZ Der Wald - Allgemeine Forst Zeitschrift für Waldwirtschaft und Umweltvorsorge, 62(7). BLV Verlagsgesellschaft, Berlin.","url":"https://dfde.efi.int/db/dfde_app.php","role":"secondary","citation":"AFZ Der Wald (2007). \"Tabellarischer Überblick: Waldschutz 2006/2007\". AFZ Der Wald - Allgemeine Forst Zeitschrift für Waldwirtschaft und Umweltvorsorge, 62(7). BLV Verlagsgesellschaft, Berlin. (accessed through the DFDE, EFI; original not seen)"},{"title":"AFZ Der Wald (2008). \"Überblick: Waldschutz 2007/2008\". AFZ Der Wald - Allgemeine Forst Zeitschrift für Waldwirtschaft und Umweltvorsorge, 63(7). BLV Verlagsgesellschaft, Berlin.","url":"https://dfde.efi.int/db/dfde_app.php","role":"secondary","citation":"AFZ Der Wald (2008). \"Überblick: Waldschutz 2007/2008\". AFZ Der Wald - Allgemeine Forst Zeitschrift für Waldwirtschaft und Umweltvorsorge, 63(7). BLV Verlagsgesellschaft, Berlin. (accessed through the DFDE, EFI; original not seen)"},{"title":"AFZ Der Wald (2009). \"Überblick: Waldschutz 2008/2009\". AFZ Der Wald - Allgemeine Forst Zeitschrift für Waldwirtschaft und Umweltvorsorge, 64(7). BLV Verlagsgesellschaft, Berlin.","url":"https://dfde.efi.int/db/dfde_app.php","role":"secondary","citation":"AFZ Der Wald (2009). \"Überblick: Waldschutz 2008/2009\". AFZ Der Wald - Allgemeine Forst Zeitschrift für Waldwirtschaft und Umweltvorsorge, 64(7). BLV Verlagsgesellschaft, Berlin. (accessed through the DFDE, EFI; original not seen)"},{"title":"AFZ Der Wald (2010). \"Tabellarischer Überblick: Waldschutz 2009/2010\". AFZ Der Wald - Allgemeine Forst Zeitschrift für Waldwirtschaft und Umweltvorsorge, 65(7). BLV Verlagsgesellschaft, Berlin.","url":"https://dfde.efi.int/db/dfde_app.php","role":"secondary","citation":"AFZ Der Wald (2010). \"Tabellarischer Überblick: Waldschutz 2009/2010\". AFZ Der Wald - Allgemeine Forst Zeitschrift für Waldwirtschaft und Umweltvorsorge, 65(7). BLV Verlagsgesellschaft, Berlin. (accessed through the DFDE, EFI; original not seen)"}]},{"id":"BB0211","name":"Rheinland-Pfalz spruce bark-beetle outbreak (2011–2013), state reports via DFDE","country":"Germany","y0":2011,"y1":2013,"lat":49.954,"lon":7.274,"sev":"high","conf":"medium","note":"Reported Ips typographus damage in Rheinland-Pfalz of at least 50,000 m3 per year in 2011-2013, peak 235,000 m3 in 2011 (Waldschutz reports, AFZ Der Wald).","sources":[{"title":"AFZ Der Wald (2012). \"Waldschutzsituation: Tabellarischer Überblick 2011/2012\". AFZ Der Wald - Allgemeine Forst Zeitschrift für Waldwirtschaft und Umweltvorsorge, 67(7). BLV Verlagsgesellschaft, Berlin.","url":"https://dfde.efi.int/db/dfde_app.php","role":"secondary","citation":"AFZ Der Wald (2012). \"Waldschutzsituation: Tabellarischer Überblick 2011/2012\". AFZ Der Wald - Allgemeine Forst Zeitschrift für Waldwirtschaft und Umweltvorsorge, 67(7). BLV Verlagsgesellschaft, Berlin. (accessed through the DFDE, EFI; original not seen)"},{"title":"AFZ Der Wald (2013). \"Tabellarischer Überblick: Waldschutzsituation 2012/2013\". AFZ Der Wald - Allgemeine Forst Zeitschrift für Waldwirtschaft und Umweltvorsorge, 68(7). BLV Verlagsgesellschaft, Berlin.","url":"https://dfde.efi.int/db/dfde_app.php","role":"secondary","citation":"AFZ Der Wald (2013). \"Tabellarischer Überblick: Waldschutzsituation 2012/2013\". AFZ Der Wald - Allgemeine Forst Zeitschrift für Waldwirtschaft und Umweltvorsorge, 68(7). BLV Verlagsgesellschaft, Berlin. (accessed through the DFDE, EFI; original not seen)"},{"title":"AFZ Der Wald (2014). \"Waldschutzsituation: Tabellarischer Überblick\". AFZ Der Wald - Allgemeine Forst Zeitschrift für Waldwirtschaft und Umweltvorsorge, 69(7). BLV Verlagsgesellschaft, Berlin.","url":"https://dfde.efi.int/db/dfde_app.php","role":"secondary","citation":"AFZ Der Wald (2014). \"Waldschutzsituation: Tabellarischer Überblick\". AFZ Der Wald - Allgemeine Forst Zeitschrift für Waldwirtschaft und Umweltvorsorge, 69(7). BLV Verlagsgesellschaft, Berlin. (accessed through the DFDE, EFI; original not seen)"}]},{"id":"BB0212","name":"Sachsen spruce bark-beetle outbreak (2002–2004), state reports via DFDE","country":"Germany","y0":2002,"y1":2004,"lat":50.928,"lon":13.522,"sev":"moderate","conf":"medium","note":"Reported Ips typographus damage in Sachsen of at least 50,000 m3 per year in 2002-2004, peak 90,000 m3 in 2003 (Waldschutz reports, AFZ Der Wald).","sources":[{"title":"AFZ Der Wald (2003). \"Tabellarischer Überblick: Waldschutz 2002/2003\". AFZ Der Wald - Allgemeine Forst Zeitschrift für Waldwirtschaft und Umweltvorsorge, 58(7). BLV Verlagsgesellschaft, Berlin.","url":"https://dfde.efi.int/db/dfde_app.php","role":"secondary","citation":"AFZ Der Wald (2003). \"Tabellarischer Überblick: Waldschutz 2002/2003\". AFZ Der Wald - Allgemeine Forst Zeitschrift für Waldwirtschaft und Umweltvorsorge, 58(7). BLV Verlagsgesellschaft, Berlin. (accessed through the DFDE, EFI; original not seen)"},{"title":"AFZ Der Wald (2004). \"Waldschutzsituation 2003/2004: Tabellarischer Überblick\". AFZ Der Wald - Allgemeine Forst Zeitschrift für Waldwirtschaft und Umweltvorsorge, 59(7). BLV Verlagsgesellschaft, Berlin.","url":"https://dfde.efi.int/db/dfde_app.php","role":"secondary","citation":"AFZ Der Wald (2004). \"Waldschutzsituation 2003/2004: Tabellarischer Überblick\". AFZ Der Wald - Allgemeine Forst Zeitschrift für Waldwirtschaft und Umweltvorsorge, 59(7). BLV Verlagsgesellschaft, Berlin. (accessed through the DFDE, EFI; original not seen)"},{"title":"AFZ Der Wald (2005). \"Tabellarischer Überblick: Waldschutz 2004/2005\". AFZ Der Wald - Allgemeine Forst Zeitschrift für Waldwirtschaft und Umweltvorsorge, 60(7). BLV Verlagsgesellschaft, Berlin.","url":"https://dfde.efi.int/db/dfde_app.php","role":"secondary","citation":"AFZ Der Wald (2005). \"Tabellarischer Überblick: Waldschutz 2004/2005\". AFZ Der Wald - Allgemeine Forst Zeitschrift für Waldwirtschaft und Umweltvorsorge, 60(7). BLV Verlagsgesellschaft, Berlin. (accessed through the DFDE, EFI; original not seen)"}]},{"id":"BB0213","name":"Sachsen spruce bark-beetle outbreak (2006), state reports via DFDE","country":"Germany","y0":2006,"y1":2006,"lat":50.928,"lon":13.522,"sev":"moderate","conf":"medium","note":"Reported Ips typographus damage in Sachsen of at least 50,000 m3 per year in 2006-2006, peak 72,000 m3 in 2006 (Waldschutz reports, AFZ Der Wald).","sources":[{"title":"AFZ Der Wald (2007). \"Tabellarischer Überblick: Waldschutz 2006/2007\". AFZ Der Wald - Allgemeine Forst Zeitschrift für Waldwirtschaft und Umweltvorsorge, 62(7). BLV Verlagsgesellschaft, Berlin.","url":"https://dfde.efi.int/db/dfde_app.php","role":"secondary","citation":"AFZ Der Wald (2007). \"Tabellarischer Überblick: Waldschutz 2006/2007\". AFZ Der Wald - Allgemeine Forst Zeitschrift für Waldwirtschaft und Umweltvorsorge, 62(7). BLV Verlagsgesellschaft, Berlin. (accessed through the DFDE, EFI; original not seen)"}]},{"id":"BB0214","name":"Sachsen spruce bark-beetle outbreak (2008), state reports via DFDE","country":"Germany","y0":2008,"y1":2008,"lat":50.928,"lon":13.522,"sev":"moderate","conf":"medium","note":"Reported Ips typographus damage in Sachsen of at least 50,000 m3 per year in 2008-2008, peak 119,000 m3 in 2008 (Waldschutz reports, AFZ Der Wald).","sources":[{"title":"AFZ Der Wald (2009). \"Überblick: Waldschutz 2008/2009\". AFZ Der Wald - Allgemeine Forst Zeitschrift für Waldwirtschaft und Umweltvorsorge, 64(7). BLV Verlagsgesellschaft, Berlin.","url":"https://dfde.efi.int/db/dfde_app.php","role":"secondary","citation":"AFZ Der Wald (2009). \"Überblick: Waldschutz 2008/2009\". AFZ Der Wald - Allgemeine Forst Zeitschrift für Waldwirtschaft und Umweltvorsorge, 64(7). BLV Verlagsgesellschaft, Berlin. (accessed through the DFDE, EFI; original not seen)"}]},{"id":"BB0215","name":"Sachsen-Anhalt spruce bark-beetle outbreak (2002–2006), state reports via DFDE","country":"Germany","y0":2002,"y1":2006,"lat":51.99,"lon":11.559,"sev":"moderate","conf":"medium","note":"Reported Ips typographus damage in Sachsen-Anhalt of at least 50,000 m3 per year in 2002-2006, peak 106,392 m3 in 2006 (Waldschutz reports, AFZ Der Wald).","sources":[{"title":"AFZ Der Wald (2003). \"Tabellarischer Überblick: Waldschutz 2002/2003\". AFZ Der Wald - Allgemeine Forst Zeitschrift für Waldwirtschaft und Umweltvorsorge, 58(7). BLV Verlagsgesellschaft, Berlin.","url":"https://dfde.efi.int/db/dfde_app.php","role":"secondary","citation":"AFZ Der Wald (2003). \"Tabellarischer Überblick: Waldschutz 2002/2003\". AFZ Der Wald - Allgemeine Forst Zeitschrift für Waldwirtschaft und Umweltvorsorge, 58(7). BLV Verlagsgesellschaft, Berlin. (accessed through the DFDE, EFI; original not seen)"},{"title":"AFZ Der Wald (2004). \"Waldschutzsituation 2003/2004: Tabellarischer Überblick\". AFZ Der Wald - Allgemeine Forst Zeitschrift für Waldwirtschaft und Umweltvorsorge, 59(7). BLV Verlagsgesellschaft, Berlin.","url":"https://dfde.efi.int/db/dfde_app.php","role":"secondary","citation":"AFZ Der Wald (2004). \"Waldschutzsituation 2003/2004: Tabellarischer Überblick\". AFZ Der Wald - Allgemeine Forst Zeitschrift für Waldwirtschaft und Umweltvorsorge, 59(7). BLV Verlagsgesellschaft, Berlin. (accessed through the DFDE, EFI; original not seen)"},{"title":"AFZ Der Wald (2005). \"Tabellarischer Überblick: Waldschutz 2004/2005\". AFZ Der Wald - Allgemeine Forst Zeitschrift für Waldwirtschaft und Umweltvorsorge, 60(7). BLV Verlagsgesellschaft, Berlin.","url":"https://dfde.efi.int/db/dfde_app.php","role":"secondary","citation":"AFZ Der Wald (2005). \"Tabellarischer Überblick: Waldschutz 2004/2005\". AFZ Der Wald - Allgemeine Forst Zeitschrift für Waldwirtschaft und Umweltvorsorge, 60(7). BLV Verlagsgesellschaft, Berlin. (accessed through the DFDE, EFI; original not seen)"},{"title":"AFZ Der Wald (2006). \"Tabellarischer Überblick: Waldschutz 2005/2006\". AFZ Der Wald - Allgemeine Forst Zeitschrift für Waldwirtschaft und Umweltvorsorge, 61(7). BLV Verlagsgesellschaft, Berlin.","url":"https://dfde.efi.int/db/dfde_app.php","role":"secondary","citation":"AFZ Der Wald (2006). \"Tabellarischer Überblick: Waldschutz 2005/2006\". AFZ Der Wald - Allgemeine Forst Zeitschrift für Waldwirtschaft und Umweltvorsorge, 61(7). BLV Verlagsgesellschaft, Berlin. (accessed through the DFDE, EFI; original not seen)"},{"title":"AFZ Der Wald (2007). \"Tabellarischer Überblick: Waldschutz 2006/2007\". AFZ Der Wald - Allgemeine Forst Zeitschrift für Waldwirtschaft und Umweltvorsorge, 62(7). BLV Verlagsgesellschaft, Berlin.","url":"https://dfde.efi.int/db/dfde_app.php","role":"secondary","citation":"AFZ Der Wald (2007). \"Tabellarischer Überblick: Waldschutz 2006/2007\". AFZ Der Wald - Allgemeine Forst Zeitschrift für Waldwirtschaft und Umweltvorsorge, 62(7). BLV Verlagsgesellschaft, Berlin. (accessed through the DFDE, EFI; original not seen)"}]},{"id":"BB0216","name":"Schleswig-Holstein spruce bark-beetle outbreak (2019–2021), state reports via DFDE","country":"Germany","y0":2019,"y1":2021,"lat":54.208,"lon":9.956,"sev":"moderate","conf":"medium","note":"Reported Ips typographus damage in Schleswig-Holstein of at least 50,000 m3 per year in 2019-2021, peak 158,050 m3 in 2020 (Waldschutz reports, AFZ Der Wald).","sources":[{"title":"AFZ Der Wald (2020). \"Übersicht Waldschutz 2019/2020\". AFZ Der Wald - Allgemeine Forst Zeitschrift für Waldwirtschaft und Umweltvorsorge, 75(11). BLV Verlagsgesellschaft, Berlin.","url":"https://dfde.efi.int/db/dfde_app.php","role":"secondary","citation":"AFZ Der Wald (2020). \"Übersicht Waldschutz 2019/2020\". AFZ Der Wald - Allgemeine Forst Zeitschrift für Waldwirtschaft und Umweltvorsorge, 75(11). BLV Verlagsgesellschaft, Berlin. (accessed through the DFDE, EFI; original not seen)"},{"title":"AFZ Der Wald (2021). \"Übersicht Waldschutz 2020/2021\". AFZ Der Wald - Allgemeine Forst Zeitschrift für Waldwirtschaft und Umweltvorsorge, 76(9). BLV Verlagsgesellschaft, Berlin.","url":"https://dfde.efi.int/db/dfde_app.php","role":"secondary","citation":"AFZ Der Wald (2021). \"Übersicht Waldschutz 2020/2021\". AFZ Der Wald - Allgemeine Forst Zeitschrift für Waldwirtschaft und Umweltvorsorge, 76(9). BLV Verlagsgesellschaft, Berlin. (accessed through the DFDE, EFI; original not seen)"},{"title":"AFZ Der Wald (2022)","url":"https://dfde.efi.int/db/dfde_app.php","role":"secondary","citation":"AFZ Der Wald (2022) (accessed through the DFDE, EFI; original not seen)"}]},{"id":"BB0217","name":"Thüringen spruce bark-beetle outbreak (2003–2009), state reports via DFDE","country":"Germany","y0":2003,"y1":2009,"lat":50.927,"lon":11.315,"sev":"high","conf":"medium","note":"Reported Ips typographus damage in Thüringen of at least 50,000 m3 per year in 2003-2009, peak 349,600 m3 in 2008 (Waldschutz reports, AFZ Der Wald).","sources":[{"title":"AFZ Der Wald (2004). \"Waldschutzsituation 2003/2004: Tabellarischer Überblick\". AFZ Der Wald - Allgemeine Forst Zeitschrift für Waldwirtschaft und Umweltvorsorge, 59(7). BLV Verlagsgesellschaft, Berlin.","url":"https://dfde.efi.int/db/dfde_app.php","role":"secondary","citation":"AFZ Der Wald (2004). \"Waldschutzsituation 2003/2004: Tabellarischer Überblick\". AFZ Der Wald - Allgemeine Forst Zeitschrift für Waldwirtschaft und Umweltvorsorge, 59(7). BLV Verlagsgesellschaft, Berlin. (accessed through the DFDE, EFI; original not seen)"},{"title":"AFZ Der Wald (2005). \"Tabellarischer Überblick: Waldschutz 2004/2005\". AFZ Der Wald - Allgemeine Forst Zeitschrift für Waldwirtschaft und Umweltvorsorge, 60(7). BLV Verlagsgesellschaft, Berlin.","url":"https://dfde.efi.int/db/dfde_app.php","role":"secondary","citation":"AFZ Der Wald (2005). \"Tabellarischer Überblick: Waldschutz 2004/2005\". AFZ Der Wald - Allgemeine Forst Zeitschrift für Waldwirtschaft und Umweltvorsorge, 60(7). BLV Verlagsgesellschaft, Berlin. (accessed through the DFDE, EFI; original not seen)"},{"title":"AFZ Der Wald (2006). \"Tabellarischer Überblick: Waldschutz 2005/2006\". AFZ Der Wald - Allgemeine Forst Zeitschrift für Waldwirtschaft und Umweltvorsorge, 61(7). BLV Verlagsgesellschaft, Berlin.","url":"https://dfde.efi.int/db/dfde_app.php","role":"secondary","citation":"AFZ Der Wald (2006). \"Tabellarischer Überblick: Waldschutz 2005/2006\". AFZ Der Wald - Allgemeine Forst Zeitschrift für Waldwirtschaft und Umweltvorsorge, 61(7). BLV Verlagsgesellschaft, Berlin. (accessed through the DFDE, EFI; original not seen)"},{"title":"AFZ Der Wald (2007). \"Tabellarischer Überblick: Waldschutz 2006/2007\". AFZ Der Wald - Allgemeine Forst Zeitschrift für Waldwirtschaft und Umweltvorsorge, 62(7). BLV Verlagsgesellschaft, Berlin.","url":"https://dfde.efi.int/db/dfde_app.php","role":"secondary","citation":"AFZ Der Wald (2007). \"Tabellarischer Überblick: Waldschutz 2006/2007\". AFZ Der Wald - Allgemeine Forst Zeitschrift für Waldwirtschaft und Umweltvorsorge, 62(7). BLV Verlagsgesellschaft, Berlin. (accessed through the DFDE, EFI; original not seen)"},{"title":"AFZ Der Wald (2008). \"Überblick: Waldschutz 2007/2008\". AFZ Der Wald - Allgemeine Forst Zeitschrift für Waldwirtschaft und Umweltvorsorge, 63(7). BLV Verlagsgesellschaft, Berlin.","url":"https://dfde.efi.int/db/dfde_app.php","role":"secondary","citation":"AFZ Der Wald (2008). \"Überblick: Waldschutz 2007/2008\". AFZ Der Wald - Allgemeine Forst Zeitschrift für Waldwirtschaft und Umweltvorsorge, 63(7). BLV Verlagsgesellschaft, Berlin. (accessed through the DFDE, EFI; original not seen)"},{"title":"AFZ Der Wald (2009). \"Überblick: Waldschutz 2008/2009\". AFZ Der Wald - Allgemeine Forst Zeitschrift für Waldwirtschaft und Umweltvorsorge, 64(7). BLV Verlagsgesellschaft, Berlin.","url":"https://dfde.efi.int/db/dfde_app.php","role":"secondary","citation":"AFZ Der Wald (2009). \"Überblick: Waldschutz 2008/2009\". AFZ Der Wald - Allgemeine Forst Zeitschrift für Waldwirtschaft und Umweltvorsorge, 64(7). BLV Verlagsgesellschaft, Berlin. (accessed through the DFDE, EFI; original not seen)"},{"title":"AFZ Der Wald (2010). \"Tabellarischer Überblick: Waldschutz 2009/2010\". AFZ Der Wald - Allgemeine Forst Zeitschrift für Waldwirtschaft und Umweltvorsorge, 65(7). BLV Verlagsgesellschaft, Berlin.","url":"https://dfde.efi.int/db/dfde_app.php","role":"secondary","citation":"AFZ Der Wald (2010). \"Tabellarischer Überblick: Waldschutz 2009/2010\". AFZ Der Wald - Allgemeine Forst Zeitschrift für Waldwirtschaft und Umweltvorsorge, 65(7). BLV Verlagsgesellschaft, Berlin. (accessed through the DFDE, EFI; original not seen)"}]},{"id":"BB0218","name":"Gdańsk regional directorate (Poland) spruce bark-beetle and drought salvage (2006–2010), State Forests via DFDE","country":"Poland","y0":2006,"y1":2010,"lat":54.397,"lon":18.22,"sev":"moderate","conf":"medium","note":"Spruce salvage felling by bark beetles and drought in the Gdańsk directorate peaked at 198,724 m3 in 2007; largest districts in the run: Kartuzy (178k m3), Gdańsk (149k m3), Kolbudy (112k m3), Wejherowo (100k m3).","sources":[{"title":"Forest damage in Poland 2005-2014, Information System of the State Forests, Directorate General of State Forest, Poland, 2020","url":"https://dfde.efi.int/db/dfde_app.php","role":"secondary","citation":"Forest damage in Poland 2005-2014, Information System of the State Forests, Directorate General of State Forest, Poland, 2020 (accessed through the DFDE, EFI; original not seen)"}]},{"id":"BB0219","name":"Gdańsk regional directorate (Poland) spruce bark-beetle and drought salvage (2013), State Forests via DFDE","country":"Poland","y0":2013,"y1":2013,"lat":54.342,"lon":18.242,"sev":"moderate","conf":"medium","note":"Spruce salvage felling by bark beetles and drought in the Gdańsk directorate peaked at 81,263 m3 in 2013; largest districts in the run: Kartuzy (26k m3), Kolbudy (16k m3), Gdańsk (14k m3), Wejherowo (8k m3).","sources":[{"title":"Forest damage in Poland 2005-2014, Information System of the State Forests, Directorate General of State Forest, Poland, 2020","url":"https://dfde.efi.int/db/dfde_app.php","role":"secondary","citation":"Forest damage in Poland 2005-2014, Information System of the State Forests, Directorate General of State Forest, Poland, 2020 (accessed through the DFDE, EFI; original not seen)"}]},{"id":"BB0220","name":"Katowice regional directorate (Poland) spruce bark-beetle and drought salvage (2006–2009), State Forests via DFDE","country":"Poland","y0":2006,"y1":2009,"lat":49.711,"lon":18.935,"sev":"extreme","conf":"medium","note":"Spruce salvage felling by bark beetles and drought in the Katowice directorate peaked at 1,104,155 m3 in 2007; largest districts in the run: Ujsoły (804k m3), Węgierska Górka (762k m3), Wisła (590k m3), Ustroń (301k m3).","sources":[{"title":"Forest damage in Poland 2005-2014, Information System of the State Forests, Directorate General of State Forest, Poland, 2020","url":"https://dfde.efi.int/db/dfde_app.php","role":"secondary","citation":"Forest damage in Poland 2005-2014, Information System of the State Forests, Directorate General of State Forest, Poland, 2020 (accessed through the DFDE, EFI; original not seen)"}]},{"id":"BB0221","name":"Kraków regional directorate (Poland) spruce bark-beetle and drought salvage (2006–2009), State Forests via DFDE","country":"Poland","y0":2006,"y1":2009,"lat":49.537,"lon":20.202,"sev":"moderate","conf":"medium","note":"Spruce salvage felling by bark beetles and drought in the Kraków directorate peaked at 79,583 m3 in 2006; largest districts in the run: Nowy Targ (95k m3), Myślenice (57k m3), Krościenko (55k m3), Piwniczna (34k m3).","sources":[{"title":"Forest damage in Poland 2005-2014, Information System of the State Forests, Directorate General of State Forest, Poland, 2020","url":"https://dfde.efi.int/db/dfde_app.php","role":"secondary","citation":"Forest damage in Poland 2005-2014, Information System of the State Forests, Directorate General of State Forest, Poland, 2020 (accessed through the DFDE, EFI; original not seen)"}]},{"id":"BB0222","name":"Olsztyn regional directorate (Poland) spruce bark-beetle and drought salvage (2005–2009), State Forests via DFDE","country":"Poland","y0":2005,"y1":2009,"lat":53.724,"lon":20.531,"sev":"moderate","conf":"medium","note":"Spruce salvage felling by bark beetles and drought in the Olsztyn directorate peaked at 102,988 m3 in 2008; largest districts in the run: Wipsowo (29k m3), Nidzica (23k m3), Strzałowo (22k m3), Lidzbark (21k m3).","sources":[{"title":"Forest damage in Poland 2005-2014, Information System of the State Forests, Directorate General of State Forest, Poland, 2020","url":"https://dfde.efi.int/db/dfde_app.php","role":"secondary","citation":"Forest damage in Poland 2005-2014, Information System of the State Forests, Directorate General of State Forest, Poland, 2020 (accessed through the DFDE, EFI; original not seen)"}]},{"id":"BB0223","name":"Olsztyn regional directorate (Poland) spruce bark-beetle and drought salvage (2013–2019), State Forests via DFDE","country":"Poland","y0":2013,"y1":2019,"lat":53.729,"lon":20.517,"sev":"moderate","conf":"medium","note":"Spruce salvage felling by bark beetles and drought in the Olsztyn directorate peaked at 66,699 m3 in 2019; largest districts in the run: Wipsowo (37k m3), Spychowo (25k m3), Mrągowo (24k m3), Susz (24k m3).","sources":[{"title":"Forest damage in Poland 2005-2014, Information System of the State Forests, Directorate General of State Forest, Poland, 2020","url":"https://dfde.efi.int/db/dfde_app.php","role":"secondary","citation":"Forest damage in Poland 2005-2014, Information System of the State Forests, Directorate General of State Forest, Poland, 2020 (accessed through the DFDE, EFI; original not seen)"},{"title":"Forest damage in Poland 2015-2019, Information System of the State Forests, Directorate General of State Forest, Poland, 2019 ","url":"https://dfde.efi.int/db/dfde_app.php","role":"secondary","citation":"Forest damage in Poland 2015-2019, Information System of the State Forests, Directorate General of State Forest, Poland, 2019  (accessed through the DFDE, EFI; original not seen)"}]},{"id":"BB0224","name":"Poznań regional directorate (Poland) spruce bark-beetle and drought salvage (2019), State Forests via DFDE","country":"Poland","y0":2019,"y1":2019,"lat":52.092,"lon":17.06,"sev":"moderate","conf":"medium","note":"Spruce salvage felling by bark beetles and drought in the Poznań directorate peaked at 71,339 m3 in 2019; largest districts in the run: Kościan (11k m3), Piaski (6k m3), Grodzisk (6k m3), Czerniejewo (6k m3).","sources":[{"title":"Forest damage in Poland 2015-2019, Information System of the State Forests, Directorate General of State Forest, Poland, 2019 ","url":"https://dfde.efi.int/db/dfde_app.php","role":"secondary","citation":"Forest damage in Poland 2015-2019, Information System of the State Forests, Directorate General of State Forest, Poland, 2019  (accessed through the DFDE, EFI; original not seen)"}]},{"id":"BB0225","name":"Szczecin regional directorate (Poland) spruce bark-beetle and drought salvage (2019), State Forests via DFDE","country":"Poland","y0":2019,"y1":2019,"lat":53.251,"lon":15.11,"sev":"moderate","conf":"medium","note":"Spruce salvage felling by bark beetles and drought in the Szczecin directorate peaked at 129,063 m3 in 2019; largest districts in the run: Resko (23k m3), Mieszkowice (15k m3), Dobrzany (9k m3), Łobez (8k m3).","sources":[{"title":"Forest damage in Poland 2015-2019, Information System of the State Forests, Directorate General of State Forest, Poland, 2019 ","url":"https://dfde.efi.int/db/dfde_app.php","role":"secondary","citation":"Forest damage in Poland 2015-2019, Information System of the State Forests, Directorate General of State Forest, Poland, 2019  (accessed through the DFDE, EFI; original not seen)"}]},{"id":"BB0226","name":"Szczecinek regional directorate (Poland) spruce bark-beetle and drought salvage (2007–2008), State Forests via DFDE","country":"Poland","y0":2007,"y1":2008,"lat":54.269,"lon":17.041,"sev":"moderate","conf":"medium","note":"Spruce salvage felling by bark beetles and drought in the Szczecinek directorate peaked at 156,579 m3 in 2007; largest districts in the run: Leśny Dwór (103k m3), Bytów (38k m3), Damnica (27k m3), Sławno (18k m3).","sources":[{"title":"Forest damage in Poland 2005-2014, Information System of the State Forests, Directorate General of State Forest, Poland, 2020","url":"https://dfde.efi.int/db/dfde_app.php","role":"secondary","citation":"Forest damage in Poland 2005-2014, Information System of the State Forests, Directorate General of State Forest, Poland, 2020 (accessed through the DFDE, EFI; original not seen)"}]},{"id":"BB0227","name":"Szczecinek regional directorate (Poland) spruce bark-beetle and drought salvage (2013–2019), State Forests via DFDE","country":"Poland","y0":2013,"y1":2019,"lat":54.155,"lon":16.809,"sev":"high","conf":"medium","note":"Spruce salvage felling by bark beetles and drought in the Szczecinek directorate peaked at 223,196 m3 in 2015; largest districts in the run: Sławno (204k m3), Leśny Dwór (125k m3), Polanów (97k m3), Warcino (84k m3).","sources":[{"title":"Forest damage in Poland 2005-2014, Information System of the State Forests, Directorate General of State Forest, Poland, 2020","url":"https://dfde.efi.int/db/dfde_app.php","role":"secondary","citation":"Forest damage in Poland 2005-2014, Information System of the State Forests, Directorate General of State Forest, Poland, 2020 (accessed through the DFDE, EFI; original not seen)"},{"title":"Forest damage in Poland 2015-2019, Information System of the State Forests, Directorate General of State Forest, Poland, 2019 ","url":"https://dfde.efi.int/db/dfde_app.php","role":"secondary","citation":"Forest damage in Poland 2015-2019, Information System of the State Forests, Directorate General of State Forest, Poland, 2019  (accessed through the DFDE, EFI; original not seen)"}]},{"id":"BB0228","name":"Wrocław regional directorate (Poland) spruce bark-beetle and drought salvage (2008), State Forests via DFDE","country":"Poland","y0":2008,"y1":2008,"lat":50.801,"lon":16.306,"sev":"high","conf":"medium","note":"Spruce salvage felling by bark beetles and drought in the Wrocław directorate peaked at 396,765 m3 in 2008; largest districts in the run: Wałbrzych (62k m3), Jugów (34k m3), Kamienna Góra (34k m3), Świdnica (28k m3).","sources":[{"title":"Forest damage in Poland 2005-2014, Information System of the State Forests, Directorate General of State Forest, Poland, 2020","url":"https://dfde.efi.int/db/dfde_app.php","role":"secondary","citation":"Forest damage in Poland 2005-2014, Information System of the State Forests, Directorate General of State Forest, Poland, 2020 (accessed through the DFDE, EFI; original not seen)"}]},{"id":"BB0229","name":"Zielona Góra regional directorate (Poland) spruce bark-beetle and drought salvage (2019), State Forests via DFDE","country":"Poland","y0":2019,"y1":2019,"lat":51.999,"lon":15.27,"sev":"moderate","conf":"medium","note":"Spruce salvage felling by bark beetles and drought in the Zielona Góra directorate peaked at 55,844 m3 in 2019; largest districts in the run: Lubsko (9k m3), Sulechów (7k m3), Świebodzin (6k m3), Torzym (5k m3).","sources":[{"title":"Forest damage in Poland 2015-2019, Information System of the State Forests, Directorate General of State Forest, Poland, 2019 ","url":"https://dfde.efi.int/db/dfde_app.php","role":"secondary","citation":"Forest damage in Poland 2015-2019, Information System of the State Forests, Directorate General of State Forest, Poland, 2019  (accessed through the DFDE, EFI; original not seen)"}]},{"id":"BB0230","name":"Slovakia national bark-beetle damage (2007–2013), JRC series via DFDE","country":"Slovakia","y0":2007,"y1":2013,"lat":48.7,"lon":19.5,"sev":"extreme","conf":"medium","note":"National yearly bark-beetle damage 2007-2013 (thousand m3): 2007: 2,265 thousand m3; 2008: 3,102 thousand m3; 2009: 3,530 thousand m3; 2010: 3,215 thousand m3; 2011: 2,179 thousand m3; 2012: 2,130 thousand m3; 2013: 1,705 thousand m3 (Camia","sources":[{"title":"Camia, A., Giuntoli, J., Jonsson, K., Robert, N., Cazzaniga, N., Jasinevičius, G., Avitabile, V., Grassi, G., Barredo Cano, J.I. and Mubareka, S., The use of woody biomass for energy production in the EU, EUR 30548 EN, Publications Office of the European U","url":"https://dfde.efi.int/db/dfde_app.php","role":"secondary","citation":"Camia, A., Giuntoli, J., Jonsson, K., Robert, N., Cazzaniga, N., Jasinevičius, G., Avitabile, V., Grassi, G., Barredo Cano, J.I. and Mubareka, S., The use of woody biomass for energy production in the EU, EUR 30548 EN, Publications Office of the European U (accessed through the DFDE, EFI; original not seen)"}]},{"id":"BB0231","name":"Slovakia national bark-beetle damage (2016–2019), JRC series via DFDE","country":"Slovakia","y0":2016,"y1":2019,"lat":48.7,"lon":19.5,"sev":"extreme","conf":"medium","note":"National yearly bark-beetle damage 2016-2019 (thousand m3): 2016: 3,263 thousand m3; 2017: 3,793 thousand m3; 2018: 4,067 thousand m3; 2019: 3,628 thousand m3 (Camia et al. 2021, JRC; compiled in DFDE).","sources":[{"title":"Camia, A., Giuntoli, J., Jonsson, K., Robert, N., Cazzaniga, N., Jasinevičius, G., Avitabile, V., Grassi, G., Barredo Cano, J.I. and Mubareka, S., The use of woody biomass for energy production in the EU, EUR 30548 EN, Publications Office of the European U","url":"https://dfde.efi.int/db/dfde_app.php","role":"secondary","citation":"Camia, A., Giuntoli, J., Jonsson, K., Robert, N., Cazzaniga, N., Jasinevičius, G., Avitabile, V., Grassi, G., Barredo Cano, J.I. and Mubareka, S., The use of woody biomass for energy production in the EU, EUR 30548 EN, Publications Office of the European U (accessed through the DFDE, EFI; original not seen)"}]},{"id":"BB0232","name":"Czechia national bark-beetle damage (2006–2011), JRC series via DFDE","country":"Czechia","y0":2006,"y1":2011,"lat":49.8,"lon":15.5,"sev":"extreme","conf":"medium","note":"National yearly bark-beetle damage 2006-2011 (thousand m3): 2006: 1,139 thousand m3; 2007: 1,556 thousand m3; 2008: 2,315 thousand m3; 2009: 2,624 thousand m3; 2010: 1,788 thousand m3; 2011: 1,054 thousand m3 (Camia et al. 2021, JRC; compil","sources":[{"title":"Camia, A., Giuntoli, J., Jonsson, K., Robert, N., Cazzaniga, N., Jasinevičius, G., Avitabile, V., Grassi, G., Barredo Cano, J.I. and Mubareka, S., The use of woody biomass for energy production in the EU, EUR 30548 EN, Publications Office of the European U","url":"https://dfde.efi.int/db/dfde_app.php","role":"secondary","citation":"Camia, A., Giuntoli, J., Jonsson, K., Robert, N., Cazzaniga, N., Jasinevičius, G., Avitabile, V., Grassi, G., Barredo Cano, J.I. and Mubareka, S., The use of woody biomass for energy production in the EU, EUR 30548 EN, Publications Office of the European U (accessed through the DFDE, EFI; original not seen)"}]},{"id":"BB0233","name":"Czechia national bark-beetle damage (2015–2019), JRC series via DFDE","country":"Czechia","y0":2015,"y1":2019,"lat":49.8,"lon":15.5,"sev":"extreme","conf":"medium","note":"National yearly bark-beetle damage 2015-2019 (thousand m3): 2015: 2,309 thousand m3; 2016: 4,420 thousand m3; 2017: 5,853 thousand m3; 2018: 13,059 thousand m3; 2019: 22,780 thousand m3 (Camia et al. 2021, JRC; compiled in DFDE).","sources":[{"title":"Camia, A., Giuntoli, J., Jonsson, K., Robert, N., Cazzaniga, N., Jasinevičius, G., Avitabile, V., Grassi, G., Barredo Cano, J.I. and Mubareka, S., The use of woody biomass for energy production in the EU, EUR 30548 EN, Publications Office of the European U","url":"https://dfde.efi.int/db/dfde_app.php","role":"secondary","citation":"Camia, A., Giuntoli, J., Jonsson, K., Robert, N., Cazzaniga, N., Jasinevičius, G., Avitabile, V., Grassi, G., Barredo Cano, J.I. and Mubareka, S., The use of woody biomass for energy production in the EU, EUR 30548 EN, Publications Office of the European U (accessed through the DFDE, EFI; original not seen)"}]},{"id":"BB0234","name":"Germany national bark-beetle damage (2006–2009), JRC series via DFDE","country":"Germany","y0":2006,"y1":2009,"lat":51.2,"lon":10.4,"sev":"extreme","conf":"medium","note":"National yearly bark-beetle damage 2006-2009 (thousand m3): 2006: 5,022 thousand m3; 2007: 2,497 thousand m3; 2008: 4,559 thousand m3; 2009: 3,892 thousand m3 (Camia et al. 2021, JRC; compiled in DFDE).","sources":[{"title":"Camia, A., Giuntoli, J., Jonsson, K., Robert, N., Cazzaniga, N., Jasinevičius, G., Avitabile, V., Grassi, G., Barredo Cano, J.I. and Mubareka, S., The use of woody biomass for energy production in the EU, EUR 30548 EN, Publications Office of the European U","url":"https://dfde.efi.int/db/dfde_app.php","role":"secondary","citation":"Camia, A., Giuntoli, J., Jonsson, K., Robert, N., Cazzaniga, N., Jasinevičius, G., Avitabile, V., Grassi, G., Barredo Cano, J.I. and Mubareka, S., The use of woody biomass for energy production in the EU, EUR 30548 EN, Publications Office of the European U (accessed through the DFDE, EFI; original not seen)"}]},{"id":"BB0235","name":"Germany national bark-beetle damage (2015–2019), JRC series via DFDE","country":"Germany","y0":2015,"y1":2019,"lat":51.2,"lon":10.4,"sev":"extreme","conf":"medium","note":"National yearly bark-beetle damage 2015-2019 (thousand m3): 2015: 3,300 thousand m3; 2016: 4,600 thousand m3; 2017: 6,000 thousand m3; 2018: 11,300 thousand m3; 2019: 31,702 thousand m3 (Camia et al. 2021, JRC; compiled in DFDE).","sources":[{"title":"Camia, A., Giuntoli, J., Jonsson, K., Robert, N., Cazzaniga, N., Jasinevičius, G., Avitabile, V., Grassi, G., Barredo Cano, J.I. and Mubareka, S., The use of woody biomass for energy production in the EU, EUR 30548 EN, Publications Office of the European U","url":"https://dfde.efi.int/db/dfde_app.php","role":"secondary","citation":"Camia, A., Giuntoli, J., Jonsson, K., Robert, N., Cazzaniga, N., Jasinevičius, G., Avitabile, V., Grassi, G., Barredo Cano, J.I. and Mubareka, S., The use of woody biomass for energy production in the EU, EUR 30548 EN, Publications Office of the European U (accessed through the DFDE, EFI; original not seen)"}]},{"id":"BB0236","name":"Austria national bark-beetle damage (2015–2019), JRC series via DFDE","country":"Austria","y0":2015,"y1":2019,"lat":47.3,"lon":13.5,"sev":"extreme","conf":"medium","note":"National yearly bark-beetle damage 2015-2019 (thousand m3): 2015: 2,291 thousand m3; 2016: 3,005 thousand m3; 2017: 3,338 thousand m3; 2018: 4,379 thousand m3; 2019: 4,261 thousand m3 (Camia et al. 2021, JRC; compiled in DFDE).","sources":[{"title":"Camia, A., Giuntoli, J., Jonsson, K., Robert, N., Cazzaniga, N., Jasinevičius, G., Avitabile, V., Grassi, G., Barredo Cano, J.I. and Mubareka, S., The use of woody biomass for energy production in the EU, EUR 30548 EN, Publications Office of the European U","url":"https://dfde.efi.int/db/dfde_app.php","role":"secondary","citation":"Camia, A., Giuntoli, J., Jonsson, K., Robert, N., Cazzaniga, N., Jasinevičius, G., Avitabile, V., Grassi, G., Barredo Cano, J.I. and Mubareka, S., The use of woody biomass for energy production in the EU, EUR 30548 EN, Publications Office of the European U (accessed through the DFDE, EFI; original not seen)"}]},{"id":"BB0237","name":"Slovenia national bark-beetle damage (2015–2019), JRC series via DFDE","country":"Slovenia","y0":2015,"y1":2019,"lat":46.1,"lon":14.8,"sev":"extreme","conf":"medium","note":"National yearly bark-beetle damage 2015-2019 (thousand m3): 2015: 1,818 thousand m3; 2016: 2,321 thousand m3; 2017: 1,802 thousand m3; 2018: 738 thousand m3; 2019: 1,341 thousand m3 (Camia et al. 2021, JRC; compiled in DFDE).","sources":[{"title":"Camia, A., Giuntoli, J., Jonsson, K., Robert, N., Cazzaniga, N., Jasinevičius, G., Avitabile, V., Grassi, G., Barredo Cano, J.I. and Mubareka, S., The use of woody biomass for energy production in the EU, EUR 30548 EN, Publications Office of the European U","url":"https://dfde.efi.int/db/dfde_app.php","role":"secondary","citation":"Camia, A., Giuntoli, J., Jonsson, K., Robert, N., Cazzaniga, N., Jasinevičius, G., Avitabile, V., Grassi, G., Barredo Cano, J.I. and Mubareka, S., The use of woody biomass for energy production in the EU, EUR 30548 EN, Publications Office of the European U (accessed through the DFDE, EFI; original not seen)"}]},{"id":"BB0238","name":"Finland national bark-beetle damage (2012–2017), JRC series via DFDE","country":"Finland","y0":2012,"y1":2017,"lat":61.5,"lon":26.0,"sev":"high","conf":"medium","note":"National yearly bark-beetle damage 2012-2017 (thousand m3): 2012: 282 thousand m3; 2013: 633 thousand m3; 2014: 432 thousand m3; 2015: 364 thousand m3; 2016: 264 thousand m3; 2017: 123 thousand m3 (Camia et al. 2021, JRC; compiled in DFDE).","sources":[{"title":"Camia, A., Giuntoli, J., Jonsson, K., Robert, N., Cazzaniga, N., Jasinevičius, G., Avitabile, V., Grassi, G., Barredo Cano, J.I. and Mubareka, S., The use of woody biomass for energy production in the EU, EUR 30548 EN, Publications Office of the European U","url":"https://dfde.efi.int/db/dfde_app.php","role":"secondary","citation":"Camia, A., Giuntoli, J., Jonsson, K., Robert, N., Cazzaniga, N., Jasinevičius, G., Avitabile, V., Grassi, G., Barredo Cano, J.I. and Mubareka, S., The use of woody biomass for energy production in the EU, EUR 30548 EN, Publications Office of the European U (accessed through the DFDE, EFI; original not seen)"}]},{"id":"BB0239","name":"Lithuania Ips typographus damage (2004), via DFDE","country":"Lithuania","y0":2004,"y1":2004,"lat":55.3,"lon":23.9,"sev":"high","conf":"medium","note":"225,000 m3 in 2004, against 44,000-75,000 m3 a year in 2005-2014 (Camia et al. 2021, JRC; compiled in DFDE).","sources":[{"title":"Camia, A., Giuntoli, J., Jonsson, K., Robert, N., Cazzaniga, N., Jasinevičius, G., Avitabile, V., Grassi, G., Barredo Cano, J.I. and Mubareka, S., The use of woody biomass for energy production in the EU, EUR 30548 EN, Publications Office of the European U","url":"https://dfde.efi.int/db/dfde_app.php","role":"secondary","citation":"Camia, A., Giuntoli, J., Jonsson, K., Robert, N., Cazzaniga, N., Jasinevičius, G., Avitabile, V., Grassi, G., Barredo Cano, J.I. and Mubareka, S., The use of woody biomass for energy production in the EU, EUR 30548 EN, Publications Office of the European U (accessed through the DFDE, EFI; original not seen)"}]},{"id":"BB0240","name":"Czechia national bark-beetle damage (1983–1988), Patacca et al. 2023 series","country":"Czechia","y0":1983,"y1":1988,"lat":49.8,"lon":15.5,"sev":"extreme","conf":"medium","note":"Reported national bark-beetle damage rose from 400,000 m3 in 1982 to 1.14 million m3 in 1983 and 1.58 million m3 in 1984, staying above 800,000 m3 until 1988.","sources":[{"title":"Significant increase in natural disturbance impacts on European forests since 1950","url":"https://github.com/MarcoPatacca/DFDE_manuscript_Patacca_etal2022","role":"secondary","citation":"Patacca M. et al. (2023). Significant increase in natural disturbance impacts on European forests since 1950. Global Change Biology 29 (doi 10.1111/gcb.16531); cleaned series in the public workbook disturbance_timeseries_database.xlsx."}]},{"id":"BB0241","name":"Slovakia national bark-beetle damage (1983–1989), Patacca et al. 2023 series","country":"Slovakia","y0":1983,"y1":1989,"lat":48.7,"lon":19.5,"sev":"high","conf":"medium","note":"Reported national bark-beetle damage of 586,000 m3 in 1983 and 736,000 m3 in 1984, 290,000-500,000 m3 a year until 1989.","sources":[{"title":"Significant increase in natural disturbance impacts on European forests since 1950","url":"https://github.com/MarcoPatacca/DFDE_manuscript_Patacca_etal2022","role":"secondary","citation":"Patacca M. et al. (2023). Significant increase in natural disturbance impacts on European forests since 1950. Global Change Biology 29 (doi 10.1111/gcb.16531); cleaned series in the public workbook disturbance_timeseries_database.xlsx."}]},{"id":"BB0242","name":"Germany national bark-beetle damage (2001–2005), Patacca et al. 2023 series","country":"Germany","y0":2001,"y1":2005,"lat":51.2,"lon":10.4,"sev":"extreme","conf":"medium","note":"Reported national bark-beetle damage of 2.85 million m3 in 2001 and 4.5-4.7 million m3 in 2003-2004, the wave after storm Lothar and the drought of 2003.","sources":[{"title":"Significant increase in natural disturbance impacts on European forests since 1950","url":"https://github.com/MarcoPatacca/DFDE_manuscript_Patacca_etal2022","role":"secondary","citation":"Patacca M. et al. (2023). Significant increase in natural disturbance impacts on European forests since 1950. Global Change Biology 29 (doi 10.1111/gcb.16531); cleaned series in the public workbook disturbance_timeseries_database.xlsx."}]},{"id":"BB0243","name":"Sudetes (Poland) bark-beetle outbreak continuing (1988–1997), Grodzki 1999 via DFDE","country":"Poland","y0":1988,"y1":1997,"lat":50.8,"lon":15.7,"sev":"high","conf":"medium","note":"Infested volume (thousand m3) per year: 1988: 921; 1989: 698; 1990: 417; 1991: 349; 1992: 272; 1993: 349; 1994: 388; 1995: 456; 1996: 349; 1997: 359; follows the 1980-1987 peak of up to 1.47 million m3.","sources":[{"title":"Grodzki, W. (1999). Phytosanitary situation in Polish mountain forests in 1998 and forecast for 1999. Methodology of forest insect and disease survey in Central Europe. In Proceedings of the Second Workshop of the IUFRO Working Party (Vol. 7, No. 10, pp. 2","url":"https://dfde.efi.int/db/dfde_app.php","role":"secondary","citation":"Grodzki, W. (1999). Phytosanitary situation in Polish mountain forests in 1998 and forecast for 1999. Methodology of forest insect and disease survey in Central Europe. In Proceedings of the Second Workshop of the IUFRO Working Party (Vol. 7, No. 10, pp. 2 (accessed through the DFDE, EFI; original not seen)"}]},{"id":"BB0244","name":"Polish Carpathians bark-beetle damage (1980–1997), Grodzki 1999 via DFDE","country":"Poland","y0":1980,"y1":1997,"lat":49.6,"lon":20.5,"sev":"high","conf":"medium","note":"DFDE gives two Carpathian series (thousand m3 per year, larger/smaller value): 1980: 281/136; 1981: 272/116; 1982: 582/272; 1983: 417/262; 1984: 708/524; 1985: 621/485; 1986: 543/465; 1987: 436/378; 1988: 427/310; 1989: 368/291; 1990: 330/2","sources":[{"title":"Grodzki, W. (1999). Phytosanitary situation in Polish mountain forests in 1998 and forecast for 1999. Methodology of forest insect and disease survey in Central Europe. In Proceedings of the Second Workshop of the IUFRO Working Party (Vol. 7, No. 10, pp. 2","url":"https://dfde.efi.int/db/dfde_app.php","role":"secondary","citation":"Grodzki, W. (1999). Phytosanitary situation in Polish mountain forests in 1998 and forecast for 1999. Methodology of forest insect and disease survey in Central Europe. In Proceedings of the Second Workshop of the IUFRO Working Party (Vol. 7, No. 10, pp. 2 (accessed through the DFDE, EFI; original not seen)"}]}];
    const DOM = {
      citationBox: document.getElementById("citationBox"),
      citationForLayer: document.getElementById("citationForLayer"),
      overlayCardStrip: document.getElementById("overlayCardStrip"),
      overlayDescription: document.getElementById("overlayDescription"),
      hexControlPanel: document.getElementById("hexControlPanel"),
      hexControlTitle: document.getElementById("hexControlTitle"),
      phenipsPanel: document.getElementById("phenipsPanel"),
      phenipsScenSelector: document.getElementById("phenipsScenSelector"),
      phenisPeriodSelector: document.getElementById("phenisPeriodSelector"),
      phenipsVarSelector: document.getElementById("phenipsVarSelector"),
      phenipsViewSelector: document.getElementById("phenipsViewSelector"),
      phenipsExpSelector: document.getElementById("phenipsExpSelector"),
      phenipsExpRow: document.getElementById("phenipsExpRow"),
      phenipsYearRow: document.getElementById("phenipsYearRow"),
      phenipsYearSlider: document.getElementById("phenipsYearSlider"),
      phenipsYearVal: document.getElementById("phenipsYearVal"),
      phenipsNote: document.getElementById("phenipsNote"),
      forecastPanel: document.getElementById("forecastPanel"),
      forecastModelSelector: document.getElementById("forecastModelSelector"),
      forecastPeriodSelector: document.getElementById("forecastPeriodSelector"),
      forecastViewSelector: document.getElementById("forecastViewSelector"),
      forecastViewDeltaOption: document.getElementById("forecastViewDeltaOption"),
      forecastViewShareOption: document.getElementById("forecastViewShareOption"),
      forecastViewExtrapOption: document.getElementById("forecastViewExtrapOption"),
      legendFlagMarker: document.getElementById("legendFlagMarker"),
      legendFlagNote: document.getElementById("legendFlagNote"),
      forecastThresholdSelector: document.getElementById("forecastThresholdSelector"),
      forecastMetricsGrid: document.getElementById("forecastMetricsGrid"),
      forecastConfusionMatrix: document.getElementById("forecastConfusionMatrix"),
      forecastMetricsNote: document.getElementById("forecastMetricsNote"),
      forecastModelNote: document.getElementById("forecastModelNote"),
      windPanel: document.getElementById("windPanel"),
      windPeriodSelector: document.getElementById("windPeriodSelector"),
      windGcmSelector: document.getElementById("windGcmSelector"),
      droughtFuturePanel: document.getElementById("droughtFuturePanel"),
      droughtFuturePeriodSelector: document.getElementById("droughtFuturePeriodSelector"),
      droughtRefPanel: document.getElementById("droughtRefPanel"),
      droughtRefMetricSelector: document.getElementById("droughtRefMetricSelector"),
      fireWeatherPanel: document.getElementById("fireWeatherPanel"),
      fireViewSelector: document.getElementById("fireViewSelector"),
      fireScenSelector: document.getElementById("fireScenSelector"),
      fireInfoDesc: document.getElementById("fireInfoDesc"),
      fireInfoValues: document.getElementById("fireInfoValues"),
      droughtRefInfoDesc: document.getElementById("droughtRefInfoDesc"),
      droughtRefInfoValues: document.getElementById("droughtRefInfoValues"),
      droughtRefInfoLegend: document.getElementById("droughtRefInfoLegend"),
      riskSlider: document.getElementById("riskSlider"),
      riskValue: document.getElementById("riskValue"),
      legendThreshMarker: document.getElementById("legendThreshMarker"),
      shareBtn: document.getElementById("shareBtn"),
      splitBtn: document.getElementById("splitBtn"),
      detailModeBtn: document.getElementById("detailModeBtn"),
      forecastModelField: document.getElementById("forecastModelField"),
      forecastMetricsBox: document.getElementById("forecastMetricsBox"),
      windGcmField: document.getElementById("windGcmField"),
      pilotActionsBtn: document.getElementById("pilotActionsBtn"),
      dssPillarPilotActions: document.getElementById("dssPillarPilotActions"),
      dssPillarPilotActionsStatus: document.getElementById("dssPillarPilotActionsStatus"),
      dssPillarNatura2000: document.getElementById("dssPillarNatura2000"),
      dssPillarNatura2000Status: document.getElementById("dssPillarNatura2000Status"),
      dssPillarDetail: document.getElementById("dssPillarDetail"),
      dssPillarRestoration: document.getElementById("dssPillarRestoration"),
      dssPillarRestorationStatus: document.getElementById("dssPillarRestorationStatus"),
      dssPillarRiskDrivers: document.getElementById("dssPillarRiskDrivers"),
      dssPillarRiskDriversStatus: document.getElementById("dssPillarRiskDriversStatus"),
      layerRibbon: document.querySelector(".layer-ribbon"),
      mapCaptionLeft: document.getElementById("mapCaptionLeft"),
      mapCaptionRight: document.getElementById("mapCaptionRight"),
      yearSlider: document.getElementById("yearSlider"),
      yearValue: document.getElementById("yearValue"),
      toggleOverlayBtn: document.getElementById("toggleOverlayBtn"),
      fitOverlayBtn: document.getElementById("fitOverlayBtn"),
      legendBar: document.querySelector(".legend-bar"),
      legendForLayer: document.getElementById("legendForLayer"),
      legendMin: document.getElementById("legendMin"),
      legendMax: document.getElementById("legendMax"),
      legendSwatches: document.getElementById("legendSwatches"),
      rcpSelector: document.getElementById("rcpSelector")
    };


    // ====================================================================
    // OVERLAY REGISTRY HELPERS (visibility, panes, image layers)
    // ====================================================================
    function getOverlayIds() {
      return Object.keys(APP.layerInstances);
    }

    function getTopVisibleOverlayId() {
      if (!APP.visibleOverlays.length) {
        return APP.focusedOverlay || "bbtl";
      }
      return APP.visibleOverlays[APP.visibleOverlays.length - 1];
    }

    function isOverlayVisible(overlayId) {
      return APP.visibleOverlays.includes(overlayId);
    }

    // Same as getTopVisibleOverlayId, but for an arbitrary state bundle (the
    // live APP for the active split-screen side, or a parked APP.panels[...]
    // snapshot for the inactive side) -- lets click-info be computed for
    // either side without switching which one is "active".
    function getTopVisibleOverlayIdFor(state) {
      if (!state.visibleOverlays || !state.visibleOverlays.length) {
        return state.focusedOverlay || "bbtl";
      }
      return state.visibleOverlays[state.visibleOverlays.length - 1];
    }

    function syncOverlayVisibilityFlag() {
      APP.overlayVisible = APP.visibleOverlays.length > 0;
    }

    function showOverlay(overlayId) {
      APP.visibleOverlays = APP.visibleOverlays.filter(function(id) {
        return id !== overlayId;
      });
      APP.visibleOverlays.push(overlayId);
      syncOverlayVisibilityFlag();
      updateOverlayPaneOrder();
      syncCatalogueEventsFromPastDisturbance();
    }

    function hideOverlay(overlayId) {
      APP.visibleOverlays = APP.visibleOverlays.filter(function(id) {
        return id !== overlayId;
      });
      syncOverlayVisibilityFlag();
      updateOverlayPaneOrder();
      syncCatalogueEventsFromPastDisturbance();
    }

    function getVisibleOverlay() {
      return OVERLAY_CATALOG[getTopVisibleOverlayId()] || OVERLAY_CATALOG.bbtl;
    }

    function getFocusedOverlay() {
      return OVERLAY_CATALOG[APP.focusedOverlay] || OVERLAY_CATALOG.bbtl;
    }

    function getOverlayPaneName(overlayId) {
      return "pane-" + overlayId;
    }

    function getHexLayer(overlayId) {
      return APP.layerInstances[overlayId || getTopVisibleOverlayId()] || null;
    }

    function updateOverlayPaneOrder() {
      APP.visibleOverlays.forEach(function(overlayId, index) {
        const pane = APP.map && APP.map.getPane(getOverlayPaneName(overlayId));
        if (pane) {
          pane.style.zIndex = String(300 + index * 20);
        }
        const layer = APP.layerInstances[overlayId];
        if (layer && layer.bringToFront) {
          layer.bringToFront();
        }
      });
    }

    // Start a fresh layer load and make any older in-flight fetch harmless.
    function beginOverlayLoad(overlayId) {
      APP.overlayLoadTokens[overlayId] += 1;
      clearLayerInstance(overlayId);
      return APP.overlayLoadTokens[overlayId];
    }

    function isCurrentOverlayLoad(overlayId, requestToken) {
      return requestToken === APP.overlayLoadTokens[overlayId];
    }

    // Crossfade: when an image layer is replaced (period, scenario or variable change), the old
    // image stays under the new one until the new one has loaded, then the two swap by opacity.
    const FADE_MS = 400;
    const LAYER_GHOSTS = {};
    const prefersReducedMotion = function() {
      return window.matchMedia && window.matchMedia("(prefers-reduced-motion: reduce)").matches;
    };

    function dropGhosts(overlayId, fade) {
      (LAYER_GHOSTS[overlayId] || []).forEach(function(g) {
        clearTimeout(g.timer);
        const el = g.layer.getElement && g.layer.getElement();
        const remove = function() { if (g.map.hasLayer(g.layer)) g.map.removeLayer(g.layer); };
        if (fade && el) { el.style.transition = "opacity " + FADE_MS + "ms ease"; el.style.opacity = 0; setTimeout(remove, FADE_MS + 50); }
        else remove();
      });
      LAYER_GHOSTS[overlayId] = [];
    }

    function clearLayerInstance(overlayId) {
      const layer = APP.layerInstances[overlayId];
      if (!layer) return;
      APP.layerInstances[overlayId] = null;
      if (prefersReducedMotion() || !layer.getElement || !layer.getElement()) { APP.map.removeLayer(layer); return; }
      // Keep it as a ghost; if no replacement claims it within 120 ms it was a plain "hide", so fade it out.
      const g = { layer: layer, map: APP.map, timer: null };
      g.timer = setTimeout(function() {
        LAYER_GHOSTS[overlayId] = (LAYER_GHOSTS[overlayId] || []).filter(function(x) { return x !== g; });
        const el = layer.getElement();
        if (el) { el.style.transition = "opacity " + FADE_MS + "ms ease"; el.style.opacity = 0; }
        setTimeout(function() { if (g.map.hasLayer(layer)) g.map.removeLayer(layer); }, FADE_MS + 50);
      }, 120);
      (LAYER_GHOSTS[overlayId] = LAYER_GHOSTS[overlayId] || []).push(g);
    }

    // Crisp pixels when zoomed in past the image's own resolution (browsers otherwise smear them into a blur).
    function sharpenImageLayer(layer) {
      const el = layer && layer.getElement && layer.getElement(), map = layer && layer._map;
      if (!el || !map || !el.naturalWidth) return;
      const b = layer.getBounds();
      const w = map.latLngToLayerPoint(b.getNorthEast()).x - map.latLngToLayerPoint(b.getSouthWest()).x;
      el.style.imageRendering = (w / el.naturalWidth > 1.3) ? "pixelated" : "auto";
    }

    function addImageOverlayLayer(overlayId, imagePath, bounds) {
      const target = getLayerOpacity(overlayId);
      const ghosts = LAYER_GHOSTS[overlayId] || [];
      ghosts.forEach(function(g) { clearTimeout(g.timer); });
      const fade = !prefersReducedMotion() && ghosts.length > 0;
      const layer = L.imageOverlay(imagePath, [
        [bounds.south, bounds.west],
        [bounds.north, bounds.east]
      ], {
        opacity: fade ? 0 : target,
        pane: getOverlayPaneName(overlayId)
      }).addTo(APP.map);
      APP.layerInstances[overlayId] = layer;
      layer.on("load", function() { sharpenImageLayer(layer); });
      if (!APP.map._sharpenBound) {
        APP.map._sharpenBound = true;
        APP.map.on("zoomend", function() {
          Object.keys(APP.layerInstances).forEach(function(k) { if (APP.layerInstances[k] && APP.layerInstances[k].getBounds && APP.layerInstances[k]._image) sharpenImageLayer(APP.layerInstances[k]); });
        });
      }
      if (fade) {
        const el = layer.getElement();
        let done = false;
        const swap = function() {
          if (done) return; done = true;
          el.style.transition = "opacity " + FADE_MS + "ms ease";
          if (APP.layerInstances[overlayId] === layer) el.style.opacity = getLayerOpacity(overlayId);
          dropGhosts(overlayId, true);
        };
        el.addEventListener("load", swap, { once: true });
        el.addEventListener("error", swap, { once: true });
        setTimeout(swap, 6000);
      }
      return layer;
    }

    function getLayerOpacity(overlayId) {
      const opacity = Number(APP.layerOpacity[overlayId]);
      return isFinite(opacity) ? opacity : 0.85;
    }

    // No visible status line in the UI anymore -- kept as a console-only
    // hook so load/error lifecycle is still inspectable while debugging.
    function setStatus(text) {
      console.debug("[status]", text);
    }

    function formatForecastThreshold(threshold) {
      return Number(threshold).toFixed(2);
    }


    // ====================================================================
    // BARK-BEETLE FORECAST BUNDLE (manifest, entries, metrics)
    // ====================================================================
    function getForecastEntryKey() {
      return APP.currentForecastModel + "::" + APP.currentForecastPeriod;
    }

    function getForecastEntry() {
      if (!APP.forecastManifest || !APP.forecastManifest.entries) {
        return null;
      }
      return APP.forecastManifest.entries[getForecastEntryKey()] || null;
    }

    function getForecastOverlayPaths(entry) {
      const tifPath = APP.currentForecastView === "probability" ? entry.ensemble_probability_tif
        : APP.currentForecastView === "delta" ? entry.delta_from_near_tif
        : APP.currentForecastView === "share05" ? entry.share05_tif
        : APP.currentForecastView === "extrapolation" ? entry.extrapolation_tif
        : entry.threshold_visible_tifs[APP.currentForecastThreshold];

      if (!tifPath) {
        return null;
      }

      return {
        tifPath: tifPath,
        imagePath: tifPath.replace(/\.tif$/i, ".png") + "?v=" + ASSET_VERSION,
        metaPath: tifPath.replace(/\.tif$/i, "_meta.json")
      };
    }

    async function ensurePastOverlayMetaLoaded() {
      const overlay = OVERLAY_CATALOG.bbtpast;
      if (overlay.meta.bounds && overlay.meta.legend) {
        return overlay.meta;
      }

      const response = await fetch(overlay.meta.metaUrl, { cache: "no-store" });
      if (!response.ok) {
        throw new Error("HTTP " + response.status + " for " + overlay.meta.metaUrl);
      }

      const overlayMeta = await response.json();
      overlay.meta.bounds = overlayMeta.bounds;
      overlay.meta.legend = overlayMeta.legend;
      if (overlayMeta.image) {
        overlay.meta.image = "webmap_data/" + overlayMeta.image;
      }
      return overlay.meta;
    }

    async function ensureForecastManifestLoaded() {
      if (APP.forecastManifest) {
        return APP.forecastManifest;
      }

      const response = await fetch(FORECAST_MANIFEST_URL, { cache: "no-store" });
      if (!response.ok) {
        throw new Error("HTTP " + response.status + " for forecast manifest");
      }

      APP.forecastManifest = await response.json();
      return APP.forecastManifest;
    }

    async function ensureForecastMetricsLoaded(modelKey) {
      if (APP.forecastMetricsCache[modelKey]) {
        return APP.forecastMetricsCache[modelKey];
      }

      const modelMeta = APP.forecastManifest && APP.forecastManifest.models ? APP.forecastManifest.models[modelKey] : null;
      if (!modelMeta) {
        throw new Error("Forecast model metadata missing for " + modelKey);
      }

      const response = await fetch(modelMeta.metrics_file, { cache: "no-store" });
      if (!response.ok) {
        throw new Error("HTTP " + response.status + " for forecast metrics");
      }

      APP.forecastMetricsCache[modelKey] = await response.json();
      return APP.forecastMetricsCache[modelKey];
    }

    function updateForecastPanelVisibility() {
      DOM.forecastPanel.hidden = APP.focusedOverlay !== "bbtforecast";
      DOM.phenipsPanel.hidden = APP.focusedOverlay !== "phenips";
      DOM.windPanel.hidden = APP.focusedOverlay !== "wind";
      DOM.droughtFuturePanel.hidden = APP.focusedOverlay !== "droughtFuture";
      DOM.droughtRefPanel.hidden = APP.focusedOverlay !== "droughtRef";
      DOM.fireWeatherPanel.hidden = APP.focusedOverlay !== "fireWeather";
      renderDroughtRefInfo();
      renderFireInfo();
    }

    const FIRE_INFO = {
      change: {
        desc: "Relative change of the mean summer (June-August) seasonal severity rating compared with 1981-2005. The rating rises steeply with the Fire Weather Index, so it reflects how difficult a fire would be to contain.",
        values: "Warm colours = more severe fire weather than in 1981-2005 (dark red = more than doubled), blue = less severe. Pale colours = change smaller than the year-to-year noise of the model. Grey = fire weather too mild for a change to be computed (mostly northern Europe)."
      },
      level: {
        desc: "Mean summer (June-August) Fire Weather Index, in the standard fire-danger classes.",
        values: "Dark colours = more dangerous fire weather. Compare the 1981-2005 level with a future period to see where the Mediterranean (high level already) and Central Europe (large relative increase) stand."
      }
    };
    function renderFireInfo() {
      if (!DOM.fireInfoDesc) return;
      const info = FIRE_INFO[APP.currentFireView] || FIRE_INFO.change;
      DOM.fireInfoDesc.textContent = info.desc;
      DOM.fireInfoValues.textContent = info.values;
    }

    // Fire image key. "change" has no historical image, so it falls back to level.
    function fireImageKey(state) {
      const scen = state.currentFireScen;
      const view = scen === "historical" ? "level" : state.currentFireView;
      return view + "::" + scen;
    }

    // Per-metric documentation (Stefan Ebner, BFW, 2026-08-31 notes): each
    // reference-period metric gets its own "what/how/legend" text, unlike
    // the future SPEI layer where that text is fixed regardless of the
    // selected scenario/period.
    const DROUGHT_REF_INFO = {
      n_events_2yr: {
        desc: "Number of multi-year drought events detected during the historical reference period 1981–2010.",
        values: "The values indicate the number of multi-year drought events detected in each pixel during 1981–2010. A MYD event is defined as SPEI-12 ≤ −1.1 occurring in at least two consecutive years.",
        legend: "Higher values indicate a greater number of separate multi-year drought events during the historical reference period."
      },
      max_duration: {
        desc: "Maximum duration of a single multi-year drought event detected during the historical reference period 1981–2010.",
        values: "The values indicate the duration, in consecutive years, of the longest multi-year drought event detected in each pixel during 1981–2010. A MYD event is defined as SPEI-12 ≤ −1.1 occurring in at least two consecutive years.",
        legend: "Higher values indicate a longer-lasting individual multi-year drought event."
      },
      total_duration: {
        desc: "Cumulative duration of all multi-year drought events detected during the historical reference period 1981–2010.",
        values: "The values indicate the total number of years affected by all detected multi-year drought events in each pixel during 1981–2010. A MYD event is defined as SPEI-12 ≤ −1.1 occurring in at least two consecutive years.",
        legend: "Higher values indicate a greater cumulative duration of multi-year drought conditions during the historical reference period."
      }
    };

    function renderDroughtRefInfo() {
      if (!DOM.droughtRefInfoDesc) return;
      const info = DROUGHT_REF_INFO[APP.currentDroughtRefMetric] || DROUGHT_REF_INFO.n_events_2yr;
      DOM.droughtRefInfoDesc.textContent = info.desc;
      DOM.droughtRefInfoValues.textContent = info.values;
      DOM.droughtRefInfoLegend.textContent = info.legend;
    }


    // ====================================================================
    // PHENIPS & WIND RASTER LAYERS
    // ====================================================================
    const PHENIPS_V5_BOUNDS = { south: 34.4781, west: -10.6989, north: 71.4273, east: 38.6006 };
    const PHENIPS_YEARS = []; for (let y = 2024; y >= 1981; y--) PHENIPS_YEARS.push(y);

    // Observed (v5, CHELSA) views are named <var>_<exposure>_<period>; "projection" views keep the old files.
    function phenipsIsObs(state) { return (state || APP).currentPhenipsScen === "obs"; }

    function getPhenipsV5Name(s) {
      const v = s.currentPhenipsVar, e = s.currentPhenipsExp, p = s.currentPhenisPeriod;
      if (v === "share") return "phenips_share_years_ge1_generation";
      const base = "phenips_" + (v === "gdd" ? "dd" : "gen") + "_" + e;
      if (s.currentPhenipsView === "diff") return base + "_change";
      if (p === "m1") return base + "_mean_1981_2010";
      if (p === "m2") return base + "_mean_2011_2024";
      return base + "_" + p.slice(1);
    }

    function getPhenipsImagePath(periodOverride) {
      if (phenipsIsObs()) return "webmap_data/phenips_v5/" + getPhenipsV5Name(APP) + ".webp?v=" + ASSET_VERSION;
      const p    = periodOverride || APP.currentPhenisPeriod;
      const v    = APP.currentPhenipsVar;
      const view = APP.currentPhenipsView;
      const s    = APP.currentPhenipsScen;
      let file;
      if (p === "hist") {
        // Historical is scenario-independent
        file = "phenips_" + v + "_hist.png";
      } else if (view === "contrast" && v === "gen") {
        file = "phenips_gen_" + s + "_contrast.png";
      } else if (view === "diff") {
        file = "phenips_" + v + "_" + s + "_diff_" + p + ".png";
      } else {
        file = "phenips_" + v + "_" + s + "_" + p + ".png";
      }
      return "webmap_data/phenips/" + file + "?v=20260531";
    }

    // Keep the five PHENIPS selectors consistent with each other (observed vs projection options).
    function syncPhenipsControls() {
      const obs = phenipsIsObs();
      const opt = function(v, t) { return '<option value="' + v + '">' + t + '</option>'; };
      if (obs) {
        if (APP.currentPhenipsVar !== "gdd" && APP.currentPhenipsVar !== "gen" && APP.currentPhenipsVar !== "share") APP.currentPhenipsVar = "gdd";
        if (APP.currentPhenipsExp !== "sun") APP.currentPhenipsExp = "shade";
        if (APP.currentPhenipsView === "contrast") APP.currentPhenipsView = "abs";
        if (!/^(m1|m2|y(19[89]\d|20[0-2]\d))$/.test(APP.currentPhenisPeriod) || (APP.currentPhenisPeriod > "y2024" && APP.currentPhenisPeriod[0] === "y")) APP.currentPhenisPeriod = "m2";
        DOM.phenisPeriodSelector.innerHTML = opt("m1", "Mean 1981–2010") + opt("m2", "Mean 2011–2024")
          + PHENIPS_YEARS.map(function(y) { return opt("y" + y, "Single year " + y); }).join("");
      } else {
        if (APP.currentPhenipsVar === "share") APP.currentPhenipsVar = "gdd";
        if (["hist", "near", "mid", "end"].indexOf(APP.currentPhenisPeriod) < 0) APP.currentPhenisPeriod = "hist";
        DOM.phenisPeriodSelector.innerHTML = opt("hist", "Historical 1985–2005") + opt("near", "Near-future 2005–2040")
          + opt("mid", "Mid-century 2041–2070") + opt("end", "End-century 2071–2100");
      }
      DOM.phenipsVarSelector.querySelector('option[value="share"]').disabled = !obs;
      DOM.phenipsViewSelector.querySelector('option[value="contrast"]').disabled = obs;
      DOM.phenipsViewSelector.querySelector('option[value="diff"]').textContent = obs ? "Change 2011–2024 vs. 1981–2010" : "Change vs. historical";
      DOM.phenipsExpRow.hidden = !obs || APP.currentPhenipsVar === "share";
      DOM.phenipsYearRow.hidden = !obs || APP.currentPhenipsVar === "share" || APP.currentPhenipsView === "diff";
      const yr = /^y\d{4}$/.test(APP.currentPhenisPeriod) ? Number(APP.currentPhenisPeriod.slice(1)) : (APP.phenipsLastYear || 2024);
      DOM.phenipsYearSlider.value = yr;
      DOM.phenipsYearVal.textContent = /^y\d{4}$/.test(APP.currentPhenisPeriod) ? String(yr) : "(pick a year)";
      DOM.phenipsViewSelector.disabled = obs && APP.currentPhenipsVar === "share";
      DOM.phenisPeriodSelector.disabled = obs && (APP.currentPhenipsVar === "share" || APP.currentPhenipsView === "diff");
      DOM.phenipsScenSelector.value = APP.currentPhenipsScen;
      DOM.phenipsVarSelector.value = APP.currentPhenipsVar;
      DOM.phenipsExpSelector.value = APP.currentPhenipsExp;
      DOM.phenisPeriodSelector.value = APP.currentPhenisPeriod;
      DOM.phenipsViewSelector.value = APP.currentPhenipsView;
      DOM.phenipsNote.textContent = obs
        ? "Observed CHELSA climate, each year 1981–2024 · bark-temperature corrected degree-days · shaded = no radiation term, sun-exposed = with radiation · 1 km"
        : "Ensemble mean · 3 GCMs · bark-temp corrected (Baier et al. 2007) · max 3 generations (photoperiod-gated)";
    }

    async function ensurePhenipsV5Manifest() {
      if (APP.phenipsV5) return APP.phenipsV5;
      const res = await fetch("webmap_data/phenips_v5/manifest.json?v=" + ASSET_VERSION);
      if (!res.ok) throw new Error("HTTP " + res.status + " for PHENIPS v5 manifest");
      APP.phenipsV5 = await res.json();
      return APP.phenipsV5;
    }

    async function ensurePhenipsDataLoaded() {
      if (APP.phenipsData) return APP.phenipsData;
      const res = await fetch(PHENIPS_DATA_URL, { cache: "no-store" });
      if (!res.ok) throw new Error("HTTP " + res.status + " for PHENIPS data");
      APP.phenipsData = await res.json();
      return APP.phenipsData;
    }

    function samplePhenipsAt(latlng) {
      if (!APP.phenipsData) return null;
      let best = null, bestDist = Infinity;
      APP.phenipsData.forEach(function(pt) {
        const d = Math.pow(pt.lat - latlng.lat, 2) + Math.pow(pt.lng - latlng.lng, 2);
        if (d < bestDist) { bestDist = d; best = pt; }
      });
      return (best && bestDist < 0.5) ? best : null;
    }

    async function updatePhenipsLayer() {
      if (!isOverlayVisible("phenips")) return;
      clearLayerInstance("phenips");
      syncPhenipsControls();
      if (phenipsIsObs()) {
        addImageOverlayLayer("phenips", getPhenipsImagePath(), PHENIPS_V5_BOUNDS);
        ensurePhenipsV5Manifest().then(setLegendLabels).catch(function(e) { console.warn("PHENIPS manifest load failed:", e); });
      } else {
        addImageOverlayLayer("phenips", getPhenipsImagePath(), PHENIPS_BOUNDS);
        ensurePhenipsDataLoaded().catch(function(e) {
          console.warn("PHENIPS data load failed:", e);
        });
      }
      setLegendLabels();
    }

    async function ensureWindMetaLoaded() {
      const overlay = OVERLAY_CATALOG.wind;
      if (overlay.meta.images) {
        return overlay.meta;
      }
      const response = await fetch(WIND_DATA_URL, { cache: "no-store" });
      if (!response.ok) {
        throw new Error("HTTP " + response.status + " for " + WIND_DATA_URL);
      }
      const windMeta = await response.json();
      overlay.meta.images = windMeta.images;
      overlay.meta.gcms = windMeta.gcms;
      return overlay.meta;
    }

    // "mpi" only has near/far periods; historical is model-independent so it
    // always uses the default image regardless of the GCM selector.
    function currentWindImageKey() {
      const gcm = APP.currentWindPeriod === "historical" ? "default" : APP.currentWindGcm;
      return gcm === "mpi" ? APP.currentWindPeriod + "_mpi" : APP.currentWindPeriod;
    }

    async function updateWindLayer() {
      if (!isOverlayVisible("wind")) return;
      const requestToken = beginOverlayLoad("wind");
      try {
        const meta = await ensureWindMetaLoaded();
        if (!isCurrentOverlayLoad("wind", requestToken)) return;
        const key = currentWindImageKey();
        const image = meta.images[key] || meta.images[APP.currentWindPeriod] || meta.images.near;
        addRiskImageOverlay("wind", "webmap_data/" + image.file + "?v=" + ASSET_VERSION, image.bounds, requestToken);
        setLegendLabels();
        setStatus("");
      } catch (error) {
        console.error(error);
        setStatus("Wind overlay could not be loaded");
      }
    }

    // MPI-ESM1-2-HR only exists for near/far; historical always renders the
    // model-independent default, so lock the selector to match.
    function updateWindGcmAvailability() {
      if (!DOM.windGcmSelector) return;
      const isHistorical = APP.currentWindPeriod === "historical";
      DOM.windGcmSelector.disabled = isHistorical;
      if (isHistorical) APP.currentWindGcm = "default";
      DOM.windGcmSelector.value = APP.currentWindGcm;
    }

    // ====================================================================
    // DROUGHT LAYERS (Stefan Ebner, BFW) — classified, categorical rasters,
    // so no risk-threshold masking (addImageOverlayLayer, not addRiskImageOverlay).
    // ====================================================================
    async function ensureDroughtFutureMetaLoaded() {
      const overlay = OVERLAY_CATALOG.droughtFuture;
      if (overlay.meta.images) return overlay.meta;
      const response = await fetch(DROUGHT_FUTURE_DATA_URL, { cache: "no-store" });
      if (!response.ok) throw new Error("HTTP " + response.status + " for " + DROUGHT_FUTURE_DATA_URL);
      overlay.meta.images = (await response.json()).images;
      return overlay.meta;
    }

    async function ensureDroughtRefMetaLoaded() {
      const overlay = OVERLAY_CATALOG.droughtRef;
      if (overlay.meta.images) return overlay.meta;
      const response = await fetch(DROUGHT_REF_DATA_URL, { cache: "no-store" });
      if (!response.ok) throw new Error("HTTP " + response.status + " for " + DROUGHT_REF_DATA_URL);
      overlay.meta.images = (await response.json()).images;
      return overlay.meta;
    }

    async function ensureFireMetaLoaded() {
      const overlay = OVERLAY_CATALOG.fireWeather;
      if (overlay.meta.images) return overlay.meta;
      const response = await fetch(overlay.meta.metaUrl, { cache: "no-store" });
      if (!response.ok) throw new Error("HTTP " + response.status + " for " + overlay.meta.metaUrl);
      const json = await response.json();
      overlay.meta.images = json.images;
      overlay.meta.legends = json.legends;
      return overlay.meta;
    }

    async function updateFireWeatherLayer() {
      if (!isOverlayVisible("fireWeather")) return;
      const requestToken = beginOverlayLoad("fireWeather");
      setStatus("Loading fire-weather layer...");
      try {
        const meta = await ensureFireMetaLoaded();
        if (!isCurrentOverlayLoad("fireWeather", requestToken)) return;
        const image = meta.images[fireImageKey(APP)];
        if (!image) throw new Error("No fire image for " + fireImageKey(APP));
        addImageOverlayLayer("fireWeather", "webmap_data/" + image.file + "?v=" + ASSET_VERSION, image.bounds);
        updateOverlayPaneOrder();
        setLegendLabels();
        setStatus("");
      } catch (error) {
        console.error(error);
        setStatus("Fire-weather overlay could not be loaded");
      }
    }

    async function updateDroughtFutureLayer() {
      if (!isOverlayVisible("droughtFuture")) return;
      const requestToken = beginOverlayLoad("droughtFuture");
      setStatus("Loading future drought classes...");
      try {
        const meta = await ensureDroughtFutureMetaLoaded();
        if (!isCurrentOverlayLoad("droughtFuture", requestToken)) return;
        const image = meta.images[APP.currentDroughtFuturePeriod];
        if (!image) throw new Error("No drought image for " + APP.currentDroughtFuturePeriod);
        addImageOverlayLayer("droughtFuture", "webmap_data/" + image.file + "?v=" + ASSET_VERSION, image.bounds);
        updateOverlayPaneOrder();
        setLegendLabels();
        setStatus("");
      } catch (error) {
        console.error(error);
        setStatus("Drought future overlay could not be loaded");
      }
    }

    async function updateDroughtRefLayer() {
      if (!isOverlayVisible("droughtRef")) return;
      const requestToken = beginOverlayLoad("droughtRef");
      setStatus("Loading reference drought indicators...");
      try {
        const meta = await ensureDroughtRefMetaLoaded();
        if (!isCurrentOverlayLoad("droughtRef", requestToken)) return;
        const image = meta.images[APP.currentDroughtRefMetric];
        if (!image) throw new Error("No drought image for " + APP.currentDroughtRefMetric);
        addImageOverlayLayer("droughtRef", "webmap_data/" + image.file + "?v=" + ASSET_VERSION, image.bounds);
        updateOverlayPaneOrder();
        setLegendLabels();
        setStatus("");
      } catch (error) {
        console.error(error);
        setStatus("Drought reference overlay could not be loaded");
      }
    }

    function populateSelectOptions(select, values, selectedValue, labelFn) {
      select.innerHTML = "";
      values.forEach(function(value) {
        const option = document.createElement("option");
        option.value = String(value);
        option.textContent = labelFn ? labelFn(value) : String(value);
        if (String(value) === String(selectedValue)) {
          option.selected = true;
        }
        select.appendChild(option);
      });
    }

    function getMatchingForecastThresholds(manifest, metrics, entry) {
      const manifestThresholds = (manifest.thresholds || []).map(Number);
      const metricThresholds = (metrics.aggregate_by_cutoff || []).map(function(item) {
        return Number(item.cutoff);
      });
      const entryThresholds = Object.keys(entry.threshold_visible_tifs || {}).map(Number);

      return manifestThresholds.filter(function(value) {
        return metricThresholds.includes(value) && entryThresholds.includes(value);
      }).map(String);
    }

    function renderForecastMetricsGrid(items) {
      DOM.forecastMetricsGrid.innerHTML = "";
      items.forEach(function(item) {
        const cell = document.createElement("div");
        cell.className = "metric-item";
        cell.innerHTML = "<strong>" + item.label + "</strong>" + item.value;
        DOM.forecastMetricsGrid.appendChild(cell);
      });
    }

    // Some models (e.g. a single climatological-normal scenario with no
    // train/test evaluation) report "NA" instead of a real number here.
    function fmtMetricNum(v, decimals) {
      const n = Number(v);
      return isFinite(n) ? n.toFixed(decimals) : "N/A";
    }

    function renderForecastPerformance(row, cutoff, auc) {
      // AUC is threshold-independent (ranks predictions across ALL cutoffs,
      // not just the selected one) -- show it whenever the model has one,
      // even if the specific cutoff's row is missing/NA.
      const aucItem = isFinite(Number(auc)) ? [{ label: "AUC", value: fmtMetricNum(auc, 3) }] : [];

      if (!row) {
        renderForecastMetricsGrid(aucItem);
        renderForecastConfusionMatrix(null);
        DOM.forecastMetricsNote.textContent = "No matching historical evaluation row is available for cutoff " + formatForecastThreshold(cutoff) + ".";
        return;
      }

      renderForecastMetricsGrid(aucItem.concat([
        { label: "Precision", value: fmtMetricNum(row.precision, 3) },
        { label: "Recall", value: fmtMetricNum(row.recall, 3) },
        { label: "Specificity", value: fmtMetricNum(row.specificity, 3) },
        { label: "F1", value: fmtMetricNum(row.F1, 3) },
        { label: "Red fraction", value: fmtMetricNum(row.red_fraction, 3) },
        { label: "Cutoff", value: formatForecastThreshold(cutoff) }
      ]));
      renderForecastConfusionMatrix(row);
      DOM.forecastMetricsNote.textContent = isFinite(Number(row.precision))
        ? "Historical evaluation summary aggregated across available GCMs at cutoff " + formatForecastThreshold(cutoff) + "."
        : "No historical evaluation exists for this scenario (single climatological-normal snapshot, not a validated forecast).";
    }

    function renderForecastConfusionMatrix(row) {
      DOM.forecastConfusionMatrix.innerHTML = "";
      if (!row) {
        return;
      }
      // Only the older, currently-hidden CORDEX models have real TP/FP/FN/TN
      // counts -- every model actually selectable today has none, so a grid
      // of "N/A" boxes would show for every single one. Skip it cleanly
      // instead, rather than rendering cells that never have data.
      const hasCounts = ["TP", "FP", "FN", "TN"].some(function(k) { return isFinite(Number(row[k])); });
      if (!hasCounts) {
        return;
      }

      // Show each cell as a share of all classified pixels (the standard
      // way to read a confusion matrix at a glance), alongside the raw
      // count -- percentage alone would lose the actual scale.
      const total = ["TP", "FP", "FN", "TN"].reduce(function(sum, k) {
        const n = Number(row[k]);
        return sum + (isFinite(n) ? n : 0);
      }, 0);
      const cellHtml = function(label, key) {
        const n = Number(row[key]);
        if (!isFinite(n)) return "<strong>" + label + "</strong>N/A";
        const pct = total > 0 ? Math.round((n / total) * 100) + "%" : "N/A";
        return "<strong>" + label + "</strong>" + n.toLocaleString() + " <span class=\"cm-pct\">(" + pct + ")</span>";
      };

      const cells = [
        { className: "cm-corner", html: "" },
        { className: "cm-head", html: "Observed +<br>disturbance" },
        { className: "cm-head", html: "Observed -<br>no disturbance" },
        { className: "cm-side", html: "Predicted +" },
        { className: "cm-cell", html: cellHtml("TP", "TP") },
        { className: "cm-cell", html: cellHtml("FP", "FP") },
        { className: "cm-side", html: "Predicted -" },
        { className: "cm-cell", html: cellHtml("FN", "FN") },
        { className: "cm-cell", html: cellHtml("TN", "TN") }
      ];

      cells.forEach(function(item) {
        const cell = document.createElement("div");
        cell.className = item.className;
        cell.innerHTML = item.html;
        DOM.forecastConfusionMatrix.appendChild(cell);
      });
    }

    function updateOverlayCardStates() {
      if (typeof syncDriverPillar === "function") syncDriverPillar();
      if (!DOM.overlayCardStrip) {
        return;
      }

      DOM.overlayCardStrip.querySelectorAll("[data-overlay]").forEach(function(button) {
        const isActive = button.dataset.overlay === APP.focusedOverlay;
        button.classList.toggle("active", isActive);
      });

      const layersBtn = document.getElementById("mobLayersBtn");
      if (layersBtn) layersBtn.textContent = APP.visibleOverlays.length ? "Layers \u00B7 " + APP.visibleOverlays.length : "Layers";

      DOM.overlayCardStrip.querySelectorAll(".bb-accessories").forEach(function(group) {
        const toggle = DOM.overlayCardStrip.querySelector('[data-toggle-target="' + group.id + '"]');
        if (!toggle) return;
        toggle.classList.toggle("has-on", Array.from(group.querySelectorAll("[data-toggle-overlay]")).some(function(b) {
          return isOverlayVisible(b.dataset.toggleOverlay);
        }));
      });

      DOM.overlayCardStrip.querySelectorAll("[data-toggle-overlay]").forEach(function(button) {
        const overlayId = button.dataset.toggleOverlay;
        const isOn = isOverlayVisible(overlayId);
        button.classList.toggle("is-on", isOn);
        button.textContent = isOn ? "Hide" : "Display";
      });

      Object.keys(APP.layerOpacity).forEach(function(overlayId) {
        const slider = DOM.overlayCardStrip.querySelector("[data-opacity-overlay=\"" + overlayId + "\"]");
        const value = DOM.overlayCardStrip.querySelector("[data-opacity-value=\"" + overlayId + "\"]");
        if (slider) slider.value = String(getLayerOpacity(overlayId));
        if (value) value.textContent = getLayerOpacity(overlayId).toFixed(2);
      });
    }

    function applyOpacityToLayers() {
      // Restyle every currently-loaded hex layer (not just the top-visible
      // one) since more than one can be visible/opacity-adjusted at once.
      HEX_OVERLAY_IDS.forEach(function(overlayId) {
        const hexLayer = APP.layerInstances[overlayId];
        if (!hexLayer) return;
        hexLayer.setStyle(baseHexStyle(overlayId));
        if (APP.selectedHexLayer && hexLayer.hasLayer && hexLayer.hasLayer(APP.selectedHexLayer)) {
          highlightHex(APP.selectedHexLayer, overlayId);
        }
      });

      RASTER_OVERLAYS.forEach(function(overlayId) {
        const layer = APP.layerInstances[overlayId];
        if (layer && layer.setOpacity) {
          layer.setOpacity(getLayerOpacity(overlayId));
        }
      });
    }

    function refreshForecastOverlayIfVisible() {
      syncOverlayControls(isOverlayVisible("bbtforecast") ? "bbtforecast" : false);
    }

    function setLayerVisibilityFromRibbon(overlayId, shouldShow) {
      if (!OVERLAY_CATALOG[overlayId] || !OVERLAY_CATALOG[overlayId].available) {
        return;
      }

      // The risk threshold is per displayed surface; clear it (and un-hide any
      // currently visible risk layer) before the displayed set changes so the
      // previous feature's hiding never lingers onto the next one.
      resetRiskThreshold();

      if (shouldShow) {
        showOverlay(overlayId);
        APP.focusedOverlay = overlayId;
        syncOverlayControls(overlayId);
      } else if (isOverlayVisible(overlayId)) {
        hideOverlay(overlayId);
        clearOverlayLayer(overlayId);
        updateOverlayCardStates();
        setLegendLabels();
        syncCitationVisibility();
        setStatus(APP.overlayVisible ? "Overlay hidden" : "All overlays hidden");
      }
    }

    // A period-strip chip's accessory layer (droughtRef, windPressure,
    // phenips, bbtl) lives inside a collapsed bb-accessories panel -- jumping
    // to it should also reveal it, not just change map state behind a closed
    // panel.
    const ACCESSORY_PARENT_GROUP = {
      droughtRef: "droughtAccessories",
      wildfire: "fireAccessories",
      bbtpast: "bbAccessories",
      windPressure: "windAccessories",
      bbtl: "bbAccessories",
      phenips: "bbAccessories"
    };

    function expandAccessoryFor(overlayId) {
      const accessoriesId = ACCESSORY_PARENT_GROUP[overlayId];
      if (!accessoriesId) return;
      const accessories = document.getElementById(accessoriesId);
      if (!accessories || !accessories.hidden) return;
      accessories.hidden = false;
      const toggle = document.querySelector('[data-toggle-target="' + accessoriesId + '"]');
      if (toggle) {
        toggle.classList.add("open");
        toggle.setAttribute("aria-expanded", "true");
      }
    }

    // Click-to-jump for period-strip chips (2026-09 follow-up to the
    // read-only coverage strip -- "give me a link to click"). `period` is a
    // semantic slot ("historical"/"near"/"far"), not necessarily a literal
    // control value -- each branch below translates it to whatever that
    // overlay's own controls actually use, since drought's period selector
    // encodes a scenario too and the hex layers use a literal year instead.
    async function jumpToOverlayPeriod(overlayId, period) {
      if (!OVERLAY_CATALOG[overlayId] || !OVERLAY_CATALOG[overlayId].available) return;

      expandAccessoryFor(overlayId);

      if (overlayId === "droughtFuture") {
        const prefix = period === "near" ? "2041_2070" : "2071_2100";
        const opt = Array.from(DOM.droughtFuturePeriodSelector.options).find(function(o) {
          return o.value.indexOf(prefix) === 0;
        });
        if (opt) APP.currentDroughtFuturePeriod = opt.value;
      } else if (overlayId === "fireWeather") {
        if (period === "historical") { APP.currentFireScen = "historical"; APP.currentFireView = "level"; }
        else {
          const sfx = period === "near" ? "2041_2070" : "2071_2098";
          const cur = APP.currentFireScen === "historical" ? "rcp85" : APP.currentFireScen.split("_")[0];
          APP.currentFireScen = cur + "_" + sfx;
        }
      } else if (overlayId === "wind") {
        APP.currentWindPeriod = period;
      } else if (overlayId === "wildfire" || overlayId === "bbtl") {
        APP.currentYear = period === "near" ? 2030 : 2100;
      } else if (overlayId === "bbtforecast") {
        // The final model has all three slots: historic 1981-2010, near 2041-2070, later 2071-2100.
        // (2011-2040 is reachable from the period selector.) The chips always show the final model.
        try {
          const manifest = await ensureForecastManifestLoaded();
          APP.currentForecastModel = P1C4_MODEL;
          const periods = (manifest.models[P1C4_MODEL] && manifest.models[P1C4_MODEL].periods) || [];
          const want = period === "historical" ? periods[0] : period === "near" ? "2041_2070" : periods[periods.length - 1];
          if (want && periods.includes(want)) APP.currentForecastPeriod = want;
        } catch (error) {
          console.error(error);
        }
      }

      setLayerVisibilityFromRibbon(overlayId, true);
      syncControlsToActive();

      const card = DOM.overlayCardStrip && DOM.overlayCardStrip.querySelector('[data-overlay="' + overlayId + '"]');
      if (card && card.scrollIntoView) card.scrollIntoView({ behavior: "smooth", block: "nearest" });
    }

    async function hydrateForecastControls() {
      updateForecastPanelVisibility();
      if (APP.focusedOverlay !== "bbtforecast") {
        return;
      }

      const manifest = await ensureForecastManifestLoaded();
      const modelKeys = Object.keys(manifest.models || {});
      if (!modelKeys.length) {
        throw new Error("Forecast bundle has no models");
      }

      if (!modelKeys.includes(APP.currentForecastModel)) {
        APP.currentForecastModel = modelKeys[0];
      }

      const periods = (manifest.models[APP.currentForecastModel] && manifest.models[APP.currentForecastModel].periods) || manifest.periods || [];
      if (!periods.includes(APP.currentForecastPeriod)) {
        APP.currentForecastPeriod = periods[0];
      }

      const entry = getForecastEntry();
      if (!entry) {
        throw new Error("Forecast entry missing for " + getForecastEntryKey());
      }

      const metrics = await ensureForecastMetricsLoaded(APP.currentForecastModel);
      const thresholds = getMatchingForecastThresholds(manifest, metrics, entry);

      if (!thresholds.length) {
        throw new Error("No matching forecast thresholds found for " + APP.currentForecastModel);
      }

      if (!thresholds.includes(APP.currentForecastThreshold)) {
        APP.currentForecastThreshold = thresholds.includes("0.75") ? "0.75" : thresholds[0];
      }

      // Most models (all but "ensemble") have no delta_from_near_tif at all --
      // selecting "delta" view for one used to fail silently (blank map, a
      // console error, and a status-bar line the user would have to notice).
      // Disable the option up front instead, and fall back off it if it was
      // already selected (e.g. switching FROM ensemble TO a model without one).
      const modelSupportsDelta = periods.some(function(p) {
        const e = manifest.entries[APP.currentForecastModel + "::" + p];
        return e && e.delta_from_near_tif;
      });
      if (DOM.forecastViewDeltaOption) {
        DOM.forecastViewDeltaOption.disabled = !modelSupportsDelta;
      }
      if (!modelSupportsDelta && APP.currentForecastView === "delta") {
        APP.currentForecastView = "probability";
      }
      // Extra views of the final model: share of years >= 0.5 and the extrapolation share (only when delivered).
      const hasField = function(field) {
        return periods.some(function(p) { const e = manifest.entries[APP.currentForecastModel + "::" + p]; return e && e[field]; });
      };
      const finalModel = isP1c4();
      if (DOM.forecastViewShareOption) DOM.forecastViewShareOption.disabled = !hasField("share05_tif");
      if (DOM.forecastViewExtrapOption) DOM.forecastViewExtrapOption.disabled = !hasField("extrapolation_tif");
      if (DOM.forecastViewDeltaOption) DOM.forecastViewDeltaOption.textContent = finalModel ? "Change vs 1981–2010 (index points)" : "Δ % change vs. current (2011–2040)";
      const thrOpt = DOM.forecastViewSelector.querySelector('option[value="threshold"]');
      if (thrOpt) thrOpt.textContent = finalModel ? "Flagged areas (hazard 0.5 or more)" : "High-hazard areas (above cut-off)";
      const probOpt = DOM.forecastViewSelector.querySelector('option[value="probability"]');
      if (probOpt) probOpt.textContent = finalModel ? "Hazard (0–0.8 scale)" : "Hazard (0–1)";
      [["share05", DOM.forecastViewShareOption], ["extrapolation", DOM.forecastViewExtrapOption]].forEach(function(pair) {
        if (APP.currentForecastView === pair[0] && pair[1] && pair[1].disabled) APP.currentForecastView = "probability";
      });

      // The per-cutoff note (below, for the currently selected threshold) can
      // be misleadingly generic when NA -- this model-level note is the real,
      // model-specific documentation already written into metrics_confusion.json
      // (e.g. explaining that a model IS validated, just not at every cutoff).
      if (DOM.forecastModelNote) {
        DOM.forecastModelNote.textContent = metrics.notes || "";
      }

      const manifestModels = manifest.models || {};
      populateSelectOptions(DOM.forecastModelSelector, modelKeys, APP.currentForecastModel, function(value) {
        const base = manifestModels[value] && manifestModels[value].label ? manifestModels[value].label : value.toUpperCase();
        return value === P1C4_MODEL ? base : base + " (earlier model version)";
      });
      populateSelectOptions(DOM.forecastPeriodSelector, periods, APP.currentForecastPeriod, function(value) {
        const labels = { "1981_2010": "Historic (1981-2010)", "2011_2040": "Recent (2011-2040)", "2041_2070": "Near future (2041-2070)", "2071_2100": "Later future, end of century (2071-2100)", "1983_2010": "Reference (1983-2010)", "reference_1984_2010": "Reference (1984-2010)", "future_2071_2100": "End century (2071-2100)" };
        return labels[value] || value.replace("_", "-");
      });
      populateSelectOptions(DOM.forecastThresholdSelector, thresholds, APP.currentForecastThreshold, function(value) {
        return formatForecastThreshold(value);
      });
      DOM.forecastViewSelector.value = APP.currentForecastView;
      DOM.forecastThresholdSelector.disabled = APP.currentForecastView !== "threshold";

      const cutoff = Number(APP.currentForecastThreshold);
      const row = (metrics.aggregate_by_cutoff || []).find(function(item) {
        return Number(item.cutoff) === cutoff;
      });

      if (metrics.grid_items) {
        renderForecastMetricsGrid(metrics.grid_items);
        renderForecastConfusionMatrix(null);
        DOM.forecastMetricsNote.textContent = metrics.validation_note || "";
      } else {
        renderForecastPerformance(row, cutoff, metrics.auc);
      }
    }

    function renderLegendSwatches(items) {
      DOM.legendSwatches.innerHTML = "";
      if (!items || !items.length) {
        return;
      }

      items.forEach(function(item) {
        const row = document.createElement("div");
        row.className = "legend-swatch";

        const color = document.createElement("span");
        color.className = "legend-swatch-color";
        color.style.background = item.color;

        const label = document.createElement("span");
        label.textContent = item.label;

        row.appendChild(color);
        row.appendChild(label);
        DOM.legendSwatches.appendChild(row);
      });
    }


    // ====================================================================
    // LEGEND
    // ====================================================================
    function setLegendLabels() {
      setLegendLabelsCore();
      updatePeriodUi();
    }

    // ====================================================================
    // PERIOD STRIP: active chip, "what is shown" badge, play-through
    // ====================================================================
    function selectedText(dom) {
      return dom && dom.selectedOptions && dom.selectedOptions[0] ? dom.selectedOptions[0].textContent.trim() : "";
    }

    function periodTextFor(id) {
      if (id === "droughtFuture") return selectedText(DOM.droughtFuturePeriodSelector);
      if (id === "droughtRef") return "1981–2010";
      if (id === "fireWeather") {
        const s = APP.currentFireScen || "";
        if (s === "historical") return "1981–2005";
        const m = s.match(/^(rcp\d+)_(\d{4})_(\d{4})$/);
        return m ? m[1].toUpperCase().replace("RCP", "RCP ") + " · " + m[2] + "–" + m[3] : s;
      }
      if (id === "wind") return selectedText(DOM.windPeriodSelector);
      if (id === "bbtforecast") return String(APP.currentForecastPeriod || "").replace(/_/g, "–");
      if (id === "bbtpast") return "1985–2024";
      if (id === "phenips") return (phenipsIsObs() ? "observed " : String(APP.currentPhenipsScen).toUpperCase().replace("RCP", "RCP ") + " · ") + selectedText(DOM.phenisPeriodSelector).replace("Single year ", "");
      if (id === "bbtl" || id === "wildfire") return String(APP.currentRcp || "").toUpperCase().replace("RCP", "RCP ") + " · " + APP.currentYear;
      return "";
    }

    function chipIsActive(chip) {
      const id = chip.dataset.jumpOverlay, p = chip.dataset.jumpPeriod;
      if (!isOverlayVisible(id)) return false;
      if (id === "droughtFuture") return String(APP.currentDroughtFuturePeriod).indexOf(p === "near" ? "2041_2070" : "2071_2100") === 0;
      if (id === "fireWeather") {
        const s = APP.currentFireScen || "";
        return p === "historical" ? s === "historical" : s.endsWith(p === "near" ? "2041_2070" : "2071_2098");
      }
      if (id === "wind") return APP.currentWindPeriod === p;
      if (id === "bbtl") return APP.currentYear === 2030;
      if (id === "wildfire") return APP.currentYear === (p === "near" ? 2030 : 2100);
      if (id === "bbtforecast") {
        const m = APP.forecastManifest, e = m && m.models && m.models[P1C4_MODEL];
        const per = (e && e.periods) || [];
        if (!per.length || !isP1c4()) return false;
        return APP.currentForecastPeriod === (p === "historical" ? per[0] : p === "near" ? "2041_2070" : per[per.length - 1]);
      }
      return true; // single-slot chips (Historic of droughtRef / bbtpast): shown = active
    }

    function updatePeriodUi() {
      if (typeof syncPillarStatuses === "function" && DOM.dssPillarDetail) syncPillarStatuses();
      document.querySelectorAll(".period-chip[data-jump-overlay]").forEach(function(chip) {
        const on = chipIsActive(chip);
        chip.classList.toggle("is-active", on);
        chip.setAttribute("aria-pressed", on ? "true" : "false");
      });
      const badge = document.getElementById("periodBadge");
      if (!badge) return;
      const top = APP.overlayVisible ? getTopVisibleOverlayId() : null;
      const ov = top && OVERLAY_CATALOG[top];
      if (!ov) { badge.hidden = true; return; }
      const per = periodTextFor(top);
      badge.textContent = ov.label + (per ? " · " + per : "");
      badge.hidden = false;
    }

    function setLegendLabelsCore() {
      if (!APP.overlayVisible) {
        syncRiskFilterBox(null);
        renderLegendSwatches([]);
        DOM.legendBar.style.display = "none";
        DOM.legendMin.textContent = "No layer on yet";
        DOM.legendMax.textContent = "";
        if (DOM.legendForLayer) DOM.legendForLayer.textContent = "";
        return;
      }

      const topOverlayId = getTopVisibleOverlayId();
      const overlay = getVisibleOverlay();
      renderLegendSwatches([]);
      if (DOM.legendFlagMarker) DOM.legendFlagMarker.style.display = "none";
      if (DOM.legendFlagNote) DOM.legendFlagNote.hidden = true;
      DOM.legendBar.style.display = "block";
      DOM.legendBar.style.background = BARK_BEETLE_GRADIENT;
      DOM.legendMin.parentElement.style.display = "";
      syncRiskFilterBox(topOverlayId);
      updateRiskSliderRange(topOverlayId);
      updateLegendThreshMarker(topOverlayId);

      // The legend describes the top VISIBLE layer, which is deliberately
      // independent of whichever card's controls are focused above (so you
      // can browse a layer's description while a different one stays
      // displayed) -- when the two differ, say so, or the numbers below
      // silently belong to something other than what "Focused layer" says.
      if (DOM.legendForLayer) {
        const topLabel = topOverlayId && OVERLAY_CATALOG[topOverlayId] && OVERLAY_CATALOG[topOverlayId].label;
        DOM.legendForLayer.textContent = (topOverlayId && topOverlayId !== APP.focusedOverlay && topLabel)
          ? "— showing: " + topLabel
          : "";
      }
      if (topOverlayId === "phenips" && phenipsIsObs()) {
        const lay = APP.phenipsV5 && APP.phenipsV5.layers[getPhenipsV5Name(APP)];
        renderLegendSwatches([]);
        if (!lay) { DOM.legendBar.style.display = "none"; DOM.legendMin.textContent = "Loading..."; DOM.legendMax.textContent = ""; return; }
        const fmt = function(x) { return (x > 0 && lay.min < 0 ? "+" : "") + (Math.abs(x) < 10 ? String(x).replace("-", "−") : Math.round(x)); };
        const unit = APP.currentPhenipsVar === "gdd" ? " DD" : (APP.currentPhenipsVar === "gen" ? " gen" : "");
        DOM.legendBar.style.display = "block";
        DOM.legendBar.style.background = buildGradientCss(lay.colours);
        DOM.legendMin.textContent = APP.currentPhenipsVar === "share" ? "0 % of years" : fmt(lay.min) + unit;
        DOM.legendMax.textContent = APP.currentPhenipsVar === "share" ? "100 % of years" : fmt(lay.max) + unit;
        return;
      }
      if (topOverlayId === "phenips") {
        const v = APP.currentPhenipsVar;
        const view = APP.currentPhenipsView;
        const p = APP.currentPhenisPeriod;
        DOM.legendBar.style.display = "none";
        if (v === "gdd" && view === "abs") {
          DOM.legendBar.style.display = "block";
          DOM.legendBar.style.background = BARK_BEETLE_GRADIENT;
          DOM.legendMin.textContent = "0 DD";
          DOM.legendMax.textContent = "4000 DD";
        } else if (v === "gdd" && view !== "abs") {
          DOM.legendBar.style.display = "block";
          DOM.legendBar.style.background = BARK_BEETLE_GRADIENT;
          DOM.legendMin.textContent = "0 DD";
          DOM.legendMax.textContent = "+1100 DD";
        } else {
          // generations — categorical swatches, sampled off the same ramp
          renderLegendSwatches([
            { color: "#abdda4", label: "1 generation" },
            { color: "#fec776", label: "2 generations" },
            { color: "#9e0142", label: "3 generations" }
          ]);
          if (view === "contrast") {
            DOM.legendMin.textContent = "Now (hist) →";
            DOM.legendMax.textContent = "near-future gain";
          } else if (view === "diff" && p !== "hist") {
            renderLegendSwatches([
              { color: "#fef0a5", label: "+1 generation" },
              { color: "#9e0142", label: "+2 generations" }
            ]);
            DOM.legendMin.textContent = "no change";
            DOM.legendMax.textContent = "max gain";
          } else {
            DOM.legendMin.textContent = "0 gen";
            DOM.legendMax.textContent = "3 gen";
          }
        }
        return;
      }

      if (topOverlayId === "bbtpast" && overlay.meta && overlay.meta.legend && overlay.meta.legend.gradient) {
        renderLegendSwatches([]);
        DOM.legendBar.style.display = "block";
        DOM.legendBar.style.background = buildGradientCss(overlay.meta.legend.gradient);
        DOM.legendMin.textContent = "1 year";
        DOM.legendMax.textContent = "10+ years";
        return;
      }

      if (topOverlayId === "bbtpast" && overlay.meta && overlay.meta.legend) {
        DOM.legendBar.style.display = "none";
        DOM.legendMin.textContent = overlay.meta.legend.labels[0] || "min";
        DOM.legendMax.textContent = overlay.meta.legend.labels[overlay.meta.legend.labels.length - 1] || "max";
        renderLegendSwatches(overlay.meta.legend.labels.map(function(label, index) {
          return {
            color: overlay.meta.legend.colors[index],
            label: label
          };
        }));
        return;
      }

      if (topOverlayId === "bbtpast") {
        DOM.legendBar.style.display = "none";
        DOM.legendMin.textContent = "Loading...";
        DOM.legendMax.textContent = "";
        return;
      }

      if (topOverlayId === "droughtFuture" || topOverlayId === "droughtRef") {
        const key = topOverlayId === "droughtFuture" ? APP.currentDroughtFuturePeriod : APP.currentDroughtRefMetric;
        const image = overlay.meta && overlay.meta.images && overlay.meta.images[key];
        DOM.legendBar.style.display = "none";
        // The min/max summary row is redundant here (and the full label
        // text is too long for it) -- the swatch list below already shows
        // every class with its own label.
        DOM.legendMin.parentElement.style.display = "none";
        if (!image) {
          renderLegendSwatches([]);
          return;
        }
        renderLegendSwatches(image.legend.labels.map(function(label, index) {
          return { color: image.legend.colors[index], label: label };
        }));
        return;
      }

      if (topOverlayId === "fireWeather") {
        const key = fireImageKey(APP);
        const legends = overlay.meta && overlay.meta.legends;
        DOM.legendBar.style.display = "none";
        DOM.legendMin.parentElement.style.display = "none";
        if (!legends) { renderLegendSwatches([]); return; }
        const lg = key.indexOf("change::") === 0 ? legends.change : legends.level;
        renderLegendSwatches(lg.labels.map(function(label, index) {
          return { color: lg.colors[index], label: label };
        }));
        return;
      }

      if (topOverlayId === "wind") {
        DOM.legendBar.style.display = "block";
        DOM.legendBar.style.background = WIND_GRADIENT;
        DOM.legendMin.textContent = "0 — low";
        DOM.legendMax.textContent = "≥ 0.5 — high";
        return;
      }

      if (topOverlayId === "bbtforecast" && isP1c4()) {
        const v = APP.currentForecastView;
        DOM.legendBar.style.display = "block";
        if (v === "delta") {
          DOM.legendBar.style.background = DELTA_GRADIENT;
          DOM.legendMin.textContent = "−0.3 (fall)"; DOM.legendMax.textContent = "+0.3 (rise)";
        } else if (v === "share05") {
          DOM.legendBar.style.background = P1C4_GRADIENT;
          DOM.legendMin.textContent = "0 — no year"; DOM.legendMax.textContent = "1 — every year";
        } else if (v === "extrapolation") {
          DOM.legendBar.style.background = EXTRAP_GRADIENT;
          DOM.legendMin.textContent = "0 — inside training range"; DOM.legendMax.textContent = "1 — outside";
        } else {
          DOM.legendBar.style.background = P1C4_GRADIENT;
          DOM.legendMin.textContent = "0"; DOM.legendMax.textContent = P1C4_CAP.toFixed(1);
          if (DOM.legendFlagMarker) { DOM.legendFlagMarker.style.left = (P1C4_FLAG / P1C4_CAP * 100) + "%"; DOM.legendFlagMarker.style.display = "block"; }
          if (DOM.legendFlagNote) DOM.legendFlagNote.hidden = false;
        }
        return;
      }

      if (topOverlayId === "bbtforecast" && APP.currentForecastView === "delta") {
        DOM.legendBar.style.display = "block";
        DOM.legendBar.style.background = DELTA_GRADIENT;
        DOM.legendMin.textContent = "-" + DELTA_CAP + "%";
        DOM.legendMax.textContent = "+" + DELTA_CAP + "%";
        return;
      }

      if (topOverlayId === "bbtforecast") {
        DOM.legendBar.style.display = "block";
        // Both views now share the same ramp (continuous 0-1 vs. the
        // classified view's 8 discrete bins over 0.55-1.0+) so switching
        // views no longer changes the color language, only the granularity.
        DOM.legendBar.style.background = BARK_BEETLE_GRADIENT;
        if (APP.currentForecastView === "threshold") {
          DOM.legendMin.textContent = formatForecastThreshold(APP.currentForecastThreshold);
          DOM.legendMax.textContent = "1.00";
        } else if (APP.currentForecastRasterRange) {
          DOM.legendMin.textContent = Number(APP.currentForecastRasterRange.min).toFixed(2);
          DOM.legendMax.textContent = Number(APP.currentForecastRasterRange.max).toFixed(2);
        } else {
          DOM.legendMin.textContent = "0.00";
          DOM.legendMax.textContent = "1.00";
        }
        return;
      }

      if (!overlay.available || !overlay.meta) {
        DOM.legendMin.textContent = "NA";
        DOM.legendMax.textContent = "NA";
        return;
      }

      DOM.legendMin.textContent = Number(overlay.meta.valueRange.zmin).toFixed(3);
      DOM.legendMax.textContent = Number(overlay.meta.valueRange.zmax).toFixed(3);
    }

    function syncCitationVisibility() {
      const topId = getTopVisibleOverlayId();
      const showCitation = HEX_OVERLAY_IDS.includes(topId) && APP.overlayVisible;
      DOM.citationBox.style.display = showCitation ? "block" : "none";
      if (showCitation && DOM.citationForLayer) {
        DOM.citationForLayer.textContent = "— " + OVERLAY_CATALOG[topId].label;
      }
    }


    // ====================================================================
    // HEX / COUNTRY STYLING, CLICK MARKER & READOUT
    // ====================================================================
    function getColor(overlayId, value) {
      if (!isFinite(value)) {
        return "#cccccc";
      }

      const overlay = OVERLAY_CATALOG[overlayId];
      if (!overlay || !overlay.available || !overlay.meta) {
        return "#cccccc";
      }

      const zmin = Number(overlay.meta.valueRange.zmin);
      const zmax = Number(overlay.meta.valueRange.zmax);

      if (zmax <= zmin) {
        return BARK_BEETLE_RAMP[0];
      }

      const t = (value - zmin) / (zmax - zmin);
      return rampColor(BARK_BEETLE_RAMP, t);
    }

    // Returns a Leaflet `style` callback for the given hex overlay -- pass
    // the result to L.geoJSON's `style` option or to `.setStyle(...)`.
    function baseHexStyle(overlayId) {
      return function(feature) {
        const value = Number(feature.properties.mean);
        let fillOpacity = getLayerOpacity(overlayId);
        // Dissolve the low-risk end: scale opacity down across the bottom of the
        // scale so the darkest-blue hexes fade into the basemap.
        const meta = OVERLAY_CATALOG[overlayId].meta;
        if (meta && meta.valueRange && isFinite(value)) {
          const zmax = Number(meta.valueRange.zmax);
          const t = zmax > 0 ? value / zmax : 1;
          const fade = RISK_FADE_FLOOR + (1 - RISK_FADE_FLOOR) * Math.max(0, Math.min(1, t / RISK_FADE_KNEE));
          fillOpacity *= fade;
        }
        // Hide hexes whose value is below the slider's absolute risk threshold,
        // but only while this overlay is the top (displayed) layer the slider describes.
        if (APP.riskThreshold > 0 && getTopVisibleOverlayId() === overlayId &&
            isFinite(value) && value < APP.riskThreshold) {
          fillOpacity = 0;
        }
        // Keep border tied to fill opacity so slider behavior feels predictable.
        const borderOpacity = fillOpacity <= 0 ? 0 : Math.min(1, fillOpacity * 0.6);

        return {
          color: "#222222",
          weight: 0.9,
          opacity: borderOpacity,
          fillColor: getColor(overlayId, value),
          fillOpacity: fillOpacity
        };
      };
    }

    function highlightHex(layer, overlayId) {
      const baseOpacity = getLayerOpacity(overlayId);
      const selectedOpacity = Math.min(baseOpacity + 0.15, 0.5);

      layer.setStyle({
        color: "#111111",
        weight: 2.0,
        opacity: 1,
        fillColor: "#ffffff",
        fillOpacity: selectedOpacity
      });

      if (layer.bringToFront) {
        layer.bringToFront();
      }
    }

    // Not told which hex overlay owns the current selection (there can be
    // several visible at once, see HEX_OVERLAY_IDS) -- so ask each layer
    // instance directly rather than guessing from the top-visible id.
    function resetSelectedHex() {
      if (!APP.selectedHexLayer) return;
      HEX_OVERLAY_IDS.forEach(function(id) {
        const layer = APP.layerInstances[id];
        if (layer && layer.hasLayer && layer.hasLayer(APP.selectedHexLayer)) {
          layer.resetStyle(APP.selectedHexLayer);
        }
      });
      APP.selectedHexLayer = null;
    }

    function resetSelectedCountry() {
      if (APP.selectedCountryLayer && APP.countryLayer) {
        APP.countryLayer.resetStyle(APP.selectedCountryLayer);
        APP.selectedCountryLayer = null;
      }
    }

    function updateMarker(latlng) {
      if (APP.clickMarker) {
        APP.clickMarker.setLatLng(latlng);
      } else {
        APP.clickMarker = L.circleMarker(latlng, {
          radius: 5,
          color: "#666666",
          weight: 1,
          fillColor: "#ffffff",
          fillOpacity: 0.7
        }).addTo(APP.map);
      }
    }

    // In split screen, drop the same marker on the OTHER map too, so both
    // sides point at the same spot for comparison. `side` is whichever map
    // was actually clicked (already the active side by the time this runs).
    function mirrorClickMarker(latlng, side) {
      if (!APP.split) return;
      const otherSide = side === "left" ? "right" : "left";
      const panel = APP.panels[otherSide];
      if (!panel || !panel.map) return;
      if (panel.clickMarker) {
        panel.clickMarker.setLatLng(latlng);
      } else {
        panel.clickMarker = L.circleMarker(latlng, {
          radius: 5,
          color: "#666666",
          weight: 1,
          fillColor: "#ffffff",
          fillOpacity: 0.7
        }).addTo(panel.map);
      }
    }

    function getCurrentFile(overlayId) {
      const overlay = OVERLAY_CATALOG[overlayId];
      if (!overlay || !overlay.available || !overlay.meta) {
        return null;
      }
      return overlay.meta.files[APP.currentRcp + "_" + APP.currentYear];
    }

    function getCountryName(feature) {
      if (!feature || !feature.properties) {
        return "Unknown";
      }

      return feature.properties.CNTR_NAME
        || feature.properties.NAME_ENGL
        || feature.properties.name
        || feature.properties.NAME_EN
        || feature.properties.ADMIN
        || feature.properties.SOVEREIGNT
        || "Unknown";
    }

    // --- Small geometry helpers (replace the 600 KB Turf.js; only these two
    // questions were ever asked of it) ----------------------------------------
    // GeoJSON Polygon / MultiPolygon -> list of polygons, each a list of rings.
    function geometryPolygons(geometry) {
      if (!geometry) return [];
      if (geometry.type === "Polygon") return [geometry.coordinates];
      if (geometry.type === "MultiPolygon") return geometry.coordinates;
      return [];
    }

    // Even-odd ray casting, so holes (rings after the first) count as outside.
    function pointInRing(x, y, ring) {
      let inside = false;
      for (let i = 0, j = ring.length - 1; i < ring.length; j = i++) {
        const xi = ring[i][0], yi = ring[i][1], xj = ring[j][0], yj = ring[j][1];
        if ((yi > y) !== (yj > y) && x < (xj - xi) * (y - yi) / (yj - yi) + xi) inside = !inside;
      }
      return inside;
    }

    function pointInFeature(lng, lat, feature) {
      return geometryPolygons(feature && feature.geometry).some(function(rings) {
        let n = 0;
        for (let k = 0; k < rings.length; k++) if (pointInRing(lng, lat, rings[k])) n++;
        return n % 2 === 1;
      });
    }

    // True if any polygon part of the feature has a bounding box that overlaps
    // the window [xmin, ymin, xmax, ymax] (a cheap test; good enough to drop
    // countries far outside the map extent).
    function featureTouchesWindow(feature, xmin, ymin, xmax, ymax) {
      return geometryPolygons(feature && feature.geometry).some(function(rings) {
        let x0 = Infinity, y0 = Infinity, x1 = -Infinity, y1 = -Infinity;
        rings[0].forEach(function(c) {
          if (c[0] < x0) x0 = c[0]; if (c[0] > x1) x1 = c[0];
          if (c[1] < y0) y0 = c[1]; if (c[1] > y1) y1 = c[1];
        });
        return x0 <= xmax && x1 >= xmin && y0 <= ymax && y1 >= ymin;
      });
    }

    function getMapExtentCountryGeojson(geojson) {
      if (!geojson || !geojson.features) {
        return geojson;
      }

      const b = MAP_META.bbox;

      return {
        type: "FeatureCollection",
        features: geojson.features.filter(function(feature) {
          return featureTouchesWindow(feature, b.xmin, b.ymin, b.xmax, b.ymax);
        })
      };
    }

    function fitDataBounds() {
      const topVisibleOverlayId = getTopVisibleOverlayId();
      const overlay = getVisibleOverlay();
      const meta = (overlay.meta && overlay.meta.bbox) ? overlay.meta : MAP_META;
      const hexLayer = getHexLayer(topVisibleOverlayId);
      const topLayer = APP.layerInstances[topVisibleOverlayId];

      if (HEX_OVERLAY_IDS.includes(topVisibleOverlayId) && hexLayer) {
        APP.map.fitBounds(hexLayer.getBounds());
      } else if (topLayer && topLayer.getBounds) {
        APP.map.fitBounds(topLayer.getBounds());
      } else {
        APP.map.fitBounds([
          [meta.bbox.ymin, meta.bbox.xmin],
          [meta.bbox.ymax, meta.bbox.xmax]
        ]);
      }
    }

    function resetLayerCache(overlayId) {
      if (HEX_OVERLAY_IDS.includes(overlayId)) {
        APP.hexGeojson[overlayId] = null;
        resetSelectedHex();
      } else if (overlayId === "bbtforecast") {
        APP.currentForecastRasterRange = null;
      }
    }

    function clearOverlayLayer(overlayId) {
      const targets = overlayId ? [overlayId] : getOverlayIds();

      targets.forEach(function(id) {
        APP.overlayLoadTokens[id] += 1;
        clearLayerInstance(id);
        resetLayerCache(id);
      });

      HEX_OVERLAY_IDS.forEach(function(id) {
        if (!isOverlayVisible(id)) APP.hexGeojson[id] = null;
      });
      resetSelectedCountry();
    }


    // ====================================================================
    // LAYER LOADERS & CONTROL-PANEL SYNC
    // ====================================================================
    async function updateOverlayLayer(overlayId) {
      if (!isOverlayVisible(overlayId)) {
        return;
      }

      const updateLayer = OVERLAY_LOADERS[overlayId];
      if (updateLayer) {
        await updateLayer();
      }

      updateOverlayPaneOrder();
    }

    async function updateClassifiedImageOverlay() {
      if (!isOverlayVisible("bbtpast")) return;

      const requestToken = beginOverlayLoad("bbtpast");
      const overlay = OVERLAY_CATALOG.bbtpast;
      setStatus("Loading historical bark beetle overlay...");

      try {
        await ensurePastOverlayMetaLoaded();

        if (!isCurrentOverlayLoad("bbtpast", requestToken)) {
          return;
        }

        addImageOverlayLayer("bbtpast", overlay.meta.image, overlay.meta.bounds);

        setLegendLabels();
        setStatus("");
      } catch (error) {
        console.error(error);
        setStatus("Historical bark beetle overlay could not be loaded");
      }
    }

    async function updateForecastOverlay() {
      if (!isOverlayVisible("bbtforecast")) return;

      const requestToken = beginOverlayLoad("bbtforecast");

      APP.currentForecastRasterRange = null;
      setStatus("Loading bark beetle forecast...");

      try {
        await ensureForecastManifestLoaded(); // a deep link can draw the layer before any panel has asked for the manifest
        await hydrateForecastControls();
        const entry = getForecastEntry();
        if (!entry) {
          throw new Error("Forecast entry missing for " + getForecastEntryKey());
        }

        const overlayPaths = getForecastOverlayPaths(entry);
        if (!overlayPaths) {
          setStatus(APP.currentForecastView === "delta" ? "No change map for the reference period: it is the baseline." : "This view is not available for the selected period.");
          clearOverlayLayer("bbtforecast");
          return;
        }

        const response = await fetch(overlayPaths.metaPath, { cache: "no-store" });
        if (!response.ok) {
          throw new Error("HTTP " + response.status + " for " + overlayPaths.metaPath);
        }

        const overlayMeta = await response.json();
        if (!isCurrentOverlayLoad("bbtforecast", requestToken)) {
          return;
        }

        APP.currentForecastRasterRange = {
          min: overlayMeta.value_range && isFinite(Number(overlayMeta.value_range.min))
            ? Number(overlayMeta.value_range.min)
            : (APP.currentForecastView === "threshold" ? Number(APP.currentForecastThreshold) : 0),
          max: overlayMeta.value_range && isFinite(Number(overlayMeta.value_range.max))
            ? Number(overlayMeta.value_range.max)
            : 1
        };
        addRiskImageOverlay("bbtforecast", overlayPaths.imagePath, overlayMeta.bounds, requestToken);

        if (!isCurrentOverlayLoad("bbtforecast", requestToken)) {
          return;
        }

        setLegendLabels();
        setStatus("");
      } catch (error) {
        console.error(error);
        setStatus("Forecast overlay could not be loaded for " + APP.currentForecastModel.toUpperCase() + " " + APP.currentForecastPeriod.replace("_", "-"));
      }
    }

    const OVERLAY_LOADERS = {
      bbtl: function() { return updateHexLayer("bbtl"); },
      wildfire: function() { return updateHexLayer("wildfire"); },
      windPressure: function() { return updateHexLayer("windPressure"); },
      phenips: updatePhenipsLayer,
      bbtpast: updateClassifiedImageOverlay,
      bbtforecast: updateForecastOverlay,
      wind: updateWindLayer,
      droughtFuture: updateDroughtFutureLayer,
      droughtRef: updateDroughtRefLayer,
      fireWeather: updateFireWeatherLayer
    };

    // Show only the first sentence; the rest sits behind a "More" link.
    function setShortDescription(el, text) {
      el.textContent = "";
      const m = /^(.+?[.!?])\s+(.+)$/s.exec(text || "");
      if (!m) { el.textContent = text || ""; return; }
      el.appendChild(document.createTextNode(m[1] + " "));
      const rest = document.createElement("span");
      rest.textContent = m[2]; rest.hidden = true;
      const link = document.createElement("a");
      link.href = "#"; link.textContent = "More";
      link.addEventListener("click", function(e) {
        e.preventDefault(); rest.hidden = !rest.hidden; link.textContent = rest.hidden ? "More" : "Less";
      });
      el.appendChild(link); el.appendChild(document.createTextNode(" ")); el.appendChild(rest);
    }

    function syncOverlayControls(refreshTarget) {
      const focusedOverlay = getFocusedOverlay();
      const showHexControls = HEX_OVERLAY_IDS.includes(APP.focusedOverlay);

      // Selection and display are deliberately separate:
      // card click changes settings, while Display/Hide changes map visibility.
      setShortDescription(DOM.overlayDescription, focusedOverlay.description);
      updateOverlayCardStates();
      updateForecastPanelVisibility();
      DOM.hexControlPanel.hidden = !showHexControls;
      if (showHexControls) DOM.hexControlTitle.textContent = focusedOverlay.label + " Controls";
      DOM.rcpSelector.disabled = !showHexControls;
      DOM.yearSlider.disabled = !showHexControls;
      DOM.yearValue.style.opacity = showHexControls ? "1" : "0.6";
      DOM.toggleOverlayBtn.textContent = isOverlayVisible(APP.focusedOverlay) ? "Hide overlay" : "Show overlay";
      syncCitationVisibility();

      if (!focusedOverlay.available) {
        setStatus(focusedOverlay.label + " overlay is not available yet");
        return;
      }

      if (!showHexControls) {
        DOM.yearValue.textContent = "NA";
        setLegendLabels();
        if (refreshTarget && APP.overlayVisible) {
          updateOverlayLayer(refreshTarget);
        }
        return;
      }

      const years = focusedOverlay.meta.years;
      DOM.yearSlider.max = years.length - 1;
      if (!years.includes(APP.currentYear)) {
        APP.currentYear = years[0];
      }
      DOM.yearSlider.value = String(years.indexOf(APP.currentYear));
      DOM.yearValue.textContent = String(APP.currentYear);
      setLegendLabels();
      if (refreshTarget && APP.overlayVisible) {
        updateOverlayLayer(refreshTarget);
      }
    }

    function findContainingFeature(geojson, latlng) {
      if (!geojson || !geojson.features) return null;

      for (const feature of geojson.features) {
        if (pointInFeature(latlng.lng, latlng.lat, feature)) return feature;
      }

      return null;
    }

    function findLayerForFeature(targetFeature, overlayId) {
      const hexLayer = getHexLayer(overlayId);
      if (!hexLayer || !targetFeature) return null;

      let foundLayer = null;

      hexLayer.eachLayer(function(layer) {
        if (layer.feature === targetFeature) {
          foundLayer = layer;
        }
      });

      return foundLayer;
    }

    // ====================================================================
    // CLICK INSPECTION & POPUP
    // ====================================================================
    // Drives the map-side effects of a click (marker, hex highlight); no
    // longer touches a side-panel readout or a country border highlight --
    // both were removed as UI clutter. setStatus() calls remain as a
    // console-only lifecycle trace (see setStatus).
    function updatePanelFromClick(latlng) {
      updateMarker(latlng);

      const topOverlayId = getTopVisibleOverlayId();
      const hexFeature = findContainingFeature(APP.hexGeojson[topOverlayId], latlng);

      resetSelectedHex();
      resetSelectedCountry();

      if (topOverlayId === "phenips") {
        if (phenipsIsObs()) { setStatus("PHENIPS observed: read the value from the legend colour scale"); return; }
        const pt = samplePhenipsAt(latlng);
        setStatus(pt ? "PHENIPS sampled at nearest grid point" : "No PHENIPS data at clicked point");
        return;
      }

      // Pre-rendered images, not queryable data -- nothing to sample.
      if (topOverlayId === "bbtforecast") return;
      if (topOverlayId === "bbtpast") return;
      if (topOverlayId === "wind") return;
      if (topOverlayId === "droughtFuture" || topOverlayId === "droughtRef" || topOverlayId === "fireWeather") return;

      if (hexFeature) {
        const layer = findLayerForFeature(hexFeature, topOverlayId);
        if (layer) {
          APP.selectedHexLayer = layer;
          highlightHex(layer, topOverlayId);
        }
        setStatus("Hex found");
      } else {
        setStatus("No hex at clicked point");
      }
    }

    // Computes the label + sampled value for a click, against an arbitrary
    // state bundle -- mirrors the per-overlay branches in updatePanelFromClick
    // but reads no DOM and mutates nothing, so it's safe to call for the
    // inactive split-screen side too (with its own parked overlay/period).
    function computeClickValueInfo(state, latlng, overlayIdOverride) {
      const overlayId = overlayIdOverride || getTopVisibleOverlayIdFor(state);
      const overlay = OVERLAY_CATALOG[overlayId] || {};
      const label = overlay.label || "Layer";

      if (overlayId === "phenips" && state.currentPhenipsScen === "obs") return { label: label, valueHtml: null };
      if (overlayId === "phenips") {
        const pt = samplePhenipsAt(latlng);
        if (!pt) return { label: label, valueHtml: "—" };
        const v    = state.currentPhenipsVar;
        const p    = state.currentPhenisPeriod;
        const unit = v === "gdd" ? " DD" : " gen";
        const prec = v === "gdd" ? 0 : 2;
        const hist = v === "gdd" ? pt.gdd_hist : pt.gen_hist;
        const scens = ["rcp85", "rcp26", "rcp45"];
        const labels = { rcp85: "RCP8.5", rcp26: "RCP2.6", rcp45: "RCP4.5" };
        let lines = [];
        if (p === "hist") {
          lines.push("Hist: " + (hist != null ? hist.toFixed(prec) + unit : "—"));
        } else {
          scens.forEach(function(s) {
            const key = (v === "gdd" ? "gdd_" : "gen_") + s + "_" + p;
            const val = pt[key];
            if (val == null) { lines.push(labels[s] + ": —"); return; }
            const d = val - hist;
            lines.push(labels[s] + ": " + val.toFixed(prec) + unit +
              " (" + (d >= 0 ? "+" : "") + d.toFixed(prec) + ")");
          });
        }
        return { label: label, valueHtml: lines.join("<br>") };
      }

      // These are pre-rendered images, not queryable data -- no sampled
      // value to show (valueHtml: null tells the popup to skip that line
      // entirely rather than print an implementation detail like "PNG overlay").
      if (overlayId === "bbtforecast") return { label: label, valueHtml: null };
      if (overlayId === "bbtpast")     return { label: label, valueHtml: null };
      if (overlayId === "wind")        return { label: label, valueHtml: null };
      if (overlayId === "droughtFuture" || overlayId === "droughtRef" || overlayId === "fireWeather") return { label: label, valueHtml: null };

      const hexFeature = findContainingFeature(state.hexGeojson && state.hexGeojson[overlayId], latlng);
      if (hexFeature) {
        const val = Number(hexFeature.properties.mean);
        return { label: label, valueHtml: isFinite(val) ? "Pressure index " + val.toFixed(3) : "NA" };
      }
      return { label: label, valueHtml: "—" };
    }

    // Lightweight on-map popup mirroring the side-panel readout, so a click
    // gives an immediate value bubble at the point (like Esri map popups).
    // `state`/`map` default to the active side, but split-screen clicks also
    // call this with the parked panel + the OTHER map so both sides get their
    // own popup instead of only the actively-clicked one.
    // Every Natura 2000 site under the click (SAC and SPA often overlap).
    function natura2000SitesAt(latlng) {
      if (!APP.showNatura2000 || !APP.natura2000Data) return [];
      return APP.natura2000Data.features.filter(function(f) {
        return pointInFeature(latlng.lng, latlng.lat, f);
      });
    }

    function openClickPopup(latlng, state, map) {
      track("map-click");
      APP.lastClick = latlng; if (APP.showSeed4Forest) updateSeedPlace();
      state = state || APP;
      map = map || APP.map;
      const sites = natura2000SitesAt(latlng);
      // With nothing switched on the popup still opens: place, copy-for-Seed4Forest and the links
      const countryFeature = findContainingFeature(APP.countryGeojson, latlng);
      const countryName = countryFeature ? getCountryName(countryFeature) : "—";
      // One row per layer that is switched on, so a click on one pixel reads
      // them all together (top-most first).
      const ids = (state.visibleOverlays && state.visibleOverlays.length
        ? state.visibleOverlays.slice().reverse() : []);
      const asyncRows = [];
      let rows = "";
      let phenipsChart = false;
      ids.forEach(function(id) {
        const info = computeClickValueInfo(state, latlng, id);
        if (id === "phenips" && state.currentPhenipsScen === "obs" && state.currentPhenipsVar !== "share") {
          rows += '<span class="map-popup-layer">' + info.label + ' · by year</span>'
            + '<span class="map-popup-val phenips-chart" data-phenips-chart>Reading 44 years…</span>';
          phenipsChart = true;
          return;
        }
        const wantsAsync = info.valueHtml == null && CLICK_VALUE_OVERLAYS.includes(id);
        if (!wantsAsync && info.valueHtml == null) return;
        rows += '<span class="map-popup-layer">' + info.label + '</span>'
          + (wantsAsync ? '<span class="map-popup-val" data-async-val="' + id + '">Reading value…</span>'
                        : '<span class="map-popup-val">' + info.valueHtml + '</span>');
        if (wantsAsync) asyncRows.push(id);
      });
      if (sites.length) {
        rows += '<span class="map-popup-layer">Natura 2000</span>'
          + sites.slice(0, 4).map(function(f) {
              const q = f.properties;
              return '<span class="map-popup-val">' + q.name + '<small> · ' + q.designation + ' · ' + q.code
                + ' · <a href="https://biodiversity.europa.eu/sites/natura2000/' + encodeURIComponent(q.code)
                + '" target="_blank" rel="noopener noreferrer" title="Official EEA site page (Biodiversity Information System for Europe)">EEA site page ↗</a></small></span>';
            }).join("");
      }
      const html = '<div class="map-popup">'
        + '<div class="map-popup-head">' + (countryFeature ? countryName : "Map point") + '</div>'
        + '<div class="map-popup-body">'
        + 'Lat ' + latlng.lat.toFixed(3) + ', Lon ' + latlng.lng.toFixed(3) + '<br>'
        + rows
        + '<button type="button" class="risk-radar-btn">Hazard profile ⬡</button>'
        + placeLinksHtml(latlng)
        + '</div></div>';
      const popup = L.popup({
        autoPan: true,
        autoPanPaddingTopLeft: [12, isMobileLayout() ? 64 : 20],
        autoPanPaddingBottomRight: [12, isMobileLayout() ? 130 : 20],
        closeButton: true,
        className: "mini-popup",
        maxWidth: isMobileLayout() ? Math.min(340, window.innerWidth - 28) : 260,
        offset: [0, -4]
      }).setLatLng(latlng).setContent(html).openOn(map);

      // Snapshot THIS side's config now, so the button always samples
      // whichever side/settings were live when it was actually opened --
      // both maps in split screen behave independently this way.
      const radarState = snapshotRadarState(state);
      const popupEl = popup.getElement();
      asyncRows.forEach(function(id) {
        describeClickValue(id, latlng, radarState).then(function(text) {
          const slot = popupEl && popupEl.querySelector('[data-async-val="' + id + '"]');
          if (slot) slot.textContent = text || "No data at this point";
        });
      });
      const btn = popupEl && popupEl.querySelector(".risk-radar-btn");
      if (btn) {
        btn.addEventListener("click", function() { openRadarChart(latlng, radarState); });
      }
      if (phenipsChart) {
        const slot = popupEl && popupEl.querySelector("[data-phenips-chart]");
        phenipsSeriesAt(latlng, state).then(function(res) {
          if (slot) slot.innerHTML = res ? phenipsChartSvg(res) : "No data at this point";
        }).catch(function(err) { console.warn(err); if (slot) slot.textContent = "Could not read the yearly maps"; });
      }
    }

    // Links that open the same place in other tools. Only the satellite link carries the coordinates;
    // the others open the tool itself (none of them documents a coordinate link).
    // Seed4Forest's location box reads "longitude, latitude" (tested: it offers that reading first, and
    // "latitude, longitude" would be taken the wrong way round); then Enter and "Confirm Location".
    function seedCoordText(ll) { return ll.lng.toFixed(5) + ", " + ll.lat.toFixed(5); }
    function copyText(txt, done) {
      (navigator.clipboard ? navigator.clipboard.writeText(txt) : Promise.reject()).then(done, function() { window.prompt("Copy this place:", txt); });
    }
    document.addEventListener("click", function(e) {
      const b = e.target.closest && e.target.closest(".seed-copy-btn");
      if (!b) return;
      e.preventDefault();
      track("seed4forest-copy");
      copyText(b.dataset.coords, function() {
        const old = b.dataset.label || b.textContent; b.dataset.label = old; b.textContent = "Copied ✓";
        setTimeout(function() { b.textContent = old; }, 1800);
      });
    }, true);  // capture: Leaflet popups stop clicks from bubbling

    function placeLinksHtml(latlng) {
      const sat = "https://www.google.com/maps/@" + latlng.lat.toFixed(5) + "," + latlng.lng.toFixed(5) + ",15z/data=!3m1!1e3";
      const ch = latlng.lat > 45.8 && latlng.lat < 47.9 && latlng.lng > 5.9 && latlng.lng < 10.6;
      const a = function(href, text, tip) { return '<a href="' + href + '" target="_blank" rel="noopener noreferrer" title="' + tip + '">' + text + ' ↗</a>'; };
      return '<div class="map-popup-seed"><button type="button" class="seed-copy-btn" data-coords="' + seedCoordText(latlng) + '" title="Copies ' + seedCoordText(latlng) + '">Copy place for Seed4Forest</button>'
        + '<span>Paste in its location box, press Enter, then Confirm Location.</span></div>'
        + '<div class="map-popup-links">Open this place elsewhere: '
        + a(sat, "satellite view", "Satellite imagery centred on this point")
        + ' · ' + a("https://app.seed4forest.org/", "Seed4Forest", "BFW tool for tree species and seed sources (opens the tool; set the place there)")
        + (ch ? ' · ' + a("https://www.tree-app.ch", "Tree app (CH)", "Swiss tree-species recommendation tool") : "")
        + ' · ' + a("https://forest-fire.emergency.copernicus.eu", "EFFIS", "European Forest Fire Information System")
        + ' · ' + a("https://drought.emergency.copernicus.eu", "drought observatory", "European Drought Observatory")
        + '</div>';
    }

    // Yearly PHENIPS values at a point: read the 44 yearly map images at that pixel and invert the colour ramp.
    async function phenipsSeriesAt(latlng, state) {
      const man = await ensurePhenipsV5Manifest();
      const v = state.currentPhenipsVar, e = state.currentPhenipsExp;
      const base = "phenips_" + (v === "gdd" ? "dd" : "gen") + "_" + e;
      const lay = man.layers[base + "_2000"];
      const fn = makeColorToRiskFraction(buildRiskLut(lay.colours));
      const years = []; for (let y = 1981; y <= 2024; y++) years.push(y);
      const vals = await Promise.all(years.map(function(y) {
        return cachedFetchImage("webmap_data/phenips_v5/" + base + "_" + y + ".webp?v=" + ASSET_VERSION).then(function(img) {
          return sampleImageAt(img, PHENIPS_V5_BOUNDS, latlng, fn);
        }).then(function(t) { return t == null ? null : lay.min + t * (lay.max - lay.min); });
      }));
      if (vals.every(function(x) { return x == null; })) return null;
      return { years: years, vals: vals, unit: v === "gdd" ? "degree-days" : "generations", max: lay.max, dd: v === "gdd" };
    }

    function phenipsChartSvg(r) {
      const W = 230, H = 96, L = 28, R = 4, T = 6, B = 16;
      const ok = r.vals.filter(function(x) { return x != null; });
      const hi = Math.max(r.max * 0.4, Math.max.apply(null, ok)), lo = 0;
      const x = function(y) { return L + (y - 1981) / 43 * (W - L - R); };
      const yy = function(val) { return T + (1 - (val - lo) / (hi - lo)) * (H - T - B); };
      let d = "", pen = false;
      r.vals.forEach(function(val, i) {
        if (val == null) { pen = false; return; }
        d += (pen ? "L" : "M") + x(r.years[i]).toFixed(1) + " " + yy(val).toFixed(1); pen = true;
      });
      const base = r.vals.slice(0, 30).filter(function(q) { return q != null; });
      const mean = base.length ? base.reduce(function(a, b) { return a + b; }, 0) / base.length : null;
      const f = function(q) { return r.dd ? Math.round(q) : q.toFixed(1); };
      const dot = function(yr) {
        const q = r.vals[yr - 1981];
        return q == null ? "" : '<circle cx="' + x(yr).toFixed(1) + '" cy="' + yy(q).toFixed(1) + '" r="2.6" fill="#b0452a"/><text x="' + x(yr).toFixed(1) + '" y="' + (yy(q) - 4).toFixed(1) + '" font-size="8" text-anchor="middle" fill="#b0452a">' + yr + '</text>';
      };
      return '<svg viewBox="0 0 ' + W + ' ' + H + '" width="100%" role="img" aria-label="Yearly ' + r.unit + ' at this point, 1981 to 2024">'
        + '<line x1="' + L + '" y1="' + yy(0) + '" x2="' + (W - R) + '" y2="' + yy(0) + '" stroke="#c9d8e8"/>'
        + (mean != null ? '<line x1="' + L + '" y1="' + yy(mean).toFixed(1) + '" x2="' + (W - R) + '" y2="' + yy(mean).toFixed(1) + '" stroke="#7a8da3" stroke-dasharray="3 3"/>' : '')
        + '<path d="' + d + '" fill="none" stroke="#244166" stroke-width="1.5"/>' + dot(2003) + dot(2018)
        + '<text x="' + (L - 3) + '" y="' + (T + 6) + '" font-size="8" text-anchor="end" fill="#5a6b80">' + f(hi) + '</text>'
        + '<text x="' + (L - 3) + '" y="' + (yy(0) + 3) + '" font-size="8" text-anchor="end" fill="#5a6b80">0</text>'
        + [1981, 2000, 2024].map(function(yr) { return '<text x="' + x(yr) + '" y="' + (H - 3) + '" font-size="8" text-anchor="' + (yr === 2024 ? "end" : "middle") + '" fill="#5a6b80">' + yr + '</text>'; }).join("")
        + '</svg><small>' + r.unit + ' per year (' + (APP.currentPhenipsExp === "sun" ? "sun-exposed" : "shaded") + '); dashed = mean 1981–2010. Approximate: read from the map colours.</small>';
    }

    // ====================================================================
    // RISK RADAR CHART -- 5-hazard profile at a clicked point
    // ====================================================================
    // Bark beetle, wind, and drought are the hazards with a real model/
    // classification right now (see OVERLAY_CATALOG: wildfire is still
    // `available: false`, ash dieback has no layer at all -- it only
    // appears as a pilot-action hazard tag). Those two remaining axes are a
    // fixed 0.2 placeholder until a model exists; the other three are
    // sampled for real, INCLUDING when hidden -- a hidden layer is fetched
    // off-screen using whichever config this side currently holds (its
    // default if never touched, or whatever was last activated if the user
    // configured it before hiding it: hiding never resets
    // currentRcp/currentWindPeriod/currentDroughtFuturePeriod/etc, so
    // reading them already gives the right answer either way). Drought is a
    // classified severity tier (Extreme..Low), not a modelled probability
    // like the other two -- colorToDroughtSeverity maps it onto the same
    // 0..1 axis for comparability, it isn't literally "risk %".
    const RADAR_PLACEHOLDER = 0.2;
    const BBT_OVERLAY_IDS = ["bbtl", "bbtforecast", "bbtpast"];

    // Per-session cache so repeated radar samples don't re-fetch the same
    // off-screen image/geojson/manifest.
    const RADAR_FETCH_CACHE = new Map();
    function cachedFetchJson(url) {
      if (!RADAR_FETCH_CACHE.has(url)) {
        RADAR_FETCH_CACHE.set(url, fetch(url, { cache: "no-store" }).then(function(r) {
          if (!r.ok) throw new Error("HTTP " + r.status + " for " + url);
          return r.json();
        }));
      }
      return RADAR_FETCH_CACHE.get(url);
    }
    function cachedFetchImage(url) {
      if (!RADAR_FETCH_CACHE.has(url)) {
        RADAR_FETCH_CACHE.set(url, new Promise(function(resolve) {
          const img = new Image();
          img.onload = function() { resolve(img); };
          img.onerror = function() { resolve(null); };
          img.src = url;
        }));
      }
      return RADAR_FETCH_CACHE.get(url);
    }

    // Reverse the shared Spectral ramp LUT to recover a risk fraction (0..1)
    // from a raster's actual pixel color -- the same LUT the threshold
    // slider uses to mask pixels, just read instead of written.
    // Some of these rasters are coarse (the bark-beetle forecast PNG is only
    // 279px wide across all of Europe -- roughly 17 km/pixel), so an exact
    // click can land one pixel off a patchy mask's edge even though the
    // point is clearly "in" visible data at the zoom level a user is
    // looking at. Rather than failing on that single exact pixel, search a
    // small neighborhood and use the nearest pixel that actually has data.
    const SAMPLE_SEARCH_RADIUS = 3;

    // `.complete`/`.naturalWidth` can report true/set before the browser has
    // actually finished decoding the bitmap for drawing -- rare, but when it
    // happens drawImage silently paints nothing, which reads as "no data" at
    // a point that very much has data. img.decode() is the reliable signal;
    // await it (falling back to drawing anyway if it rejects) before
    // reading any pixels.
    async function sampleImageAt(img, bounds, latlng, colorFn) {
      if (!img || !img.naturalWidth) return null;
      const fracX = (latlng.lng - bounds.west) / (bounds.east - bounds.west);
      // The PNGs are Web-Mercator images stretched between their lat/lon
      // bounds, so the vertical position must be computed in Mercator space.
      const mercY = function(lat) { return Math.asinh(Math.tan(lat * Math.PI / 180)); };
      const fracY = (mercY(bounds.north) - mercY(latlng.lat)) / (mercY(bounds.north) - mercY(bounds.south));
      if (!(fracX >= 0 && fracX <= 1 && fracY >= 0 && fracY <= 1)) return null;
      const w = img.naturalWidth, h = img.naturalHeight;
      const px = Math.min(w - 1, Math.max(0, Math.floor(fracX * w)));
      const py = Math.min(h - 1, Math.max(0, Math.floor(fracY * h)));
      if (img.decode) { try { await img.decode(); } catch (e) { /* fall through and try anyway */ } }
      try {
        const c = document.createElement("canvas");
        c.width = w; c.height = h;
        const ctx = c.getContext("2d");
        ctx.drawImage(img, 0, 0);
        const x0 = Math.max(0, px - SAMPLE_SEARCH_RADIUS), x1 = Math.min(w - 1, px + SAMPLE_SEARCH_RADIUS);
        const y0 = Math.max(0, py - SAMPLE_SEARCH_RADIUS), y1 = Math.min(h - 1, py + SAMPLE_SEARCH_RADIUS);
        const region = ctx.getImageData(x0, y0, x1 - x0 + 1, y1 - y0 + 1).data;
        const rw = x1 - x0 + 1;
        let best = null, bestDist = Infinity;
        for (let yy = y0; yy <= y1; yy++) {
          for (let xx = x0; xx <= x1; xx++) {
            const idx = ((yy - y0) * rw + (xx - x0)) * 4;
            if (region[idx + 3] === 0) continue;
            const dist = (xx - px) * (xx - px) + (yy - py) * (yy - py);
            if (dist < bestDist) {
              bestDist = dist;
              best = { r: region[idx], g: region[idx + 1], b: region[idx + 2] };
            }
          }
        }
        if (!best) return null;
        return (colorFn || colorToRiskFraction)(best.r, best.g, best.b);
      } catch (e) {
        return null; // tainted canvas or similar -- treat as "no data"
      }
    }

    // Sample a PNG raster overlay's value at a point for a given side's
    // `state` (live APP for the active side, a parked panel snapshot for
    // the inactive one) -- whether or not that overlay is actually visible
    // right now. Prefers the real displayed <img> (cheap) when the layer
    // happens to already be loaded; otherwise fetches off-screen using
    // this side's own current model/period/threshold selection.
    async function sampleRasterOverlayFor(overlayId, latlng, state) {
      // droughtFuture's severity is already a native 0..1 reading (see
      // colorToDroughtSeverity above) -- unlike the other overlays here, it
      // isn't a colour-ramp fraction that still needs an overlayScaleCap.
      // wind has its own palette (WIND_RAMP) since it was split from bark
      // beetle's -- must use colorToWindFraction, not the default LUT.
      const isDrought = overlayId === "droughtFuture";
      const colorFn = isDrought ? colorToDroughtSeverity
        : overlayId === "wind" ? colorToWindFraction
        : colorToRiskFraction;
      const liveLayer = state.layerInstances && state.layerInstances[overlayId];
      if (liveLayer && liveLayer.getElement && liveLayer.getBounds) {
        const img = liveLayer.getElement();
        if (img && img.complete && img.naturalWidth) {
          const b = liveLayer.getBounds();
          const frac = await sampleImageAt(img, { west: b.getWest(), east: b.getEast(), north: b.getNorth(), south: b.getSouth() }, latlng, colorFn);
          if (frac != null) return isDrought ? frac : frac * overlayScaleCap(overlayId);
        }
      }

      if (overlayId === "wind") {
        let meta;
        try { meta = await ensureWindMetaLoaded(); } catch (e) { return null; }
        const gcm = state.currentWindPeriod === "historical" ? "default" : state.currentWindGcm;
        const key = gcm === "mpi" ? state.currentWindPeriod + "_mpi" : state.currentWindPeriod;
        const image = meta.images[key] || meta.images[state.currentWindPeriod] || meta.images.near;
        if (!image) return null;
        const img = await cachedFetchImage("webmap_data/" + image.file + "?v=" + ASSET_VERSION);
        const frac = await sampleImageAt(img, image.bounds, latlng, colorToWindFraction);
        return frac == null ? null : frac * overlayScaleCap("wind");
      }

      if (overlayId === "bbtforecast") {
        try { await ensureForecastManifestLoaded(); } catch (e) { return null; }
        const manifest = APP.forecastManifest;
        if (!manifest || !manifest.entries) return null;
        const entry = manifest.entries[state.currentForecastModel + "::" + state.currentForecastPeriod];
        if (!entry) return null;
        // The radar chart always wants an absolute risk reading (0..1), never
        // a delta -- a negative change doesn't fit its 0-100% axis semantic
        // -- so only the "threshold" view diverts to a classified raster.
        const tifPath = state.currentForecastView === "threshold"
          ? (entry.threshold_visible_tifs || {})[state.currentForecastThreshold]
          : entry.ensemble_probability_tif;
        if (!tifPath) return null;
        const imagePath = tifPath.replace(/\.tif$/i, ".png") + "?v=" + ASSET_VERSION;
        const metaPath = tifPath.replace(/\.tif$/i, "_meta.json");
        let overlayMeta;
        try { overlayMeta = await cachedFetchJson(metaPath); } catch (e) { return null; }
        const img = await cachedFetchImage(imagePath);
        if (state.currentForecastModel === P1C4_MODEL) {
          const f = await sampleImageAt(img, overlayMeta.bounds, latlng, colorToP1c4Fraction);
          return f == null ? null : f * P1C4_CAP;
        }
        return await sampleImageAt(img, overlayMeta.bounds, latlng); // already 0..1
      }

      if (overlayId === "bbtpast") {
        let meta;
        try { meta = await ensurePastOverlayMetaLoaded(); } catch (e) { return null; }
        if (!meta.image || !meta.bounds) return null;
        const img = await cachedFetchImage(meta.image + "?v=" + ASSET_VERSION);
        const ramp = meta.legend && meta.legend.gradient;
        return await sampleImageAt(img, meta.bounds, latlng, ramp ? makeColorToRiskFraction(buildRiskLut(ramp)) : undefined); // 0..1 of the 0-10 years scale
      }

      if (overlayId === "droughtFuture") {
        let meta;
        try { meta = await ensureDroughtFutureMetaLoaded(); } catch (e) { return null; }
        const image = meta.images[state.currentDroughtFuturePeriod];
        if (!image) return null;
        const img = await cachedFetchImage("webmap_data/" + image.file + "?v=" + ASSET_VERSION);
        return await sampleImageAt(img, image.bounds, latlng, colorToDroughtSeverity); // already 0..1 severity
      }

      return null;
    }

    // Reading for a point from any vector hex overlay (see HEX_OVERLAY_IDS),
    // for a given side's state: its own cached geojson if already loaded,
    // otherwise fetched off-screen (cached per-URL) using that side's own
    // current RCP/year selection.
    async function sampleHexOverlayFor(overlayId, latlng, state) {
      const cached = state.hexGeojson && state.hexGeojson[overlayId];
      if (cached) {
        const f = findContainingFeature(cached, latlng);
        if (!f) return null;
        const v = Number(f.properties.mean);
        return isFinite(v) ? v : null;
      }
      const meta = OVERLAY_CATALOG[overlayId].meta;
      const fileName = meta.files[state.currentRcp + "_" + state.currentYear];
      if (!fileName) return null;
      let geojson;
      try { geojson = await cachedFetchJson(fileName); } catch (e) { return null; }
      const f = findContainingFeature(geojson, latlng);
      if (!f) return null;
      const v = Number(f.properties.mean);
      return isFinite(v) ? v : null;
    }

    // Reading of a classified (legend-coloured) PNG at a point: the pixel
    // colour is matched to the layer's own legend, so what a user reads in a
    // popup is exactly the class shown in the legend.
    const CLICK_VALUE_OVERLAYS = ["droughtFuture", "droughtRef", "fireWeather", "wind", "bbtforecast"];

    function hexToRgb(hex) {
      const h = hex.replace("#", "");
      return [parseInt(h.slice(0, 2), 16), parseInt(h.slice(2, 4), 16), parseInt(h.slice(4, 6), 16)];
    }

    async function classifiedSpecFor(overlayId, state, forRadar) {
      if (overlayId === "droughtFuture") {
        const m = await ensureDroughtFutureMetaLoaded();
        const im = m.images[state.currentDroughtFuturePeriod];
        return im ? { file: im.file, bounds: im.bounds, legend: im.legend } : null;
      }
      if (overlayId === "droughtRef") {
        const m = await ensureDroughtRefMetaLoaded();
        const im = m.images[state.currentDroughtRefMetric];
        return im ? { file: im.file, bounds: im.bounds, legend: im.legend } : null;
      }
      if (overlayId === "fireWeather") {
        const m = await ensureFireMetaLoaded();
        const key = forRadar ? "level::" + state.currentFireScen : fireImageKey(state);
        const im = m.images[key];
        if (!im) return null;
        return { file: im.file, bounds: im.bounds, legend: key.indexOf("change::") === 0 ? m.legends.change : m.legends.level };
      }
      return null;
    }

    async function readLegendClass(overlayId, latlng, state, forRadar) {
      let spec;
      try { spec = await classifiedSpecFor(overlayId, state, forRadar); } catch (e) { return null; }
      if (!spec) return null;
      const rgb = spec.legend.colors.map(hexToRgb);
      const img = await cachedFetchImage("webmap_data/" + spec.file + "?v=" + ASSET_VERSION);
      const idx = await sampleImageAt(img, spec.bounds, latlng, function(r, g, b) {
        let best = 0, bestD = Infinity;
        rgb.forEach(function(c, i) {
          const d = (c[0] - r) * (c[0] - r) + (c[1] - g) * (c[1] - g) + (c[2] - b) * (c[2] - b);
          if (d < bestD) { bestD = d; best = i; }
        });
        return best;
      });
      if (idx == null) return null;
      return { index: idx, count: spec.legend.labels.length, label: spec.legend.labels[idx] };
    }

    const shortLabel = function(label) { return String(label).replace(/\s*\(.*$/, "").replace(/ drought hazard$/i, ""); };
    const BB_CLASSES = [[0.25, "low"], [0.5, "moderate"], [0.75, "high"], [Infinity, "very high"]];
    const WIND_CLASSES = [[0.25, "low"], [0.5, "moderate"], [Infinity, "high"]];
    const classOf = function(table, v) { return table.find(function(row) { return v < row[0]; })[1]; };

    // Share-valued extra layers of the final model at a point (0..1): "extrapolation" (greys) or "share05" (hazard ramp).
    async function sampleForecastExtra(latlng, state, kind) {
      try { await ensureForecastManifestLoaded(); } catch (e) { return null; }
      const entry = APP.forecastManifest.entries[state.currentForecastModel + "::" + state.currentForecastPeriod];
      const tif = entry && entry[kind + "_tif"];
      if (!tif) return null;
      try {
        const meta = await cachedFetchJson(tif.replace(/\.tif$/i, "_meta.json"));
        const img = await cachedFetchImage(tif.replace(/\.tif$/i, ".png") + "?v=" + ASSET_VERSION);
        return await sampleImageAt(img, meta.bounds, latlng, kind === "extrapolation" ? colorToExtrapFraction : colorToP1c4Fraction);
      } catch (e) { return null; }
    }

    async function describeClickValue(overlayId, latlng, state) {
      if (overlayId === "wind") {
        const v = await sampleRasterOverlayFor("wind", latlng, state);
        return v == null ? null : "Damage probability " + v.toFixed(2) + (v >= 0.5 - 1e-6 ? " or more" : "") + " (" + classOf(WIND_CLASSES, v) + ")";
      }
      if (overlayId === "bbtforecast") {
        const v = await sampleRasterOverlayFor("bbtforecast", latlng, state);
        if (v == null) return null;
        if (state.currentForecastModel !== P1C4_MODEL) return "Hazard " + v.toFixed(2) + " (" + classOf(BB_CLASSES, v) + ")";
        // Final model: the extrapolation share and the share of years at or above 0.5 sit beside the hazard when delivered
        const [ex, sh] = await Promise.all([sampleForecastExtra(latlng, state, "extrapolation"), sampleForecastExtra(latlng, state, "share05")]);
        return "Hazard " + v.toFixed(2) + " (" + classOf(P1C4_CLASSES, v) + ")"
          + (sh == null ? "" : "; 0.5 or more in " + Math.round(sh * 100) + "% of years")
          + (ex == null ? "" : "; " + Math.round(ex * 100) + "% of climate predictors outside the training range");
      }
      const c = await readLegendClass(overlayId, latlng, state, false);
      return c ? c.label : null;
    }

    // Fields relevant to radar sampling, snapshotted at click/popup time so
    // a button's behavior stays tied to whichever side (and whatever it was
    // configured to) it was actually opened from, even if the user keeps
    // interacting with the map(s) afterward -- split screen included.
    const RADAR_STATE_FIELDS = [
      "focusedOverlay", "layerInstances", "hexGeojson",
      "currentRcp", "currentYear",
      "currentForecastModel", "currentForecastPeriod", "currentForecastView", "currentForecastThreshold",
      "currentWindPeriod", "currentWindGcm", "currentDroughtFuturePeriod", "currentDroughtRefMetric",
      "currentFireScen", "currentFireView"
    ];
    function snapshotRadarState(state) {
      const snap = {};
      RADAR_STATE_FIELDS.forEach(function(f) { snap[f] = state[f]; });
      return snap;
    }

    // One axis per driver that has a real layer. Each axis shows the class or
    // index of that driver's own layer, rescaled to 0..1 only to draw the
    // shape -- no invented placeholder values, no single integrated number.
    async function computeRadarValues(latlng, radarState) {
      const state = radarState || snapshotRadarState(APP);
      const [bb, wind, drought, fire] = await Promise.all([
        sampleRasterOverlayFor("bbtforecast", latlng, state),
        sampleRasterOverlayFor("wind", latlng, state),
        readLegendClass("droughtFuture", latlng, state, true),
        readLegendClass("fireWeather", latlng, state, true)
      ]);
      const fromClass = function(c) { return c ? (c.count - c.index) / c.count : null; };
      return [
        { id: "bbtl", label: "Bark beetle", value: bb == null ? null : Math.min(1, bb / (state.currentForecastModel === P1C4_MODEL ? P1C4_CAP : 1)),
          text: bb == null ? null : "hazard " + bb.toFixed(2) + " (" + classOf(state.currentForecastModel === P1C4_MODEL ? P1C4_CLASSES : BB_CLASSES, bb) + ")" },
        { id: "wind", label: "Wind", value: wind == null ? null : Math.min(1, wind / 0.5),
          text: wind == null ? null : "damage probability " + wind.toFixed(2) + " (" + classOf(WIND_CLASSES, wind) + ")" },
        { id: "fire", label: "Fire", value: fire ? (fire.index + 1) / fire.count : null,
          text: fire ? "summer fire weather " + shortLabel(fire.label).toLowerCase() : null },
        { id: "drought", label: "Drought", value: fromClass(drought),
          text: drought ? shortLabel(drought.label).toLowerCase() + " hazard" : null }
      ];
    }

    function radarPolarPoint(cx, cy, r, angleDeg) {
      const a = (angleDeg - 90) * Math.PI / 180;
      return [cx + r * Math.cos(a), cy + r * Math.sin(a)];
    }

    function renderRadarSvg(values) {
      const cx = 120, cy = 112, R = 82;
      const n = values.length;
      const step = 360 / n;
      const rings = [0.2, 0.4, 0.6, 0.8, 1.0];

      let svg = '<svg viewBox="-34 0 308 224" width="100%" height="auto" class="radar-svg" style="max-width:260px">';

      rings.forEach(function(frac) {
        const pts = values.map(function(_, i) {
          const p = radarPolarPoint(cx, cy, R * frac, i * step);
          return p[0].toFixed(1) + "," + p[1].toFixed(1);
        }).join(" ");
        svg += '<polygon points="' + pts + '" class="radar-grid-ring" />';
      });

      values.forEach(function(v, i) {
        const p = radarPolarPoint(cx, cy, R, i * step);
        svg += '<line x1="' + cx + '" y1="' + cy + '" x2="' + p[0].toFixed(1)
          + '" y2="' + p[1].toFixed(1) + '" class="radar-spoke" />';
        const lp = radarPolarPoint(cx, cy, R + 20, i * step);
        svg += '<text x="' + lp[0].toFixed(1) + '" y="' + lp[1].toFixed(1)
          + '" class="radar-axis-label" text-anchor="middle">' + v.label + '</text>';
      });

      const frac = function(v) { return v.value == null ? 0 : Math.max(0, Math.min(1, v.value)); };
      const points = values.map(function(v, i) { return radarPolarPoint(cx, cy, R * frac(v), i * step); });
      const poly = points.map(function(p) { return p[0].toFixed(1) + "," + p[1].toFixed(1); }).join(" ");
      // Fill is one shape (fine as a rough gestalt), but the OUTLINE is drawn
      // as separate segments below: a solid line only where both endpoints
      // are real measurements, dashed/muted wherever it touches a
      // placeholder (no model exists, e.g. fire/ash) or missing (a model
      // exists but this point has no data) value -- a single solid outline
      // through fabricated and real numbers alike would misrepresent
      // placeholder guesses as equally-measured.
      svg += '<polygon points="' + poly + '" class="radar-value-poly" />';
      const isConfident = function(v) { return !v.placeholder && v.value != null; };
      for (let i = 0; i < n; i++) {
        const j = (i + 1) % n;
        const segCls = (isConfident(values[i]) && isConfident(values[j])) ? "radar-seg-confident" : "radar-seg-unconfident";
        svg += '<line x1="' + points[i][0].toFixed(1) + '" y1="' + points[i][1].toFixed(1)
          + '" x2="' + points[j][0].toFixed(1) + '" y2="' + points[j][1].toFixed(1)
          + '" class="' + segCls + '" />';
      }

      values.forEach(function(v, i) {
        const p = radarPolarPoint(cx, cy, R * frac(v), i * step);
        const cls = v.placeholder ? "radar-dot-placeholder" : (v.value == null ? "radar-dot-missing" : "radar-dot-real");
        svg += '<circle cx="' + p[0].toFixed(1) + '" cy="' + p[1].toFixed(1) + '" r="4" class="radar-dot ' + cls + '" />';
      });

      svg += '</svg>';
      return svg;
    }

    function ensureRadarPanel() {
      if (APP.radarPanelEl) return APP.radarPanelEl;
      const el = document.createElement("div");
      el.id = "riskRadarPanel";
      el.className = "risk-radar-panel";
      document.body.appendChild(el);
      APP.radarPanelEl = el;
      return el;
    }

    function closeRadarChart() {
      if (APP.radarPanelEl) APP.radarPanelEl.classList.remove("is-open");
    }

    function radarHeaderHtml(latlng) {
      return '<div class="risk-radar-header"><strong>Hazard profile</strong>'
        + '<span class="risk-radar-coords">' + latlng.lat.toFixed(3) + ', ' + latlng.lng.toFixed(3) + '</span>'
        + '<button type="button" class="risk-radar-close" aria-label="Close">&times;</button></div>';
    }

    // `radarState` is a snapshot from snapshotRadarState (see attachMapClick)
    // identifying which side's config to sample hidden layers with; falls
    // back to the live active side if called without one.
    //
    // Sampling is async (off-screen fetches for hidden layers can take a
    // moment), so clicking a second point before the first finishes must
    // not let the slower, stale request clobber the newer one once both
    // resolve -- guarded the same way the app's other async overlay loads
    // already are, with a monotonic request token checked before render.
    async function openRadarChart(latlng, radarState) {
      const token = ++APP.radarRequestToken;
      const el = ensureRadarPanel();
      el.innerHTML = radarHeaderHtml(latlng) + '<div class="risk-radar-body"><div class="risk-radar-loading">Sampling hazard layers…</div></div>';
      el.classList.add("is-open");
      el.querySelector(".risk-radar-close").addEventListener("click", closeRadarChart);

      const values = await computeRadarValues(latlng, radarState);
      if (token !== APP.radarRequestToken) return; // a newer click superseded this one
      if (!el.classList.contains("is-open")) return; // closed while sampling

      const rows = values.map(function(v) {
        return '<div class="radar-legend-row"><span class="radar-legend-swatch"></span>' + v.label
          + ': <strong>' + (v.text || "no data at this point") + '</strong></div>';
      }).join("") + '<div class="radar-legend-row"><span class="radar-legend-swatch is-placeholder"></span>Ash dieback: <strong>no layer yet</strong></div>';
      const hasUnconfident = values.some(function(v) { return v.value == null; });
      const edgeNote = '<div class="radar-legend-note">Each axis shows the class or index of that driver’s own layer (settings as currently chosen), so the shape compares drivers loosely and is not a combined risk value.'
        + (hasUnconfident ? ' Dashed edges touch a driver with no data at this point.' : '') + '</div>';
      el.innerHTML = radarHeaderHtml(latlng) + '<div class="risk-radar-body">' + renderRadarSvg(values)
        + '<div class="risk-radar-legend">' + rows + edgeNote + '</div></div>';
      el.querySelector(".risk-radar-close").addEventListener("click", closeRadarChart);
    }

    // ====================================================================
    // BASEMAPS
    // ====================================================================
    // Basemap options. `build` returns a fresh tile layer (or null for the
    // plain border-only style); `bg` sets the map background behind the tiles.
    const CARTO_ATTRIB = "&copy; OpenStreetMap contributors &copy; CARTO";
    const BASEMAPS = {
      osm: {
        label: "Streets",
        bg: "#aadaff",
        build: function() {
          return L.tileLayer("https://{s}.tile.openstreetmap.org/{z}/{x}/{y}.png", {
            maxZoom: 19,
            attribution: "&copy; OpenStreetMap contributors"
          });
        }
      },
      light: {
        label: "Light",
        bg: "#f3f3f1",
        build: function() {
          return L.tileLayer("https://{s}.basemaps.cartocdn.com/light_all/{z}/{x}/{y}{r}.png", {
            subdomains: "abcd", maxZoom: 20, attribution: CARTO_ATTRIB
          });
        }
      },
      dark: {
        label: "Dark",
        bg: "#1b1b1b",
        build: function() {
          return L.tileLayer("https://{s}.basemaps.cartocdn.com/dark_all/{z}/{x}/{y}{r}.png", {
            subdomains: "abcd", maxZoom: 20, attribution: CARTO_ATTRIB
          });
        }
      },
      plain: {
        label: "Plain (borders only)",
        bg: "#d9dde0",
        build: function() { return null; }
      }
    };

    function setBasemap(key) {
      const def = BASEMAPS[key] || BASEMAPS.osm;
      APP.currentBasemap = BASEMAPS[key] ? key : "osm";

      // Active map.
      if (APP.basemapLayer) { APP.map.removeLayer(APP.basemapLayer); APP.basemapLayer = null; }
      const mapEl = document.getElementById("map");
      if (mapEl) mapEl.style.background = def.bg;
      const layer = def.build();
      if (layer) { layer.addTo(APP.map); if (layer.bringToBack) layer.bringToBack(); APP.basemapLayer = layer; }

      // Inactive (parked) map, if split.
      if (APP.split) {
        const inactive = APP.activeSide === "left" ? "right" : "left";
        const pr = APP.panels[inactive];
        if (pr && pr.map) {
          if (pr.basemapLayer) pr.map.removeLayer(pr.basemapLayer);
          const l2 = def.build();
          if (l2) { l2.addTo(pr.map); if (l2.bringToBack) l2.bringToBack(); }
          pr.basemapLayer = l2;
        }
        const map2El = document.getElementById("map2");
        if (map2El) map2El.style.background = def.bg;
      }

      if (APP.basemapControlEl) {
        APP.basemapControlEl.querySelectorAll("button[data-basemap]").forEach(function(btn) {
          btn.classList.toggle("is-on", btn.dataset.basemap === APP.currentBasemap);
        });
      }
    }

    function addBasemapControl() {
      const control = L.control({ position: "topright" });
      control.onAdd = function() {
        const el = L.DomUtil.create("div", "basemap-control");
        let html = '<span class="basemap-title">Base map</span>';
        Object.keys(BASEMAPS).forEach(function(key) {
          html += '<button type="button" data-basemap="' + key + '">' + BASEMAPS[key].label + '</button>';
        });
        el.innerHTML = html;
        L.DomEvent.disableClickPropagation(el);
        L.DomEvent.disableScrollPropagation(el);
        el.querySelectorAll("button[data-basemap]").forEach(function(btn) {
          btn.addEventListener("click", function() {
            setBasemap(btn.dataset.basemap);
          });
        });
        APP.basemapControlEl = el;
        return el;
      };
      control.addTo(APP.map);
    }

    // Split-screen toggle belongs next to the thing it splits (the map),
    // not in the header row of unrelated project/info links. Moves the
    // existing #splitBtn node (with all its listeners already wired) into
    // a Leaflet control instead of recreating it.
    function refreshMapSize() {
      setTimeout(function() {
        if (APP.map) APP.map.invalidateSize();
        const right = APP.split && APP.panels && (APP.panels.right || APP.panels.left);
        if (right && right.map && right.map.invalidateSize) right.map.invalidateSize();
      }, 60);
    }
    // Three independent folds, each with its own arrow in its own gutter:
    //   top panel (banner + DSS components) -> fold-top
    //   driver panel (layer ribbon)         -> fold-drivers
    //   side panel (map inspector)          -> fold-side
    const FOLDS = {
      "fold-top": { btn: "foldTopBtn", what: "top panel", gutter: "gutterTop", name: "RE-ENFORCE DSS" },
      "fold-drivers": { btn: "expandMapBtn", what: "driver panel", gutter: "gutterMap", name: "Disturbance drivers" },
      "fold-side": { btn: "foldSideBtn", what: "side panel", gutter: "gutterSide", side: true, name: "Map inspector" }
    };
    // Phone layout (<= 700 px): panels are bottom sheets, so the desktop fold arrows do not apply.
    const MOBILE_MQ = window.matchMedia("(max-width: 700px)");
    function isMobileLayout() { return MOBILE_MQ.matches; }
    function setFold(cls, on) {
      if (on && isMobileLayout()) return;
      document.body.classList.toggle(cls, on);
      const f = FOLDS[cls];
      const btn = document.getElementById(f.btn);
      if (btn) {
        btn.setAttribute("aria-pressed", String(!on));
        btn.title = (on ? "Show " : "Hide ") + f.what;
        btn.setAttribute("aria-label", btn.title);
      }
      refreshMapSize();
    }
    function setMapExpanded(on) { Object.keys(FOLDS).forEach(function(c) { setFold(c, on); }); }
    const FOLD_ARROW = '<svg viewBox="0 0 18 14" aria-hidden="true"><path d="M3 9.5 9 3.5l6 6" fill="none" stroke="currentColor" stroke-width="2.4" stroke-linecap="round" stroke-linejoin="round"/></svg>';
    function addMenuToggleTab() {
      Object.keys(FOLDS).forEach(function(cls) {
        const f = FOLDS[cls];
        const host = document.getElementById(f.gutter);
        if (!host) return;
        const tab = document.createElement("button");
        tab.type = "button";
        tab.id = f.btn;
        tab.className = "menu-toggle-tab" + (f.side ? " is-side" : "");
        tab.setAttribute("aria-pressed", "true");
        tab.innerHTML = FOLD_ARROW;
        tab.title = "Hide " + f.what;
        tab.setAttribute("aria-label", tab.title);
        tab.addEventListener("click", function(e) {
          e.stopPropagation();
          setFold(cls, !document.body.classList.contains(cls));
        });
        const label = document.createElement("span");
        label.className = "fold-label";
        label.textContent = f.name;
        label.title = "Show / hide: " + f.name;
        label.addEventListener("click", function() { tab.click(); });
        if (f.side) { host.appendChild(tab); host.appendChild(label); }
        else { host.appendChild(label); host.appendChild(tab); }
      });
    }

    document.addEventListener("keydown", function(e) {
      if (e.key === "Escape") {
        document.body.classList.remove("m-menu", "m-layers", "m-legend");
        if (APP.map) APP.map.closePopup();
        if (APP.panels && APP.panels.right && APP.panels.right.map) APP.panels.right.map.closePopup();
        closeRadarChart();
      }
      if (e.key === "Escape" && Object.keys(FOLDS).some(function(c) { return document.body.classList.contains(c); })) setMapExpanded(false);
    });

    function initMobileLayout() {
      const body = document.body;
      const menuBtn = document.getElementById("mobMenuBtn");
      const sheet = function(name) {
        body.classList.toggle("m-layers", name === "layers");
        body.classList.toggle("m-legend", name === "legend");
        body.classList.remove("m-menu");
        if (menuBtn) menuBtn.setAttribute("aria-expanded", "false");
        if (name === "legend") requestAnimationFrame(function() {
          const p = document.getElementById("panel"), box = p && p.querySelector(".legend-box"), card = box && box.closest(".panel-card");
          const head = p && p.querySelector(".m-sheet-head");
          if (card) p.scrollTop += card.getBoundingClientRect().top - p.getBoundingClientRect().top - (head ? head.offsetHeight : 0) - 6;
        });
      };
      [[".layer-ribbon", "Layers", "layers"], ["#panel", "Map inspector and legend", "legend"]]
        .forEach(function(def) {
          const host = document.querySelector(def[0]);
          if (!host) return;
          const head = document.createElement("div");
          head.className = "m-sheet-head";
          head.innerHTML = "<span></span><button type=\"button\" aria-label=\"Close\">&times;</button>";
          head.firstChild.textContent = def[1];
          head.lastChild.addEventListener("click", function() { sheet(null); });
          host.insertBefore(head, host.firstChild);
        });
      if (menuBtn) menuBtn.addEventListener("click", function() {
        const open = !body.classList.contains("m-menu");
        sheet(null);
        body.classList.toggle("m-menu", open);
        menuBtn.setAttribute("aria-expanded", String(open));
      });
      APP.openLayersSheet = function() { sheet("layers"); };
      document.getElementById("mobLayersBtn").addEventListener("click", function() { sheet("layers"); });
      document.getElementById("mobLegendBtn").addEventListener("click", function() { sheet("legend"); });
      if (DOM.shareBtn) DOM.shareBtn.addEventListener("click", function() { body.classList.remove("m-menu"); });
      // Leaving the phone layout: forget phone-only state and any desktop folds that were skipped.
      MOBILE_MQ.addEventListener("change", function() {
        body.classList.remove("m-menu", "m-layers", "m-legend");
        Object.keys(FOLDS).forEach(function(c) { body.classList.remove(c); });
        refreshMapSize();
      });
    }

    // "Show my location": asks the browser for the position (it stays in the browser, nothing is sent), drops a pin,
    // zooms in and opens the usual readout for that spot. The pin stays until removed with the small x button.
    function addLocateControl() {
      const control = L.control({ position: "topleft" });
      control.onAdd = function() {
        const bar = L.DomUtil.create("div", "leaflet-bar locate-control");
        L.DomEvent.disableClickPropagation(bar);
        const go = L.DomUtil.create("a", "locate-go", bar);
        go.href = "#"; go.setAttribute("role", "button"); go.setAttribute("aria-label", "Show my location");
        go.title = "Show my location (your position stays in your browser)";
        go.innerHTML = '<svg width="18" height="18" viewBox="0 0 24 24" fill="none" stroke="#173F69" stroke-width="2" stroke-linecap="round"><circle cx="12" cy="12" r="6"/><circle cx="12" cy="12" r="1.8" fill="#173F69"/><path d="M12 1v4M12 19v4M1 12h4M19 12h4"/></svg>';
        const off = L.DomUtil.create("a", "locate-off", bar);
        off.href = "#"; off.setAttribute("role", "button"); off.setAttribute("aria-label", "Remove my location pin"); off.title = "Remove my location pin"; off.textContent = "×"; off.hidden = true;
        L.DomEvent.on(go, "click", function(e) { L.DomEvent.preventDefault(e); locateMe(off); });
        L.DomEvent.on(off, "click", function(e) { L.DomEvent.preventDefault(e); clearMyLocation(off); });
        return bar;
      };
      control.addTo(APP.map);
    }

    function clearMyLocation(offBtn) {
      ["myLocationPin", "myLocationRing"].forEach(function(k) { if (APP[k]) { APP.map.removeLayer(APP[k]); APP[k] = null; } });
      APP.myLatLng = null;
      if (offBtn) offBtn.hidden = true;
    }

    // Short message on the map (setStatus only logs to the console)
    function locateNote(text) {
      let el = document.getElementById("locateNote");
      if (!el) { el = document.createElement("div"); el.id = "locateNote"; el.setAttribute("role", "status"); document.getElementById("map").appendChild(el); }
      el.textContent = text; el.hidden = !text;
      clearTimeout(locateNote.t); if (text) locateNote.t = setTimeout(function() { el.hidden = true; }, 7000);
    }

    function locateMe(offBtn, auto) {
      if (!navigator.geolocation) { if (!auto) locateNote("This browser cannot give a location."); return; }
      track(auto ? "locate-auto" : "locate");
      if (!auto) locateNote("Asking for your location…");
      navigator.geolocation.getCurrentPosition(function(pos) {
        const latlng = L.latLng(pos.coords.latitude, pos.coords.longitude);
        locateNote("");
        if (!L.latLngBounds(MAP_SOFT_BOUNDS).contains(latlng)) { if (!auto) locateNote("Your location is outside the area this map covers."); return; }
        clearMyLocation(offBtn);
        APP.myLatLng = latlng;
        const acc = pos.coords.accuracy;
        if (isFinite(acc) && acc < 5000) APP.myLocationRing = L.circle(latlng, { radius: acc, color: "#1a73e8", weight: 1, fillColor: "#1a73e8", fillOpacity: 0.1, interactive: false }).addTo(APP.map);
        APP.myLocationPin = L.circleMarker(latlng, { radius: 8, color: "#ffffff", weight: 3, fillColor: "#1a73e8", fillOpacity: 1, className: "my-location-dot" }).addTo(APP.map);
        APP.myLocationPin.on("click", function() { openClickPopup(latlng, APP, APP.map); });
        if (offBtn) offBtn.hidden = false;
        if (!auto) APP.map.panTo(latlng);   // never zoom
        if (auto) { locateNote("Blue dot: your position. × removes it."); return; }
        updatePanelFromClick(latlng);
        openClickPopup(latlng, APP, APP.map);
      }, function(err) {
        if (auto) return;   // a refusal at start is silent; the button stays available
        locateNote(err.code === 1 ? "Location permission was denied. Allow it in your browser's site settings to use this button."
          : err.code === 3 ? "Finding your location took too long. Try again." : "Your location is not available right now.");
      }, { enableHighAccuracy: false, timeout: 12000, maximumAge: 60000 });
    }

    // At start: show where you are. Asks once per browser; if the browser already allows it, marks the position on every visit.
    // A shared link keeps its own view and is never moved.
    function autoLocateOnce() {
      if (location.hash.length > 1 || !navigator.geolocation) return;
      const off = document.querySelector(".locate-off");
      let asked = false;
      try { asked = localStorage.getItem("locateAsked") === "1"; } catch (e) { /* storage blocked: treat as not asked */ }
      const ask = function() {
        try { localStorage.setItem("locateAsked", "1"); } catch (e) { /* ignore */ }
        locateMe(off, true);
      };
      if (navigator.permissions && navigator.permissions.query) {
        navigator.permissions.query({ name: "geolocation" }).then(function(r) {
          if (r.state === "granted") locateMe(off, true);
          else if (r.state === "prompt" && !asked) ask();
        }, function() { if (!asked) ask(); });
      } else if (!asked) ask();
    }

    function addSplitScreenControl() {
      const control = L.control({ position: "topright" });
      control.onAdd = function() {
        const el = L.DomUtil.create("div", "split-screen-control");
        L.DomEvent.disableClickPropagation(el);
        L.DomEvent.disableScrollPropagation(el);
        if (DOM.splitBtn) {
          DOM.splitBtn.hidden = false;
          el.appendChild(DOM.splitBtn);
        }
        return el;
      };
      control.addTo(APP.map);
    }

    // --- Split screen: two real maps; the panel drives the ACTIVE one -------
    // The active side's state lives directly on APP (map, layerInstances, all
    // settings). Switching side parks the current state in `panels` and loads
    // the other side's state, so the whole existing pipeline only ever touches
    // the active map.

    // ====================================================================
    // SPLIT SCREEN — two synced maps; the panel drives the ACTIVE one
    // ====================================================================
    const PANEL_FIELDS = [
      "map", "basemapLayer", "clickMarker", "countryHaloLayer", "countryLayer",
      "selectedHexLayer", "selectedCountryLayer", "currentYear", "currentRcp",
      "riskThreshold", "visibleOverlays", "focusedOverlay", "layerOpacity",
      "currentPhenipsScen", "currentPhenisPeriod", "currentPhenipsVar",
      "currentPhenipsView", "currentPhenipsExp", "currentForecastModel", "currentForecastPeriod",
      "currentForecastView", "currentForecastThreshold", "currentWindPeriod", "currentWindGcm",
      "currentDroughtFuturePeriod", "currentDroughtRefMetric", "currentFireView", "currentFireScen",
      "currentForecastRasterRange", "overlayLoadTokens", "overlayVisible",
      "currentHexGeojson", "layerInstances"
    ];

    function savePanel(side) {
      const o = {};
      PANEL_FIELDS.forEach(function(f) { o[f] = APP[f]; });
      APP.panels[side] = o;
    }
    function loadPanel(side) {
      const o = APP.panels[side];
      if (!o) return;
      PANEL_FIELDS.forEach(function(f) { APP[f] = o[f]; });
    }

    // Soft edge: the map can be dragged past the Europe frame, but it pulls
    // back (0 = free, 1 = hard wall) and springs back when released.
    const MAP_SOFT_BOUNDS = [[24.0, -40.0], [80.0, 60.0]];
    const MAP_BOUNDS_VISCOSITY = 0.4;

    function createMapPanes(map) {
      getOverlayIds().forEach(function(overlayId, index) {
        const paneName = getOverlayPaneName(overlayId);
        map.createPane(paneName);
        map.getPane(paneName).style.zIndex = String(300 + index * 20);
      });
      map.createPane("countryHaloPane");
      map.getPane("countryHaloPane").style.zIndex = 700;
      map.createPane("countryPane");
      map.getPane("countryPane").style.zIndex = 710;
      // Natura 2000 polygons hit the same country-border click-theft bug as
      // everything else interactive on this map -- own pane above it.
      map.createPane("natura2000Pane");
      map.getPane("natura2000Pane").style.zIndex = 712;
      // Catalogue-event dots must sit above the (interactive) country-border
      // pane, or the border layer's SVG intercepts their clicks -- this bit
      // hardest near dense multi-country borders (e.g. Germany), where many
      // short border segments cover a wide swath of the map.
      map.createPane("catalogueEventsPane");
      map.getPane("catalogueEventsPane").style.zIndex = 715;
      // Pilot-action pins hit the exact same country-border click-theft bug
      // (they used Leaflet's default markerPane, z-index 600 -- well below
      // countryPane's 710) -- same fix, own pane above it.
      map.createPane("pilotActionsPane");
      map.getPane("pilotActionsPane").style.zIndex = 716;
      // Keep click popups above the country borders.
      if (map.getPane("popupPane")) map.getPane("popupPane").style.zIndex = 800;
    }

    function buildCountryLayersOn(map) {
      if (!APP.countryGeojson) return { countryHaloLayer: null, countryLayer: null };
      const halo = L.geoJSON(APP.countryGeojson, {
        pane: "countryHaloPane", interactive: false, smoothFactor: 0, noClip: true,
        style: function() { return { color: "#ffffff", weight: 3.4, opacity: 0.85, fillOpacity: 0, lineCap: "round", lineJoin: "round" }; }
      }).addTo(map);
      const line = L.geoJSON(APP.countryGeojson, {
        pane: "countryPane", interactive: true, smoothFactor: 0, noClip: true,
        style: function() { return { color: "#29465b", weight: 1.05, opacity: 0.95, fillOpacity: 0, lineCap: "round", lineJoin: "round" }; }
      }).addTo(map);
      return { countryHaloLayer: halo, countryLayer: line };
    }

    // Build the pilot-action pins for exactly one map instance and return
    // the layer group. Each side (left/right) gets its OWN group tied to
    // its OWN map -- never cross-added, so split screen never duplicates.
    function buildPilotLayer(map) {
      const group = L.layerGroup();
      PILOT_ACTIONS.forEach(function(p) {
        const titleText = p.pa + " " + p.label + " (" + p.country + "): " + p.type;
        const icon = L.divIcon({
          className: "",
          html: '<div class="pilot-marker-hit" title="' + titleText + '">'
            + '<div class="pilot-marker-icon"><img src="webmap_data/pilot_icons/' + p.hazard + '.png?v=' + ASSET_VERSION + '" alt="' + p.pa + '" /></div>'
            + '<div class="pilot-marker-label">' + p.pa + '</div>'
            + '</div>',
          iconSize: [36, 36],
          iconAnchor: [18, 18],
          popupAnchor: [0, -8]
        });
        const coordNote = p.approx
          ? "Approx. location (exact coordinates not published)"
          : "Coordinates: " + p.lat.toFixed(4) + "°N, " + p.lng.toFixed(4) + "°E";
        const popupHtml = '<strong>' + p.label + '</strong> — ' + p.country + '<br>'
          + '<strong>' + p.pa + '</strong>: ' + p.type + '<br>'
          + '<span style="color:#666;font-size:11px;">' + coordNote + '</span><br>'
          + '<a href="' + p.link + '" target="_blank" rel="noopener noreferrer">About this pilot action ↗</a>';
        L.marker([p.lat, p.lng], { icon: icon, keyboard: false, pane: "pilotActionsPane" })
          .bindPopup(popupHtml, { className: "pilot-popup", closeButton: true })
          .addTo(group);
      });
      if (APP.showPilotActions) group.addTo(map);
      return group;
    }

    // Apply APP.showPilotActions to whichever pilot layer groups currently
    // exist. Shared by the toggle button and share-link restoration (which
    // sets the flag directly rather than flipping it).
    function syncPilotActionsVisibility() {
      [[APP.pilotLayerLeft, APP.leftMap], [APP.pilotLayerRight, APP.panels.right && APP.panels.right.map]]
        .forEach(function(pair) {
          const group = pair[0], m = pair[1];
          if (!group || !m) return;
          if (APP.showPilotActions) { if (!m.hasLayer(group)) group.addTo(m); }
          else if (m.hasLayer(group)) { m.removeLayer(group); }
        });
      if (DOM.pilotActionsBtn) DOM.pilotActionsBtn.classList.toggle("is-on", APP.showPilotActions);
      syncPillarStatuses();
    }

    function togglePilotActions() {
      APP.showPilotActions = !APP.showPilotActions;
      syncPilotActionsVisibility();
    }

    // ====================================================================
    // NATURA 2000 (EEA, 2024 release) -- ALL 11,672 SAC/SPA sites of the 7
    // project countries (was: only the 440 within 50km of a Pilot Action).
    // Fetched from the EEA ArcGIS REST API with server-side generalisation
    // (maxAllowableOffset 0.01 deg ~ 1 km, 3-decimal coordinates): ~6.4 MB
    // raw, ~1.3 MB gzipped. Sites that are both SAC and SPA appear once.
    // Lazy-loaded on first toggle-on, not on page load, and rendered with
    // Leaflet's canvas renderer (not SVG).
    // ====================================================================
    const NATURA2000_URL = "webmap_data/natura2000/natura2000_7countries.geojson";

    function buildNatura2000LayerOn(map) {
      return L.geoJSON(APP.natura2000Data, {
        pane: "natura2000Pane",
        interactive: false,
        renderer: L.canvas({ pane: "natura2000Pane" }),
        style: natura2000Style
      });
    }

    function natura2000Style() {
      const k = APP.natura2000Opacity;
      return { color: "#1b5e20", weight: 1, opacity: 0.55 * k, fillColor: "#66bb6a", fillOpacity: 0.12 * k };
    }

    function syncNatura2000Visibility() {
      [["natura2000LayerLeft", APP.leftMap], ["natura2000LayerRight", APP.panels.right && APP.panels.right.map]]
        .forEach(function(pair) {
          const key = pair[0], m = pair[1];
          if (!m) return;
          if (APP.showNatura2000 && APP.natura2000Data) {
            if (!APP[key]) APP[key] = buildNatura2000LayerOn(m);
            if (!m.hasLayer(APP[key])) APP[key].addTo(m);
          } else if (APP[key] && m.hasLayer(APP[key])) {
            m.removeLayer(APP[key]);
          }
        });
      const naturaBox = document.getElementById("naturaBox");
      if (naturaBox) naturaBox.hidden = !APP.showNatura2000;
      syncPillarStatuses();
    }

    async function toggleNatura2000() {
      APP.showNatura2000 = !APP.showNatura2000;
      if (APP.showNatura2000 && !APP.natura2000Data) {
        const nb = document.getElementById("tabNaturaBtn"); if (nb) nb.textContent = "Loading…";
        try {
          const response = await fetch(NATURA2000_URL, { cache: "no-store" });
          if (!response.ok) throw new Error("HTTP " + response.status + " for " + NATURA2000_URL);
          APP.natura2000Data = await response.json();
        } catch (error) {
          console.error(error);
          APP.showNatura2000 = false;
          const nb2 = document.getElementById("tabNaturaBtn"); if (nb2) nb2.textContent = "Failed to load";
          return;
        }
      }
      syncNatura2000Visibility();
    }

    // ====================================================================
    // DSS TABS: the four numbered boxes switch the menu below them. The map keeps whatever is on;
    // the "On the map now" strip lists every layer with a x to remove it, so nothing is hidden silently.
    // ====================================================================
    const TAB_PANELS = { drivers: "overlayCardStrip", ecosystem: "tabEcosystem", restoration: "tabRestoration", pilot: "tabPilot" };
    const TAB_PILLARS = { drivers: "dssPillarRiskDrivers", ecosystem: "dssPillarNatura2000", restoration: "dssPillarRestoration", pilot: "dssPillarPilotActions" };
    const DRIVER_MAIN_LAYERS = ["droughtFuture", "fireWeather", "wind", "bbtforecast"];
    const DRIVER_ALL_LAYERS = DRIVER_MAIN_LAYERS.concat(["droughtRef", "wildfire", "windPressure", "bbtpast", "phenips", "bbtl"]);
    const ON_MAP_NAMES = {
      droughtFuture: "Future drought", droughtRef: "Drought (past)", fireWeather: "Fire weather", wildfire: "Fire pressure",
      wind: "Wind damage", windPressure: "Wind pressure", bbtforecast: "Beetle hazard", bbtl: "Beetle pressure",
      bbtpast: "Past disturbance", phenips: "PHENIPS"
    };

    function setActiveTab(key) {
      if (!TAB_PANELS[key]) key = "drivers";
      APP.activeTab = key;
      Object.keys(TAB_PANELS).forEach(function(k) {
        const panel = document.getElementById(TAB_PANELS[k]), pillar = document.getElementById(TAB_PILLARS[k]);
        if (panel) panel.hidden = (k !== key);
        if (pillar) { pillar.classList.toggle("is-tab-active", k === key); pillar.setAttribute("aria-selected", k === key ? "true" : "false"); }
      });
      if (key === "pilot") buildPilotCards();
      const t = document.querySelector(".layer-ribbon .m-sheet-head span");
      if (t) t.textContent = { drivers: "Disturbance drivers", ecosystem: "Ecosystem services", restoration: "Restoration priorities", pilot: "Pilot actions" }[key];
    }

    function tabLayerCounts() {
      return {
        drivers: DRIVER_ALL_LAYERS.filter(isOverlayVisible).length,
        ecosystem: APP.showNatura2000 ? 1 : 0,
        restoration: APP.showSeed4Forest ? 1 : 0,
        pilot: APP.showPilotActions ? 1 : 0
      };
    }

    function syncPillarStatuses() {
      const counts = tabLayerCounts();
      Object.keys(TAB_PILLARS).forEach(function(k) {
        const pillar = document.getElementById(TAB_PILLARS[k]);
        if (!pillar) return;
        pillar.classList.toggle("has-layers", counts[k] > 0);
        const st = pillar.querySelector(".dss-pillar-status");
        if (!st) return;
        st.textContent = counts[k] ? counts[k] + " on map" : (k === "restoration" ? "Preview" : "None on");
      });
      const set = function(id, on) { const b = document.getElementById(id); if (b) { b.textContent = on ? "Hide" : "Display"; b.classList.toggle("is-on", on); } };
      set("tabNaturaBtn", APP.showNatura2000); set("tabSeedBtn", APP.showSeed4Forest);
      const pins = document.getElementById("tabPilotAll"); if (pins) { pins.textContent = APP.showPilotActions ? "Hide all pins" : "Show all pins"; pins.classList.toggle("is-on", APP.showPilotActions); }
      updateOnMapStrip();
    }

    function updateOnMapStrip() {
      const strip = DOM.dssPillarDetail;
      if (!strip) return;
      const items = [];
      APP.visibleOverlays.slice().reverse().forEach(function(id) {
        items.push({ label: ON_MAP_NAMES[id] || (OVERLAY_CATALOG[id] && OVERLAY_CATALOG[id].label) || id, title: periodTextFor(id), kind: "overlay", id: id });
      });
      if (APP.showNatura2000) items.push({ label: "Natura 2000", kind: "natura" });
      if (APP.showPilotActions) items.push({ label: "Pilot pins", kind: "pins" });
      if (APP.showSeed4Forest) items.push({ label: "Seed4Forest", kind: "seed" });
      strip.innerHTML = '<span class="on-map-label">On the map now:</span>'
        + (items.length ? items.map(function(it, i) {
            return '<span class="on-map-chip"' + (it.title ? ' title="' + it.title + '"' : "") + '>' + it.label
              + '<button type="button" data-onmap-remove="' + i + '" aria-label="Remove ' + it.label + '">&times;</button></span>';
          }).join("") + (items.length > 1 ? '<button type="button" class="on-map-clear" data-onmap-clear>clear all</button>' : "")
          : '<span class="on-map-none">nothing — pick a layer in a tab below</span>');
      strip._items = items;
    }

    function removeOnMapItem(it) {
      if (it.kind === "overlay") setLayerVisibilityFromRibbon(it.id, false);
      else if (it.kind === "natura") toggleNatura2000();
      else if (it.kind === "pins") togglePilotActions();
      else if (it.kind === "seed") setSeedPane(false);
    }

    // Seed4Forest 2.0 (BFW) beside the map. The tool opens at its own default view: it takes no
    // coordinates from a link, so the place is shown here to copy or to click again over there.
    function seedPlaceText() {
      const c = APP.lastClick || APP.myLatLng || (APP.map && APP.map.getCenter());
      return c ? seedCoordText(c) : "—";
    }
    function updateSeedPlace() {
      const el = document.getElementById("seedPlace");
      if (el) el.textContent = seedPlaceText() + (APP.lastClick ? "" : APP.myLatLng ? " (your position)" : " (map centre)");
    }
    function setSeedPane(on) {
      APP.showSeed4Forest = !!on;
      const pane = document.getElementById("seedPane"), frame = document.getElementById("seedFrame"), app = document.getElementById("app");
      if (!pane) return;
      if (on) {
        if (!frame.getAttribute("src")) frame.setAttribute("src", "https://app.seed4forest.org/");
        updateSeedPlace();
        track("seed4forest-open");
      }
      pane.hidden = !on;
      app.classList.toggle("seed-on", !!on);
      syncPillarStatuses();
      refreshMapSize();
      if (APP.map) setTimeout(function() { APP.map.invalidateSize(); }, 50);
    }

    function buildPilotCards() {
      const host = document.getElementById("tabPilotCards");
      if (!host || host.dataset.built) return;
      host.dataset.built = "1";
      const groups = {};
      PILOT_ACTIONS.forEach(function(p) { (groups[p.pa] = groups[p.pa] || []).push(p); });
      host.innerHTML = '<div class="tab-card"><span class="tab-card-eyebrow">All pins</span><span class="tab-card-title">Pilot Actions</span>'
        + '<span class="tab-card-copy">PA-1 to PA-6, ' + PILOT_ACTIONS.length + ' pins.</span>'
        + '<button class="overlay-tool-btn" type="button" id="tabPilotAll">Show all pins</button></div>'
        + Object.keys(groups).sort().map(function(pa) {
            const g = groups[pa], countries = Array.from(new Set(g.map(function(p) { return p.country; }))).join(", ");
            return '<div class="tab-card is-clickable" role="button" tabindex="0" data-pa="' + pa + '"><span class="tab-card-eyebrow">' + pa + '</span>'
              + '<span class="tab-card-title">' + g[0].type.replace(/\s*\(.*$/, "") + '</span>'
              + '<span class="tab-card-copy">' + countries + ' · ' + g.length + (g.length > 1 ? ' pins' : ' pin') + '</span></div>';
          }).join("");
      const focusPa = function(pa) {
        const g = groups[pa]; if (!g) return;
        if (!APP.showPilotActions) togglePilotActions();
        const b = L.latLngBounds(g.map(function(p) { return [p.lat, p.lng]; }));
        if (g.length > 1 && b.getSouthWest().distanceTo(b.getNorthEast()) > 500) APP.map.fitBounds(b.pad(0.4), { maxZoom: 9 });
        else APP.map.setView(b.getCenter(), 9);
      };
      host.addEventListener("click", function(e) {
        const card = e.target.closest("[data-pa]"); if (card) focusPa(card.dataset.pa);
        if (e.target.id === "tabPilotAll") togglePilotActions();
      });
      host.addEventListener("keydown", function(e) {
        const card = e.target.closest && e.target.closest("[data-pa]");
        if (card && (e.key === "Enter" || e.key === " ")) { e.preventDefault(); focusPa(card.dataset.pa); }
      });
      syncPillarStatuses();
    }

    function syncDriverPillar() { syncPillarStatuses(); }

    function initDssTabs() {
      document.querySelectorAll(".dss-pillar[data-tab]").forEach(function(pillar) {
        const pick = function() {
          setActiveTab(pillar.dataset.tab);
          if (isMobileLayout() && APP.openLayersSheet) APP.openLayersSheet();
        };
        pillar.addEventListener("click", pick);
        pillar.addEventListener("keydown", function(e) { if (e.key === "Enter" || e.key === " ") { e.preventDefault(); pick(); } });
      });
      DOM.dssPillarDetail.addEventListener("click", function(e) {
        const rm = e.target.closest("[data-onmap-remove]");
        if (rm) { removeOnMapItem(DOM.dssPillarDetail._items[Number(rm.dataset.onmapRemove)]); return; }
        if (e.target.closest("[data-onmap-clear]")) DOM.dssPillarDetail._items.slice().forEach(removeOnMapItem);
      });
      const tn = document.getElementById("tabNaturaBtn"); if (tn) tn.addEventListener("click", toggleNatura2000);
      const ts = document.getElementById("tabSeedBtn"); if (ts) ts.addEventListener("click", function() { setSeedPane(!APP.showSeed4Forest); });
      const sc = document.getElementById("seedClose"); if (sc) sc.addEventListener("click", function() { setSeedPane(false); });
      const cp = document.getElementById("seedCopy");
      if (cp) cp.addEventListener("click", function() {
        copyText(seedPlaceText(), function() { cp.textContent = "Copied ✓"; });
        setTimeout(function() { cp.textContent = "Copy"; }, 1800);
      });
      APP.map.on("moveend", function() { if (APP.showSeed4Forest && !APP.lastClick) updateSeedPlace(); });
      setActiveTab(APP.activeTab);
      syncPillarStatuses();
    }

    // Bark-beetle event-catalogue markers -- same simple pattern as
    // buildPilotLayer above (independent layer group, not registered in
    // OVERLAY_CATALOG/PANEL_FIELDS, present on both map instances at all
    // times, toggled by a single boolean). Small neutral-grey circles so they
    // never get mistaken for the risk colour ramp used by the raster/hex
    // layers; each is a real, individually-sourced historical event, not a
    // model prediction.
    // Severity palette matches the standalone barkbeetle-catalogue site's own
    // explorer.html legend exactly, so the two tools read consistently.
    const CATALOGUE_SEVERITY_COLOR = {
      minor: "#6b7280", moderate: "#b5891b", high: "#b0452a", extreme: "#8a1f2b"
    };

    function buildCatalogueLayer(map) {
      const group = L.layerGroup();
      CATALOGUE_EVENTS.forEach(function(e) {
        const yearLabel = e.y0 === e.y1 ? String(e.y0) : (e.y0 + '–' + e.y1);
        const sources = e.sources || [];
        const sourcesHtml = sources.length
          ? '<div style="margin-top:4px;font-size:11px;">'
            + sources.map(function(s) {
                const roleTag = '<span style="color:#888;">[' + s.role + ']</span> ';
                return s.url
                  ? roleTag + '<a href="' + s.url + '" target="_blank" rel="noopener noreferrer">' + s.title + ' ↗</a>'
                  : roleTag + s.title + ' <span style="color:#888;">(no link)</span>';
              }).join('<br>')
            + '</div>'
          : '<span style="color:#888;font-size:11px;">no linked source</span>';
        const popupHtml = '<strong>' + e.name + '</strong><br>'
          + e.country + ' · ' + yearLabel + '<br>'
          + '<span style="color:#666;font-size:11px;">severity: ' + e.sev + ' · confidence: ' + e.conf + '</span><br>'
          + '<span style="font-size:12px;">' + e.note + '</span>'
          + sourcesHtml;
        L.circleMarker([e.lat, e.lon], {
          pane: "catalogueEventsPane",
          radius: 5,
          color: "#ffffff",
          weight: 1,
          fillColor: CATALOGUE_SEVERITY_COLOR[e.sev] || "#2a2a2a",
          fillOpacity: 0.85,
          keyboard: false,
          // Path layers bubble click events to the map by default (unlike
          // L.marker, which doesn't) -- without this, clicking a dot also
          // fires the map's own click-inspection popup, which then replaces
          // this popup before the user ever sees it.
          bubblingMouseEvents: false
        })
          .bindPopup(popupHtml, { className: "pilot-popup", closeButton: true })
          // Belt-and-suspenders on top of bubblingMouseEvents:false -- this
          // app's own attachMapClick() listens on the map with map.on("click"),
          // and without an explicit stop here that handler's click-inspection
          // popup wins the race and immediately replaces this one.
          .on("click", function(evt) { L.DomEvent.stopPropagation(evt); })
          .addTo(group);
      });
      if (APP.showCatalogueEvents) group.addTo(map);
      return group;
    }

    // Mirrors syncPilotActionsVisibility's addTo/removeLayer mechanics, but
    // there is no independent toggle button -- catalogue-event markers are
    // an accessory of the "Past disturbance" (bbtpast) layer card, not a
    // standalone layer, so their visibility is derived from bbtpast's own
    // display state rather than stored/toggled separately.
    function syncCatalogueEventsVisibility() {
      [[APP.catalogueLayerLeft, APP.leftMap], [APP.catalogueLayerRight, APP.panels.right && APP.panels.right.map]]
        .forEach(function(pair) {
          const group = pair[0], m = pair[1];
          if (!group || !m) return;
          if (APP.showCatalogueEvents) { if (!m.hasLayer(group)) group.addTo(m); }
          else if (m.hasLayer(group)) { m.removeLayer(group); }
        });
    }

    // Called whenever any overlay is shown/hidden (see showOverlay/
    // hideOverlay) and after share-URL overlay-visibility restoration --
    // recomputes catalogue-marker visibility from bbtpast's current state
    // rather than requiring a dedicated event for it.
    function syncCatalogueEventsFromPastDisturbance() {
      const shouldShow = isOverlayVisible("bbtpast");
      if (shouldShow === APP.showCatalogueEvents) return;
      APP.showCatalogueEvents = shouldShow;
      syncCatalogueEventsVisibility();
    }

    // Closing the click popup also removes the little marker at the clicked
    // point. Deferred one tick: opening a new popup closes the old one first,
    // and that must not take the fresh marker with it.
    function attachMapClick(map, side) {
      map.on("popupclose", function() {
        setTimeout(function() {
          if (map._popup && map.hasLayer(map._popup)) return;   // a new popup took over
          [APP, APP.panels && APP.panels.left, APP.panels && APP.panels.right].forEach(function(st) {
            if (st && st.map === map && st.clickMarker) {
              map.removeLayer(st.clickMarker);
              st.clickMarker = null;
            }
          });
        }, 0);
      });
      map.on("click", function(e) {
        if (APP.split && side !== APP.activeSide) setActiveSide(side);
        updatePanelFromClick(e.latlng);
        openClickPopup(e.latlng, APP, APP.map);
        mirrorClickMarker(e.latlng, side);

        // Unify the info popup across both maps in split screen: the other
        // side keeps its own parked overlay/period settings, so sample and
        // show its value too, not just the side that was actually clicked.
        if (APP.split) {
          const otherSide = side === "left" ? "right" : "left";
          const otherPanel = APP.panels[otherSide];
          if (otherPanel && otherPanel.map) {
            openClickPopup(e.latlng, otherPanel, otherPanel.map);
          }
        }
      });
    }

    function newMapBasemapLayer(map) {
      const def = BASEMAPS[APP.currentBasemap] || BASEMAPS.osm;
      const layer = def.build();
      if (layer) { layer.addTo(map); if (layer.bringToBack) layer.bringToBack(); }
      return layer;
    }

    // Reflect the active config in the control panel (no map re-render).
    function syncControlsToActive() {
      if (DOM.rcpSelector) DOM.rcpSelector.value = APP.currentRcp;
      syncPhenipsControls();
      if (DOM.windPeriodSelector) DOM.windPeriodSelector.value = APP.currentWindPeriod;
      updateWindGcmAvailability();
      if (DOM.droughtFuturePeriodSelector) DOM.droughtFuturePeriodSelector.value = APP.currentDroughtFuturePeriod;
      if (DOM.droughtRefMetricSelector) DOM.droughtRefMetricSelector.value = APP.currentDroughtRefMetric;
      if (DOM.fireViewSelector) DOM.fireViewSelector.value = APP.currentFireView;
      if (DOM.fireScenSelector) DOM.fireScenSelector.value = APP.currentFireScen;
      if (DOM.riskSlider) DOM.riskSlider.value = Math.round(APP.riskThreshold * 100);
      if (DOM.riskValue) DOM.riskValue.textContent = Math.round(APP.riskThreshold * 100) + "%";
      if (APP.focusedOverlay === "bbtforecast") {
        hydrateForecastControls().catch(function(e) { console.error(e); });
      }
      updateOverlayCardStates();
      syncOverlayControls(false);
    }

    function updateSplitHighlight() {
      const a = document.getElementById("map");
      const b = document.getElementById("map2");
      if (a) a.classList.toggle("compare-active", APP.split && APP.activeSide === "left");
      if (b) b.classList.toggle("compare-active", APP.split && APP.activeSide === "right");
      if (DOM.mapCaptionLeft) DOM.mapCaptionLeft.textContent = (APP.activeSide === "left" ? "● Active — Left" : "Left");
      if (DOM.mapCaptionRight) DOM.mapCaptionRight.textContent = (APP.activeSide === "right" ? "● Active — Right" : "Right");
      if (DOM.splitBtn) {
        DOM.splitBtn.textContent = APP.split ? "Single map" : "Split screen";
        DOM.splitBtn.classList.toggle("is-on", APP.split);
      }
    }

    function setActiveSide(side) {
      if (!APP.split || side === APP.activeSide) return;
      savePanel(APP.activeSide);
      APP.activeSide = side;
      loadPanel(side);
      syncControlsToActive();
      updateSplitHighlight();
      setStatus("Editing the " + side + " map.");
    }

    async function activateSplit() {
      if (APP.split) return;
      APP.split = true;
      document.getElementById("app").classList.add("compare-on");
      APP.map.invalidateSize();

      const mapR = L.map("map2", {
        center: APP.map.getCenter(), zoom: APP.map.getZoom(),
        renderer: L.svg(), zoomControl: false, attributionControl: false,
        maxBounds: MAP_SOFT_BOUNDS, maxBoundsViscosity: MAP_BOUNDS_VISCOSITY
      });
      createMapPanes(mapR);
      const rb = buildCountryLayersOn(mapR);
      const rBasemap = newMapBasemapLayer(mapR);
      attachMapClick(mapR, "right");
      APP.pilotLayerRight = buildPilotLayer(mapR);
      APP.catalogueLayerRight = buildCatalogueLayer(mapR);
      document.getElementById("map2").style.background = (BASEMAPS[APP.currentBasemap] || BASEMAPS.osm).bg;

      // The right panel starts as a clone of the left config with its own map.
      APP.panels.right = {
        map: mapR, basemapLayer: rBasemap, clickMarker: null,
        countryHaloLayer: rb.countryHaloLayer, countryLayer: rb.countryLayer,
        selectedHexLayer: null, selectedCountryLayer: null,
        currentYear: APP.currentYear, currentRcp: APP.currentRcp, riskThreshold: APP.riskThreshold,
        visibleOverlays: APP.visibleOverlays.slice(), focusedOverlay: APP.focusedOverlay,
        layerOpacity: Object.assign({}, APP.layerOpacity),
        currentPhenipsScen: APP.currentPhenipsScen, currentPhenisPeriod: APP.currentPhenisPeriod,
        currentPhenipsVar: APP.currentPhenipsVar, currentPhenipsView: APP.currentPhenipsView, currentPhenipsExp: APP.currentPhenipsExp,
        currentForecastModel: APP.currentForecastModel, currentForecastPeriod: APP.currentForecastPeriod,
        currentForecastView: APP.currentForecastView, currentForecastThreshold: APP.currentForecastThreshold,
        currentWindPeriod: APP.currentWindPeriod, currentWindGcm: APP.currentWindGcm, currentForecastRasterRange: null,
        currentDroughtFuturePeriod: APP.currentDroughtFuturePeriod, currentDroughtRefMetric: APP.currentDroughtRefMetric,
        currentFireView: APP.currentFireView, currentFireScen: APP.currentFireScen,
        overlayLoadTokens: { bbtl: 0, bbtforecast: 0, phenips: 0, bbtpast: 0, wind: 0, droughtFuture: 0, droughtRef: 0, fireWeather: 0, wildfire: 0, windPressure: 0 },
        overlayVisible: APP.overlayVisible, hexGeojson: { bbtl: null, wildfire: null, windPressure: null },
        layerInstances: { bbtl: null, bbtforecast: null, phenips: null, bbtpast: null, wind: null, droughtFuture: null, droughtRef: null, fireWeather: null, wildfire: null, windPressure: null }
      };

      if (APP.map.sync) { APP.map.sync(mapR); mapR.sync(APP.map); }
      syncNatura2000Visibility();

      // Render the right map: become right, render its overlays, return to left.
      savePanel("left");
      APP.activeSide = "right";
      loadPanel("right");
      mapR.invalidateSize();
      try {
        await Promise.all(APP.visibleOverlays.map(function(id) { return updateOverlayLayer(id); }));
      } catch (e) { console.error(e); }
      savePanel("right");
      APP.activeSide = "left";
      loadPanel("left");

      updateSplitHighlight();
      setStatus("Split on — click a map to make it active; the panel edits the active map.");
    }

    function deactivateSplit() {
      if (!APP.split) return;
      if (APP.activeSide === "right") { savePanel("right"); APP.activeSide = "left"; loadPanel("left"); }
      const right = APP.panels.right;
      if (right && right.map) {
        if (APP.map.unsync) { try { APP.map.unsync(right.map); right.map.unsync(APP.map); } catch (e) {} }
        right.map.remove();
      }
      APP.panels.right = null;
      APP.pilotLayerRight = null;
      APP.catalogueLayerRight = null;
      APP.natura2000LayerRight = null;
      APP.split = false;
      document.getElementById("app").classList.remove("compare-on");
      updateSplitHighlight();
      syncControlsToActive();
      setTimeout(function() { if (APP.map) APP.map.invalidateSize(); }, 0);
    }

    function toggleSplit() {
      if (APP.split) deactivateSplit(); else activateSplit();
    }


    // ====================================================================
    // INIT, SHAREABLE LINK & CONTROL WIRING
    // ====================================================================
    function initMap() {
      const defaultOverlay = getVisibleOverlay();
      if (defaultOverlay.meta && Array.isArray(defaultOverlay.meta.years) && defaultOverlay.meta.years.length) {
        APP.currentYear = defaultOverlay.meta.years[0];
      }

      APP.map = L.map("map", {
        center: [50.2, 12.5],
        zoom: 4.7,
        renderer: L.svg(),
        maxBounds: MAP_SOFT_BOUNDS,
        maxBoundsViscosity: MAP_BOUNDS_VISCOSITY
      });
      APP.leftMap = APP.map;
      createMapPanes(APP.map);
      APP.pilotLayerLeft = buildPilotLayer(APP.leftMap);
      APP.catalogueLayerLeft = buildCatalogueLayer(APP.leftMap);

      addBasemapControl();
      addSplitScreenControl();
      setBasemap(APP.currentBasemap);

      // Place search (Nominatim) — type a town/region and fly there.
      if (L.Control.geocoder) {
        L.Control.geocoder({
          defaultMarkGeocode: true,
          placeholder: "Search place…",
          position: "topleft"
        }).addTo(APP.map);
      }

      addLocateControl();
      attachMapClick(APP.map, "left");
      setTimeout(autoLocateOnce, 600);   // after the map has settled

      setLegendLabels();
      initControls();
      loadCountryBorders();
      applyStateFromUrl();
      // Without a shared link, draw the default layer(s) explicitly.
      if (!location.hash) APP.visibleOverlays.forEach(function(id) { updateOverlayLayer(id); });
      addMenuToggleTab();
    }

    // --- Shareable view state (URL hash) -----------------------------------
    // Encode the meaningful view so a link reproduces what the user sees.
    // Fixed key order for the compact "lo" (layer opacity) param -- must match
    // the order restored in applyStateFromUrl.
    const LAYER_OPACITY_KEYS = ["bbtl", "bbtforecast", "phenips", "bbtpast", "wind", "droughtFuture", "droughtRef", "wildfire", "windPressure", "fireWeather"];

    // The plain-text settings of a shared link: [URL key, APP field]. Used by BOTH
    // buildShareUrl and applyStateFromUrl, so a setting can never be saved without being restored.
    const SHARE_FIELDS = [
      ["fm", "currentForecastModel"], ["fp", "currentForecastPeriod"], ["fv", "currentForecastView"],
      ["ft", "currentForecastThreshold"], ["ps", "currentPhenipsScen"], ["pp", "currentPhenisPeriod"],
      ["pv", "currentPhenipsVar"], ["pw", "currentPhenipsView"], ["pe", "currentPhenipsExp"], ["tb", "activeTab"], ["wp", "currentWindPeriod"],
      ["wg", "currentWindGcm"], ["dp", "currentDroughtFuturePeriod"], ["dm", "currentDroughtRefMetric"],
      ["fw", "currentFireView"], ["fs", "currentFireScen"], ["rc", "currentRcp"]
    ];

    function buildShareUrl() {
      const p = new URLSearchParams();
      p.set("bm", APP.currentBasemap);
      p.set("ov", APP.visibleOverlays.join(","));
      p.set("fo", APP.focusedOverlay);
      SHARE_FIELDS.forEach(function(f) { p.set(f[0], APP[f[1]]); });
      if (APP.currentYear != null) p.set("yr", String(APP.currentYear));
      p.set("th", String(Math.round(APP.riskThreshold * 100)));
      // Map view: center + zoom, so a shared link opens on the same place.
      if (APP.map) {
        const c = APP.map.getCenter();
        p.set("lat", c.lat.toFixed(5));
        p.set("lng", c.lng.toFixed(5));
        p.set("z", APP.map.getZoom().toFixed(2));
      }
      // Per-layer opacity sliders.
      p.set("lo", LAYER_OPACITY_KEYS.map(function(k) { return APP.layerOpacity[k]; }).join(","));
      p.set("sp", APP.split ? "1" : "0");
      p.set("sv", APP.simpleView ? "1" : "0");
      p.set("pa", APP.showPilotActions ? "1" : "0");
      p.set("n2", APP.showNatura2000 ? "1" : "0");
      const folded = Object.keys(FOLDS).filter(function(c) { return document.body.classList.contains(c); });
      if (folded.length) p.set("fd", folded.join(","));
      // No "ce" param -- catalogue-event visibility is derived from "ov"
      // (whether bbtpast is in it), not stored independently.
      return location.origin + location.pathname + "#" + p.toString();
    }

    function applyStateFromUrl() {
      const hash = location.hash.replace(/^#/, "");
      if (!hash) return;
      const p = new URLSearchParams(hash);
      if (!p.has("ov") && !p.has("bm")) return;

      if (p.has("bm")) setBasemap(p.get("bm"));

      SHARE_FIELDS.forEach(function(f) { if (p.has(f[0])) APP[f[1]] = p.get(f[0]); });
      if (p.has("yr")) APP.currentYear = Number(p.get("yr"));
      if (p.has("th")) APP.riskThreshold = Math.max(0, Math.min(100, Number(p.get("th")) || 0)) / 100;
      if (p.has("sv")) APP.simpleView = p.get("sv") !== "0";
      if (p.has("fd")) p.get("fd").split(",").forEach(function(c) { if (FOLDS[c]) setFold(c, true); });
      if (p.get("n2") === "1" && !APP.showNatura2000) toggleNatura2000();
      if (p.has("pa")) { APP.showPilotActions = p.get("pa") === "1"; syncPilotActionsVisibility(); }

      if (p.has("lo")) {
        const vals = p.get("lo").split(",");
        LAYER_OPACITY_KEYS.forEach(function(key, i) {
          const v = Number(vals[i]);
          if (isFinite(v)) APP.layerOpacity[key] = Math.max(0, Math.min(1, v));
        });
      }
      if (p.has("lat") && p.has("lng") && p.has("z")) {
        const lat = Number(p.get("lat")), lng = Number(p.get("lng")), z = Number(p.get("z"));
        if (isFinite(lat) && isFinite(lng) && isFinite(z)) APP.map.setView([lat, lng], z);
      }

      const setSel = function(dom, val) { if (dom && val != null) dom.value = val; };
      syncPhenipsControls();
      setSel(DOM.windPeriodSelector, APP.currentWindPeriod);
      updateWindGcmAvailability();
      setSel(DOM.droughtFuturePeriodSelector, APP.currentDroughtFuturePeriod);
      setSel(DOM.droughtRefMetricSelector, APP.currentDroughtRefMetric);
      setSel(DOM.fireViewSelector, APP.currentFireView);
      setSel(DOM.fireScenSelector, APP.currentFireScen);
      setSel(DOM.rcpSelector, APP.currentRcp);
      if (DOM.riskSlider) DOM.riskSlider.value = Math.round(APP.riskThreshold * 100);
      if (DOM.riskValue) DOM.riskValue.textContent = Math.round(APP.riskThreshold * 100) + "%";

      if (p.has("ov")) {
        const asked = p.get("ov").split(",").filter(Boolean);
        APP.visibleOverlays = asked.filter(function(id) {
          return OVERLAY_CATALOG[id] && OVERLAY_CATALOG[id].available;
        });
        // A link naming only layers that no longer exist opens the default (empty) map.
        if (asked.length && !APP.visibleOverlays.length) APP.visibleOverlays = [];
        syncOverlayVisibilityFlag();
        syncCatalogueEventsFromPastDisturbance();
      }
      if (p.has("fo") && OVERLAY_CATALOG[p.get("fo")]) APP.focusedOverlay = p.get("fo");

      getOverlayIds().forEach(function(id) { clearLayerInstance(id); });
      updateOverlayPaneOrder();
      syncOverlayControls(false);
      applySimpleView();
      APP.visibleOverlays.forEach(function(id) { updateOverlayLayer(id); });
      updateOverlayCardStates();
      setLegendLabels();
      syncCitationVisibility();

      // Restore split screen last, so the right map clones the already-
      // restored view/overlay/opacity state from the left one.
      if (p.get("sp") === "1" && !APP.split) activateSplit();
    }

    function bindShareButton() {
      if (!DOM.shareBtn) return;
      DOM.shareBtn.addEventListener("click", function() {
        const url = buildShareUrl();
        history.replaceState(null, "", url);
        const done = function() {
          DOM.shareBtn.textContent = "Link copied!";
          setTimeout(function() { DOM.shareBtn.textContent = "Share view"; }, 1600);
        };
        if (navigator.clipboard && navigator.clipboard.writeText) {
          navigator.clipboard.writeText(url).then(done, function() { window.prompt("Copy this link:", url); });
        } else {
          window.prompt("Copy this link:", url);
        }
      });
    }

    function stopCardControlClick(event) {
      event.stopPropagation();
    }

    function bindForecastSelector(select, stateKey) {
      select.addEventListener("change", function() {
        APP[stateKey] = select.value;
        refreshForecastOverlayIfVisible();
      });
    }

    // Show/hide the "advanced" controls (climate-model/GCM pickers,
    // validation detail) per APP.simpleView. Every time period, including
    // end-of-century, stays available in both modes -- only pinned to a
    // model/GCM ensemble is hidden. Reversible and idempotent -- safe to
    // call any time, not just when the toggle itself is clicked.
    function applySimpleView() {
      const simple = APP.simpleView;
      if (DOM.detailModeBtn) {
        DOM.detailModeBtn.textContent = simple ? "Show full detail" : "Simple view";
        DOM.detailModeBtn.setAttribute("aria-pressed", String(!simple));
      }

      // Bark beetle forecast: hide the Model picker + validation metrics,
      // pinned to the one model with a real held-out evaluation and a
      // construction-guaranteed (not just empirically observed) non-
      // decreasing future response -- see CHELSA_NATIVE_TRAINING project docs.
      if (DOM.forecastModelField) DOM.forecastModelField.hidden = simple;
      if (DOM.forecastMetricsBox) DOM.forecastMetricsBox.hidden = simple;
      if (simple && APP.currentForecastModel !== SIMPLE_VIEW_DEFAULT_FORECAST_MODEL) {
        APP.currentForecastModel = SIMPLE_VIEW_DEFAULT_FORECAST_MODEL;
        refreshForecastOverlayIfVisible();
      }

      // Wind: hide only the GCM picker. Every period option (including far
      // future) stays available in both modes -- not something to hide.
      if (DOM.windGcmField) DOM.windGcmField.hidden = simple;
      if (simple && APP.currentWindGcm !== "mpi") {
        APP.currentWindGcm = "mpi";
        if (DOM.windGcmSelector) DOM.windGcmSelector.value = "mpi";
        if (isOverlayVisible("wind")) updateWindLayer();
      }
    }

    function initControls() {
      const defaultOverlay = getFocusedOverlay();
      if (defaultOverlay.meta && Array.isArray(defaultOverlay.meta.years) && defaultOverlay.meta.years.length) {
        DOM.yearSlider.max = defaultOverlay.meta.years.length - 1;
        DOM.yearSlider.value = 0;
        DOM.yearValue.textContent = String(APP.currentYear);
      } else {
        DOM.yearValue.textContent = "NA";
      }

      DOM.overlayCardStrip.querySelectorAll("[data-overlay]").forEach(function(card) {
        card.addEventListener("click", function() {
          if (card.classList.contains("coming-soon")) {
            return;
          }
          APP.focusedOverlay = card.dataset.overlay;
          syncOverlayControls(false);
        });
      });

      DOM.overlayCardStrip.querySelectorAll("[data-toggle-overlay]").forEach(function(toggleButton) {
        toggleButton.addEventListener("pointerdown", stopCardControlClick);
        toggleButton.addEventListener("click", function(event) {
          stopCardControlClick(event);
          const overlayId = toggleButton.dataset.toggleOverlay;
          const isOn = toggleButton.classList.contains("is-on");
          setLayerVisibilityFromRibbon(overlayId, !isOn);
        });
      });

      DOM.overlayCardStrip.querySelectorAll("[data-opacity-overlay]").forEach(function(slider) {
        slider.addEventListener("pointerdown", stopCardControlClick);
        slider.addEventListener("click", stopCardControlClick);
        slider.addEventListener("input", function(event) {
          stopCardControlClick(event);
          const overlayId = slider.dataset.opacityOverlay;
          APP.layerOpacity[overlayId] = Number(slider.value);
          updateOverlayCardStates();
          applyOpacityToLayers();
        });
      });

      bindForecastSelector(DOM.forecastModelSelector, "currentForecastModel");
      bindForecastSelector(DOM.forecastPeriodSelector, "currentForecastPeriod");
      bindForecastSelector(DOM.forecastViewSelector, "currentForecastView");
      bindForecastSelector(DOM.forecastThresholdSelector, "currentForecastThreshold");

      let yearTimer = null;
      DOM.phenipsYearSlider.addEventListener("input", function() {
        APP.phenipsLastYear = Number(DOM.phenipsYearSlider.value);
        APP.currentPhenisPeriod = "y" + APP.phenipsLastYear;
        DOM.phenipsYearVal.textContent = String(APP.phenipsLastYear);
        DOM.phenisPeriodSelector.value = APP.currentPhenisPeriod;
        clearTimeout(yearTimer);
        yearTimer = setTimeout(function() { if (isOverlayVisible("phenips")) updatePhenipsLayer(); else setLegendLabels(); }, 150);
      });
      [
        [DOM.phenipsScenSelector, "currentPhenipsScen"], [DOM.phenisPeriodSelector, "currentPhenisPeriod"],
        [DOM.phenipsVarSelector, "currentPhenipsVar"], [DOM.phenipsViewSelector, "currentPhenipsView"],
        [DOM.phenipsExpSelector, "currentPhenipsExp"]
      ].forEach(function(pair) {
        pair[0].addEventListener("change", function() {
          APP[pair[1]] = pair[0].value;
          syncPhenipsControls();
          if (isOverlayVisible("phenips")) updatePhenipsLayer();
          else setLegendLabels();
        });
      });

      DOM.windPeriodSelector.addEventListener("change", function() {
        APP.currentWindPeriod = DOM.windPeriodSelector.value;
        updateWindGcmAvailability();
        if (isOverlayVisible("wind")) updateWindLayer();
        else setLegendLabels();
      });

      DOM.windGcmSelector.addEventListener("change", function() {
        APP.currentWindGcm = DOM.windGcmSelector.value;
        if (isOverlayVisible("wind")) updateWindLayer();
        else setLegendLabels();
      });

      DOM.droughtFuturePeriodSelector.addEventListener("change", function() {
        APP.currentDroughtFuturePeriod = DOM.droughtFuturePeriodSelector.value;
        if (isOverlayVisible("droughtFuture")) updateDroughtFutureLayer();
        else setLegendLabels();
      });

      DOM.droughtRefMetricSelector.addEventListener("change", function() {
        APP.currentDroughtRefMetric = DOM.droughtRefMetricSelector.value;
        renderDroughtRefInfo();
        if (isOverlayVisible("droughtRef")) updateDroughtRefLayer();
        else setLegendLabels();
      });

      DOM.fireViewSelector.addEventListener("change", function() {
        APP.currentFireView = DOM.fireViewSelector.value;
        renderFireInfo();
        if (isOverlayVisible("fireWeather")) updateFireWeatherLayer();
        else setLegendLabels();
      });
      DOM.fireScenSelector.addEventListener("change", function() {
        APP.currentFireScen = DOM.fireScenSelector.value;
        if (isOverlayVisible("fireWeather")) updateFireWeatherLayer();
        else setLegendLabels();
      });

      // Live while dragging (at most one redraw per ~100 ms), final redraw on release.
      let riskTimer = null;
      const applyRiskSlider = function() {
        riskTimer = null;
        APP.riskThreshold = Number(DOM.riskSlider.value) / 100;
        DOM.riskValue.textContent = DOM.riskSlider.value + "%";
        applyRiskThreshold();
      };
      DOM.riskSlider.addEventListener("input", function() {
        DOM.riskValue.textContent = DOM.riskSlider.value + "%";
        if (!riskTimer) riskTimer = setTimeout(applyRiskSlider, 100);
      });
      DOM.riskSlider.addEventListener("change", function() {
        clearTimeout(riskTimer);
        applyRiskSlider();
        setLegendLabels();
      });

      bindShareButton();
      if (DOM.detailModeBtn) DOM.detailModeBtn.addEventListener("click", function() {
        APP.simpleView = !APP.simpleView;
        applySimpleView();
      });
      applySimpleView();
      if (DOM.splitBtn) DOM.splitBtn.addEventListener("click", toggleSplit);
      if (DOM.pilotActionsBtn) DOM.pilotActionsBtn.addEventListener("click", togglePilotActions);
      initMobileLayout();
      initDssTabs();
      initSiteSettings();

      // Any driver box's accessory layers stay collapsed until its header
      // strip is clicked (currently Bark beetle +3, Drought +1).
      document.querySelectorAll(".bb-toggle[data-toggle-target]").forEach(function(toggle) {
        const accessories = document.getElementById(toggle.dataset.toggleTarget);
        if (!accessories) return;
        const toggleGroup = function() {
          accessories.hidden = !accessories.hidden;
          toggle.classList.toggle("open", !accessories.hidden);
          toggle.setAttribute("aria-expanded", String(!accessories.hidden));
        };
        toggle.addEventListener("click", toggleGroup);
        toggle.addEventListener("keydown", function(e) {
          if (e.key === "Enter" || e.key === " ") { e.preventDefault(); toggleGroup(); }
        });
      });

      // Period-strip chips with real coverage are click-to-jump links;
      // empty chips (no data-jump-overlay) stay plain and non-interactive.
      document.querySelectorAll(".period-chip[data-jump-overlay]").forEach(function(chip) {
        chip.setAttribute("role", "button");
        chip.setAttribute("tabindex", "0");
        const activate = function() {
          jumpToOverlayPeriod(chip.dataset.jumpOverlay, chip.dataset.jumpPeriod).catch(function(error) {
            console.error(error);
          });
        };
        chip.addEventListener("click", activate);
        chip.addEventListener("keydown", function(e) {
          if (e.key === "Enter" || e.key === " ") { e.preventDefault(); activate(); }
        });
      });

      DOM.rcpSelector.value = APP.currentRcp;
      DOM.rcpSelector.addEventListener("change", function() {
        APP.currentRcp = DOM.rcpSelector.value;
        HEX_OVERLAY_IDS.forEach(function(id) {
          if (isOverlayVisible(id)) updateHexLayer(id);
        });
      });

      DOM.yearSlider.addEventListener("input", function() {
        const showHexControls = HEX_OVERLAY_IDS.includes(APP.focusedOverlay);
        if (!showHexControls) {
          return;
        }
        const yearIndex = Number(DOM.yearSlider.value);
        APP.currentYear = getFocusedOverlay().meta.years[yearIndex];
        DOM.yearValue.textContent = String(APP.currentYear);
        HEX_OVERLAY_IDS.forEach(function(id) {
          if (isOverlayVisible(id)) updateHexLayer(id);
        });
      });

      DOM.toggleOverlayBtn.addEventListener("click", function() {
        const focusedVisible = isOverlayVisible(APP.focusedOverlay);
        setLayerVisibilityFromRibbon(APP.focusedOverlay, !focusedVisible);
      });

      DOM.fitOverlayBtn.addEventListener("click", function() {
        fitDataBounds();
      });

      syncOverlayControls("bbtforecast");
    }

    async function loadCountryBorders() {
      try {
        const response = await fetch(COUNTRY_GEOJSON_URL);
        if (!response.ok) {
          throw new Error("HTTP " + response.status);
        }

        APP.countryGeojson = getMapExtentCountryGeojson(await response.json());

        const layers = buildCountryLayersOn(APP.map);
        APP.countryHaloLayer = layers.countryHaloLayer;
        APP.countryLayer = layers.countryLayer;
      } catch (error) {
        console.error(error);
        setStatus("Country borders could not be loaded");
      }
    }

    async function updateHexLayer(overlayId) {
      if (!isOverlayVisible(overlayId)) return;

      const requestToken = beginOverlayLoad(overlayId);

      const fileName = getCurrentFile(overlayId);

      if (!fileName) {
        setStatus("No file configured for this year");
        return;
      }

      APP.hexGeojson[overlayId] = null;
      resetSelectedHex();
      resetSelectedCountry();

      try {
        setStatus("Loading " + fileName + "...");

        const response = await fetch(fileName);
        if (!response.ok) {
          throw new Error("HTTP " + response.status + " for " + fileName);
        }

        const nextGeojson = await response.json();
        if (!isCurrentOverlayLoad(overlayId, requestToken)) {
          return;
        }

        APP.hexGeojson[overlayId] = nextGeojson;
        APP.layerInstances[overlayId] = L.geoJSON(nextGeojson, {
          pane: getOverlayPaneName(overlayId),
          style: baseHexStyle(overlayId),
          onEachFeature: function(feature, layer) {
            const value = Number(feature.properties.mean);

            layer.bindTooltip(
              isFinite(value) ? ("mean: " + value.toFixed(6)) : "mean: NA",
              {
                sticky: true,
                direction: "top",
                className: "hex-tooltip"
              }
            );

            layer.on("mouseover", function() {
              if (APP.selectedHexLayer !== layer) {
                layer.setStyle({
                  color: "#111111",
                  weight: 1.4,
                  opacity: 1,
                  fillOpacity: Math.min(getLayerOpacity(overlayId) + 0.10, 1)
                });
              }
            });

            layer.on("mouseout", function() {
              if (APP.selectedHexLayer !== layer && APP.layerInstances[overlayId]) {
                APP.layerInstances[overlayId].resetStyle(layer);
              }
            });
          }
        }).addTo(APP.map);

        if (!isCurrentOverlayLoad(overlayId, requestToken)) {
          return;
        }
        setStatus("");

      } catch (error) {
        console.error(error);
        setStatus("Failed loading " + fileName + " :: " + error.message);
      }
    }

    initMap();
