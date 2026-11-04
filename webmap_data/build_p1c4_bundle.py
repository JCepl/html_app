#!/usr/bin/env python3
"""Final bark-beetle hazard model (P1C4, SSP5-8.5) -> forecast_bundle/p1c4_ssp585.

Reads the delivered rasters (EPSG:4326, Float32, nodata -9999), reprojects each to Web Mercator
(so the Leaflet imageOverlay lines up with the basemap), colours it and writes, per period:
  probability_ensemble_mean.png/_meta.json   hazard mean, sequential ramp, 0-0.8
  risk_visible_t50.png/_meta.json            the same, only cells >= 0.5 (the flag threshold)
  delta_from_near.png/_meta.json             change vs 1981-2010, diverging, -0.3..+0.3 (future periods)
  share05.png/_meta.json                     share of years with hazard >= 0.5, 0-1   (if delivered)
  extrapolation.png/_meta.json               share of predictors outside the training range, greys (if delivered)
and registers the model in forecast_bundle/manifest.json. The .tif names in the manifest are only the stems the app
swaps for .png/_meta.json; no GeoTIFF is copied into the bundle.

Usage: python3 webmap_data/build_p1c4_bundle.py <folder with the rasters> [--max-width 3200]
Files are found by name: 'hazard_mean' + period, 'hazard_change' + period, 'share05' + period, 'extrapolation' + period."""
import argparse, glob, json, os, sys
import numpy as np
from osgeo import gdal, osr
from PIL import Image

gdal.UseExceptions()
HERE = os.path.dirname(os.path.abspath(__file__))
BUNDLE = os.path.join(HERE, "forecast_bundle")
MODEL = "p1c4_ssp585"
PERIODS = [("1981_2010", "obs_1981-2010"), ("2011_2040", "2011-2040"), ("2041_2070", "2041-2070"), ("2071_2100", "2071-2100")]
CAP = 0.8           # top of the hazard scale
CHANGE_CAP = 0.3    # +/- of the change scale
# Sequential light -> dark (ColorBrewer YlOrRd, 9 classes); must equal P1C4_RAMP in app.js
SEQ = ["#ffffcc", "#ffeda0", "#fed976", "#feb24c", "#fd8d3c", "#fc4e2a", "#e31a1c", "#bd0026", "#800026"]
# Diverging blue (fall) - white - red (rise); same stops as DELTA_GRADIENT in app.js
DIV = [(0.00, "#053061"), (0.18, "#2166ac"), (0.32, "#4393c3"), (0.42, "#92c5de"), (0.48, "#d1e5f0"), (0.50, "#f7f7f7"),
       (0.52, "#fddbc7"), (0.58, "#f4a582"), (0.68, "#d6604d"), (0.82, "#b2182b"), (1.00, "#67001f")]
GREY = [(0.0, "#f0f0f0"), (1.0, "#252525")]
log = lambda *a: print(*a, flush=True)


def hexrgb(h): return [int(h[i:i + 2], 16) for i in (1, 3, 5)]


def ramp(stops, t):
    """Linear sRGB interpolation of stops [(pos, hex)] at t in [0,1] -> (n,3) uint8."""
    pos = np.array([s[0] for s in stops]); cols = np.array([hexrgb(s[1]) for s in stops], float)
    return np.stack([np.interp(t, pos, cols[:, k]) for k in range(3)], -1).round().astype(np.uint8)


SEQ_STOPS = [(i / (len(SEQ) - 1), c) for i, c in enumerate(SEQ)]


def find(src, kind, tag):
    hits = sorted(glob.glob(os.path.join(src, "**", f"*{kind}*{tag}*.tif"), recursive=True))
    return hits[0] if hits else None


def spruce_domain(ds, mask_path, min_share=0.25):
    """Boolean array on the raster's own grid: cells where >= min_share of the 1 km spruce >= 0.20 mask is 1.
    This is the model's domain ('cells with spruce >= 0.20', as in 14_projection_summary.py); the maps are shown only there."""
    gt = ds.GetGeoTransform(); w, h = ds.RasterXSize, ds.RasterYSize
    m = gdal.Warp("", gdal.Open(mask_path), format="MEM", dstSRS=ds.GetProjection() or "EPSG:4326", resampleAlg="average",
                  outputBounds=(gt[0], gt[3] + gt[5] * h, gt[0] + gt[1] * w, gt[3]), width=w, height=h, outputType=gdal.GDT_Float32)
    return m.GetRasterBand(1).ReadAsArray() >= min_share


def prepare(path, mask_path, fill_px=3):
    """Source raster as an in-memory dataset with NaN outside the spruce domain; gaps inside it (a few cells of the
    approximate build where the colour bar cannot be inverted) are filled from the nearest valid cell within fill_px."""
    from scipy import ndimage as ndi
    ds = gdal.Open(path); b = ds.GetRasterBand(1); nd = b.GetNoDataValue()
    a = b.ReadAsArray().astype("f4")
    if nd is not None: a[a == nd] = np.nan
    if mask_path:
        dom = spruce_domain(ds, mask_path)
        gap = dom & ~np.isfinite(a)
        if gap.any() and np.isfinite(a).any():
            d, (ri, ci) = ndi.distance_transform_edt(~np.isfinite(a), return_indices=True)
            fill = gap & (d <= fill_px); a[fill] = a[ri[fill], ci[fill]]
        a[~dom] = np.nan
    out = gdal.GetDriverByName("MEM").Create("", ds.RasterXSize, ds.RasterYSize, 1, gdal.GDT_Float32)
    out.SetGeoTransform(ds.GetGeoTransform()); out.SetProjection(ds.GetProjection() or "EPSG:4326")
    a[~np.isfinite(a)] = -9999; out.GetRasterBand(1).WriteArray(a); out.GetRasterBand(1).SetNoDataValue(-9999)
    return out


def to_mercator(path, max_width, mask_path=None):
    """Read, mask to the spruce domain, reproject to EPSG:3857 (bilinear, nodata kept); return array (NaN = no data) and lat/lon bounds."""
    ds = prepare(path, mask_path); b = ds.GetRasterBand(1); nd = b.GetNoDataValue()
    gt = ds.GetGeoTransform(); w = ds.RasterXSize
    merc_w = 2 * 20037508.34 * (abs(gt[1]) * w) / 360.0
    res = merc_w / min(max_width, int(w * 1.5))
    out = gdal.Warp("", ds, format="MEM", dstSRS="EPSG:3857", resampleAlg="bilinear", xRes=res, yRes=res,
                    srcNodata=nd, dstNodata=-9999, outputType=gdal.GDT_Float32)
    a = out.GetRasterBand(1).ReadAsArray().astype("f4"); a[a == -9999] = np.nan
    g = out.GetGeoTransform(); x0, y1 = g[0], g[3]; x1, y0 = x0 + g[1] * out.RasterXSize, y1 + g[5] * out.RasterYSize
    s3857 = osr.SpatialReference(); s3857.ImportFromEPSG(3857); s4326 = osr.SpatialReference(); s4326.ImportFromEPSG(4326)
    s4326.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER); s3857.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    ct = osr.CoordinateTransformation(s3857, s4326)
    west, south, _ = ct.TransformPoint(x0, y0); east, north, _ = ct.TransformPoint(x1, y1)
    return a, {"south": south, "west": west, "north": north, "east": east}


def write(outdir, stem, rgba, bounds, source, vmin, vmax):
    os.makedirs(outdir, exist_ok=True)
    # 256-colour palette with alpha: about a quarter of the size of RGBA, and the colours stay on the ramp
    Image.fromarray(rgba, "RGBA").quantize(colors=256, method=Image.FASTOCTREE, dither=Image.NONE).save(
        os.path.join(outdir, stem + ".png"), optimize=True)
    json.dump({"image": stem + ".png", "bounds": bounds, "dimensions": {"width": rgba.shape[1], "height": rgba.shape[0]},
               "source": source, "value_range": {"min": vmin, "max": vmax}},
              open(os.path.join(outdir, stem + "_meta.json"), "w"), indent=2)


def colour(a, stops, lo, hi, fade_floor=70, fade_knee=0.15, mask=None):
    """Colour a float array on lo..hi; alpha 0 on NaN, fading in over the lowest `fade_knee` of the scale."""
    ok = np.isfinite(a)
    t = np.clip((np.nan_to_num(a) - lo) / (hi - lo), 0, 1)
    rgba = np.zeros(a.shape + (4,), np.uint8)
    rgba[..., :3] = ramp(stops, t)
    rgba[..., 3] = np.where(ok, (fade_floor + (255 - fade_floor) * np.clip(t / fade_knee, 0, 1)).round(), 0)
    if mask is not None: rgba[..., 3] = np.where(mask & ok, rgba[..., 3], 0)
    return rgba


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("src"); ap.add_argument("--max-width", type=int, default=3200)
    ap.add_argument("--label-suffix", default="")
    ap.add_argument("--mask", default="/Volumes/Extreme SSD/RE-ENFORCE/MODEL_RUN_v6_2026-10-29/maps/native/spruce_ge020_mask.tif",
                    help="1 km mask of cells with spruce >= 0.20 (model domain); 'none' shows all land")
    args = ap.parse_args(); mask = None if args.mask == "none" else args.mask
    entries, periods, stats = {}, [], {}
    for key, tag in PERIODS:
        hz = find(args.src, "hazard_mean", tag)
        if not hz: log("missing hazard_mean for", key); continue
        outdir = os.path.join(BUNDLE, MODEL, key)
        stem = f"webmap_data/forecast_bundle/{MODEL}/{key}/"
        a, bounds = to_mercator(hz, args.max_width, mask)
        src_note = "Final model P1C4 (16 predictors), SSP5-8.5, CHELSA-CMIP6 MPI-ESM1-2-HR, year replay" + args.label_suffix
        write(outdir, "probability_ensemble_mean", colour(a, SEQ_STOPS, 0, CAP), bounds, src_note, 0, CAP)
        write(outdir, "risk_visible_t50", colour(a, SEQ_STOPS, 0, CAP, fade_floor=255, mask=np.nan_to_num(a) >= 0.5), bounds, src_note, 0, CAP)
        land = np.isfinite(a)
        stats[key] = {"mean": float(np.nanmean(a)), "flag_share": float((a[land] >= 0.5).mean())}
        entry = {"model": MODEL, "period": key, "scenario": "ssp585", "gcms": ["CHELSA-CMIP6 MPI-ESM1-2-HR"],
                 "note": "Hazard for spruce forest from bark-beetle and wind disturbance: a ranking calibrated at the sample prevalence, not an annual probability.",
                 "ensemble_probability_tif": stem + "probability_ensemble_mean.tif",
                 "threshold_visible_tifs": {"0.5": stem + "risk_visible_t50.tif"}}
        ch = find(args.src, "hazard_change", tag)
        if ch and key != "1981_2010":
            d, db = to_mercator(ch, args.max_width, mask)
            write(outdir, "delta_from_near", _delta_rgba(d), db, src_note + "; period minus 1981-2010", -CHANGE_CAP, CHANGE_CAP)
            entry["delta_from_near_tif"] = stem + "delta_from_near.tif"
            stats[key]["change_mean"] = float(np.nanmean(d)); stats[key]["rise_gt_0.1"] = float((d[np.isfinite(d)] > 0.1).mean())
        for kind, stemname, stops, field in (("share05", "share05", SEQ_STOPS, "share05_tif"), ("extrapolation", "extrapolation", GREY, "extrapolation_tif")):
            f = find(args.src, kind, tag)
            if f:
                s, sb = to_mercator(f, args.max_width, mask)
                write(outdir, stemname, colour(s, stops, 0, 1), sb, src_note, 0, 1); entry[field] = stem + stemname + ".tif"
        entries[f"{MODEL}::{key}"] = entry; periods.append(key); log("built", key, stats[key])
    if not periods: sys.exit("no hazard_mean rasters found in " + args.src)
    # metrics file: the validation numbers of the final model, in the shape the app's panel reads
    json.dump({"model": MODEL, "auc": 0.87,
               "notes": "Final bark-beetle model: boosted trees on 16 predictors (terrain, climate normals and 4-year climate memory), learnt from 223 documented outbreak events overlaid with the European Forest Disturbance Atlas, 1985-2024. AUC 0.87 with a region and a period both held out; 0.82 trained to 2014 and tested on 2017-2024. The layer is a susceptibility index, a ranking calibrated at the sample prevalence, not an annual attack probability and not a forecast of single years. One climate model (MPI-ESM1-2-HR), SSP5-8.5; SSP1-2.6 and SSP3-7.0 to follow.",
               "aggregate_by_cutoff": [{"cutoff": "0.5", "precision": "NA", "recall": "NA", "specificity": "NA", "F1": "NA", "red_fraction": "NA"}],
               "grid_items": [{"label": "AUC, region and period held out", "value": "0.87"}, {"label": "AUC, trained to 2014, tested 2017–2024", "value": "0.82"}, {"label": "Flag threshold", "value": "0.5"}],
               "validation_note": "Case-control samples: precision and recall do not describe real prevalence, so they are not shown. The index ranks where disturbance is more likely over the long term; it does not track single years."},
              open(os.path.join(BUNDLE, MODEL, "metrics_confusion.json"), "w"), indent=2, ensure_ascii=False)
    mp = os.path.join(BUNDLE, "manifest.json"); m = json.load(open(mp))
    m["models"].pop(MODEL, None); [m["entries"].pop(k) for k in [k for k in m["entries"] if k.startswith(MODEL + "::")]]
    m["models"] = {MODEL: {"metrics_file": f"webmap_data/forecast_bundle/{MODEL}/metrics_confusion.json",
                           "label": "Final model (P1C4), SSP5-8.5", "scenarios": ["ssp585"], "periods": periods,
                           "value_max": CAP, "change_cap": CHANGE_CAP, "flag_threshold": 0.5, "default": True}, **m["models"]}
    m["entries"].update(entries)
    json.dump(m, open(mp, "w"), indent=2)
    json.dump(stats, open(os.path.join(BUNDLE, MODEL, "stats.json"), "w"), indent=2)
    log("manifest updated:", MODEL, periods)


def _delta_rgba(d):
    """Diverging colours centred on 0; alpha fades to a faint floor where the change is ~0."""
    ok = np.isfinite(d); t = (np.clip(np.nan_to_num(d), -CHANGE_CAP, CHANGE_CAP) + CHANGE_CAP) / (2 * CHANGE_CAP)
    rgba = np.zeros(d.shape + (4,), np.uint8); rgba[..., :3] = ramp(DIV, t)
    mag = np.clip(np.abs(np.nan_to_num(d)) / 0.05, 0, 1)
    rgba[..., 3] = np.where(ok, (40 + 215 * mag).round(), 0)
    return rgba


if __name__ == "__main__":
    main()
