#!/usr/bin/env python3
"""Re-render the drought class maps at 2400 px (export_drought_overlay.R used 900 px with modal aggregation,
i.e. ~9x coarser than the 1 km source). Same colours/classes/band order; palette PNG (lossless, small).
Usage: python3 build_drought_overlays_hires.py [dir_with_the_two_tifs]"""
import json, math, os, subprocess, sys, tempfile
import numpy as np
from PIL import Image
SRC = sys.argv[1] if len(sys.argv) > 1 else "/Volumes/Extreme SSD/RE-ENFORCE_D211_DELIVERY_2026-10-31/3_Drought_SPEI_Stefan"
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = HERE + "/drought"; W = 2400
REF = {"n_events_2yr": (2, {0: "#F7F7F7", 1: "#fee391", 2: "#fe9929", 3: "#8c2d04"}),
       "max_duration": (3, {0: "#F7F7F7", 2: "#fee391", 3: "#fe9929", 4: "#cc4c02", 5: "#8c2d04"}),
       "total_duration": (4, {0: "#F7F7F7", 2: "#fee391", 3: "#fe9929", 4: "#ec7014", 5: "#cc4c02", 6: "#8c2d04", 7: "#8c2d04"})}
FUT = ["2041_2070_ssp370", "2041_2070_ssp585", "2071_2100_ssp370", "2071_2100_ssp585"]
FCOL = {1: "#B2182B", 2: "#EF8A62", 3: "#FDDBC7", 4: "#F7F7F7"}
from osgeo import gdal
gdal.UseExceptions()
R = 6378137.0
def render(tif, band, colors, out_png):
    one = gdal.Translate("", tif, format="VRT", bandList=[band])
    ds = gdal.Warp("", one, format="MEM", dstSRS="EPSG:3857", resampleAlg="mode", width=W, height=0,
                   dstNodata=-9999, outputType=gdal.GDT_Float32)
    a = ds.GetRasterBand(1).ReadAsArray(); gt = ds.GetGeoTransform()
    idx = np.zeros(a.shape, np.uint8); pal = [0, 0, 0]
    for i, (v, hexc) in enumerate(colors.items(), start=1):
        idx[np.abs(a - v) < 1e-6] = i
        pal += [int(hexc[1:3], 16), int(hexc[3:5], 16), int(hexc[5:7], 16)]
    im = Image.fromarray(idx, "P"); im.putpalette(pal + [0] * (768 - len(pal)))
    im.save(out_png, optimize=True, transparency=0)
    x0, y1 = gt[0], gt[3]; x1, y0 = x0 + gt[1] * ds.RasterXSize, y1 + gt[5] * ds.RasterYSize
    ll = lambda y: math.degrees(math.atan(math.sinh(y / R)))
    return {"south": ll(y0), "west": math.degrees(x0 / R), "north": ll(y1), "east": math.degrees(x1 / R)}

for key, (band, cols) in REF.items():
    b = render(SRC + "/myd_SPEI-12_reference_DSS.tif", band, cols, f"{OUT}/ref_{key}.png")
    m = json.load(open(f"{OUT}/drought_ref_meta.json")); m["images"][key]["bounds"] = b
    json.dump(m, open(f"{OUT}/drought_ref_meta.json", "w"), indent=2); print(key, os.path.getsize(f"{OUT}/ref_{key}.png") // 1000, "KB")
for i, key in enumerate(FUT, start=1):
    b = render(SRC + "/SPEI-equivalent-12_future_classes.tif", i, FCOL, f"{OUT}/future_{key}.png")
    m = json.load(open(f"{OUT}/drought_future_meta.json")); m["images"][key]["bounds"] = b
    json.dump(m, open(f"{OUT}/drought_future_meta.json", "w"), indent=2); print(key, os.path.getsize(f"{OUT}/future_{key}.png") // 1000, "KB")
