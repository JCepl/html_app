#!/usr/bin/env python3
"""Build the local country-border file for the DSS from GISCO (EU) CNTR_RG_<res>_2024_4326.geojson.

Usage: python3 build_country_borders.py <downloaded geojson> <output geojson>
Keeps only polygon parts that touch the map window (drops overseas territories), keeps the three
name properties the app reads, and rounds coordinates to 3 decimals (~100 m).
Source: GISCO, (C) EuroGeographics for the administrative boundaries.
"""
import json, sys

WINDOW = (-30.0, 30.0, 50.0, 75.0)   # xmin, ymin, xmax, ymax: map extent plus margin
KEEP = ("CNTR_ID", "CNTR_NAME", "NAME_ENGL")

def ring(r):
    out = []
    for x, y in r:
        p = [round(x, 3), round(y, 3)]
        if not out or out[-1] != p:
            out.append(p)
    return out if len(out) >= 4 else None

def touches(poly):
    xs = [c[0] for c in poly[0]]; ys = [c[1] for c in poly[0]]
    return min(xs) <= WINDOW[2] and max(xs) >= WINDOW[0] and min(ys) <= WINDOW[3] and max(ys) >= WINDOW[1]

src, dst = sys.argv[1], sys.argv[2]
data = json.load(open(src))
feats = []
for f in data["features"]:
    g = f["geometry"]
    polys = [g["coordinates"]] if g["type"] == "Polygon" else g["coordinates"] if g["type"] == "MultiPolygon" else []
    keep = []
    for poly in polys:
        if not touches(poly): continue
        rings = [ring(r) for r in poly]
        if rings[0]: keep.append([r for r in rings if r])
    if not keep: continue
    feats.append({"type": "Feature", "properties": {k: f["properties"].get(k) for k in KEEP},
                  "geometry": {"type": "MultiPolygon", "coordinates": keep}})
json.dump({"type": "FeatureCollection", "features": feats}, open(dst, "w"), separators=(",", ":"), ensure_ascii=False)
print(len(feats), "countries ->", dst)
