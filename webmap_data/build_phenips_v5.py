#!/usr/bin/env python3
"""Stage the v5 PHENIPS (degree-days, generations; sun/shade) map images as WebP + one manifest.
Usage: python3 build_phenips_v5.py [dss_staging_dir]"""
import json, os, sys, glob
from PIL import Image
src = sys.argv[1] if len(sys.argv) > 1 else "/Volumes/Extreme SSD/RE-ENFORCE/MODEL_RUN_v5_2026-10-28/maps/dss_staging"
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "phenips_v5")
os.makedirs(out, exist_ok=True)
man = {"bounds": None, "layers": {}}
n = 0
for d in sorted(glob.glob(src + "/*/")):
    name = os.path.basename(d.rstrip("/"))
    if not name.startswith("phenips_"):
        continue
    m = json.load(open(f"{d}{name}_meta.json"))
    Image.open(f"{d}{name}.png").save(f"{out}/{name}.webp", quality=80, method=6)
    man["bounds"] = man["bounds"] or m["bounds"]
    man["layers"][name] = {"min": m["value_range"]["min"], "max": m["value_range"]["max"], "colours": m["colours"], "note": m["note"]}
    n += 1
json.dump(man, open(out + "/manifest.json", "w"), separators=(",", ":"))
print(n, "layers,", sum(os.path.getsize(f) for f in glob.glob(out + "/*")) // 1000000, "MB")
