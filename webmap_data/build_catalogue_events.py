#!/usr/bin/env python3
"""Rebuild the CATALOGUE_EVENTS line of app.js from the enhanced event catalogue
(base catalogue + DFDE additions). Usage: python3 build_catalogue_events.py [catalogue_enhanced_dir] [v2_data_dir]"""
import csv, json, re, sys, os
here = os.path.dirname(os.path.abspath(__file__))
enh = sys.argv[1] if len(sys.argv) > 1 else "/Volumes/Extreme SSD/RE-ENFORCE_D211_DELIVERY_2026-10-31/4_BarkBeetle_Jaroslav/catalogue_enhanced"
base = sys.argv[2] if len(sys.argv) > 2 else "/Volumes/Extreme SSD/RE-ENFORCE_BB_literature_v2_2026-10-30/data"
rd = lambda p: list(csv.DictReader(open(p, encoding="utf-8-sig")))
srcs = {s["source_id"]: s for s in rd(base + "/sources.csv") + rd(enh + "/sources_dfde.csv")}
links = rd(base + "/event_sources.csv") + rd(enh + "/event_sources_dfde.csv")
by_ev = {}
for l in links:
    by_ev.setdefault(l["event_id"], {})[l["source_id"]] = l["source_role"]
out, skipped = [], []
for e in rd(enh + "/events_enhanced_with_dfde.csv"):
    try:
        lat, lon = float(e["latitude"]), float(e["longitude"])
    except ValueError:
        skipped.append(e["event_id"]); continue
    sl = []
    for sid, role in by_ev.get(e["event_id"], {}).items():
        s = srcs.get(sid)
        if s: sl.append({"title": s["title"], "url": s["url"], "role": role, "citation": s["citation"]})
    sl.sort(key=lambda x: x["role"] != "primary")
    out.append({"id": e["event_id"], "name": e["event_name"], "country": e["country"],
                "y0": int(e["year_start"]), "y1": int(e["year_end"] or e["year_start"]),
                "lat": lat, "lon": lon, "sev": e["severity"], "conf": e["confidence"],
                "note": e["short_comment"], "sources": sl})
js = "const CATALOGUE_EVENTS = " + json.dumps(out, ensure_ascii=False, separators=(",", ":")) + ";"
p = os.path.join(here, "..", "app.js")
t = open(p, encoding="utf-8").read()
t2, n = re.subn(r"^    const CATALOGUE_EVENTS = .*;$", lambda m: "    " + js, t, count=1, flags=re.M)
assert n == 1
open(p, "w", encoding="utf-8").write(t2)
print(len(out), "events,", sum(len(x["sources"]) for x in out), "source links; skipped (no coordinates):", skipped)
print("severity:", sorted({x["sev"] for x in out}), "countries:", len({x["country"] for x in out}))
