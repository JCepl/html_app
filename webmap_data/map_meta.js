const MAP_META = {
  bbox: {
    xmin: -10.75788667,
    ymin: 34.83306535,
    xmax: 31.77414497,
    ymax: 71.24297619
  },
  valueRange: { zmin: 0, zmax: 0.455476857 },
  years: [2030, 2040, 2050, 2060, 2070, 2080, 2090, 2100],
  rcps: ["rcp26", "rcp45", "rcp85"],
  files: {
    "rcp26_2030": "webmap_data/geojson/bbtl_rcp26_2030.geojson",
    "rcp26_2040": "webmap_data/geojson/bbtl_rcp26_2040.geojson",
    "rcp26_2050": "webmap_data/geojson/bbtl_rcp26_2050.geojson",
    "rcp26_2060": "webmap_data/geojson/bbtl_rcp26_2060.geojson",
    "rcp26_2070": "webmap_data/geojson/bbtl_rcp26_2070.geojson",
    "rcp26_2080": "webmap_data/geojson/bbtl_rcp26_2080.geojson",
    "rcp26_2090": "webmap_data/geojson/bbtl_rcp26_2090.geojson",
    "rcp26_2100": "webmap_data/geojson/bbtl_rcp26_2100.geojson",
    "rcp45_2030": "webmap_data/geojson/bbtl_rcp45_2030.geojson",
    "rcp45_2040": "webmap_data/geojson/bbtl_rcp45_2040.geojson",
    "rcp45_2050": "webmap_data/geojson/bbtl_rcp45_2050.geojson",
    "rcp45_2060": "webmap_data/geojson/bbtl_rcp45_2060.geojson",
    "rcp45_2070": "webmap_data/geojson/bbtl_rcp45_2070.geojson",
    "rcp45_2080": "webmap_data/geojson/bbtl_rcp45_2080.geojson",
    "rcp45_2090": "webmap_data/geojson/bbtl_rcp45_2090.geojson",
    "rcp45_2100": "webmap_data/geojson/bbtl_rcp45_2100.geojson",
    "rcp85_2030": "webmap_data/geojson/bbtl_rcp85_2030.geojson",
    "rcp85_2040": "webmap_data/geojson/bbtl_rcp85_2040.geojson",
    "rcp85_2050": "webmap_data/geojson/bbtl_rcp85_2050.geojson",
    "rcp85_2060": "webmap_data/geojson/bbtl_rcp85_2060.geojson",
    "rcp85_2070": "webmap_data/geojson/bbtl_rcp85_2070.geojson",
    "rcp85_2080": "webmap_data/geojson/bbtl_rcp85_2080.geojson",
    "rcp85_2090": "webmap_data/geojson/bbtl_rcp85_2090.geojson",
    "rcp85_2100": "webmap_data/geojson/bbtl_rcp85_2100.geojson"
  }
};

// Same source (Grunig et al. 2026, Science; SVD disturbance simulations,
// processed_data/dist_rates_50km_shinyapp/) and same 50km reference grid as
// MAP_META (bark beetle) above -- just the fire and wind agents instead of
// bbtl. Built the same way: agent+metric='rate' filter, null (non-forest)
// hexes dropped, tiny negative numerical-smoothing artifacts clamped to 0,
// split by rcp x year, reprojected EPSG:3035 -> EPSG:4326.
const MAP_META_FIRE = {
  bbox: MAP_META.bbox,
  valueRange: { zmin: 0, zmax: 0.5 },
  years: MAP_META.years,
  rcps: MAP_META.rcps,
  files: {
    "rcp26_2030": "webmap_data/geojson/fire_rcp26_2030.geojson",
    "rcp26_2040": "webmap_data/geojson/fire_rcp26_2040.geojson",
    "rcp26_2050": "webmap_data/geojson/fire_rcp26_2050.geojson",
    "rcp26_2060": "webmap_data/geojson/fire_rcp26_2060.geojson",
    "rcp26_2070": "webmap_data/geojson/fire_rcp26_2070.geojson",
    "rcp26_2080": "webmap_data/geojson/fire_rcp26_2080.geojson",
    "rcp26_2090": "webmap_data/geojson/fire_rcp26_2090.geojson",
    "rcp26_2100": "webmap_data/geojson/fire_rcp26_2100.geojson",
    "rcp45_2030": "webmap_data/geojson/fire_rcp45_2030.geojson",
    "rcp45_2040": "webmap_data/geojson/fire_rcp45_2040.geojson",
    "rcp45_2050": "webmap_data/geojson/fire_rcp45_2050.geojson",
    "rcp45_2060": "webmap_data/geojson/fire_rcp45_2060.geojson",
    "rcp45_2070": "webmap_data/geojson/fire_rcp45_2070.geojson",
    "rcp45_2080": "webmap_data/geojson/fire_rcp45_2080.geojson",
    "rcp45_2090": "webmap_data/geojson/fire_rcp45_2090.geojson",
    "rcp45_2100": "webmap_data/geojson/fire_rcp45_2100.geojson",
    "rcp85_2030": "webmap_data/geojson/fire_rcp85_2030.geojson",
    "rcp85_2040": "webmap_data/geojson/fire_rcp85_2040.geojson",
    "rcp85_2050": "webmap_data/geojson/fire_rcp85_2050.geojson",
    "rcp85_2060": "webmap_data/geojson/fire_rcp85_2060.geojson",
    "rcp85_2070": "webmap_data/geojson/fire_rcp85_2070.geojson",
    "rcp85_2080": "webmap_data/geojson/fire_rcp85_2080.geojson",
    "rcp85_2090": "webmap_data/geojson/fire_rcp85_2090.geojson",
    "rcp85_2100": "webmap_data/geojson/fire_rcp85_2100.geojson"
  }
};

const MAP_META_WIND = {
  bbox: MAP_META.bbox,
  valueRange: { zmin: 0, zmax: 0.5 },
  years: MAP_META.years,
  rcps: MAP_META.rcps,
  files: {
    "rcp26_2030": "webmap_data/geojson/wind_rcp26_2030.geojson",
    "rcp26_2040": "webmap_data/geojson/wind_rcp26_2040.geojson",
    "rcp26_2050": "webmap_data/geojson/wind_rcp26_2050.geojson",
    "rcp26_2060": "webmap_data/geojson/wind_rcp26_2060.geojson",
    "rcp26_2070": "webmap_data/geojson/wind_rcp26_2070.geojson",
    "rcp26_2080": "webmap_data/geojson/wind_rcp26_2080.geojson",
    "rcp26_2090": "webmap_data/geojson/wind_rcp26_2090.geojson",
    "rcp26_2100": "webmap_data/geojson/wind_rcp26_2100.geojson",
    "rcp45_2030": "webmap_data/geojson/wind_rcp45_2030.geojson",
    "rcp45_2040": "webmap_data/geojson/wind_rcp45_2040.geojson",
    "rcp45_2050": "webmap_data/geojson/wind_rcp45_2050.geojson",
    "rcp45_2060": "webmap_data/geojson/wind_rcp45_2060.geojson",
    "rcp45_2070": "webmap_data/geojson/wind_rcp45_2070.geojson",
    "rcp45_2080": "webmap_data/geojson/wind_rcp45_2080.geojson",
    "rcp45_2090": "webmap_data/geojson/wind_rcp45_2090.geojson",
    "rcp45_2100": "webmap_data/geojson/wind_rcp45_2100.geojson",
    "rcp85_2030": "webmap_data/geojson/wind_rcp85_2030.geojson",
    "rcp85_2040": "webmap_data/geojson/wind_rcp85_2040.geojson",
    "rcp85_2050": "webmap_data/geojson/wind_rcp85_2050.geojson",
    "rcp85_2060": "webmap_data/geojson/wind_rcp85_2060.geojson",
    "rcp85_2070": "webmap_data/geojson/wind_rcp85_2070.geojson",
    "rcp85_2080": "webmap_data/geojson/wind_rcp85_2080.geojson",
    "rcp85_2090": "webmap_data/geojson/wind_rcp85_2090.geojson",
    "rcp85_2100": "webmap_data/geojson/wind_rcp85_2100.geojson"
  }
};
