# Martian Map — Methodology

## 1. Datasets

| # | Exact dataset | Resolution | Source URL | Use | Live or static |
|---|---|---|---|---|---|
| 1 | Mars MGS MOLA Global Color Shaded Relief 463m (Trek layer `Mars_MGS_MOLA_ClrShade_merge_global_463m`) | 128 px/deg, ~463 m/px at equator (shown at WMTS zoom 3) | https://astrogeology.usgs.gov/search/map/mars_mgs_mola_global_color_shaded_relief_463m — served by https://trek.nasa.gov/tiles/apidoc/trekAPI.html?body=mars | Default basemap (NASA imagery) | Live tiles; offline fallback = illustrative map |
| 2 | Viking MDIM 2.1 Colorized Global Mosaic (Trek layer `Viking_MDIM21_ClrMosaic_global_232m`) | ~232 m/px | https://trek.nasa.gov/tiles/Mars/EQ/Viking_MDIM21_ClrMosaic_global_232m/1.0.0/WMTSCapabilities.xml | Alternative basemap | Live tiles |
| 3 | Mars MGS MOLA DEM 463m (from MOLA MEGDR data; USGS Astrogeology / NASA GSFC) | 128 px/deg, ~463 m/px, 16-bit; vertical uncertainty at least ±3 m; heights relative to the Mars areoid (global range ~ −8,200 to +21,229 m) | https://astrogeology.usgs.gov/search/map/mars_mgs_mola_dem_463m ; https://pds-geosciences.wustl.edu/missions/mgs/megdr.html | Numerical site elevation | Elevation order: (1) **backend raster sampling** of the local file, (2) **live** image-service query, (3) **static approximate reference** only if both fail |
| 4 | Heuristic risk scores, hazard overlays, route | n/a | n/a (project-authored) | Demonstration comparison | Static, built into the page (not NASA data) |

Elevation fallback chain (first success wins; the Region Inspector names the source that answered):
1. NASA Mars Trek MOLA DEM image service.
2. NASA Mars Trek MOLA + MEX HRSC blended DEM (200 m/px) — https://astrogeology.usgs.gov/search/map/mars_mgs_mola_mex_hrsc_blended_dem_global_200m
3. Esri-hosted copy of the USGS MOLA + HRSC blended DEM.
4. Static approximate reference values, labelled "NASA-derived reference · approximate" (hand-entered, not retrieved): Jezero ≈ −2,600 m, Gale ≈ −4,500 m, Olympus Mons ≈ +21,200 m, Valles Marineris ≈ −4,000 m, Arcadia ≈ −3,500 m.

## 2. How coordinates are matched to elevation data
Each site keeps its fixed latitude/longitude (decimal degrees, east-positive longitude). The point is sent to the DEM image service's `identify` operation, which returns the single DEM pixel (~463 m across) covering that point; there is no averaging or interpolation. The app tries longitude as −180…180 first, then 0…360, because services differ. Values outside −9,000…+22,000 m are rejected. One cell does not describe local slope or relief and may differ from published landing-site elevations.

## 3. Risk formulas (unchanged)
- Site terrain risk (0–100) and the slope, rock and thermal scores are preset demonstration values. Labels: Low < 40, Moderate < 70, High ≥ 70; CAUTION ≥ 75.
- Route distance is great-circle on a 3,389.5 km sphere.
- Route risk = mean of the two endpoint risk scores − 8 (if "Avoid steep terrain") − 5 (if "Minimize rocky zones") + 0.08 × latitude difference in degrees, rounded and limited to 8–99. Status: REVIEW ≥ 70, DRAFT ≥ 40, otherwise LOWER RISK.

## 4. How NASA data influences exploration decisions
- **Imagery** gives visual context for where a candidate sits (basins, highlands, canyons).
- **Elevation** is shown beside the estimated risk so a person can spot extreme depth or height (e.g. Olympus Mons) and request closer study.
- NASA data is **not an input to the risk formulas**. The output is a shortlist to validate with high-resolution products (HiRISE, CTX), not an operational decision.

## 5. Data labels
NASA IMAGERY · NASA-DERIVED ELEVATION · ESTIMATED · ILLUSTRATIVE. The map also shows a basemap status ("NASA imagery loaded" / "Illustrative fallback active").

## 6. Backend: elevation, DEM slope and terrain routing
- **Raster sampling:** pixel column = floor((lon + 180) x 128), row = floor((90 - lat) x 128) on the 46080 x 23040 USGS MOLA DEM; the covering pixel's value is returned (no interpolation). Planetocentric latitude, east-positive longitude, heights above the Mars areoid.
- **DEM slope (NASA-derived):** Horn's method on the 3 x 3 neighbouring native cells using the true cell size in metres at that latitude. It is shown next to, and does not replace, the estimated "Slope exposure" score; no risk score was changed.
- **Terrain route (NASA-derived path, estimated risk):** block-averaged DEM corridor; 8-connected A* with step cost = length x (1 + (slope_along / ref)^2 + 0.5 (slope_cell / ref)^2), ref = 3 deg with "Avoid steep terrain" (8 deg otherwise); steps over 20 deg (30 deg) or into no-data are blocked; the path is then shortened into straight legs where that costs no more; the direct great-circle line is scored with the same model and the cheaper is returned. The route risk index keeps its original heuristic formula.
- **What it does not use:** no rock-abundance, thermal, dust or lighting data. "Minimize rocky zones" changes only the estimated risk index.
- **Limits:** grid-scale slopes (several km per cell) understate local slopes; cost units are a relative index, not energy or time; routes over about 175 degrees of arc are rejected.
- **Safeguards:** start-up sanity check (Olympus Mons > +15,000 m, Hellas < -5,000 m); tests use synthetic rasters that the app labels "SYNTHETIC TEST DATA".
