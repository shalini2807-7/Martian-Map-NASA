# Martian Map backend (elevation, slope, terrain routing)

A small Flask service that reads the **NASA MGS MOLA DEM** raster and answers three questions:

| Endpoint | Purpose |
|---|---|
| `GET /api/health` | Is the DEM available? Dataset details and a sanity-check result. |
| `GET /api/elevation?lat=18.38&lon=77.58` | Elevation of the pixel covering the point + slope from the 8 neighbouring cells (Horn 3x3). |
| `GET /api/route?from_lat=&from_lon=&to_lat=&to_lon=&avoid_steep=1` | Terrain-aware route (A* least-cost + smoothing), compared with the direct line. |

## Dataset
**Mars MGS MOLA DEM 463m** (USGS Astrogeology, from NASA MOLA MEGDR)
- page: https://astrogeology.usgs.gov/search/map/mars_mgs_mola_dem_463m
- file: https://planetarymaps.usgs.gov/mosaic/Mars_MGS_MOLA_DEM_mosaic_global_463m.tif (about 2 GB)
- 46080 x 23040 px, 16-bit, 128 px/degree (~463 m/px), simple cylindrical, planetocentric latitude,
  longitude positive east (-180..180), heights in metres above the Mars areoid, vertical uncertainty >= +/-3 m.

## Run it
```bash
cd backend
python -m venv .venv && source .venv/bin/activate      # Windows: .venv\Scripts\activate
pip install -r requirements.txt
python download_dem.py        # one time, about 2 GB, resumable -> backend/data/
python app.py                 # then open http://127.0.0.1:5000
```
`index.html` is served at `/`. Opening `index.html` directly from disk also works while the backend runs.
Without the download, the backend reads the USGS file over HTTP range requests: point elevation works,
but routes read large windows and are slow. Set `MOLA_DEM_PATH=/path/to/file.tif` to use another copy.

## How it works
* **Coordinate to pixel:** `col = floor((lon + 180) * 128)`, `row = floor((90 - lat) * 128)`; the pixel value is returned (no interpolation). Longitude wraps at +/-180.
* **Slope:** Horn's method on the 3x3 native cells, with the real cell size in metres at that latitude.
* **Route:** the corridor is read at reduced resolution (block average, <= ~300 cells per side), then 8-connected A* with step cost `length * (1 + (slope_along/ref)^2 + 0.5 (slope_cell/ref)^2)` (ref 3 deg with "avoid steep", else 8 deg; steps over 20/30 deg or into no-data are blocked), then shortened into straight legs wherever that costs no more. The direct great-circle line is scored with the same model and the cheaper result is returned. Routes may cross the +/-180 seam.
* **Sanity check:** at start-up Olympus Mons must read above +15,000 m and Hellas below -5,000 m, otherwise `/api/health` reports a warning (wrong units or orientation).

## Limits (also shown in the app)
* Slopes are grid-scale averages over several kilometres; rocks, sand, local slopes and thermal conditions are not captured. "Minimize rocky zones" only affects the estimated risk index, not the path.
* Cost units are a relative planning index, not rover energy or time. Not for operational use.
* Routes over about 175 degrees of arc are rejected.

## Tests
`python -m unittest discover -s tests` uses small **synthetic** rasters (tagged, and shown in the app as "SYNTHETIC TEST DATA"), never NASA data. They check pixel lookup, slope against an analytic ramp, antimeridian wrap, detours around a ridge, the API and the sanity check.
