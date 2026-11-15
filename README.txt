MARTIAN MAP — NASA TREK + MOLA BACKEND UPDATE

QUICK START (full features)
  cd backend
  pip install -r requirements.txt
  python download_dem.py      # once, ~2 GB NASA/USGS MOLA DEM
  python app.py               # open http://127.0.0.1:5000
See backend/README.md. Without the backend, index.html still works on its own (open it directly):
NASA Trek imagery, live image-service elevation, reference fallback, straight-line route draft.

WHAT'S IN THE APP
- Basemap: NASA Mars Trek tiles (MOLA colorized elevation default, Viking alternative, original
  illustrative map as offline fallback) with a status pill ("NASA imagery loaded" / "Illustrative fallback active").
- Elevation (Region Inspector), in this order:
    1. Backend raster sampling of the NASA MOLA DEM   -> "NASA-DERIVED · RASTER" (+ DEM slope from neighbouring cells)
    2. Live NASA Mars Trek / USGS image-service query -> "NASA-DERIVED · LIVE"
    3. Static approximate reference                   -> "NASA-derived reference · approximate"
  A details panel shows data source, coordinates, retrieval method and DEM slope.
- Terrain-based route (needs backend): least-cost path around steep terrain, with max grid slope,
  elevation gain, comparison with the straight line and an elevation profile.
  Without the backend: the original straight-line draft, clearly labelled.
- Labels: NASA IMAGERY / NASA-DERIVED (elevation, DEM slope, terrain-based path) / ESTIMATED (risk
  scores, route risk index) / ILLUSTRATIVE (hazard overlays). Existing risk scores and formulas are unchanged.
- Methodology: in-app card, METHODOLOGY.md, JSON export provenance.

LIMITS
Risk, rock and thermal scores are heuristic demonstration values, not NASA hazard products.
Terrain routes use grid-scale slopes only. Not for operational use.
