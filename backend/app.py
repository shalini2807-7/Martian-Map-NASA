"""
Martian Map backend: coordinate -> elevation/slope from the NASA MOLA DEM, and terrain-aware routes.

    python app.py            # then open http://127.0.0.1:5000

Configuration (environment variables):
    MOLA_DEM_PATH   path or URL of the DEM GeoTIFF
                    (default: backend/data/Mars_MGS_MOLA_DEM_mosaic_global_463m.tif if present,
                     otherwise the USGS file over HTTP range requests)
    PORT            default 5000
"""
import functools
import math
import os
import threading
from pathlib import Path

from flask import Flask, jsonify, request, send_from_directory

from mola_dem import DEFAULT_URL, DemError, MolaDEM, wrap_lon
from routing import plan_route

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent  # folder containing index.html
LOCAL_DEM = HERE / "data" / "Mars_MGS_MOLA_DEM_mosaic_global_463m.tif"

app = Flask(__name__)
_dem, _dem_error, _dem_lock = None, None, threading.Lock()


def dem_path():
    env = os.environ.get("MOLA_DEM_PATH")
    if env:
        return env
    return str(LOCAL_DEM) if LOCAL_DEM.exists() else "/vsicurl/" + DEFAULT_URL


def get_dem():
    """Open the DEM lazily so the server still starts (and reports a clear error) if it is missing."""
    global _dem, _dem_error
    with _dem_lock:
        if _dem is None:
            try:
                dem = MolaDEM(dem_path())
                dem.sanity = dem.sanity_check()
                _dem, _dem_error = dem, None
            except DemError as exc:
                _dem_error = str(exc)
                raise
        return _dem


@app.after_request
def add_cors(resp):
    resp.headers["Access-Control-Allow-Origin"] = "*"  # lets index.html opened from disk call the API
    return resp


@app.errorhandler(DemError)
def handle_dem_error(exc):
    return jsonify({"error": str(exc)}), 422


@app.errorhandler(ValueError)
def handle_bad_input(exc):
    return jsonify({"error": str(exc)}), 400


def num(name, lo=None, hi=None):
    raw = request.args.get(name)
    if raw is None:
        raise ValueError(f"Missing query parameter '{name}'.")
    try:
        val = float(raw)
    except ValueError:
        raise ValueError(f"Parameter '{name}' must be a number.") from None
    if not math.isfinite(val) or (lo is not None and val < lo) or (hi is not None and val > hi):
        raise ValueError(f"Parameter '{name}' is out of range.")
    return val


@app.get("/")
def index():
    return send_from_directory(ROOT, "index.html")


@app.get("/api/health")
def health():
    try:
        dem = get_dem()
        return jsonify({"status": "ok", "dem_available": True, "dataset": dem.info()})
    except DemError:
        return jsonify({"status": "degraded", "dem_available": False, "error": _dem_error,
                        "hint": "Run `python download_dem.py` or set MOLA_DEM_PATH."})


@functools.lru_cache(maxsize=512)
def _sample_cached(lat, lon):
    return get_dem().sample(lat, lon)


@app.get("/api/elevation")
def elevation():
    lat, lon = num("lat", -90, 90), wrap_lon(num("lon"))
    dem = get_dem()
    s = _sample_cached(round(lat, 5), round(lon, 5))
    slope = None if s["slope_deg"] is None else round(s["slope_deg"], 2)
    return jsonify({
        "dataset": dem.info(),
        "query": {"lat": lat, "lon": lon},
        "elevation_m": round(s["elevation_m"], 1),
        "slope_deg": slope,
        "slope_method": "Horn 3x3 window on native cells (about 463 m)",
        "pixel": {"row": s["row"], "col": s["col"], "center_lat": round(s["pixel_center_lat"], 5),
                  "center_lon": round(s["pixel_center_lon"], 5), "size_m": [round(v) for v in s["pixel_size_m"]]},
        "retrieval": "raster sampling (pixel covering the point; no interpolation)",
    })


@app.get("/api/route")
def route():
    a = (num("from_lat", -90, 90), num("from_lon"))
    b = (num("to_lat", -90, 90), num("to_lon"))
    steep = request.args.get("avoid_steep", "1") not in ("0", "false", "False")
    result = plan_route(get_dem(), a, b, avoid_steep=steep, max_dim=150)
    result["dataset"] = get_dem().info()
    return jsonify(result)


if __name__ == "__main__":
    port = int(os.environ.get("PORT", 5000))
    print(f"DEM source: {dem_path()}")
    app.run(host="127.0.0.1", port=port, threaded=True)
