"""
Reader for the USGS/NASA "Mars MGS MOLA DEM 463m" global elevation raster.

Dataset facts (USGS Astropedia, https://astrogeology.usgs.gov/search/map/mars_mgs_mola_dem_463m):
  * 46080 x 23040 pixels, 16-bit, 128 pixels/degree, ~463 m/pixel at the equator
  * Simple cylindrical, planetocentric latitude, longitude positive east, -180..180
  * Heights in metres above the Mars areoid; sphere radius 3,396,190 m (IAU 2000)
  * Vertical uncertainty is at least +/-3 m (areoid error)

The reader assumes a *global, north-up* raster and maps latitude/longitude to pixels
with the global grid (it works whether the file is stored in degrees or in metres).
"""
import math
import threading

import numpy as np
import rasterio
from rasterio.enums import Resampling
from rasterio.windows import Window

MARS_RADIUS_M = 3396190.0
DEFAULT_URL = "https://planetarymaps.usgs.gov/mosaic/Mars_MGS_MOLA_DEM_mosaic_global_463m.tif"


class DemError(Exception):
    """Raised for data problems (missing file, nodata at a point, bad geometry)."""


def wrap_lon(lon):
    """Wrap a longitude to -180 <= lon < 180."""
    return ((lon + 180.0) % 360.0) - 180.0


def haversine_m(lat1, lon1, lat2, lon2):
    p1, p2 = math.radians(lat1), math.radians(lat2)
    dphi, dlmb = p2 - p1, math.radians(lon2 - lon1)
    a = math.sin(dphi / 2) ** 2 + math.cos(p1) * math.cos(p2) * math.sin(dlmb / 2) ** 2
    return 2 * MARS_RADIUS_M * math.asin(min(1.0, math.sqrt(a)))


def horn_slope_deg(win, dx_m, dy_m):
    """Slope (degrees) at the centre of a 3x3 window using Horn's method (as in GDAL)."""
    if win.shape != (3, 3) or not np.isfinite(win).all():
        return None
    a, b, c, d, _, f, g, h, i = win.ravel()
    dzdx = ((c + 2 * f + i) - (a + 2 * d + g)) / (8.0 * dx_m)
    dzdy = ((g + 2 * h + i) - (a + 2 * b + c)) / (8.0 * dy_m)
    return math.degrees(math.atan(math.hypot(dzdx, dzdy)))


class MolaDEM:
    def __init__(self, path):
        self.path = path
        try:
            self.ds = rasterio.open(path)
        except Exception as exc:  # rasterio raises many types
            raise DemError(f"Cannot open DEM '{path}': {exc}") from exc
        self.width, self.height = self.ds.width, self.ds.height
        if abs(self.width / self.height - 2.0) > 0.01:
            raise DemError(f"DEM is not a global 2:1 raster ({self.width}x{self.height}).")
        if self.ds.transform.e >= 0:
            raise DemError("DEM must be north-up.")
        self.lon0 = -180.0 if self.ds.bounds.left < -1 else 0.0  # -180..180 or 0..360 domain
        self.ppd = self.width / 360.0  # pixels per degree
        self.nodata = self.ds.nodata
        self.scale = self.ds.scales[0] or 1.0
        self.offset = self.ds.offsets[0] or 0.0
        self.tags = self.ds.tags()
        self.synthetic = str(self.tags.get("MARTIAN_MAP_SYNTHETIC", "0")) == "1"
        self.m_per_px_equator = 2 * math.pi * MARS_RADIUS_M / self.width
        self._lock = threading.Lock()
        self.sanity = None

    # ------------------------------------------------------------------ metadata
    def info(self):
        return {
            "name": "Synthetic test raster (NOT NASA data)" if self.synthetic else "Mars MGS MOLA DEM 463m",
            "publisher": "Test fixture" if self.synthetic else "USGS Astrogeology Science Center (NASA MGS MOLA Team, MEGDR)",
            "source_url": "n/a" if self.synthetic else "https://astrogeology.usgs.gov/search/map/mars_mgs_mola_dem_463m",
            "file": str(self.path),
            "pixels": [self.width, self.height],
            "pixels_per_degree": round(self.ppd, 3),
            "resolution_m_per_px": round(self.m_per_px_equator, 1),
            "vertical_datum": "Mars areoid (GMM-2B); heights in metres",
            "latitude_type": "planetocentric",
            "longitude": "positive east",
            "vertical_uncertainty_m": "at least +/-3 m (areoid error)",
            "synthetic": self.synthetic,
            "sanity_check": self.sanity,
        }

    def sanity_check(self):
        """
        Guard against wrong units/orientation: Olympus Mons must read very high and the Hellas basin
        very low (MOLA: about +21,000 m and about -8,000 m). Skipped for synthetic test rasters.
        """
        if self.synthetic:
            return {"performed": False, "reason": "synthetic test raster"}
        def extreme(lat, lon, fn):
            col, row = self.pixel_of(lat, lon)
            half = int(0.3 * self.ppd)  # +/-0.3 degrees tolerates small position differences
            win = self._read_native(col - half, row - half, 2 * half + 1, 2 * half + 1)
            return float(fn(win)) if np.isfinite(win).any() else float("nan")
        hi = extreme(18.65, -133.8, np.nanmax)   # Olympus Mons summit
        lo = extreme(-42.4, 70.5, np.nanmin)     # Hellas Planitia floor
        return {"performed": True, "passed": bool(hi > 15000 and lo < -5000),
                "olympus_max_m": round(hi), "hellas_min_m": round(lo),
                "expected": "Olympus > +15,000 m and Hellas < -5,000 m (metres above the areoid)"}

    # ------------------------------------------------------------------ low level reads
    def _rd(self, col, row, w, h, out_shape=None):
        with self._lock:
            arr = self.ds.read(1, window=Window(col, row, w, h), out_shape=out_shape,
                               resampling=Resampling.average, masked=True)
        data = np.ma.filled(arr.astype("float32"), np.nan)
        return data * np.float32(self.scale) + np.float32(self.offset)

    def _read_native(self, col0, row0, w, h):
        """Native-resolution window; columns wrap around the globe, rows clamp at the poles."""
        rows = np.clip(np.arange(row0, row0 + h), 0, self.height - 1)
        r_lo, r_hi = int(rows.min()), int(rows.max())
        col0 %= self.width
        if col0 + w <= self.width:
            block = self._rd(col0, r_lo, w, r_hi - r_lo + 1)
        else:
            w1 = self.width - col0
            block = np.hstack([self._rd(col0, r_lo, w1, r_hi - r_lo + 1),
                               self._rd(0, r_lo, w - w1, r_hi - r_lo + 1)])
        return block[rows - r_lo, :]

    def pixel_of(self, lat, lon):
        lon = self.lon0 + ((lon - self.lon0) % 360.0)
        col = int(math.floor((lon - self.lon0) * self.ppd)) % self.width
        row = min(self.height - 1, max(0, int(math.floor((90.0 - lat) * self.ppd))))
        return col, row

    def pixel_center(self, col, row):
        lat = 90.0 - (row + 0.5) / self.ppd
        lon = wrap_lon(self.lon0 + (col + 0.5) / self.ppd)
        return lat, lon

    def pixel_size_m(self, lat):
        dy = MARS_RADIUS_M * math.radians(1.0 / self.ppd)
        dx = MARS_RADIUS_M * max(0.02, math.cos(math.radians(lat))) * math.radians(1.0 / self.ppd)
        return dx, dy

    # ------------------------------------------------------------------ point query
    def sample(self, lat, lon):
        """Elevation of the pixel covering (lat, lon) plus Horn slope from its 8 neighbours."""
        if not (-90.0 <= lat <= 90.0) or not math.isfinite(lon):
            raise DemError("Latitude must be within -90..90 and longitude finite.")
        col, row = self.pixel_of(lat, lon)
        win = self._read_native(col - 1, row - 1, 3, 3)
        z = float(win[1, 1])
        if not math.isfinite(z):
            raise DemError("The DEM has no data at this location.")
        clat, clon = self.pixel_center(col, row)
        dx, dy = self.pixel_size_m(clat)
        slope = horn_slope_deg(win, dx, dy)
        return {
            "elevation_m": z,
            "slope_deg": slope,
            "col": col, "row": row,
            "pixel_center_lat": clat, "pixel_center_lon": clon,
            "pixel_size_m": [dx, dy],
        }

    # ------------------------------------------------------------------ block read for routing
    def read_block(self, lat_min, lat_max, lon_lo, lon_hi, max_dim=300):
        """
        Read a (possibly decimated) block. lon_lo < lon_hi are *unwrapped* degrees, so a block that
        crosses the +/-180 seam is read as two pieces and joined. Decimation uses block averaging.
        Returns dict(elev[ny,nx] float32 NaN=nodata, lats[ny], lons[nx] unwrapped, cell_deg, factor).
        """
        lat_min, lat_max = max(lat_min, -88.0), min(lat_max, 88.0)
        r0 = max(0, int(math.floor((90.0 - lat_max) * self.ppd)))
        r1 = min(self.height, int(math.ceil((90.0 - lat_min) * self.ppd)))
        c0 = int(math.floor((lon_lo - self.lon0) * self.ppd))
        c1 = int(math.ceil((lon_hi - self.lon0) * self.ppd))
        ncols, nrows = c1 - c0, r1 - r0
        if ncols <= 0 or nrows <= 0:
            raise DemError("Empty block.")
        if ncols > self.width:
            raise DemError("Route corridor is wider than the planet; choose closer points.")
        f = max(1, math.ceil(max(ncols, nrows) / float(max_dim)))
        nrows_f = nrows // f
        if nrows_f < 3:
            raise DemError("Route corridor too small.")
        pieces, lon_centres, cur = [], [], c0
        while cur < c1:
            base = cur % self.width
            take = min(c1 - cur, self.width - base)
            nf = take // f
            if nf > 0:
                pieces.append(self._rd(base, r0, nf * f, nrows_f * f, out_shape=(nrows_f, nf)))
                lon_centres.append(self.lon0 + (cur + (np.arange(nf) + 0.5) * f) / self.ppd)
            cur += take
        elev = np.hstack(pieces)
        lons = np.concatenate(lon_centres)
        lats = 90.0 - (r0 + (np.arange(nrows_f) + 0.5) * f) / self.ppd
        return {"elev": elev, "lats": lats, "lons": lons, "cell_deg": f / self.ppd, "factor": f}
