"""Synthetic global rasters used ONLY for automated tests (tagged MARTIAN_MAP_SYNTHETIC=1)."""
import numpy as np
import rasterio
from rasterio.transform import from_origin

PPD = 8  # pixels per degree (the real DEM has 128)


def write(path, func, ppd=PPD):
    W, H = 360 * ppd, 180 * ppd
    lon = -180 + (np.arange(W) + 0.5) / ppd
    lat = 90 - (np.arange(H) + 0.5) / ppd
    LON, LAT = np.meshgrid(lon, lat)
    z = func(LAT, LON).astype("float32")
    with rasterio.open(path, "w", driver="GTiff", width=W, height=H, count=1, dtype="float32",
                       transform=from_origin(-180, 90, 1 / ppd, 1 / ppd), nodata=-32768.0) as ds:
        ds.write(z, 1)
        ds.update_tags(MARTIAN_MAP_SYNTHETIC="1")
    return path


def ramp(lat, lon):          # 2000 m per degree of latitude (within +/-30 deg)
    return 2000.0 * np.clip(lat, -30, 30)


def seam(lat, lon):          # smooth across the +/-180 seam
    return 3000.0 * np.cos(np.radians(lon)) + 0 * lat


def ridge(lat, lon):         # steep N-S ridge at lon=0 with a low pass near lat=+4
    gap = np.exp(-(((lat - 4) / 1.5) ** 2))
    return 6000.0 * np.exp(-((lon / 0.35) ** 2)) * (1 - gap)


def index_pattern(lat, lon):  # unique value per pixel: row + col/10000
    H, W = lat.shape
    rows, cols = np.mgrid[0:H, 0:W]
    return rows + cols / 10000.0
