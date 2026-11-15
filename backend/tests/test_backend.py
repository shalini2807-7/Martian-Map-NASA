import math
import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.dirname(__file__))

import synthetic
from mola_dem import MARS_RADIUS_M as R, DemError, MolaDEM, haversine_m
from routing import plan_route

TMP = tempfile.mkdtemp()
P = {n: synthetic.write(os.path.join(TMP, n + ".tif"), f)
     for n, f in [("ramp", synthetic.ramp), ("seam", synthetic.seam), ("ridge", synthetic.ridge),
                  ("index", synthetic.index_pattern)]}


class TestSampling(unittest.TestCase):
    def test_exact_pixel_value_and_coordinates(self):
        dem = MolaDEM(P["index"])
        s = dem.sample(18.38, 77.58)
        col, row = int((77.58 + 180) * 8), int((90 - 18.38) * 8)
        self.assertEqual((s["col"], s["row"]), (col, row))
        self.assertAlmostEqual(s["elevation_m"], row + col / 10000.0, places=2)
        self.assertAlmostEqual(s["pixel_center_lat"], 90 - (row + .5) / 8, places=6)

    def test_slope_on_known_ramp(self):
        dem = MolaDEM(P["ramp"])
        s = dem.sample(5.0, 40.0)
        expected = math.degrees(math.atan(250.0 / (R * math.radians(1 / 8))))  # 250 m per pixel
        self.assertAlmostEqual(s["slope_deg"], expected, places=3)
        self.assertAlmostEqual(s["elevation_m"], 10000.0, delta=150)

    def test_seam_wrap(self):
        dem = MolaDEM(P["seam"])
        a, b = dem.sample(0, 179.99), dem.sample(0, -179.99)
        self.assertAlmostEqual(a["elevation_m"], b["elevation_m"], delta=3)
        self.assertIsNotNone(a["slope_deg"])
        self.assertLess(a["slope_deg"], 0.2)

    def test_pole_and_bad_input(self):
        dem = MolaDEM(P["ramp"])
        self.assertIsNotNone(dem.sample(89.99, 10)["elevation_m"])
        with self.assertRaises(DemError):
            dem.sample(95, 0)

    def test_missing_file(self):
        with self.assertRaises(DemError):
            MolaDEM("/nonexistent/file.tif")

    def test_synthetic_flag(self):
        self.assertTrue(MolaDEM(P["ramp"]).info()["synthetic"])


class TestSanityCheck(unittest.TestCase):
    def _make(self, name, func):
        import rasterio
        from rasterio.transform import from_origin
        import numpy as np
        ppd = 8
        W, H = 360 * ppd, 180 * ppd
        lon = -180 + (np.arange(W) + .5) / ppd
        lat = 90 - (np.arange(H) + .5) / ppd
        LON, LAT = np.meshgrid(lon, lat)
        path = os.path.join(TMP, name + ".tif")
        with rasterio.open(path, "w", driver="GTiff", width=W, height=H, count=1, dtype="float32",
                           transform=from_origin(-180, 90, 1 / ppd, 1 / ppd)) as ds:
            ds.write(func(LAT, LON).astype("float32"), 1)  # NOT tagged synthetic: behaves like real data
        return MolaDEM(path)

    def test_plausible_mars_raster_passes(self):
        import numpy as np
        f = lambda la, lo: 21000 * np.exp(-(((la - 18.65) ** 2 + (lo + 133.8) ** 2) / 4)) - 8000 * np.exp(-(((la + 42.4) ** 2 + (lo - 70.5) ** 2) / 9))
        self.assertTrue(self._make("plausible", f).sanity_check()["passed"])

    def test_wrong_units_or_flipped_raster_fails(self):
        import numpy as np
        km = lambda la, lo: 21 * np.exp(-(((la - 18.65) ** 2 + (lo + 133.8) ** 2) / 4)) - 8 * np.exp(-(((la + 42.4) ** 2 + (lo - 70.5) ** 2) / 9))
        self.assertFalse(self._make("km_units", km).sanity_check()["passed"])
        flipped = lambda la, lo: 21000 * np.exp(-(((-la - 18.65) ** 2 + (lo + 133.8) ** 2) / 4))
        self.assertFalse(self._make("flipped", flipped).sanity_check()["passed"])


class TestRouting(unittest.TestCase):
    def test_route_goes_through_the_pass_not_over_the_ridge(self):
        dem = MolaDEM(P["ridge"])
        r = plan_route(dem, (0.0, -15.0), (0.0, 15.0), avoid_steep=True)
        t, s = r["terrain_route"], r["straight_line"]
        self.assertGreater(s["max_elevation_m"], 2500)          # straight line climbs the ridge
        self.assertLess(t["max_elevation_m"], 1500)             # terrain route uses the low pass
        self.assertLess(t["cost_index_km"], s["cost_index_km"])
        self.assertLess(t["max_slope_deg"], s["max_slope_deg"])
        self.assertGreater(t["distance_km"], s["distance_km"])  # longer but easier
        lats = [w["latitude"] for w in r["waypoints"]]
        self.assertGreater(max(lats), 2)                        # detoured north toward the pass
        self.assertLess(len(r["waypoints"]), 12)                # smoothed into a few legs
        self.assertIn("detour", r["result_type"])
        self.assertGreater(len(r["path"]), 10)
        self.assertTrue(all(-180 <= p[1] < 180 for p in r["path"]))
        self.assertEqual(r["waypoints"][0]["latitude"], 0.0)
        self.assertEqual(r["waypoints"][-1]["longitude"], 15.0)

    def test_flat_terrain_is_nearly_straight(self):
        dem = MolaDEM(P["ramp"])
        r = plan_route(dem, (0.0, 10.0), (0.0, 30.0), avoid_steep=False)
        self.assertLess(r["terrain_route"]["distance_km"] / r["straight_line"]["distance_km"], 1.05)
        self.assertLessEqual(r["terrain_route"]["cost_index_km"], r["straight_line"]["cost_index_km"])  # never worse than direct

    def test_route_across_antimeridian(self):
        dem = MolaDEM(P["seam"])
        r = plan_route(dem, (0.0, 170.0), (0.0, -170.0))
        expected = haversine_m(0, 170, 0, -170) / 1000
        self.assertAlmostEqual(r["terrain_route"]["distance_km"], expected, delta=expected * 0.06)
        lons = [w["longitude"] for w in r["waypoints"]]
        self.assertTrue(any(abs(abs(x) - 180) < 12 for x in lons))

    def test_same_point_rejected(self):
        with self.assertRaises(DemError):
            plan_route(MolaDEM(P["ramp"]), (0, 0), (0, 0))


class TestApi(unittest.TestCase):
    def client(self, name):
        os.environ["MOLA_DEM_PATH"] = P[name]
        import app as appmod
        appmod._dem, appmod._dem_error = None, None
        appmod._sample_cached.cache_clear()
        return appmod.app.test_client()

    def test_health_elevation_route(self):
        c = self.client("ridge")
        h = c.get("/api/health").get_json()
        self.assertTrue(h["dem_available"])
        e = c.get("/api/elevation?lat=0&lon=0").get_json()
        self.assertTrue(e["dataset"]["synthetic"])
        self.assertAlmostEqual(e["elevation_m"], 6000.0, delta=400)
        self.assertIn("raster sampling", e["retrieval"])
        r = c.get("/api/route?from_lat=0&from_lon=-15&to_lat=0&to_lon=15&avoid_steep=1").get_json()
        self.assertIn("waypoints", r)
        self.assertLess(r["terrain_route"]["cost_index_km"], r["straight_line"]["cost_index_km"])

    def test_errors_are_json(self):
        c = self.client("ramp")
        self.assertEqual(c.get("/api/elevation?lat=95&lon=0").status_code, 400)
        self.assertEqual(c.get("/api/elevation?lat=abc&lon=0").status_code, 400)
        self.assertEqual(c.get("/api/elevation?lat=1").status_code, 400)
        self.assertEqual(c.get("/api/route?from_lat=0&from_lon=0&to_lat=0&to_lon=0").status_code, 422)

    def test_missing_dem_reports_clearly(self):
        os.environ["MOLA_DEM_PATH"] = "/nonexistent.tif"
        import app as appmod
        appmod._dem, appmod._dem_error = None, None
        c = appmod.app.test_client()
        h = c.get("/api/health").get_json()
        self.assertFalse(h["dem_available"])
        self.assertIn("Cannot open DEM", h["error"])
        self.assertEqual(c.get("/api/elevation?lat=0&lon=0").status_code, 422)


if __name__ == "__main__":
    unittest.main(verbosity=2)
