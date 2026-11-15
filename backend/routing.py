"""
Terrain-aware route planning on the (block-averaged) MOLA DEM.

Method: 8-connected A* least-cost path. Cost of a step = length * (1 + (theta_dir/theta_ref)^2
+ 0.5*(theta_local/theta_ref)^2), where theta_dir is the slope along the step (from the elevation
difference of the two cell centres) and theta_local is the Horn slope of the destination cell.
Steps steeper than theta_max, or into no-data cells, are not allowed. The grid path is then smoothed
into straight legs wherever that does not cost more (any-angle shortcuts), and is compared with the
direct great-circle line using the very same cost model.

Limits: slopes are computed at grid scale (several km per cell), so they understate local slopes,
rocks, sand traps and other hazards. This is a planning aid for *regional* corridors only.
"""
import heapq
import math

import numpy as np

from mola_dem import DemError, MARS_RADIUS_M as R, haversine_m, wrap_lon

PROFILES = {
    True: {"theta_ref_deg": 3.0, "theta_max_deg": 20.0, "label": "avoid steep terrain"},
    False: {"theta_ref_deg": 8.0, "theta_max_deg": 30.0, "label": "standard"},
}
OFFSETS = [(-1, -1), (-1, 0), (-1, 1), (0, -1), (0, 1), (1, -1), (1, 0), (1, 1)]


def slope_grid(elev, lats, cell_deg):
    """Horn slope (degrees) for every cell; NaN neighbours give NaN (treated as 0 later)."""
    z = np.pad(elev, 1, mode="edge")
    a, b, c = z[:-2, :-2], z[:-2, 1:-1], z[:-2, 2:]
    d, f = z[1:-1, :-2], z[1:-1, 2:]
    g, h, i = z[2:, :-2], z[2:, 1:-1], z[2:, 2:]
    dy = R * math.radians(cell_deg)
    dx = R * np.maximum(0.02, np.cos(np.radians(lats))) * math.radians(cell_deg)
    dzdx = ((c + 2 * f + i) - (a + 2 * d + g)) / (8.0 * dx[:, None])
    dzdy = ((g + 2 * h + i) - (a + 2 * b + c)) / (8.0 * dy)
    return np.degrees(np.arctan(np.hypot(dzdx, dzdy)))


def great_circle_points(lat1, lon1, lat2, lon2, n):
    """n points along the great circle; longitudes returned *unwrapped* (continuous)."""
    p1, l1, p2, l2 = map(math.radians, (lat1, lon1, lat2, lon2))
    v1 = np.array([math.cos(p1) * math.cos(l1), math.cos(p1) * math.sin(l1), math.sin(p1)])
    v2 = np.array([math.cos(p2) * math.cos(l2), math.cos(p2) * math.sin(l2), math.sin(p2)])
    omega = math.acos(max(-1.0, min(1.0, float(v1 @ v2))))
    if omega < 1e-9:
        return np.array([lat1, lat2]), np.array([lon1, lon2])
    t = np.linspace(0.0, 1.0, max(2, n))[:, None]
    v = (np.sin((1 - t) * omega) * v1 + np.sin(t * omega) * v2) / math.sin(omega)
    lat = np.degrees(np.arcsin(np.clip(v[:, 2], -1, 1)))
    lon = np.degrees(np.unwrap(np.arctan2(v[:, 1], v[:, 0])))
    return lat, lon


def _astar(zl, sl, dx_row, dy, shape, start, goal, theta_ref, theta_max):
    ny, nx = shape
    n = ny * nx
    g = [math.inf] * n
    prev = [-1] * n
    closed = bytearray(n)
    s, t = start[0] * nx + start[1], goal[0] * nx + goal[1]
    gj, gi = goal

    def heur(idx):
        j, i = divmod(idx, nx)
        return math.hypot((i - gi) * dx_row[j], (j - gj) * dy)

    g[s] = 0.0
    heap = [(heur(s), s)]
    while heap:
        _, u = heapq.heappop(heap)
        if closed[u]:
            continue
        closed[u] = 1
        if u == t:
            break
        j, i = divmod(u, nx)
        zu, gu = zl[u], g[u]
        for dj, di in OFFSETS:
            nj, ni = j + dj, i + di
            if nj < 0 or nj >= ny or ni < 0 or ni >= nx:
                continue
            v = nj * nx + ni
            if closed[v]:
                continue
            zv = zl[v]
            if zv != zv:  # NaN = no data
                continue
            length = math.hypot(0.5 * (dx_row[j] + dx_row[nj]) * di, dy * dj)
            theta = math.degrees(math.atan2(abs(zv - zu), length))
            if theta > theta_max:
                continue
            ng = gu + length * (1.0 + (theta / theta_ref) ** 2 + 0.5 * (sl[v] / theta_ref) ** 2)
            if ng < g[v]:
                g[v] = ng
                prev[v] = u
                heapq.heappush(heap, (ng + heur(v), v))
    if not math.isfinite(g[t]):
        raise DemError("No traversable path found in the corridor (all routes blocked by steep or no-data cells).")
    path, u = [], t
    while u != -1:
        path.append(u)
        u = prev[u]
    path.reverse()
    return [divmod(p, nx) for p in path]


def _rdp(points, tol):
    """Ramer-Douglas-Peucker simplification (iterative). points: Nx2 array in cell units."""
    n = len(points)
    if n < 3:
        return list(range(n))
    keep = np.zeros(n, bool)
    keep[0] = keep[-1] = True
    stack = [(0, n - 1)]
    while stack:
        lo, hi = stack.pop()
        if hi <= lo + 1:
            continue
        p, q, seg = points[lo], points[hi], points[hi] - points[lo]
        norm = np.hypot(*seg)
        rel = points[lo + 1:hi] - p
        dist = np.hypot(*rel.T) if norm == 0 else np.abs(seg[0] * rel[:, 1] - seg[1] * rel[:, 0]) / norm
        k = int(np.argmax(dist))
        if dist[k] > tol:
            m = lo + 1 + k
            keep[m] = True
            stack.extend([(lo, m), (m, hi)])
    return list(np.flatnonzero(keep))


def hav_vec(lat1, lon1, lat2, lon2):
    p1, p2 = np.radians(lat1), np.radians(lat2)
    a = np.sin((p2 - p1) / 2) ** 2 + np.cos(p1) * np.cos(p2) * np.sin(np.radians(lon2 - lon1) / 2) ** 2
    return 2 * R * np.arcsin(np.minimum(1.0, np.sqrt(a)))


def nearest_idx(sorted_arr, vals):
    """Index of the nearest element of an ascending array for each value."""
    k = np.clip(np.searchsorted(sorted_arr, vals), 1, len(sorted_arr) - 1)
    left, right = sorted_arr[k - 1], sorted_arr[k]
    return np.where(np.abs(vals - left) <= np.abs(right - vals), k - 1, k)


class Grid:
    """Block-averaged DEM plus helpers to sample it along lines."""

    def __init__(self, blk, tref, tmax):
        self.z, self.lats, self.lons, self.cell = blk["elev"], blk["lats"], blk["lons"], blk["cell_deg"]
        self.ny, self.nx = self.z.shape
        self.tref, self.tmax = tref, tmax
        self.dy = R * math.radians(self.cell)
        self.sloc = np.nan_to_num(slope_grid(self.z, self.lats, self.cell), nan=0.0)
        self._neg_lats = -self.lats  # ascending

    def cell_of(self, lat, lon_u):
        return int(nearest_idx(self._neg_lats, np.array([-lat]))[0]), int(nearest_idx(self.lons, np.array([lon_u]))[0])

    def sample(self, lat, lon):
        j, i = nearest_idx(self._neg_lats, -np.asarray(lat)), nearest_idx(self.lons, np.asarray(lon))
        return self.z[j, i], self.sloc[j, i]

    def line(self, lat_a, lon_a, lat_b, lon_b):
        """Points along a straight (lat/lon-linear) leg at about one-cell spacing, with z and cell slope."""
        n = int(max(2, math.ceil(max(abs(lon_b - lon_a), abs(lat_b - lat_a)) / self.cell) + 1))
        lat, lon = np.linspace(lat_a, lat_b, n), np.linspace(lon_a, lon_b, n)
        z, sl = self.sample(lat, lon)
        return lat, lon, z, sl

    def cost(self, lat, lon, z, sl):
        """Cost of a sampled polyline with the same model as A*; inf if it hits no-data or a too-steep step."""
        if len(lat) < 2:
            return 0.0
        seg = hav_vec(lat[:-1], lon[:-1], lat[1:], lon[1:])
        theta = np.degrees(np.arctan2(np.abs(np.diff(z)), np.maximum(seg, 1e-6)))
        if not np.isfinite(z).all() or (theta > self.tmax).any():
            return math.inf
        return float(np.sum(seg * (1 + (theta / self.tref) ** 2 + 0.5 * (sl[1:] / self.tref) ** 2)))


def smooth_path(grid, plats, plons, step_cum, max_tries=120, window=400):
    """Greedy any-angle shortcutting: jump to the farthest vertex whose straight leg is not costlier."""
    n, out, i = len(plats), [0], 0
    while i < n - 1:
        best = i + 1
        for tries, k in enumerate(range(min(n - 1, i + window), i + 1, -1)):
            if tries >= max_tries:
                break
            leg = grid.cost(*grid.line(plats[i], plons[i], plats[k], plons[k]))
            if leg <= step_cum[k] - step_cum[i] + 1e-6:
                best = k
                break
        out.append(best)
        i = best
    return out


def path_stats(lat, lon, z, sl, theta_ref):
    seg = hav_vec(lat[:-1], lon[:-1], lat[1:], lon[1:])
    dz = np.diff(z)
    theta = np.degrees(np.arctan2(np.abs(dz), np.maximum(seg, 1e-6)))
    cost = float(np.sum(seg * (1 + (theta / theta_ref) ** 2 + 0.5 * (np.asarray(sl)[1:] / theta_ref) ** 2)))
    total = float(seg.sum())
    return {
        "distance_km": round(total / 1000.0, 1),
        "max_slope_deg": round(float(theta.max()), 2) if len(theta) else 0.0,
        "mean_slope_deg": round(float(np.sum(theta * seg) / total), 2) if total else 0.0,
        "elevation_gain_m": round(float(dz[dz > 0].sum())),
        "elevation_loss_m": round(float(-dz[dz < 0].sum())),
        "min_elevation_m": round(float(z.min())),
        "max_elevation_m": round(float(z.max())),
        "cost_index_km": round(cost / 1000.0, 1),
    }


def profile(lat, lon, z, points=60):
    seg = hav_vec(lat[:-1], lon[:-1], lat[1:], lon[1:])
    d = np.concatenate([[0.0], np.cumsum(seg)]) / 1000.0
    pick = np.unique(np.linspace(0, len(d) - 1, points).astype(int))
    return [[round(float(d[k]), 1), round(float(z[k]))] for k in pick]


def plan_route(dem, a, b, avoid_steep=True, max_dim=300):
    (lat1, lon1), (lat2, lon2) = a, b
    lon1, lon2 = wrap_lon(lon1), wrap_lon(lon2)
    if haversine_m(lat1, lon1, lat2, lon2) < 1000:
        raise DemError("Start and destination are the same place.")
    prof = PROFILES[bool(avoid_steep)]
    tref, tmax = prof["theta_ref_deg"], prof["theta_max_deg"]

    gc_lat, gc_lon = great_circle_points(lat1, lon1, lat2, lon2, 400)
    span = max(gc_lat.max() - gc_lat.min(), gc_lon.max() - gc_lon.min())
    if span > 175:
        raise DemError("Route is too long for one request (over 175 degrees).")
    pad = min(max(2.5, 0.25 * span), max(0.5, (350.0 - span) / 2.0))
    blk = dem.read_block(gc_lat.min() - pad, gc_lat.max() + pad, gc_lon.min() - pad, gc_lon.max() + pad, max_dim)
    grid = Grid(blk, tref, tmax)
    z, lats, lons, ny, nx = grid.z, grid.lats, grid.lons, grid.ny, grid.nx
    dx_row = (R * np.maximum(0.02, np.cos(np.radians(lats))) * math.radians(grid.cell)).tolist()

    start, goal = grid.cell_of(lat1, gc_lon[0]), grid.cell_of(lat2, gc_lon[-1])
    if not (np.isfinite(z[start]) and np.isfinite(z[goal])):
        raise DemError("Start or destination cell has no DEM data.")

    cells = _astar(z.ravel().tolist(), grid.sloc.ravel().tolist(), dx_row, grid.dy, (ny, nx), start, goal, tref, tmax)
    pj, pi = np.array([c[0] for c in cells]), np.array([c[1] for c in cells])
    plats, plons = lats[pj].copy(), lons[pi].copy()
    plats[0], plons[0], plats[-1], plons[-1] = lat1, gc_lon[0], lat2, gc_lon[-1]
    # cumulative A* step cost (same formula as the search) for the shortcut test
    seg = hav_vec(plats[:-1], plons[:-1], plats[1:], plons[1:])
    th = np.degrees(np.arctan2(np.abs(np.diff(z[pj, pi])), np.maximum(seg, 1e-6)))
    step = seg * (1 + (th / tref) ** 2 + 0.5 * (grid.sloc[pj, pi][1:] / tref) ** 2)
    cum = np.concatenate([[0.0], np.cumsum(step)])

    anchors = smooth_path(grid, plats, plons, cum)
    legs = [grid.line(plats[anchors[k]], plons[anchors[k]], plats[anchors[k + 1]], plons[anchors[k + 1]])
            for k in range(len(anchors) - 1)]
    dense = [np.concatenate(parts) for parts in zip(*legs)]  # lat, lon, z, sl

    # direct great-circle line, sampled at the same spacing and judged by the same cost model
    n_line = int(min(3000, max(2, haversine_m(lat1, lon1, lat2, lon2) / grid.dy)))
    sl_lat, sl_lon = great_circle_points(lat1, lon1, lat2, lon2, n_line)
    sz, ss = grid.sample(sl_lat, sl_lon)
    straight_cost = grid.cost(sl_lat, sl_lon, sz, ss)
    ok = np.isfinite(sz)
    straight = path_stats(sl_lat[ok], sl_lon[ok], sz[ok], ss[ok], tref)
    straight["blocked_by_terrain"] = not math.isfinite(straight_cost)
    terrain_cost = path_stats(*dense, tref)["cost_index_km"]
    if math.isfinite(straight_cost) and straight["cost_index_km"] <= terrain_cost:
        # terrain is benign: the direct line is already the cheapest path under this model
        dense, anchors_ll = (sl_lat[ok], sl_lon[ok], sz[ok], ss[ok]), [(lat1, gc_lon[0]), (lat2, gc_lon[-1])]
        mode = "direct line is already optimal (benign terrain)"
    else:
        anchors_ll = [(float(plats[k]), float(plons[k])) for k in anchors]
        mode = "detour around costly terrain" if terrain_cost < straight["cost_index_km"] else "least-cost path"
    terrain = path_stats(*dense, tref)
    pick = np.unique(np.linspace(0, len(dense[0]) - 1, min(200, len(dense[0]))).astype(int))
    path = [[round(float(dense[0][k]), 4), round(float(wrap_lon(dense[1][k])), 4)] for k in pick]

    waypoints = []
    for lat, lon in anchors_ll:
        zz, _ = grid.sample(np.array([lat]), np.array([lon]))
        waypoints.append({"latitude": round(float(lat), 4), "longitude": round(float(wrap_lon(lon)), 4),
                          "elevation_m": round(float(zz[0])) if np.isfinite(zz[0]) else None})
    saving = 100.0 * (1.0 - terrain["cost_index_km"] / straight["cost_index_km"]) if straight["cost_index_km"] else 0.0
    return {
        "method": "A* least-cost search on block-averaged NASA MOLA DEM, then any-angle smoothing",
        "result_type": mode,
        "parameters": {"profile": prof["label"], "theta_ref_deg": tref, "theta_max_deg": tmax,
                       "step_cost": "length * (1 + (slope_along/theta_ref)^2 + 0.5*(slope_cell/theta_ref)^2)"},
        "grid": {"cell_size_m": round(grid.dy), "cells": [ny, nx], "decimation_factor": blk["factor"]},
        "waypoints": waypoints,
        "path": path,
        "terrain_route": terrain,
        "straight_line": straight,
        "profile": {"terrain_route": profile(dense[0], dense[1], dense[2]), "straight_line": profile(sl_lat[ok], sl_lon[ok], sz[ok])},
        "improvement": {"cost_reduction_pct": round(saving, 1),
                        "max_slope_change_deg": round(terrain["max_slope_deg"] - straight["max_slope_deg"], 2)},
        "notes": [
            "Slopes are grid-scale averages over several kilometres; local slopes, rocks and sand are not captured.",
            "No rock-abundance or thermal dataset is used: 'Minimize rocky zones' does not change this path.",
            "Heights are block averages of the 463 m DEM, so peaks and pits appear smoother than in the raw data.",
            "Cost units are a relative planning index, not rover energy or time.",
        ],
    }
