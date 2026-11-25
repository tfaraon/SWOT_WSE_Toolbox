"""Extract one station from SWOT L2 HR Raster granules, retaining a QC audit."""
from concurrent.futures import ProcessPoolExecutor
from functools import partial
from pathlib import Path
import re
import json
import warnings
from contextlib import contextmanager

from netCDF4 import Dataset
import numpy as np
import pandas as pd
from pyproj import CRS, Transformer, Geod
from .optical import sample_mask

PRODUCT = re.compile(r"SWOT_L2_HR_Raster_(\d+)m_.*?_(\d{3})_(\d{3})_([^_]+)_(\d{8}T\d{6})_")
GROUP = ["pass", "resolution", "tile"]


# Match only the terminal processing suffix. All footprint fields and both
# acquisition times remain in the key.
REVISION = re.compile(r"_(P[A-Z]{2})([0-9]+)_([0-9]+)$")
COPY = re.compile(r"(?: \(\d+\)| copy(?: \d+)?)$", re.IGNORECASE)


def select_acquisitions(items, strict=False):
    """Select one local granule per full product identity, before any QC."""
    groups = {}
    for item in items:
        path = item[0]
        stem = COPY.sub("", path.stem)
        rev = REVISION.search(stem)
        key = stem[:rev.start()] if rev else stem
        rank = (rev[1], int(rev[2]), int(rev[3])) if rev else ("", 0, 0)
        groups.setdefault(key, []).append((rank, stem == path.stem, item))
    selected, ignored = [], {}
    for key, group in groups.items():
        if strict and len(group) > 1:
            raise ValueError(f"Duplicate acquisition {key}: " + ", ".join(str(g[2][0]) for g in group))
        # items arrive in sorted path order, so equal ranks have stable ties.
        ranked = sorted(group, key=lambda entry: entry[:2], reverse=True)
        winner = ranked[0][2]
        selected.append(winner)
        ignored[str(winner[0])] = [str(entry[2][0]) for entry in ranked[1:]]
    return selected, ignored


class GranuleReadError(OSError):
    """I/O failure while opening, reading or closing a SWOT NetCDF."""


@contextmanager
def open_granule(path):
    try:
        ds = Dataset(path)
    except (OSError, RuntimeError) as exc:
        raise GranuleReadError(str(exc)) from exc
    try:
        yield ds
    finally:
        try:
            ds.close()
        except (OSError, RuntimeError) as exc:
            raise GranuleReadError(str(exc)) from exc


def numbers(variable, key):
    try:
        data = variable[key]
    except (OSError, RuntimeError) as exc:
        raise GranuleReadError(str(exc)) from exc
    return np.ma.asarray(data, dtype=float).filled(np.nan)


def extract_candidates(task, **options):
    preferred, alternatives = task
    paths = [preferred[0], *map(Path, alternatives)]
    failures = []
    for path in paths:
        try:
            row = extract_granule((path, *preferred[1:]), **options)
            break
        except GranuleReadError as exc:
            failures.append(dict(filename=str(path), error=str(exc)))
    else:
        path, res, cycle, orbit, tile, stamp = preferred
        row = dict(filename=str(path), resolution=res, cycle=cycle, tile=tile,
                   date=pd.to_datetime(stamp, format="%Y%m%dT%H%M%S"), status="unreadable",
                   wse_raw=np.nan, **{"pass": orbit})
    row["duplicate_count"] = len(alternatives)
    row["skipped_files"] = json.dumps([str(p) for p in paths if str(p) != row["filename"]], ensure_ascii=False)
    row["read_errors"] = json.dumps(failures, ensure_ascii=False)
    return row


def centre(ds, lon, lat, spacing):
    """Use projected axes when available, including cells with missing WSE."""
    vs = ds.variables
    if "x" in vs and "y" in vs:
        crs = vs[getattr(vs["wse"], "grid_mapping", "crs")]
        attrs = {key: crs.getncattr(key) for key in crs.ncattrs()}
        projection = CRS.from_wkt(attrs["spatial_ref"]) if "spatial_ref" in attrs else CRS.from_cf(attrs)
        x, y = Transformer.from_crs(4326, projection, always_xy=True).transform(lon, lat)
        xs, ys = numbers(vs["x"], slice(None)), numbers(vs["y"], slice(None))
        if vs["wse"].dimensions != (vs["y"].dimensions[0], vs["x"].dimensions[0]):
            raise ValueError("Expected WSE indexed by y, x")
        col, row = np.abs(xs-x).argmin(), np.abs(ys-y).argmin()
        if max(abs(xs[col]-x), abs(ys[row]-y)) > spacing * 0.50001:
            return None
        return row, col, xs, ys, projection

    lons = numbers(vs["longitude"], slice(None))
    lats = numbers(vs["latitude"], slice(None))
    if lons.ndim == lats.ndim == 1:
        lons, lats = np.meshgrid(lons, lats)
    if lons.shape != lats.shape or lons.shape != vs["wse"].shape:
        raise ValueError("Incompatible geographic coordinates")
    dx = ((lons-lon+180) % 360-180) * np.cos(np.deg2rad(lat))
    dist = 111320 * np.hypot(dx, lats-lat)
    if not np.isfinite(dist).any() or np.nanmin(dist) > spacing * 0.72:
        return None
    return (*np.unravel_index(np.nanargmin(dist), dist.shape), lons, lats, None)


def extract_granule(item, *, lon, lat, width_m, radius, levels, quality_mode, min_fraction, estimator, optical, detail_dates):
    path, res, cycle, orbit, tile, stamp = item
    spacing = int(res)
    # Nearest odd number of pixels, round ties upwards.
    r = radius if width_m is None else int(np.floor(width_m / (2 * spacing)))
    side = 2*r + 1
    row = dict(filename=str(path), resolution=res, cycle=cycle, tile=tile,
               date=pd.to_datetime(stamp, format="%Y%m%dT%H%M%S"),
               window_px=side, width_m=side*spacing, n_valid=0, n_available=0,
               valid_fraction=0., wse_raw=np.nan, status="outside")
    row["pass"] = orbit
    with open_granule(path) as ds:
        loc = centre(ds, lon, lat, spacing)
        if loc is None:
            return row
        y, x, xs, ys, projection = loc
        wse, qual = ds["wse"], ds["wse_qual"]
        if qual.dimensions != wse.dimensions:
            raise ValueError(f"Inconsistent quality grid in {path}")
        box = (slice(max(0, y-r), y+r+1), slice(max(0, x-r), x+r+1))
        values = numbers(wse, box)
        valid = np.isfinite(values)
        row["n_available"] = int(valid.sum())
        if levels is not None:
            q = numbers(qual, box if quality_mode == "pixel" else (y, x))
            valid &= np.isin(q, levels)
        n = int(valid.sum())
        row["n_in_tile"] = values.size
        row["tile_fraction"] = values.size / side**2
        row["available_fraction_in_tile"] = row["n_available"] / values.size
        row["quality_retained_fraction"] = n / row["n_available"] if row["n_available"] else np.nan
        row["centre_valid"] = bool(valid[y-box[0].start, x-box[1].start])
        if projection is None:
            lons, lats = xs[box], ys[box]
        else:
            xx, yy = np.meshgrid(xs[box[1]], ys[box[0]])
            lons, lats = Transformer.from_crs(projection, 4326, always_xy=True).transform(xx, yy)
        # WGS84 geodesic offsets from the exact station, including subpixel shifts.
        az, _, distance = Geod(ellps="WGS84").inv(np.full(values.shape, float(lon)),
            np.full(values.shape, float(lat)), lons, lats)
        east = distance*np.sin(np.deg2rad(az))
        north = distance*np.cos(np.deg2rad(az))
        water = np.full(values.shape, np.nan)
        if optical is not None:
            obs, samples = sample_mask(optical, row["date"], np.r_[lon, lons.ravel()], np.r_[lat, lats.ravel()])
            if obs is not None:
                water = samples[1:].reshape(values.shape)
                row.update(optical_time=obs, optical_lag_h=(obs-row["date"].tz_localize("UTC")).total_seconds()/3600,
                           optical_centre_water=samples[0], optical_known_cells=int(np.isfinite(water).sum()),
                           optical_wet_fraction=np.nanmean(water) if np.isfinite(water).any() else np.nan)
        if row["date"].strftime("%Y-%m-%d") in detail_dates:
            row["_pixels"] = pd.DataFrame(dict(east_m=east.ravel(), north_m=north.ravel(),
                lon=lons.ravel(), lat=lats.ravel(), used=valid.ravel(), optical_water=water.ravel(),
                pixel_wse=values.ravel()))
        row.update(n_valid=n, valid_fraction=n/side**2)
        if not row["n_available"]:
            row["status"] = "no_wse"
        elif not n:
            row["status"] = "quality"
        else:
            retained = values[valid]
            median = float(np.median(retained))
            row.update(wse_mean_m=float(retained.mean()), wse_median_m=median,
                       wse_mad_m=float(np.median(np.abs(retained-median))),
                       nearest_valid_m=float(distance[valid].min()),
                       centroid_distance_m=float(np.hypot(east[valid].mean(), north[valid].mean())))
            row["wse_raw"] = median if estimator == "median" else row["wse_mean_m"]
            row["status"] = "accepted" if row["valid_fraction"] >= min_fraction else "coverage"
    return row


def station(directory, lon, lat, *, width_m=None, radius=0, levels=(0, 1),
            quality_mode="pixel", min_fraction=0., bounds=None, iqr="flag",
            iqr_factor=1.5, resolution=None, orbit=None, tile=None, workers=1,
            estimator="mean", optical=None, detail_dates=(), historical_pool=False, iqr_fit_period=None, duplicate_policy="preferred"):
    """Return a row per input granule. Dates are UTC granule start times."""
    if estimator not in ("mean", "median"):
        raise ValueError("estimator must be mean or median")
    if not Path(directory).is_dir():
        raise ValueError(f"Missing granule directory: {directory}")
    if not (-180 <= lon <= 180 and -90 <= lat <= 90):
        raise ValueError("Invalid station coordinates")
    if width_m is not None and (not np.isfinite(width_m) or width_m <= 0 or radius):
        raise ValueError("Supply a positive width_m or a pixel radius")
    if not isinstance(radius, (int, np.integer)) or radius < 0 or workers < 1:
        raise ValueError("Invalid radius or worker count")
    if quality_mode not in ("centre", "pixel") or iqr not in ("off", "flag", "reject"):
        raise ValueError("Invalid quality_mode or iqr action")
    if not 0 <= min_fraction <= 1 or not np.isfinite(iqr_factor) or iqr_factor <= 0:
        raise ValueError("Invalid coverage or IQR factor")
    if bounds is not None and (len(bounds) != 2 or not np.isfinite(bounds).all() or bounds[0] >= bounds[1]):
        raise ValueError("Expected ordered lower and upper WSE bounds in metres")
    if levels is not None:
        levels = np.atleast_1d(levels)
        if not len(levels) or not np.isin(levels, [0, 1, 2, 3]).all():
            raise ValueError("Quality classes must be drawn from 0, 1, 2, 3")
    if duplicate_policy not in ("preferred", "error"):
        raise ValueError("duplicate_policy must be preferred or error")
    items = []
    for path in sorted(Path(directory).rglob("SWOT_L2_HR_Raster_*.nc")):
        match = PRODUCT.match(path.name)
        if not match:
            raise ValueError(f"Unrecognised Raster granule name: {path.name}")
        res, cycle, p, t, stamp = match.groups()
        if resolution is not None and res != str(resolution):
            continue
        if orbit is not None and p != str(orbit).zfill(3) or tile is not None and t != str(tile):
            continue
        items.append((path, *match.groups()))
    items, ignored = select_acquisitions(items, strict=duplicate_policy == "error")
    columns = ["date", "filename", "cycle", *GROUP, "window_px", "width_m",
               "n_valid", "n_available", "valid_fraction", "wse_raw", "status"]
    tasks = [(item, ignored[str(item[0])]) for item in items]
    extract = partial(extract_candidates, lon=lon, lat=lat, width_m=width_m, radius=radius,
                      levels=levels, quality_mode=quality_mode, min_fraction=min_fraction, estimator=estimator,
                      optical=optical, detail_dates=detail_dates)
    if workers == 1 or len(items) < 2:
        rows = list(map(extract, tasks))
    else:
        with ProcessPoolExecutor(max_workers=min(workers, len(items))) as pool:
            rows = list(pool.map(extract, tasks))
    details = []
    for row in rows:
        pixels = row.pop("_pixels", None)
        if pixels is not None:
            details.append(pixels.assign(filename=row["filename"], date=row["date"], cycle=row["cycle"],
                                         resolution=row["resolution"], tile=row["tile"], **{"pass":row["pass"]}))
    columns += ["read_errors", "duplicate_count", "skipped_files", "n_in_tile", "tile_fraction", "available_fraction_in_tile", "quality_retained_fraction",
                "centre_valid", "nearest_valid_m", "centroid_distance_m", "wse_mean_m", "wse_median_m", "wse_mad_m"]
    if optical is not None:
        columns += ["optical_time","optical_lag_h","optical_centre_water","optical_known_cells","optical_wet_fraction"]
    df = pd.DataFrame(rows, columns=columns).sort_values(["date", "filename"]).reset_index(drop=True)
    failed_reads = df[df.read_errors.ne("[]")]
    if not failed_reads.empty:
        first = json.loads(failed_reads.iloc[0].read_errors)[0]
        warnings.warn(f"NetCDF read errors for {len(failed_reads)} acquisition(s); "
                      f"{df.status.eq('unreadable').sum()} without a readable copy. "
                      f"First failed file: {first['filename']}. See read_errors for details.",
                      RuntimeWarning, stacklevel=2)
    df["iqr_flag"] = False
    if bounds is not None and (not historical_pool or df.status.eq("accepted").sum() >= 4):
        outside = ~df.wse_raw.between(*bounds)
        df.loc[df.status.eq("accepted") & outside, "status"] = "bounds"
    if iqr != "off":
        accepted = df[df.status.eq("accepted")]
        groups = [(None, accepted)] if historical_pool else accepted.groupby(GROUP)
        for _, g in groups:
            fit = g
            if iqr_fit_period is not None:
                start, end = pd.to_datetime(iqr_fit_period, utc=True)
                if not start < end:
                    raise ValueError("Invalid IQR fitting period")
                utc = pd.to_datetime(g.date, utc=True)
                fit = g[(utc >= start) & (utc < end)]
            if len(fit) < 4:
                continue
            lo, hi = fit.wse_raw.quantile([.25, .75])
            margin = iqr_factor * (hi-lo)
            flagged = g.index[~g.wse_raw.between(lo-margin, hi+margin)]
            df.loc[flagged, "iqr_flag"] = True
            if iqr == "reject":
                df.loc[flagged, "status"] = "iqr"
    df["wse"] = df.wse_raw.where(df.status.eq("accepted"))
    if details:
        df.attrs["pixel_details"] = pd.concat(details, ignore_index=True)
    return df
