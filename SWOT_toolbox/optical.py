"""Independent wet/dry evidence from user-classified georeferenced imagery."""
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio
from rasterio.warp import transform


def load_masks(spec, base):
    """Prepare a small time index. Mask values are 1 water, 0 dry; others unknown."""
    table = pd.read_csv(base / spec["index"])
    times = pd.to_datetime(table["date"], utc=True, format="mixed")
    if times.duplicated().any() or times.isna().any():
        raise ValueError("Provide one optical mask per UTC timestamp")
    tolerance = pd.Timedelta(spec["tolerance"])
    if pd.isna(tolerance) or tolerance < pd.Timedelta(0):
        raise ValueError("Invalid optical matching tolerance")
    folder = (base / spec["index"]).parent
    return sorted(zip(times, (folder / p for p in table.file))), tolerance


def sample_mask(index, date, lons, lats):
    """Sample cell centres and the station in one operation."""
    images, tolerance = index
    values = np.full(len(lons), np.nan)
    if not images:
        return None, values
    time = pd.to_datetime(date, utc=True)
    observed, filename = min(images, key=lambda pair: abs(pair[0] - time))
    if abs(observed - time) > tolerance:
        return None, values
    with rasterio.open(filename) as src:
        if src.crs is None or src.nodata in (0, 1):
            raise ValueError(f"Mask needs a CRS and distinct water/dry/nodata codes: {filename}")
        good = np.flatnonzero(np.isfinite(lons) & np.isfinite(lats))
        xs, ys = transform("EPSG:4326", src.crs, lons[good], lats[good])
        for j, x, y, v in zip(good, xs, ys, src.sample(zip(xs, ys), indexes=1, masked=True)):
            row, col = src.index(x, y)
            if 0 <= row < src.height and 0 <= col < src.width and not np.ma.is_masked(v[0]) and v[0] in (0, 1):
                values[j] = float(v[0])
    return observed, values
