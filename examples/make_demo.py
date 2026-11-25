"""Create a small synthetic flood sequence in a new directory."""
from datetime import datetime, timedelta
import json
from pathlib import Path
import sys
from netCDF4 import Dataset
import numpy as np
import pandas as pd
import rasterio
from rasterio.transform import from_origin
from pyproj import CRS, Transformer


def create(folder):
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=False)
    data = folder / "rasters"
    data.mkdir()
    lon, lat = 137.028098, -28.893022
    crs = CRS.from_epsg(32753)
    x0, y0 = Transformer.from_crs(4326, crs, always_xy=True).transform(lon, lat)
    gauge, lake, masks = [], [], []
    for i, level in enumerate([-14., -14., -14., -11., -14., -14., -14.]):
        date = datetime(2025, 1, 1) + timedelta(days=i*21)
        stamp = date.strftime("%Y%m%dT%H%M%S")
        gauge.append({"date": date.isoformat(), "level_m": level - .2})
        lake.append({"time_str": date.isoformat(), "lake_id": "99", "wse": level + .1,
                     "quality_f": 0, "cycle": i+1, "pass": 394})
        if i in (0, 3, 4):
            mask = np.ones((15, 15), dtype="uint8")
            if i == 4:
                mask[:, :9] = 0
            filename = f"optical_{stamp}.tif"
            with rasterio.open(folder / filename, "w", driver="GTiff", height=15, width=15,
                               count=1, dtype="uint8", crs=crs, nodata=255,
                               transform=from_origin(x0-750, y0+750, 100, 100)) as src:
                src.write(mask, 1)
            masks.append(dict(date=(date + timedelta(hours=2)).isoformat(), file=filename))
        for res in (100, 250):
            filename = f"SWOT_L2_HR_Raster_{res}m_UTM53J_demo_{i+1:03d}_394_102F_{stamp}_{stamp}_demo.nc"
            with Dataset(data / filename, "w") as ds:
                for name, origin in (("x", x0), ("y", y0)):
                    ds.createDimension(name, 15)
                    ds.createVariable(name, "f8", (name,))[:] = origin + np.arange(-7,8)*res
                ds.createVariable("crs", "i4").setncatts(crs.to_cf())
                dims = ("y","x")
                w = ds.createVariable("wse", "f4", dims, fill_value=-9999.)
                q = ds.createVariable("wse_qual", "u1", dims, fill_value=255)
                w.grid_mapping = "crs"
                w[:] = level
                q[:] = 0
                w[7,8] = level + 9
                q[7,8] = 3
                if i == 4:
                    w[:, :8 if res == 250 else 9] = w._FillValue
                    q[7,7] = 3
                if i == 6:
                    w[:] = w._FillValue
                # Geographic coordinates allow the original notebook API to run too.
                x, y = np.meshgrid(ds["x"][:], ds["y"][:])
                lons, lats = Transformer.from_crs(crs, 4326, always_xy=True).transform(x, y)
                ds.createVariable("longitude", "f8", dims)[:] = lons
                ds.createVariable("latitude", "f8", dims)[:] = lats
                ds.setncatts(dict(geospatial_lon_min=lons.min(), geospatial_lon_max=lons.max(),
                                  geospatial_lat_min=lats.min(), geospatial_lat_max=lats.max()))
    pd.DataFrame(gauge).to_csv(folder / "gauge.csv", index=False)
    pd.DataFrame(lake).to_csv(folder / "lakesp.csv", index=False)
    pd.DataFrame(masks).to_csv(folder / "optical.csv", index=False)
    refs = [dict(kind="gauge", file="gauge.csv", tolerance="1h", datum="demo"),
            dict(kind="lakesp", file="lakesp.csv", lake_id="99", tolerance="30min",
                 datum="demo", allow={"quality_f":[0]})]
    cfg = dict(granules="rasters", output="results", raster_datum="demo",
               sites=[dict(name="Synthetic station", lon=lon, lat=lat, references=refs,
                           detail_dates=["2025-01-01","2025-03-05","2025-03-26"],
                           optical=dict(index="optical.csv", tolerance="1D"))],
               methods=[dict(name="point", radius=0), dict(name="centre", width_m=1300, quality_mode="centre"),
                        dict(name="pixel", width_m=1300), dict(name="reject", width_m=1300, iqr="reject")])
    (folder / "config.json").write_text(json.dumps(cfg, indent=2))


if __name__ == "__main__":
    create(sys.argv[1])
