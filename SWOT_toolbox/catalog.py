"""Discover SWOT Raster granules at a site and audit selected local observations."""
from datetime import date
import hashlib
import json
from pathlib import Path
from tempfile import TemporaryDirectory
from urllib.parse import urlencode, urlsplit
from urllib.request import Request, urlopen

import pandas as pd

from .stations import PRODUCT, station

CMR = "https://cmr.earthdata.nasa.gov/search/granules.umm_json"
COLUMNS = ["date_utc", "filename", "collection", "cycle", "pass", "tile", "resolution",
           "size_bytes", "algorithm", "checksum", "url", "local_path", "local_state"]


def local_state(row):
    path = Path(row["local_path"])
    if not path.is_file():
        return "missing"
    if path.stat().st_size != row["size_bytes"]:
        return "invalid_size"
    if row["checksum"]:
        digest = hashlib.new(row["algorithm"].lower().replace("-", ""))
        with path.open("rb") as stream:
            for block in iter(lambda: stream.read(1024 * 1024), b""):
                digest.update(block)
        if digest.hexdigest().lower() != row["checksum"].lower():
            return "invalid_checksum"
    return "available"


def search_site(lon, lat, start, end, directory, resolution=100, collection="SWOT_L2_HR_Raster_D"):
    """Return one row per matching catalogue file. Both dates are inclusive (UTC)."""
    if not (-180 <= lon <= 180 and -90 <= lat <= 90 and resolution in (100, 250)):
        raise ValueError("Check coordinates and choose 100 or 250 m.")
    if date.fromisoformat(start) > date.fromisoformat(end):
        raise ValueError("The start date must precede or equal the end date.")
    query = urlencode({"short_name": collection, "point": f"{lon},{lat}",
                       "temporal": f"{start}T00:00:00Z,{end}T23:59:59.999Z", "page_size": 500})
    directory = Path(directory).expanduser().resolve()
    rows = {}
    headers = {}
    while True:
        with urlopen(Request(f"{CMR}?{query}", headers=headers), timeout=90) as response:
            items = json.load(response)["items"]
            cursor = response.headers.get("CMR-Search-After")
        for item in items:
            u = item["umm"]
            for info in u["DataGranule"]["ArchiveAndDistributionInformation"]:
                name = info["Name"]
                match = PRODUCT.match(name)
                if not (name.endswith(".nc") and Path(name).name == name and match
                        and int(match[1]) == resolution):
                    continue
                url = next(r["URL"] for r in u["RelatedUrls"] if r["Type"] == "GET DATA"
                           and r["URL"].startswith("https://") and urlsplit(r["URL"]).path.endswith("/" + name))
                checksum = info.get("Checksum", {})
                row = dict(date_utc=u["TemporalExtent"]["RangeDateTime"]["BeginningDateTime"],
                           filename=name, collection=collection, cycle=match[2], tile=match[4],
                           resolution=resolution, size_bytes=int(info["SizeInBytes"]), url=url,
                           checksum=checksum.get("Value", ""), algorithm=checksum.get("Algorithm", ""),
                           local_path=str(directory / name), **{"pass": match[3]})
                row["local_state"] = local_state(row)
                rows[name] = row
        if len(items) < 500:
            break
        if not cursor or cursor == headers.get("CMR-Search-After"):
            raise RuntimeError("CMR pagination interrupted. Repeat the search.")
        headers["CMR-Search-After"] = cursor
    return pd.DataFrame(rows.values(), columns=COLUMNS).sort_values(["date_utc", "filename"]).reset_index(drop=True)


def download_selected(selection):
    """Reuse verified files and install completed downloads atomically; retry by rerunning."""
    import earthaccess

    report = selection.copy()
    report["download_error"] = ""
    authenticated = False
    for index, row in report.iterrows():
        path = Path(row.local_path)
        try:
            if local_state(row) != "available":
                if not authenticated:
                    authenticated = earthaccess.login().authenticated
                    if not authenticated:
                        raise PermissionError("Earthdata authentication required.")
                path.parent.mkdir(parents=True, exist_ok=True)
                with TemporaryDirectory(dir=path.parent) as temp:
                    earthaccess.download([row.url], temp, provider="POCLOUD")
                    staged = Path(temp) / row.filename
                    if local_state({**row.to_dict(), "local_path": staged}) != "available":
                        raise OSError("Downloaded file missing or failed integrity check")
                    staged.replace(path)
            report.at[index, "local_state"] = "available"
        except PermissionError:
            raise
        except Exception as exc:
            report.at[index, "local_state"] = "download_error"
            report.at[index, "download_error"] = str(exc)
    return report


def extract_selection(selection, lon, lat, directory, **options):
    """Run station QC on available selected files and retain every catalogue row."""
    report = selection.copy()
    report["local_state"] = [local_state(row) for _, row in report.iterrows()]
    paths = report.loc[report.local_state.eq("available"), "local_path"].tolist()
    if paths:
        audit = station(directory, lon, lat, filenames=paths, **options)
        audit = audit.rename(columns={"filename": "local_path"})
        fields = ["local_path", "status", "n_valid", "valid_fraction", "wse", "iqr_flag"]
        report = report.merge(audit[fields], how="left", on="local_path", validate="one_to_one")
        report.loc[report.status.isna() & report.local_state.eq("available"), "status"] = "not_selected_by_station"
    else:
        for field in ["status", "n_valid", "valid_fraction", "wse", "iqr_flag"]:
            report[field] = None
    report["status"] = report.status.fillna("not_processed")
    return report
