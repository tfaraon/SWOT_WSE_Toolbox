# SWOT WSE Toolbox

At Kati Thanda–Lake Eyre North, floodwater reaches different parts of the lake at different times. A water level measured at Belt Bay therefore answers a local question. I developed this toolbox to follow water levels at fixed locations through filling and drying, including at sites with few field measurements.

The toolbox extracts water surface elevation (WSE, metres) from the SWOT Level 2 High Rate Raster product. Each station is defined by its longitude and latitude. The output records the contributing pixels and the observations excluded by each filter. Belt Bay and Madigan Gulf are included as examples. Their application in the observatory is described in the [KTLE Observatory repository](https://github.com/tfaraon/ktle-observatory).

## Search and download observations at a site

Start from the repository root:

```bash
python -m pip install -e '.[download]' jupyterlab matplotlib
jupyter lab SWOT_site_download.ipynb
```

Enter the site coordinates and an inclusive UTC date range. The notebook searches NASA CMR and lists the Raster granules covering the point. Select files by row index and set `download = True` to retrieve them with an Earthdata account. The default example uses the 100 m Version D product at Belt Bay in January 2026. Change the directory to reuse an existing flat NetCDF archive. Files are checked against their catalogue size and, when provided, checksum. Rerun with an extended period to update the series. Each selected granule has its own row in the output, including overlapping tiles from the same pass. A catalogue match establishes footprint coverage. The extraction supplies the retained pixel count, quality status and local elevation. Results are written under `output/` as CSV tables and, when accepted observations exist, a PDF plot.

## Installation and first run

From the repository root, create a Python environment and install the package. The small example creates two synthetic NetCDF time series at 100 and 250 m, a gauge record and a LakeSP-style CSV. It runs offline and takes a few seconds. Use a new directory for each generated example.

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install .
python examples/make_demo.py example-run
swot-stations example-run/config.json
```

Open `example-run/results/summary.csv`. For the 100 m raster and the method named `pixel`, six of seven overpasses have a usable level. The expected mean difference is +0.20 m against the synthetic gauge and −0.10 m against the synthetic LakeSP series. In `audit.csv`, compare the flood peak under the methods named `pixel` and `reject`. Both flag the peak. Only the first retains it in the series.

## Working with SWOT data

Place the original Raster NetCDF granules under `data/rasters/`, then run `swot-stations examples/ktle.json`. The input directory is searched recursively. Copies and processing revisions are resolved before filtering. The full filename up to the processing suffix identifies the granule, including the projection, footprint and both acquisition timestamps. Following the Observatory convention, processing families are ranked lexically (e.g. PID after PIC), then by numeric processing and file counters. This is a selection convention for local files. To use a single scientific collection throughout, point the tool at that collection. For ties, the original filename is preferred to a numbered copy; the first path in sorted order wins among identical names. The `filename`, `duplicate_count` and `skipped_files` columns record the decision. To stop on competing candidates instead, set `duplicate_policy` to `error`. The configuration uses coordinates in decimal degrees (WGS 84). File paths are relative to the JSON file.

If a NetCDF cannot be opened or its pixel values cannot be read, extraction tries the next local copy or revision of that acquisition. Each attempt is recorded in the audit column `read_errors`; the selected filename identifies the version actually used. An acquisition with no readable candidate receives status `unreadable` and is excluded from coverage statistics because its station coverage cannot be established. The notebook wrappers return accepted levels and attach the list of failed files and messages to `result.attrs["read_errors"]`. A readable acquisition rejected by scientific quality checks is never replaced with a lower-ranking revision.

The optional downloader uses an Earthdata account. Install it with `python -m pip install '.[download]'`. For example, the following request covers the southern lake during June 2025.

```bash
python -m SWOT_toolbox.SWOT_Downloader data/rasters 2025-06-01 2025-06-30 --bbox 136.9 -29.1 137.7 -28.7
```

### Choosing a window and filtering observations

| Setting | Effect |
| --- | --- |
| `radius` or `width_m` | A radius of 0 extracts one pixel. A radius of 6 extracts 13 × 13 pixels. A width in metres selects the nearest odd number of pixels, rounding ties upwards. A request for 1,300 m gives 1,300 m at 100 m spacing and 1,250 m at 250 m spacing. |
| `quality_mode` and `levels` | With `pixel`, the mean uses finite WSE pixels in the selected quality classes (default 0 and 1). With `centre`, the central quality flag determines whether all finite WSE pixels in the window may contribute. |
| `min_fraction` | Minimum fraction of the requested square occupied by usable pixels. The denominator includes the part of a window outside the tile. The default is 0, and every estimate requires at least one usable pixel. |
| `bounds_m` | Optional lower and upper levels for each site, in the raster vertical datum. The example uses −16 to +6 m for KTLE. Omit this field for a new site until suitable bounds are available. |
| `iqr` | Choose `off`, `flag` (default) or `reject`. For each pass, tile and resolution with at least four retained observations, values beyond 1.5 interquartile ranges from the quartiles receive a flag. Set `iqr_factor` to change the multiplier. |
| `resolution`, `orbit`, `tile` | Optional selectors for raster spacing in metres, pass number and tile identifier. Each resolution is summarised separately. |

At a retreating shoreline, a larger window may retain water pixels after the central pixel becomes dry. Inspect `n_valid`, `valid_fraction` and the accepted levels together. A granule without a usable WSE can also reflect missing coverage or poor retrieval conditions. Establish the wet or dry state with an independent observation where possible.

## Comparing with a gauge or LakeSP

Add a `references` list to the relevant site in the JSON file. Gauge CSV files use `date,level_m`. Use metres and UTC times or timestamps with explicit offsets. The Raster match time is the granule start. Set the time tolerance to suit the reference sampling interval and the rate of water-level change.

```json
"references": [
  {"name": "Station gauge", "kind": "gauge", "file": "../data/gauge.csv",
   "datum": "EGM2008", "tolerance": "3h"},
  {"name": "LakeSP Prior", "kind": "lakesp", "file": "../data/lakesp.csv",
   "lake_id": "YOUR_LAKE_ID", "datum": "EGM2008", "tolerance": "30min",
   "allow": {"quality_f": [0]}}
]
```

Use the gauge's actual vertical datum. Where a surveyed station-specific conversion is available, declare `offset_to_raster_m` and `offset_source`. The code adds that offset to the reference level. A gauge reading relative to an unknown local zero cannot establish an absolute bias. Supply `date_column` and `level_column` for other column names.

For LakeSP Prior, select one lake ID and processing version; the default fields are `lake_id,time_str,wse`. A CSV can combine passes. A vector layer can be read directly after installing `python -m pip install '.[vector]'`. Set `layer` for a multi-layer file. Set the quality field to the name in your product. If the input provides `cycle` or `pass` columns, those identifiers must match the Raster observation in addition to the time criterion. Use the same vertical reference for both products. The comparison measures agreement between two products derived from SWOT. The gauge provides the independent evaluation.

### Comparing filling and drying

To summarise performance during independently identified stages, add a `periods` list to a site, for example `[{"name":"rising water","start":"2025-05-01","end":"2025-06-01"}]`. Choose dates from your field record or imagery. Start times are inclusive and end times exclusive. The IQR flag is calculated on the full time series before splitting by period. In the result tables, select one `phase` at a time; `all` includes every period.

## Reading the outputs

| File | Contents |
| --- | --- |
| `audit.csv` | One row per site, method and granule. Includes UTC date, cycle, pass, tile, actual window width, pixel counts, WSE before temporal or coverage rejection (`wse_raw`), accepted WSE (`wse`), IQR flag and rejection reason. |
| `summary.csv` | Eligible and usable overpasses, usable fraction, number of matched overpasses, mean signed difference (Raster minus reference), median absolute difference and RMSE in metres. |
| `matches.csv` | Accepted overpasses with the nearest reference time, signed time offset in minutes, both levels and their difference. Empty reference fields indicate no match within tolerance. |

An eligible overpass is a cycle/pass with a supplied tile covering the station, including tiles with no valid WSE. Missing input files are outside this denominator. Where tiles overlap, the comparison retains the accepted tile with the largest valid fraction, breaking ties by filename. A sparse reference reading may match two overpasses within the specified tolerance. Inspect the matched timestamps before interpreting the error statistics.

## When the shoreline moves

The centre flag records whether the SWOT cell containing the station contributed to the estimate. The nearest retained cell and the centroid of all retained cell centres are located relative to the actual station coordinates, even when the station lies between grid points. Distances are in metres. The centroid uses unweighted east and north offsets derived from WGS84 geodesics. A centroid close to the station can occur when water lies on either side of dry ground. Read it together with `centre_valid` and `nearest_valid_m`.

| Audit columns | Interpretation |
| --- | --- |
| `n_in_tile`, `tile_fraction` | Cells within tile bounds, and their share of the requested square. |
| `available_fraction_in_tile` | Cells with finite WSE divided by cells within the tile. |
| `quality_retained_fraction` | Cells used after quality screening divided by cells with finite WSE. An undefined denominator gives a missing value. |
| `wse_mean_m`, `wse_median_m`, `wse_mad_m` | Mean, median and unscaled median absolute deviation within the window. Use `"estimator": "median"` to pass the median to subsequent filters. MAD describes spatial variation; gauge errors quantify vertical agreement. |
| `nearest_valid_m`, `centroid_distance_m` | Distance to the closest retained pixel centre and the offset of their centroid. Missing when no pixel contributes. |

### Comparing configurations fairly

Read `summary.csv` for availability and error on each method's full set of matched dates. Read `common_passes.csv` for error on acquisitions shared by all configured methods for a given site, resolution, phase and reference. A match must have the same cycle, pass and reference timestamp. A method with zero matches produces an empty intersection and undefined errors. `common_pass` identifies these records in `matches.csv`. Name calibration and evaluation intervals in the site's `periods`. If a method uses IQR, an optional `"iqr_fit_period": ["2025-01-01","2025-06-01"]` fits the quartiles on that interval alone (UTC, exclusive end). A group with fewer than four fitting observations stays unflagged. Select thresholds on the calibration interval, freeze them, and report performance on the subsequent evaluation interval.

### Independent water masks and a station figure

For Sentinel-2 or other imagery, supply georeferenced single-band categorical masks with 1 for water, 0 for dry and a separate nodata value for cloud or unknown. The user supplies the classification. Add `"optical": {"index":"optical.csv","tolerance":"1D"}` to the site. The index has columns `date,file` with UTC image times and paths relative to the index. The nearest image within tolerance supplies a wet/dry observation at the exact station and at the in-tile SWOT cell centres. Wet fraction is the share of known cell-centre samples classified as wet. Inspect the signed image-to-SWOT lag (`optical_lag_h`) and the known sample count. An unmatched image or unknown classification gives a missing observation. Optical evidence is separate from the WSE filter.

Set a site's `detail_dates` to a list of UTC dates to export individual pixel locations, WSE, inclusion and optical classes to `pixels.csv` for those dates only. For the synthetic example, run the following from the repository root after extraction:

```bash
python -m pip install '.[figures]'
python examples/plot_station.py example-run/results "Synthetic station" pixel 100 example-run/station.pdf 2025-01-01 2025-03-05 2025-03-26
```

The upper panels show the classified optical samples and the SWOT cells used. The lower panels track accepted WSE and the centroid offset. At the last mapped acquisition the synthetic station is dry, while retained water cells begin approximately 200 m east of the station. Inspect the optical image at native resolution when interpreting a real shoreline.

## Python and notebooks

```python
from SWOT_toolbox.stations import station

belt_bay = station("data/rasters", 137.028098, -28.893022,
                   width_m=1300, levels=[0, 1], bounds=(-16, 6))
flagged_dates = belt_bay.loc[belt_bay.iqr_flag, ["date", "wse_raw", "status"]]
```

The two original working notebooks are included with their outputs cleared. Install `python -m pip install '.[notebook]'` for the plotting functions. Set the raster and bathymetry paths in their first cells. Notebook WSE functions, the Observatory single-file adapter and the command line use the same spatial extractor. The notebook wrappers return accepted rows and retain their explicit pixel-radius, quality and optional bed-correction arguments. For equal options, the extracted dates and levels agree. To replay the historical filtering conventions use `examples/ktle_historical.json`. The named sequential and parallel configurations restore their respective bounds and pooled IQR rules; original input versions and station coordinates must also match. Comparisons for new analyses use `examples/ktle.json`, which varies only window size. Subtracting bed elevation gives a local depth only when bed and water levels share a vertical datum.

## Citation and development

Use the citation metadata in `CITATION.cff` and cite the SWOT collection version used in your analysis. The code is distributed under the MIT License. Run `python -m unittest discover -s tests -v` from the repository root to check the extraction and comparison. Contribution guidance is in `CONTRIBUTING.md`.
