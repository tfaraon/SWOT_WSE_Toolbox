"""Compare one Raster configuration with a gauge or a LakeSP Prior time series."""
from pathlib import Path
import numpy as np
import pandas as pd


def reference_table(spec):
    path = Path(spec["file"])
    kind = spec["kind"]
    if kind not in ("gauge", "lakesp"):
        raise ValueError("Reference kind must be gauge or lakesp")
    if path.suffix.lower() == ".csv":
        source = pd.read_csv(path, dtype=str)
    elif kind == "lakesp":
        import pyogrio
        source = pyogrio.read_dataframe(path, read_geometry=False, layer=spec.get("layer"))
    else:
        raise ValueError("Gauge inputs must be CSV files")
    if kind == "lakesp":
        field = spec.get("id_column", "lake_id")
        source = source[source[field].astype(str) == str(spec["lake_id"])]
        if source.empty:
            raise ValueError(f"No matching lake_id in {path}")
    for col, accepted in spec.get("allow", {}).items():
        source = source[source[col].astype(str).isin([str(v) for v in accepted])]
    date_col = spec.get("date_column", "time_str" if kind == "lakesp" else "date")
    field = spec.get("level_column", "wse" if kind == "lakesp" else "level_m")
    ref = pd.DataFrame({"reference_date": pd.to_datetime(source[date_col], utc=True, errors="coerce", format="mixed"),
                        "reference_m": pd.to_numeric(source[field], errors="coerce")})
    ref = ref[np.isfinite(ref.reference_m) & ref.reference_date.notna()].copy()
    if kind == "lakesp":
        ref = ref[ref.reference_m > -1e8].copy()
    for col in ("cycle", "pass"):
        if col in source:
            ref[col] = pd.to_numeric(source.loc[ref.index, col], errors="raise").astype(int)
    if not spec.get("datum"):
        raise ValueError("Declare the reference vertical datum")
    offset = float(spec.get("offset_to_raster_m", 0))
    if not np.isfinite(offset):
        raise ValueError("Invalid reference offset")
    ref["reference_m"] += offset
    ref["reference_date"] = ref.reference_date.astype("datetime64[ns, UTC]")
    return ref


def compare(audit, ref, spec, raster_datum):
    """One eligible overpass per cycle/pass; never average across resolutions."""
    if spec["datum"] != raster_datum and "offset_to_raster_m" not in spec:
        raise ValueError("Different vertical datums require a documented conversion")
    if spec["datum"] != raster_datum and not spec.get("offset_source"):
        raise ValueError("Describe the source of the station-specific offset")
    tol = pd.Timedelta(spec["tolerance"])
    if tol < pd.Timedelta(0) or pd.isna(tol):
        raise ValueError("Invalid time tolerance")
    if audit.resolution.nunique() > 1:
        raise ValueError("Compare one resolution at a time")
    keys = ["cycle", "pass"]
    eligible = audit[~audit.status.isin(["outside", "unreadable"])]
    n = len(eligible.drop_duplicates(keys))
    # On overlapping tiles retain the accepted estimate with greatest coverage.
    obs = (eligible[eligible.status.eq("accepted")].sort_values(
        ["valid_fraction", "filename"], ascending=[False, True])
        .drop_duplicates(keys).copy())
    obs["date"] = pd.to_datetime(obs.date, utc=True).astype("datetime64[ns, UTC]")
    by = [k for k in keys if k in ref and spec["kind"] == "lakesp"]
    for k in by:
        obs[k] = obs[k].astype(int)
    if ref.duplicated([*by, "reference_date"]).any():
        raise ValueError("Multiple reference levels at one time; select one series or product version")
    ref = ref.copy()
    ref["reference_date"] = pd.to_datetime(ref.reference_date, utc=True).astype("datetime64[ns, UTC]")
    pairs = pd.merge_asof(obs.sort_values("date"), ref[[*by, "reference_date", "reference_m"]].sort_values("reference_date"),
                          left_on="date", right_on="reference_date", by=by or None,
                          direction="nearest", tolerance=tol)
    pairs["lag_minutes"] = (pairs.reference_date-pairs.date).dt.total_seconds()/60
    pairs["difference_m"] = pairs.wse-pairs.reference_m
    diff = pairs.difference_m.dropna().astype(float)
    metrics = dict(eligible_passes=n, usable_passes=len(obs), matched_passes=len(diff),
                   usable_fraction=len(obs)/n if n else np.nan,
                   bias_m=diff.mean(), median_absolute_difference_m=diff.abs().median(),
                   rmse_m=np.sqrt(diff.pow(2).mean()))
    return pairs, metrics


def common_passes(matches, summary):
    """Score identical overpasses matched to the identical reference reading."""
    keys = ["site", "resolution", "phase", "reference", "reference_kind"]
    scores = []
    matches["common_pass"] = False
    available = summary[summary.reference_kind.notna()] if "reference_kind" in summary else summary.iloc[:0]
    for group, rows in available.groupby(keys) if len(available) else ():
        selected = pd.Series(True, index=matches.index)
        for key, value in zip(keys, group):
            selected &= matches[key].eq(value)
        valid = matches[selected & matches.difference_m.notna()].copy()
        paired = valid.groupby(["cycle", "pass", "reference_date"]).method.transform("nunique")
        subset = valid[paired.eq(rows.method.nunique())]
        matches.loc[subset.index, "common_pass"] = True
        for method in rows.method:
            errors = subset.loc[subset.method.eq(method), "difference_m"].astype(float)
            scores.append(dict(zip(keys, group), method=method, common_passes=len(errors),
                               bias_m=errors.mean(), median_absolute_difference_m=errors.abs().median(),
                               rmse_m=np.sqrt(errors.pow(2).mean())))
    return pd.DataFrame(scores, columns=[*keys, "method", "common_passes", "bias_m",
                                         "median_absolute_difference_m", "rmse_m"])
