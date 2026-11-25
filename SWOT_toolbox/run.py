"""Run the sensitivity analysis described in a JSON configuration."""
import argparse
import json
from pathlib import Path
import pandas as pd
from .stations import station
from .comparison import reference_table, compare, common_passes
from .optical import load_masks


def run(config_path):
    path = Path(config_path).resolve()
    cfg = json.loads(path.read_text())
    granules = path.parent / cfg["granules"]
    output = path.parent / cfg["output"]
    output.mkdir(parents=True, exist_ok=True)
    audits, summary, matches = [], [], []
    footprints = []
    for field in ("methods", "sites"):
        names = [item["name"] for item in cfg[field]]
        if len(names) != len(set(names)):
            raise ValueError(f"Duplicate labels in {field}")
    for site in cfg["sites"]:
        masks = load_masks(site["optical"], path.parent) if "optical" in site else None
        labels = [spec.get("name", spec["kind"]) for spec in site.get("references", [])]
        if len(labels) != len(set(labels)):
            raise ValueError("Use distinct reference names at each station")
        refs = []
        for entry in site.get("references", []):
            spec = dict(entry, file=path.parent / entry["file"])
            refs.append((spec, reference_table(spec)))
        for variant in cfg["methods"]:
            options = {k: v for k, v in variant.items() if k != "name"}
            options.setdefault("bounds", site.get("bounds_m"))
            audit = station(granules, site["lon"], site["lat"], optical=masks,
                            detail_dates=site.get("detail_dates", []), **options)
            pixels = audit.attrs.pop("pixel_details", None)
            if pixels is not None:
                footprints.append(pixels.assign(site=site["name"], method=variant["name"]))
            if audit.empty:
                raise ValueError(f"No matching Raster granules for {site['name']}, {variant['name']}")
            audit = audit.assign(site=site["name"], method=variant["name"], datum=cfg["raster_datum"])
            audits.append(audit)
            for res, group in audit.groupby("resolution"):
                periods = [("all", group)]
                dates = pd.to_datetime(group.date, utc=True)
                for period in site.get("periods", []):
                    start, end = pd.to_datetime([period["start"], period["end"]], utc=True)
                    if start >= end or period["name"] == "all":
                        raise ValueError("Use a named period with start before end")
                    periods.append((period["name"], group[(dates >= start) & (dates < end)]))
                for phase, part in periods:
                    keys = dict(site=site["name"], method=variant["name"], resolution=res, phase=phase)
                    eligible = part[~part.status.isin(["outside", "unreadable"])].drop_duplicates(["cycle", "pass"])
                    usable = part[part.status.eq("accepted")].drop_duplicates(["cycle", "pass"])
                    counts = dict(eligible_passes=len(eligible), usable_passes=len(usable),
                                  usable_fraction=len(usable)/len(eligible) if len(eligible) else float("nan"))
                    summary.append(dict(keys, reference="coverage", **counts))
                    for spec, ref in refs:
                        pairs, metrics = compare(part, ref, spec, cfg["raster_datum"])
                        label = spec.get("name", spec["kind"])
                        matches.append(pairs.assign(phase=phase, reference=label, reference_kind=spec["kind"]))
                        summary.append(dict(keys, reference=label, reference_kind=spec["kind"], **metrics))
    if not audits:
        raise ValueError("Configure at least one site and one method")
    tables = {"audit": pd.concat(audits), "summary": pd.DataFrame(summary),
              "matches": pd.concat(matches, ignore_index=True) if matches else pd.DataFrame(columns=["reference","date","difference_m"])}
    tables["common_passes"] = common_passes(tables["matches"], tables["summary"])
    tables["pixels"] = pd.concat(footprints, ignore_index=True) if footprints else pd.DataFrame(
        columns=["site","method","date","cycle","pass","tile","resolution","filename","east_m","north_m",
                 "lon","lat","used","pixel_wse","optical_water"])
    for name, frame in tables.items():
        dest = output / (name + ".csv")
        tmp = dest.with_suffix(".csv.tmp")
        frame.to_csv(tmp, index=False)
        tmp.replace(dest)
    return tables


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("config")
    run(parser.parse_args().config)


if __name__ == "__main__":
    main()
