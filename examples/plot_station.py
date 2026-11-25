"""Plot selected station footprints and the corresponding water-level series."""
import argparse
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import pandas as pd


def plot(folder, site, method, resolution, dates, output):
    folder = Path(folder)
    def read(name):
        df = pd.read_csv(folder / name)
        return df[df.site.eq(site) & df.method.eq(method) & df.resolution.eq(resolution)].copy()
    pixels, audit = read("pixels.csv"), read("audit.csv")
    if not 1 <= len(dates) <= 3:
        raise ValueError("Choose one to three acquisition dates")
    # Same overlap tie-breaker as the validation. Prefer accepted estimates.
    audit["accepted"] = audit.status.eq("accepted")
    audit = (audit[~audit.status.isin(["outside", "unreadable"])].sort_values(["accepted","valid_fraction","filename"],
             ascending=[False, False, True]).drop_duplicates(["cycle","pass"]).sort_values("date"))
    audit["date"] = pd.to_datetime(audit.date)
    fig = plt.figure(figsize=(3.3*max(2,len(dates)),5.8), layout="constrained")
    grid = fig.add_gridspec(3,len(dates), height_ratios=[2.3,1,1])
    for i, day in enumerate(dates):
        rows = audit[audit.date.dt.strftime("%Y-%m-%d").eq(day)]
        if len(rows) != 1:
            raise ValueError(f"Expected exactly one covering passage on {day}; check the audit.")
        data = pixels[pixels.filename.eq(rows.filename.iloc[0])]
        if data.empty:
            raise ValueError(f"Add {day} to detail_dates and rerun extraction")
        ax = fig.add_subplot(grid[0,i])
        colors = data.optical_water.map({1:"lightblue",0:"tan"}).fillna("lightgrey")
        ax.scatter(data.east_m, data.north_m, c=colors, marker="s", s=18)
        used = data[data.used]
        ax.scatter(used.east_m,used.north_m,c="navy",s=5)
        ax.scatter(0,0,marker="+",color="crimson",s=100)
        ax.set(aspect="equal",title=day,xlabel="East (m)")
        if i == 0:
            ax.set_ylabel("North (m)")
    series = fig.add_subplot(grid[1,:])
    dist = fig.add_subplot(grid[2,:], sharex=series)
    series.plot(audit.date,audit.wse,"o-",color="navy",ms=3)
    dist.plot(audit.date,audit.centroid_distance_m.where(audit.accepted),"o-",color="darkgreen",ms=3)
    series.set_ylabel("WSE (m)")
    dist.set_ylabel("Centroid offset (m)")
    for day in dates:
        for ax in (series,dist):
            ax.axvline(pd.Timestamp(day),ls=":",c="0.6",lw=.7)
    symbols = [("s","lightblue","Wet"),("s","tan","Dry"),("s","lightgrey","Unknown"),
               (".","navy","Used by SWOT"),("+","crimson","Station")]
    fig.legend(handles=[Line2D([],[],ls="",marker=m,color=c,label=l) for m,c,l in symbols],
               loc="outside lower center",ncol=5,frameon=False)
    fig.savefig(output,dpi=250)
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for field in ("folder","site","method"):
        parser.add_argument(field)
    parser.add_argument("resolution", type=int)
    parser.add_argument("output")
    parser.add_argument("dates", nargs="+")
    a = parser.parse_args()
    plot(a.folder, a.site, a.method, a.resolution, a.dates, a.output)
