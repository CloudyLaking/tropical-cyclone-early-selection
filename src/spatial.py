"""Render Figure S15 from saved gridded errors and the offline coastline."""

from pathlib import Path
import gzip
import hashlib
import json
import pickle
import sys
import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.cm import ScalarMappable
import cartopy.crs as ccrs
from cartopy.io import shapereader
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from shapely.geometry import box

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "results"
DATA = ROOT / "data"
QA = ROOT / "results/checks"
for target in [OUTPUT, QA]:
    target.mkdir(parents=True, exist_ok=True)
COAST = ROOT / "assets/natural_earth/ne_50m_coastline.shp"
grid = pd.read_csv(DATA / "spatial120_grid.csv")
LON_EDGES = np.arange(90, 181, 5)
LAT_EDGES = np.arange(0, 66, 5)
assert set(grid.center) == {"ECMWF", "NCEP"}
assert grid.shape[0] == 2 * (len(LON_EDGES) - 1) * (len(LAT_EDGES) - 1)
limits = {}
for variable in ["ATE", "CTE"]:
    values = grid.loc[
        grid.shown, ["full_" + variable, "selected_" + variable]
    ].to_numpy()
    maximum = float(np.nanmax(np.abs(values)))
    limits[variable] = max(100, int(np.ceil(maximum / 100)) * 100)
cmaps = {
    "ATE": LinearSegmentedColormap.from_list(
        "blue_white_orange", ["#2F70B5", "#B5D9ED", "#FFFFFF", "#FFD183", "#ED8419"]
    ),
    "CTE": LinearSegmentedColormap.from_list(
        "yellow_white_purple", ["#F2BE27", "#FFF0AA", "#FFFFFF", "#CFB8E7", "#8A4EBA"]
    ),
}
plt.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "DejaVu Sans"],
        "font.size": 6.5,
        "axes.labelsize": 6.5,
        "axes.titlesize": 7,
        "xtick.labelsize": 6,
        "ytick.labelsize": 6,
        "axes.titleweight": "normal",
        "font.weight": "normal",
        "axes.linewidth": 0.6,
        "pdf.fonttype": 42,
        "svg.fonttype": "none",
        "figure.facecolor": "white",
        "axes.facecolor": "white",
    }
)
width_mm = 180
height_mm = 112
fig = plt.figure(figsize=(width_mm / 25.4, height_mm / 25.4))
projection = ccrs.Mercator(central_longitude=135)
geographic = ccrs.PlateCarree()
corners = projection.transform_points(
    geographic, np.array([90.0, 180.0]), np.array([0.0, 0.0])
)
left, right_edge = corners[:, 0]
span = right_edge - left
north = geographic.transform_point(0, span, projection)[1]
map_window = box(90 + 1e-07, 1e-07, 180 - 1e-07, north - 1e-07)
coast_geometries = [
    g.intersection(map_window)
    for g in shapereader.Reader(COAST).geometries()
    if g.intersects(map_window)
]
axes = []
panels = [
    ("full_ATE", "full ATE", "ATE"),
    ("full_CTE", "full CTE", "CTE"),
    ("selected_ATE", "selected ATE", "ATE"),
    ("selected_CTE", "selected CTE", "CTE"),
]
for row, center in enumerate(["ECMWF", "NCEP"]):
    subset = grid[grid.center == center]
    for col, (metric, title, variable) in enumerate(panels):
        x, y, side = (15 + 43 * col, 64 - 45 * row, 32)
        ax = fig.add_axes(
            [x / width_mm, y / height_mm, side / width_mm, side / height_mm],
            projection=projection,
        )
        axes.append(ax)
        ax.set_xlim(left, right_edge)
        ax.set_ylim(0, span)
        ax.set_aspect("equal", adjustable="box")
        values = subset.pivot(index="lat_lower", columns="lon_lower", values=metric)
        counts = subset.pivot(
            index="lat_lower", columns="lon_lower", values="n_initializations"
        )
        values = values.reindex(index=LAT_EDGES[:-1], columns=LON_EDGES[:-1]).to_numpy()
        counts = counts.reindex(index=LAT_EDGES[:-1], columns=LON_EDGES[:-1]).to_numpy()
        masked = np.ma.masked_where((counts < 5) | ~np.isfinite(values), values)
        ax.pcolormesh(
            LON_EDGES,
            LAT_EDGES,
            masked,
            transform=geographic,
            cmap=cmaps[variable],
            norm=TwoSlopeNorm(vmin=-limits[variable], vcenter=0, vmax=limits[variable]),
            shading="flat",
            linewidth=0,
            edgecolors="none",
            zorder=2,
        )
        for lon in [120, 150]:
            ax.plot(
                [lon, lon],
                [0, north],
                transform=geographic,
                color="#C9CDD1",
                linestyle=(0, (2, 2)),
                linewidth=0.35,
                zorder=1,
            )
        for lat in [20, 40, 60]:
            ax.plot(
                [90, 180],
                [lat, lat],
                transform=geographic,
                color="#C9CDD1",
                linestyle=(0, (2, 2)),
                linewidth=0.35,
                zorder=1,
            )
        ax.add_geometries(
            coast_geometries,
            crs=geographic,
            facecolor="none",
            edgecolor="#242424",
            linewidth=0.42,
            zorder=3,
        )
        ax.set_xticks([90, 120, 150, 180], crs=geographic)
        ax.set_yticks([0, 20, 40, 60], crs=geographic)
        ax.xaxis.set_major_formatter(
            LongitudeFormatter(number_format=".0f", degree_symbol="°")
        )
        ax.yaxis.set_major_formatter(
            LatitudeFormatter(number_format=".0f", degree_symbol="°")
        )
        ax.tick_params(length=2, width=0.5, pad=2)
        ax.spines["geo"].set_linewidth(0.6)
        ax.spines["geo"].set_edgecolor("#333333")
        ax.text(
            0,
            1.055,
            chr(97 + row * 4 + col),
            fontsize=8.5,
            fontweight="bold",
            transform=ax.transAxes,
            va="baseline",
            ha="left",
        )
        ax.text(
            0.12,
            1.055,
            center + " · " + title,
            fontsize=6.8,
            transform=ax.transAxes,
            va="baseline",
            ha="left",
        )
for col, (metric, title, variable) in enumerate(panels):
    x, side = (15 + 43 * col, 32)
    cax = fig.add_axes([x / width_mm, 8 / height_mm, side / width_mm, 1.7 / height_mm])
    bar = fig.colorbar(
        ScalarMappable(
            norm=TwoSlopeNorm(vmin=-limits[variable], vcenter=0, vmax=limits[variable]),
            cmap=cmaps[variable],
        ),
        cax=cax,
        orientation="horizontal",
        ticks=[-limits[variable], 0, limits[variable]],
    )
    bar.ax.tick_params(labelsize=6.3, pad=1, length=2, width=0.5)
    bar.set_label(variable + " (km)", fontsize=6.5, labelpad=2)
    bar.outline.set_linewidth(0.5)
fig.canvas.draw()
fig.savefig(OUTPUT / "FigureS15.pdf", facecolor="white")
fig.savefig(OUTPUT / "FigureS15.svg", facecolor="white")
fig.savefig(OUTPUT / "FigureS15.png", dpi=600, facecolor="white")
plt.close(fig)
print("S15 redrawn from the frozen spatial grid; no aggregation or new experiment.")
