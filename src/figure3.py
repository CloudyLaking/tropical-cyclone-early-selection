"""Render Figure 3: paired Hagibis tracks and selection-displacement geometry."""

from pathlib import Path
import sys, json, hashlib
import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Circle, Rectangle
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm

B = Path(__file__).resolve().parent.parent
D = B / "data"
O = B / "results"
O.mkdir(exist_ok=True, parents=True)
COL = {
    "full": "#0095A6",
    "early": "#F18A24",
    "selected": "#F18A24",
    "path": "#9457BA",
    "terminal": "#D0A018",
    "random": "#777C86",
    "selected_recentered": "#7355B4",
}
MARK = {
    "full": "s",
    "early": "D",
    "selected": "D",
    "path": "^",
    "terminal": "v",
    "random": "o",
}
plt.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "DejaVu Sans"],
        "font.size": 8.4,
        "axes.labelsize": 8.4,
        "xtick.labelsize": 8.4,
        "ytick.labelsize": 8.4,
        "axes.titlesize": 8.4,
        "axes.linewidth": 0.65,
        "axes.edgecolor": "#262626",
        "axes.spines.top": True,
        "axes.spines.right": True,
        "pdf.fonttype": 42,
        "svg.fonttype": "none",
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "legend.frameon": False,
        "lines.linewidth": 1.05,
        "savefig.facecolor": "white",
    }
)


def read(n):
    return pd.read_csv(D / n)


def sel(t, **kw):
    m = np.ones(len(t), bool)
    for k, v in kw.items():
        m &= t[k] == v
    return t[m]


def figgrid(rows, cols, size=34, gapx=12, gapy=21, bottom=19, top=17, left=14, right=3):
    w = left + right + cols * size + (cols - 1) * gapx
    h = bottom + top + rows * size + (rows - 1) * gapy
    f = plt.figure(figsize=(w / 25.4, h / 25.4))
    axs = []
    for i in range(rows * cols):
        r, c = divmod(i, cols)
        axs.append(
            f.add_axes(
                [
                    (left + c * (size + gapx)) / w,
                    (bottom + (rows - 1 - r) * (size + gapy)) / h,
                    size / w,
                    size / h,
                ]
            )
        )
    return (f, axs)


def finish(f, axs, name, cols):
    ids = [chr(97 + i) for i in range(len(axs))]
    f.canvas.draw()
    f.savefig(O / (name + ".pdf"), dpi=600)
    f.savefig(O / (name + ".svg"), dpi=600)
    f.savefig(O / (name + ".png"), dpi=600)
    dims = {
        "width_mm": f.get_figwidth() * 25.4,
        "height_mm": f.get_figheight() * 25.4,
        "axes_mm": [
            [
                a.get_position().width * f.get_figwidth() * 25.4,
                a.get_position().height * f.get_figheight() * 25.4,
            ]
            for a in axs
        ],
    }
    (O / (name + ".dimensions.json")).write_text(json.dumps(dims, indent=2))
    plt.close(f)


def zero(a):
    a.axhline(0, color="#888888", lw=0.6, zorder=0)


def handle(m, label=None, ls="-"):
    return Line2D(
        [],
        [],
        color=COL[m],
        marker=MARK.get(m, "o"),
        lw=1.1,
        ms=3,
        ls=ls,
        label=label or m,
    )


centers = ["ECMWF", "NCEP"]
geo = read("geometry_cases.csv")
leadcols = ["#7142C4", "#BC4FA0", "#EC527A", "#F48B3A", "#CEA526", "#078C98"]


def readable_heading(a, letter, title):
    title = title.replace("aligned paths", "paths").replace(
        "120 h, all cases", "all cases, 120 h"
    )
    a.text(
        0,
        1.065,
        letter,
        transform=a.transAxes,
        ha="left",
        va="bottom",
        fontsize=9.2,
        weight="bold",
    )
    a.text(
        0.13, 1.065, title, transform=a.transAxes, ha="left", va="bottom", fontsize=8.4
    )


heading = readable_heading


def readable_legend(f, handles, y=0.01, ncol=4):
    f.legend(
        handles=handles,
        loc="lower center",
        bbox_to_anchor=(0.5, y),
        ncol=ncol,
        fontsize=8.4,
        handlelength=1.5,
        columnspacing=1.3,
    )


legend = readable_legend
f, axs = figgrid(1, 4, size=34, gapx=12, bottom=24, top=12, left=14, right=4)
examples = {}
for c, name in [
    ("ECMWF", "illustrative_track.npz"),
    ("NCEP", "illustrative_track_ncep.npz"),
]:
    z = np.load(D / name, allow_pickle=False)
    truth = z["truth_xy"]
    tracks = z["p"] + truth[:, None, :]
    origin = truth[2]
    examples[c] = (truth - origin, tracks - origin, z["selected"])
allpts = np.concatenate(
    [np.concatenate([t[2:], p[2:].reshape(-1, 2)]) for t, p, ix in examples.values()]
)
lo = allpts.min(0)
hi = allpts.max(0)
mid = (lo + hi) / 2
span = max(hi - lo) * 1.08
off = []
for r, c in enumerate(centers):
    a = axs[r * 2]
    heading(a, chr(97 + r * 2), c + " · aligned paths")
    truth, tracks, ix = examples[c]
    for j in range(tracks.shape[1]):
        a.plot(
            tracks[2:, j, 0],
            tracks[2:, j, 1],
            color="#F9B35D" if j in ix else "#B4BDC2",
            alpha=0.4,
            lw=0.35,
        )
    for arr, m in [(tracks.mean(1), "full"), (tracks[:, ix].mean(1), "early")]:
        a.plot(
            arr[2:, 0],
            arr[2:, 1],
            color=COL[m],
            lw=1.35,
            marker=MARK[m],
            markevery=[0, 4, 8],
            ms=2.7,
            zorder=5,
        )
    a.plot(
        truth[2:, 0],
        truth[2:, 1],
        color="#202327",
        lw=1.2,
        marker="+",
        markevery=[0, 4, 8],
        ms=3.7,
        zorder=6,
    )
    a.set(
        xlim=(mid[0] - span / 2, mid[0] + span / 2),
        ylim=(mid[1] - span / 2, mid[1] + span / 2),
        xlabel="Eastward distance (km)",
        ylabel="Northward distance (km)",
    )
    a.set_aspect("equal", adjustable="box")
    a.tick_params(labelsize=8.4)
    if r == 1:
        a.set_ylabel("")
    a = axs[r * 2 + 1]
    heading(a, chr(98 + r * 2), c + " · 120 h, all cases")
    g = sel(geo, center=c)
    a.add_patch(Circle((1, 0), 1, facecolor="#F8FAFD", ec="#7B8BA3", lw=0.6))
    a.axhline(0, color="#D8D8D8", lw=0.5)
    a.scatter(
        g.circle_x,
        g.circle_y,
        s=4,
        color=COL["early"],
        alpha=0.45,
        lw=0,
        rasterized=True,
    )
    a.plot(0, 0, "s", color=COL["full"], ms=3.4)
    a.plot(1, 0, "+", color="#222222", ms=5)
    a.set(
        xlim=(-1, 3),
        ylim=(-2, 2),
        xticks=[-1, 0, 1, 2, 3],
        yticks=[-2, -1, 0, 1, 2],
        xlabel="Toward-truth shift /\nfull error",
        ylabel="Crosswise shift / full error",
    )
    a.set_aspect("equal", adjustable="box")
    improved = 100 * g.loc[(1 - g.circle_x) ** 2 + g.circle_y**2 < 1, "weight"].sum()
    a.text(
        0.04,
        0.94,
        f"{improved:.1f}% improve",
        transform=a.transAxes,
        va="top",
        fontsize=8.4,
    )
    outside = (g.circle_x < -1) | (g.circle_x > 3) | (abs(g.circle_y) > 2)
    off.append(
        f"{c}: {outside.sum()} off-scale ({100 * g.loc[outside, 'weight'].sum():.1f}% storm weight)"
    )
f.text(
    0.5,
    0.96,
    "a, c: Typhoon Hagibis (2019), initialized 07 Oct 2019, 12 UTC",
    ha="center",
    fontsize=8.4,
)
legend(
    f,
    [
        handle("full", "Full mean / origin"),
        handle("early", "Selected mean / shift"),
        Line2D([], [], color="#202327", marker="+", label="Truth / target"),
    ],
    y=0.045,
    ncol=3,
)
f.text(
    0.5, 0.015, "; ".join(off) + "; full ranges: Figure S1", ha="center", fontsize=8.4
)
finish(f, axs, "Figure3", 4)
out = O / "Figure3.dimensions.json"
dimensions = json.loads(out.read_text())
dimensions.update(
    {
        "manuscript_width_mm": 175,
        "minimum_source_font_pt": 8.4,
        "minimum_placed_font_pt": 7.7368421052631575,
        "arxiv_width_mm": 159.2,
        "minimum_arxiv_font_pt": 7.038315789473684,
        "case": "Typhoon Hagibis (2019), initialized 07 Oct 2019, 12 UTC",
        "case_panels": ["a", "c"],
        "all_case_panels": ["b", "d"],
    }
)
out.write_text(json.dumps(dimensions, indent=2), encoding="utf8")
print("Figure 3: readable final-size typography, unchanged data and canvas.")
