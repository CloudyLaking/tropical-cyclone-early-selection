"""Render Figure 4: retained fractions and saved hindsight-list errors."""

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
        "font.size": 6.5,
        "axes.labelsize": 6.5,
        "axes.titlesize": 7,
        "xtick.labelsize": 6,
        "ytick.labelsize": 6,
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


def heading(a, letter, title):
    a.text(
        0,
        1.065,
        letter,
        transform=a.transAxes,
        ha="left",
        va="bottom",
        fontsize=7.8,
        weight="bold",
    )
    a.text(
        0.105, 1.065, title, transform=a.transAxes, ha="left", va="bottom", fontsize=6.8
    )


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


def legend(f, handles, y=0.01, ncol=4):
    f.legend(
        handles=handles,
        loc="lower center",
        bbox_to_anchor=(0.5, y),
        ncol=ncol,
        fontsize=6.2,
        handlelength=1.8,
        columnspacing=1.25,
    )


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
storms = read("fixed_storms.csv")
fractions = read("fraction_summary.csv")
leadcols = ["#7142C4", "#BC4FA0", "#EC527A", "#F48B3A", "#CEA526", "#078C98"]
FIGURE4_FIGSIZE = (5.708661417322835, 4.409448818897638)
FIGURE5_FIGSIZE = (5.708661417322835, 4.52755905511811)


def compact_grid(figsize, bottom, gapy):
    f = plt.figure(figsize=figsize)
    width, height = f.get_size_inches() * 25.4
    axes = []
    for r in range(2):
        for c in range(3):
            axes.append(
                f.add_axes(
                    [
                        (14 + c * 47) / width,
                        (bottom + (1 - r) * (32 + gapy)) / height,
                        32 / width,
                        32 / height,
                    ]
                )
            )
    return (f, axes)


f, axs = compact_grid(FIGURE4_FIGSIZE, bottom=15, gapy=16)
atecm = LinearSegmentedColormap.from_list(
    "ate", ["#FDBB24", "#FFF9D6", "#FFFFFF", "#BD8ACC", "#7737A5"]
)
ctecm = LinearSegmentedColormap.from_list(
    "cte", ["#338ED0", "#BFE2F0", "#FFFFFF", "#FEDC95", "#F28C2B"]
)
for r, c in enumerate(centers):
    for j, metric in enumerate(["ATE", "CTE"]):
        a = axs[r * 3 + j]
        heading(a, chr(97 + r * 3 + j), c + " · " + metric + " by count")
        t = sel(fractions, center=c, metric=metric)
        mat = (
            t.pivot(index="fraction", columns="lead", values="mean")
            .loc[range(10, 101, 10), [48, 72, 96, 120]]
            .to_numpy()
        )
        cm = atecm if metric == "ATE" else ctecm
        vmax = 180 if metric == "ATE" else 120
        im = a.imshow(
            mat,
            origin="upper",
            cmap=cm,
            norm=TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax),
            aspect="auto",
            interpolation="nearest",
        )
        a.set(
            xticks=range(4),
            xticklabels=[48, 72, 96, 120],
            yticks=[0, 2, 4, 6, 9],
            yticklabels=[10, 30, 50, 70, 100],
            xlabel="Lead (h)",
            ylabel="Retained members (%)",
        )
        for k in [4, 9]:
            a.add_patch(
                Rectangle((-0.5, k - 0.5), 4, 1, fill=False, ec="#444444", lw=0.55)
            )
        for k in range(10):
            a.text(
                3,
                k,
                f"{mat[k, 3]:.0f}",
                ha="center",
                va="center",
                fontsize=5.8,
                color="#222222",
            )
        if r == 0:
            pos = a.get_position()
            cb = f.add_axes([pos.x0, pos.y1 + 0.08, pos.width, 0.018])
            f.colorbar(im, cax=cb, orientation="horizontal", ticks=[-vmax, 0, vmax])
            cb.tick_params(labelsize=5.8, length=2, pad=1)
            cb.set_title(metric + " (km)", fontsize=6, pad=2)
    a = axs[r * 3 + 2]
    heading(a, chr(99 + r * 3), c + " · 120-h ATE")
    for k, m in enumerate(["full", "early", "path", "terminal"]):
        t = sel(storms, center=c, lead=120, method=m).sort_values("sid")
        v = t.ATE.to_numpy()
        a.boxplot(
            [v],
            positions=[k],
            widths=0.5,
            whis=(5, 95),
            showfliers=False,
            patch_artist=True,
            manage_ticks=False,
            boxprops={
                "facecolor": COL[m] + "55",
                "edgecolor": COL[m],
                "linewidth": 0.9,
            },
            medianprops={"color": COL[m], "linewidth": 1.2},
            whiskerprops={"color": COL[m], "linewidth": 0.8},
            capprops={"color": COL[m], "linewidth": 0.8},
        )
        jitter = np.array(
            [
                (
                    int(hashlib.sha256(str(s).encode()).hexdigest()[:8], 16)
                    % 1000
                    / 999
                    - 0.5
                )
                * 0.28
                for s in t.sid
            ]
        )
        a.scatter(k + jitter, v, s=3, color=COL[m], alpha=0.35, lw=0, rasterized=True)
        a.plot(k, float(np.mean(v)), "+", color="#222222", ms=4, mew=0.8, zorder=5)
    a.set(
        xticks=range(4),
        xticklabels=["Full\nensemble", "Early\nhalf", "Path\nhalf", "120-h\nhalf"],
        ylabel="Storm-mean ATE (km)",
        xlim=(-0.5, 3.5),
    )
    a.tick_params(axis="x", labelsize=5.8, pad=3)
    zero(a)
    a.plot(
        [0.55, 0.55, 0.95, 0.95],
        [-0.19, -0.22, -0.22, -0.19],
        transform=a.transAxes,
        clip_on=False,
        color="#888888",
        lw=0.5,
    )
    a.text(
        0.75,
        -0.25,
        "Hindsight",
        transform=a.transAxes,
        ha="center",
        va="top",
        fontsize=5.8,
    )
finish(f, axs, "Figure4", 3)
