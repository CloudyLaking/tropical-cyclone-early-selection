"""Render Figure 2: member ranks, group transitions, and storm error distributions."""

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
        "font.size": 7.85,
        "axes.labelsize": 7.85,
        "xtick.labelsize": 7.85,
        "ytick.labelsize": 7.85,
        "axes.titlesize": 7.85,
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
storms = read("fixed_storms.csv")
ranks = read("rank_806_summary.csv")
hist = read("member_histories.csv.gz")
leadcols = ["#7142C4", "#BC4FA0", "#EC527A", "#F48B3A", "#CEA526", "#078C98"]


def readable_heading(a, letter, title):
    title = title.replace("member ranks", "ranks").replace("later groups", "groups")
    a.text(
        0,
        1.065,
        letter,
        transform=a.transAxes,
        ha="left",
        va="bottom",
        fontsize=8.8,
        weight="bold",
    )
    a.text(
        0.13, 1.065, title, transform=a.transAxes, ha="left", va="bottom", fontsize=7.85
    )


heading = readable_heading
f, axs = figgrid(2, 4, size=29, gapx=15, gapy=19, bottom=25, top=11, left=13, right=4)
rankcols = ["#EDAE28", "#EA7B72", "#D4539A", "#955BC4", "#55319B"]
rankmarks = ["o", "v", "s", "^", "D"]
distribution_rows = []


def stormboxes(a, c, metric, methods, leads):
    for i, h in enumerate(leads):
        for j, m in enumerate(methods):
            vals = sel(storms, center=c, lead=h, method=m).sort_values("sid")
            v = vals[metric].to_numpy()
            x = i + (j - (len(methods) - 1) / 2) * 0.29
            bp = a.boxplot(
                [v],
                positions=[x],
                widths=0.23,
                whis=(5, 95),
                showfliers=False,
                patch_artist=True,
                manage_ticks=False,
                boxprops={
                    "facecolor": COL[m] + "55",
                    "edgecolor": COL[m],
                    "linewidth": 0.8,
                },
                medianprops={"color": COL[m], "linewidth": 1.2},
                whiskerprops={"color": COL[m], "linewidth": 0.7},
                capprops={"color": COL[m], "linewidth": 0.7},
            )
            jitter = np.array(
                [
                    (
                        int(hashlib.sha256(str(s).encode()).hexdigest()[:8], 16)
                        % 1000
                        / 999
                        - 0.5
                    )
                    * 0.17
                    for s in vals.sid
                ]
            )
            a.scatter(
                x + jitter,
                v,
                s=2.1,
                color=COL[m],
                alpha=0.19,
                lw=0,
                rasterized=True,
                zorder=1,
            )
            a.plot(
                x, float(np.mean(v)), "+", color="#222222", ms=3.2, mew=0.7, zorder=5
            )
            distribution_rows.append(
                {
                    "center": c,
                    "lead": h,
                    "method": m,
                    "metric": metric,
                    "storms": len(v),
                    "q05": np.quantile(v, 0.05),
                    "q25": np.quantile(v, 0.25),
                    "median": np.median(v),
                    "q75": np.quantile(v, 0.75),
                    "q95": np.quantile(v, 0.95),
                }
            )
    a.set(
        xlim=(-0.6, len(leads) - 0.4),
        xticks=range(len(leads)),
        xticklabels=leads,
        xlabel="Lead (h)",
        ylabel=f"Storm-mean {metric} (km)",
    )
    zero(a)


for r, c in enumerate(centers):
    a = axs[r * 4]
    heading(a, chr(97 + r * 4), c + " · member ranks")
    date = ["2015-05-12T12:00:00", "2015-10-18T00:00:00"][r]
    z = sel(hist, center=c, base=date, scoring="raw")
    for mid, t in z.groupby("member"):
        t = t.sort_values("lead")
        a.plot(
            t.lead,
            100 * t.rank_fraction,
            color=rankcols[min(4, int(t.early_fraction.iloc[0] * 5))],
            lw=0.6,
            ls="-" if t.selected.iloc[0] else "--",
            alpha=0.85,
        )
    a.set(
        xlim=(24, 120),
        ylim=(100, 0),
        xticks=[24, 72, 120],
        yticks=[0, 50, 100],
        xlabel="Lead (h)",
        ylabel="Member rank (%)",
    )
    a = axs[r * 4 + 1]
    heading(a, chr(98 + r * 4), c + " · later groups")
    t = sel(ranks, center=c, lead=120).set_index("metric")
    a.axhline(20, color="#777777", ls=":", lw=0.7)
    for k in range(5):
        z = t.loc[[f"m{k}{j}" for j in range(5)]]
        v = 100 * z["mean"].to_numpy()
        a.errorbar(
            range(5),
            v,
            yerr=[v - 100 * z.ci_low, 100 * z.ci_high - v],
            color=rankcols[k],
            marker=rankmarks[k],
            ms=2.7,
            lw=1.1 if k in [0, 4] else 0.65,
            capsize=1,
            elinewidth=0.5,
        )
    a.set(
        xlim=(-0.3, 4.3),
        ylim=(10, 34),
        xticks=range(5),
        xticklabels=["Best", "2", "3", "4", "Worst"],
        xlabel="120-h error group",
        ylabel="Later-group members (%)",
    )
    a.tick_params(axis="x", labelsize=7.85, pad=2)
    for j, metric in enumerate(["ATE", "CTE"]):
        a = axs[r * 4 + 2 + j]
        heading(a, chr(99 + r * 4 + j), c + " · " + metric)
        stormboxes(a, c, metric, ["full", "early"], [48, 72, 96, 120])
left_pair_center = (axs[0].get_position().x0 + axs[1].get_position().x1) / 2
right_pair_center = (axs[2].get_position().x0 + axs[3].get_position().x1) / 2
f.legend(
    handles=[
        Line2D([], [], color=col, label=label)
        for col, label in zip(rankcols, ["Best", "2", "3", "4", "Worst"])
    ],
    loc="lower center",
    bbox_to_anchor=(left_pair_center, 0.005),
    ncol=5,
    fontsize=7.85,
    handlelength=0.9,
    columnspacing=0.65,
    title="24-h error group (a, b, e, f)",
    title_fontsize=7.85,
)
f.legend(
    handles=[handle("full", "Full"), handle("early", "Early half")],
    loc="lower center",
    bbox_to_anchor=(right_pair_center, 0.005),
    ncol=2,
    fontsize=7.85,
    handlelength=1.1,
    columnspacing=0.8,
    title="Distributions (c, d, g, h)",
    title_fontsize=7.85,
)
f.text(
    left_pair_center,
    0.122,
    "a, e: solid selected; dashed excluded",
    ha="center",
    fontsize=7.85,
)
finish(f, axs, "Figure2", 4)
out = O / "Figure2.dimensions.json"
dimensions = json.loads(out.read_text())
dimensions.update(
    {
        "manuscript_width_mm": 160,
        "minimum_source_font_pt": 7.85,
        "minimum_placed_font_pt": 7.056179775280898,
        "arxiv_width_mm": 159.2,
        "minimum_arxiv_font_pt": 7.020898876404493,
        "legend_scope": {
            "left pair": "24-h rank group; solid selected/dashed excluded in a,e",
            "right pair": "full versus early-selected half distributions",
        },
        "date_labels": "Exact example initialization dates are retained in the manuscript caption.",
    }
)
out.write_text(json.dumps(dimensions, indent=2), encoding="utf8")
print("Figure 2: readable final-size typography, unchanged data and canvas.")
