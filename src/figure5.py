"""Render Figure 5: empirical member ranges, coverage, and center adjustment."""

from pathlib import Path
import json
import hashlib
import sys
import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

B = Path(__file__).resolve().parent.parent
D = B / "data"
O = B / "results"
Q = B / "results/checks"
O.mkdir(exist_ok=True, parents=True)
Q.mkdir(exist_ok=True, parents=True)
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


def sel(table, **kwargs):
    mask = np.ones(len(table), bool)
    for key, value in kwargs.items():
        mask &= table[key] == value
    return table[mask]


def heading(axis, letter, title):
    axis.text(
        0,
        1.065,
        letter,
        transform=axis.transAxes,
        ha="left",
        va="bottom",
        fontsize=7.8,
        weight="bold",
    )
    axis.text(
        0.105,
        1.065,
        title,
        transform=axis.transAxes,
        ha="left",
        va="bottom",
        fontsize=6.8,
    )


def compact_grid(figsize, bottom, gapy):
    fig = plt.figure(figsize=figsize)
    width, height = fig.get_size_inches() * 25.4
    axes = []
    for row in range(2):
        for column in range(3):
            axes.append(
                fig.add_axes(
                    [
                        (14 + column * 47) / width,
                        (bottom + (1 - row) * (32 + gapy)) / height,
                        32 / width,
                        32 / height,
                    ]
                )
            )
    return (fig, axes)


def handle(method, label=None, ls="-"):
    return Line2D(
        [],
        [],
        color=COL[method],
        marker=MARK.get(method, "o"),
        lw=1.1,
        ms=3,
        ls=ls,
        label=label or method,
    )


def finish(fig, axes, name, columns):
    ids = [chr(97 + i) for i in range(len(axes))]
    fig.canvas.draw()
    fig.savefig(O / (name + ".pdf"), dpi=600)
    fig.savefig(O / (name + ".svg"), dpi=600)
    fig.savefig(O / (name + ".png"), dpi=600)
    dimensions = {
        "width_mm": fig.get_figwidth() * 25.4,
        "height_mm": fig.get_figheight() * 25.4,
        "axes_mm": [
            [
                axis.get_position().width * fig.get_figwidth() * 25.4,
                axis.get_position().height * fig.get_figheight() * 25.4,
            ]
            for axis in axes
        ],
    }
    (O / (name + ".dimensions.json")).write_text(json.dumps(dimensions, indent=2))
    plt.close(fig)


FIGURE5_FIGSIZE = (5.708661417322835, 4.52755905511811)
cov = pd.read_csv(D / "coverage_summary.csv")
cross = pd.read_csv(D / "coverage_crossings.csv")
contrast_records = []
fig, axes = compact_grid(FIGURE5_FIGSIZE, bottom=24, gapy=15)
for row, center in enumerate(["ECMWF", "NCEP"]):
    axis = axes[row * 3]
    heading(axis, chr(97 + row * 3), center + " · center and radius")
    for j, metric in enumerate(["distance", "r70"]):
        for method in ["full", "selected"]:
            table = sel(
                cov, center=center, cohort="motion806", comparison=method, metric=metric
            ).sort_values("lead")
            axis.plot(
                table.lead,
                table["mean"],
                color=COL[method],
                ls="-" if j == 0 else "--",
                marker=MARK[method],
                mfc=COL[method] if j == 0 else "white",
                ms=3,
                lw=1.15,
            )
    axis.set(
        xlabel="Lead (h)",
        ylabel="Distance or radius (km)",
        xticks=[48, 72, 96, 120],
        ylim=(0, None),
    )
    axis = axes[row * 3 + 1]
    heading(axis, chr(98 + row * 3), center + " · observed coverage")
    for method in ["full", "selected", "random"]:
        table = sel(
            cov, center=center, cohort="motion806", comparison=method, metric="in70"
        ).sort_values("lead")
        axis.fill_between(
            table.lead,
            100 * table.ci_low,
            100 * table.ci_high,
            color=COL[method],
            alpha=0.1,
            lw=0,
        )
        axis.plot(
            table.lead,
            100 * table["mean"],
            color=COL[method],
            marker=MARK[method],
            ms=3,
            ls="--" if method == "random" else "-",
            lw=1.2,
        )
    axis.set(
        xlabel="Lead (h)",
        ylabel="Truth inside range (%)",
        xticks=[48, 72, 96, 120],
        ylim=(0, 100),
    )
    contrast = sel(
        cov,
        center=center,
        cohort="motion806",
        lead=120,
        comparison="selected_minus_random",
        metric="in70",
    )
    assert len(contrast) == 1, "Expected one saved paired contrast per center"
    item = contrast.iloc[0]
    mean, low, high = [100 * float(item[key]) for key in ["mean", "ci_low", "ci_high"]]
    assert np.isfinite([mean, low, high]).all() and low <= mean <= high

    def signed(value):
        return f"{value:+.2f}".replace("-", "−")

    label = (
        f"120 h: selected − random\n{signed(mean)} pp [{signed(low)}, {signed(high)}]"
    )
    axis.text(
        0.04,
        0.22,
        label,
        transform=axis.transAxes,
        ha="left",
        va="top",
        fontsize=5.8,
        linespacing=1.2,
        color="#262626",
    )
    contrast_records.append(
        {
            "center": center,
            "panel": chr(98 + row * 3),
            "source_csv": "data/coverage_summary.csv",
            "predicate": {
                "cohort": "motion806",
                "center": center,
                "lead": 120,
                "comparison": "selected_minus_random",
                "metric": "in70",
            },
            "saved_mean": float(item["mean"]),
            "saved_ci_low": float(item["ci_low"]),
            "saved_ci_high": float(item["ci_high"]),
            "units_source": "fraction",
            "display_scale": 100,
            "units_display": "percentage points",
            "cases": int(item["cases"]),
            "storms": int(item["storms"]),
            "annotation": label,
            "annotation_axes_position": [0.04, 0.22],
        }
    )
    axis = axes[row * 3 + 2]
    heading(axis, chr(99 + row * 3), center + " · center adjustment")
    table = sel(cross, center=center)
    xmax = 3.5
    axis.plot([0, xmax], [0, xmax], color="#CCCCCC", lw=0.55)
    axis.axvline(1, color="#888888", lw=0.6)
    axis.axhline(1, color="#888888", lw=0.6)
    for state, color, marker in [
        ("outside", "#B9C2CB", "o"),
        ("inside", "#9ABACB", "o"),
        ("enter", "#8A4EB7", "^"),
        ("exit", "#DEA42A", "v"),
    ]:
        values = table[
            (table.state == state) & table.x.between(0, xmax) & table.y.between(0, xmax)
        ]
        axis.scatter(
            values.x,
            values.y,
            s=4 if state in ["outside", "inside"] else 7,
            color=color,
            marker=marker,
            alpha=0.5 if state in ["outside", "inside"] else 0.8,
            lw=0,
            rasterized=False,
        )
    axis.set(
        xlim=(0, xmax),
        ylim=(0, xmax),
        xticks=[0, 1, 2, 3],
        yticks=[0, 1, 2, 3],
        xlabel="Before: distance / radius",
        ylabel="After: distance / radius",
    )
    axis.set_aspect("equal", adjustable="box")
    entry = 100 * table.loc[table.state == "enter", "weight"].sum()
    exit_fraction = 100 * table.loc[table.state == "exit", "weight"].sum()
    axis.text(
        0.03,
        0.95,
        f"Enter\n{entry:.1f}%\nExit\n{exit_fraction:.1f}%",
        transform=axis.transAxes,
        va="top",
        fontsize=5.8,
        linespacing=1.25,
    )
column_centers = [
    (axes[j].get_position().x0 + axes[j].get_position().x1) / 2 for j in range(3)
]
fig.legend(
    handles=[
        handle("full", "Full: center error"),
        Line2D(
            [],
            [],
            color=COL["full"],
            marker="s",
            mfc="white",
            ls="--",
            ms=3,
            lw=1.1,
            label="Full: radius",
        ),
        handle("early", "Early half: center error"),
        Line2D(
            [],
            [],
            color=COL["early"],
            marker="D",
            mfc="white",
            ls="--",
            ms=3,
            lw=1.1,
            label="Early half: radius",
        ),
    ],
    loc="lower center",
    bbox_to_anchor=(column_centers[0], 0.005),
    ncol=1,
    fontsize=5.8,
    handlelength=1.8,
    labelspacing=0.3,
    title="a, d",
    title_fontsize=5.8,
)
fig.legend(
    handles=[
        handle("full", "Full ensemble"),
        handle("early", "Early-selected half"),
        handle("random", "Random half", ls="--"),
    ],
    loc="lower center",
    bbox_to_anchor=(column_centers[1], 0.005),
    ncol=1,
    fontsize=5.8,
    handlelength=1.8,
    labelspacing=0.3,
    title="b, e",
    title_fontsize=5.8,
)
fig.legend(
    handles=[
        Line2D(
            [],
            [],
            ls="none",
            marker="^",
            color="#8A4EB7",
            ms=3,
            label="Truth enters range",
        ),
        Line2D(
            [],
            [],
            ls="none",
            marker="v",
            color="#DEA42A",
            ms=3,
            label="Truth leaves range",
        ),
        Line2D(
            [],
            [],
            ls="none",
            marker="o",
            color="#B9C2CB",
            ms=2.8,
            label="Remains outside",
        ),
        Line2D(
            [],
            [],
            ls="none",
            marker="o",
            color="#9ABACB",
            ms=2.8,
            label="Remains inside",
        ),
    ],
    loc="lower center",
    bbox_to_anchor=(column_centers[2], 0.005),
    ncol=1,
    fontsize=5.8,
    handlelength=1.8,
    labelspacing=0.3,
    title="c, f",
    title_fontsize=5.8,
)
fig.text(
    column_centers[2], 0.965, "Future-informed diagnostic", fontsize=5.8, ha="center"
)
finish(fig, axes, "Figure5", 3)
source_path = D / "coverage_summary.csv"
report = {
    "source_sha256": hashlib.sha256(source_path.read_bytes()).hexdigest(),
    "statistical_calculation": "None; read the saved paired contrast row directly.",
    "no_marginal_interval_subtraction": True,
    "no_new_bootstrap": True,
    "figsize_inches": list(FIGURE5_FIGSIZE),
    "canvas_mm": [145, 115],
    "panel_mm": [32, 32],
    "font_size_pt": 5.8,
    "source_data_unchanged": True,
    "contrasts": contrast_records,
    "other_figures_modified": False,
}
(Q / "figure5_annotation_check.json").write_text(
    json.dumps(report, indent=2, ensure_ascii=False), encoding="utf8"
)
print(json.dumps(report, indent=2, ensure_ascii=False))
