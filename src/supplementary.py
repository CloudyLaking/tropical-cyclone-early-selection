"""Render supporting Figures S1-S14 from frozen summaries; no new scoring."""

from pathlib import Path
import hashlib
import json
import sys
import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Circle, Rectangle

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
OUT = ROOT / "results"
QA = ROOT / "results/checks"
OUT.mkdir(parents=True, exist_ok=True)
QA.mkdir(parents=True, exist_ok=True)
WIDTH_MM = 180
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
    "selected_recentered": "o",
}
LINE = {
    "full": "-",
    "early": "-",
    "selected": "-",
    "path": "--",
    "terminal": "-.",
    "random": "--",
    "selected_recentered": "--",
}
LABEL = {
    "full": "Full ensemble",
    "early": "Early-selected half",
    "selected": "Early-selected half",
    "path": "Hindsight fixed path half",
    "terminal": "Hindsight endpoint half",
    "random": "Same-size random half",
    "selected_recentered": "Retrospective center adjustment",
}
LEAD_COLORS = ["#7142C4", "#BC4FA0", "#EC527A", "#F48B3A", "#CEA526", "#078C98"]
RANK_COLORS = ["#EDAE28", "#EA7B72", "#D4539A", "#955BC4", "#55319B"]
plt.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "DejaVu Sans"],
        "font.size": 6.5,
        "axes.labelsize": 6.5,
        "axes.titlesize": 6.8,
        "xtick.labelsize": 6,
        "ytick.labelsize": 6,
        "axes.spines.top": True,
        "axes.spines.right": True,
        "axes.linewidth": 0.65,
        "axes.edgecolor": "#333333",
        "font.weight": "normal",
        "axes.titleweight": "normal",
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "legend.frameon": False,
        "legend.fontsize": 6.2,
        "pdf.fonttype": 42,
        "svg.fonttype": "none",
    }
)
CENTERS = ["ECMWF", "NCEP"]
inputs = {}
checks = {}


def read(name):
    path = DATA / name
    inputs[name] = hashlib.sha256(path.read_bytes()).hexdigest()
    return pd.read_csv(path)


def choose(frame, **conditions):
    keep = np.ones(len(frame), dtype=bool)
    for key, value in conditions.items():
        keep &= frame[key] == value
    return frame[keep]


def layout(rows, cols, top=8, bottom=20):
    if cols == 4:
        side, gutter, left, row_gap = (32, 12, 13, 18)
    elif cols == 3:
        side, gutter, left, row_gap = (43, 13, 20, 16)
    else:
        side, gutter, left, row_gap = (61, 24, 20, 20)
    height = top + bottom + rows * side + (rows - 1) * row_gap
    fig = plt.figure(figsize=(WIDTH_MM / 25.4, height / 25.4))
    axes = []
    for index in range(rows * cols):
        row, col = divmod(index, cols)
        rect = [
            (left + col * (side + gutter)) / WIDTH_MM,
            (bottom + (rows - row - 1) * (side + row_gap)) / height,
            side / WIDTH_MM,
            side / height,
        ]
        ax = fig.add_axes(rect)
        ax.tick_params(length=2.2, width=0.55, pad=2)
        axes.append(ax)
    return (fig, axes)


def heading(ax, index, title):
    ax.text(
        0,
        1.065,
        chr(97 + index),
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=8,
        fontweight="bold",
    )
    ax.text(
        0.11, 1.065, title, transform=ax.transAxes, ha="left", va="bottom", fontsize=6.8
    )


def key(method, label=None, line=None):
    return Line2D(
        [],
        [],
        color=COL[method],
        marker=MARK[method],
        ls=line or LINE[method],
        lw=1.1,
        ms=3,
        label=label or LABEL[method],
    )


def legend(fig, handles, columns=2, y_mm=2.5):
    height = fig.get_figheight() * 25.4
    fig.legend(
        handles=handles,
        loc="lower center",
        ncol=columns,
        bbox_to_anchor=(0.52, y_mm / height),
        fontsize=6.2,
        columnspacing=1.25,
        handlelength=1.6,
        handletextpad=0.5,
    )


def save(fig, axes, filename, cols, si_number):
    filename = f"FigureS{si_number:02d}"
    ids = list("abcdefghijklmnop")[: len(axes)]
    fig.canvas.draw()
    fig.savefig(OUT / f"{filename}.pdf", facecolor="white")
    fig.savefig(OUT / f"{filename}.svg", facecolor="white")
    fig.savefig(OUT / f"{filename}.png", dpi=600, facecolor="white")
    dimensions = dict(
        si_number=si_number,
        filename=filename,
        width_mm=fig.get_figwidth() * 25.4,
        height_mm=fig.get_figheight() * 25.4,
        plot_areas_mm=[
            dict(
                width=a.get_position().width * WIDTH_MM,
                height=a.get_position().height * fig.get_figheight() * 25.4,
            )
            for a in axes
        ],
    )
    checks[filename] = dimensions
    (QA / (filename + ".dimensions.json")).write_text(json.dumps(dimensions, indent=2))
    plt.close(fig)


def interval(ax, row, y, method, scale=1, marker=None, open_marker=False):
    point, lower, upper = (float(row[k]) * scale for k in ["mean", "ci_low", "ci_high"])
    ax.errorbar(
        point,
        y,
        xerr=[[point - lower], [upper - point]],
        fmt=marker or MARK[method],
        color=COL[method],
        mfc="white" if open_marker else COL[method],
        ms=3,
        lw=0.8,
        capsize=1.8,
        zorder=5,
    )


def curve(ax, rows, method, x="lead"):
    rows = rows.sort_values(x)
    ax.fill_between(
        rows[x], rows.ci_low, rows.ci_high, color=COL[method], alpha=0.12, linewidth=0
    )
    ax.plot(
        rows[x],
        rows["mean"],
        color=COL[method],
        marker=MARK[method],
        ls=LINE[method],
        ms=3,
        lw=1.1,
    )


def zero(ax):
    ax.axhline(0, color="#888888", lw=0.55, zorder=0)


def horizontal_distribution(ax, values, sids, y, method):
    jitter = np.array(
        [
            (
                int(hashlib.sha256(str(sid).encode()).hexdigest()[:8], 16) % 1000 / 999
                - 0.5
            )
            * 0.25
            for sid in sids
        ]
    )
    ax.scatter(
        values,
        y + jitter,
        color=COL[method],
        alpha=0.35,
        s=5,
        linewidth=0,
        rasterized=True,
        zorder=2,
    )
    quartiles = np.quantile(values, [0.25, 0.5, 0.75])
    ax.plot(
        quartiles[[0, 2]],
        [y, y],
        color=COL[method],
        alpha=0.65,
        lw=2.8,
        solid_capstyle="butt",
    )
    ax.plot(quartiles[1], y, "|", color="#222222", ms=7, mew=0.7)


fixed = read("fixed_summary.csv")
storm = read("fixed_storms.csv")
geometry = read("geometry_cases.csv")
crossings = read("coverage_crossings.csv")
lengths = read("storm_geometry_summary.csv")
counts = read("count_curves.csv")
minima = read("count_minima.csv")
history = read("member_histories.csv.gz")
rank817 = read("rank_817_summary.csv")
geometry_summary = read("geometry_summary.csv")
coverage = read("coverage_summary.csv")
fractions = read("fraction_summary.csv")
performance = read("performance_summary.csv").rename(
    columns={"value": "mean", "lo": "ci_low", "hi": "ci_high"}
)
fig, axes = layout(1, 4, top=8, bottom=15)
for index, ax in enumerate(axes):
    center = CENTERS[index % 2]
    if index < 2:
        data = geometry[geometry.center == center]
        ax.add_patch(Circle((1, 0), 1, facecolor="white", edgecolor="#A9C4E8", lw=0.65))
        ax.axhline(0, color="#DDDDDD", lw=0.5)
        ax.scatter(
            data.circle_x,
            data.circle_y,
            color=COL["early"],
            s=4,
            alpha=0.65,
            linewidth=0,
            rasterized=True,
        )
        ax.plot(0, 0, "s", color=COL["full"], ms=3.5)
        ax.plot(1, 0, "+", color="#252A30", ms=5)
        xmin, xmax = (
            min(-0.7, data.circle_x.min() - 0.1),
            max(2.3, data.circle_x.max() + 0.1),
        )
        ymin, ymax = (
            min(-1.5, data.circle_y.min() - 0.1),
            max(1.5, data.circle_y.max() + 0.1),
        )
        span = max(xmax - xmin, ymax - ymin)
        xm, ym = ((xmax + xmin) / 2, (ymax + ymin) / 2)
        ax.set(
            xlim=(xm - span / 2, xm + span / 2),
            ylim=(ym - span / 2, ym + span / 2),
            xlabel="Toward-truth shift / full error",
            ylabel="Perpendicular shift / full error",
        )
        heading(ax, index, center + " · shifts")
        assert len(data) == (300 if center == "ECMWF" else 506)
    else:
        data = crossings[crossings.center == center]
        maximum = max(data.x.max(), data.y.max()) * 1.04
        ax.add_patch(
            Rectangle((1, 0), maximum - 1, 1, facecolor="#E9E3F5", alpha=0.18, lw=0)
        )
        ax.add_patch(
            Rectangle((0, 1), 1, maximum - 1, facecolor="#FAECCD", alpha=0.18, lw=0)
        )
        ax.plot([0, maximum], [0, maximum], color="#D9D9D9", lw=0.55)
        ax.axvline(1, color="#777777", lw=0.6)
        ax.axhline(1, color="#777777", lw=0.6)
        for state, color, marker in [
            ("outside", "#CDD4DB", "o"),
            ("inside", "#ADC6E7", "o"),
            ("enter", "#6554A6", "^"),
            ("exit", "#BD8A22", "v"),
        ]:
            part = data[data.state == state]
            ax.scatter(
                part.x,
                part.y,
                color=color,
                marker=marker,
                s=5,
                alpha=0.6,
                linewidth=0,
                rasterized=True,
            )
        ax.set(
            xlim=(0, maximum),
            ylim=(0, maximum),
            xlabel="Before: distance / radius",
            ylabel="After: distance / radius",
        )
        heading(ax, index, center + " · coverage")
        assert len(data) == (300 if center == "ECMWF" else 506)
    ax.set_aspect("equal", adjustable="box")
save(fig, axes, "S01_full_geometry", 4, 1)
fig, axes = layout(1, 2, top=8, bottom=14)
quantities = ["mean_case_shift", "mean_storm_shift", "overall_mean_shift"]
descriptions = [
    "Each initialization",
    "Within each storm",
    "Overall target-frame vector",
]
for i, (center, ax) in enumerate(zip(CENTERS, axes)):
    for j, (quantity, label) in enumerate(zip(quantities, descriptions)):
        row = (
            choose(lengths, cohort="motion806", center=center, quantity=quantity)
            .iloc[0]
            .copy()
        )
        row["mean"] = row.value
        interval(
            ax, row, 2 - j, "early", marker="o" if i == 0 else "s", open_marker=i == 1
        )
        ax.text(3, 2 - j + 0.27, label, fontsize=6.4, va="center")
    ax.set(
        xlim=(0, 105),
        ylim=(-0.45, 2.65),
        xticks=[0, 25, 50, 75, 100],
        yticks=[],
        xlabel="Mean-track shift length (km)",
    )
    heading(ax, i, center + " · averaging levels")
save(fig, axes, "S14_shift_lengths", 2, 2)
fig, axes = layout(2, 3, top=8, bottom=25)
for i, ax in enumerate(axes):
    center, column = (CENTERS[i // 3], i % 3)
    if column == 0:
        for raw in [True, False]:
            for method, suffix in [("full", "full"), ("early", "selected")]:
                data = choose(
                    performance,
                    center=center,
                    cohort="common168",
                    metric=("raw_" if raw else "translated_") + suffix,
                ).sort_values("lead")
                ax.plot(
                    data.lead,
                    data["mean"],
                    color=COL[method],
                    ls="--" if raw else "-",
                    marker=MARK[method],
                    mfc="white" if raw else COL[method],
                    ms=3,
                    lw=1.1,
                )
        ax.set(ylabel="Great-circle error (km)")
        title = "track error"
    elif column == 1:
        for metric, method in [
            ("translation_full_km", "full"),
            ("selection_after_translation_km", "early"),
        ]:
            curve(
                ax,
                choose(performance, center=center, cohort="common168", metric=metric),
                method,
            )
        zero(ax)
        ax.set(ylabel="Paired error reduction (km)")
        title = "paired gains"
    else:
        for cohort, color, marker in [
            ("available", "#3269B2", "o"),
            ("common144", "#7142C4", "s"),
            ("common168", "#D97706", "^"),
        ]:
            data = choose(
                performance,
                center=center,
                cohort=cohort,
                metric="selection_after_translation_km",
            ).sort_values("lead")
            ax.errorbar(
                data.lead,
                data["mean"],
                yerr=[data["mean"] - data.ci_low, data.ci_high - data["mean"]],
                color=color,
                marker=marker,
                lw=1,
                ms=3,
                capsize=1.6,
            )
        zero(ax)
        ax.set(ylabel="Added-selection gain (km)")
        title = "sample comparison"
    ax.set(xlabel="Lead (h)", xticks=[48, 96, 144, 168])
    heading(ax, i, center + " · " + title)
legend(
    fig,
    [
        key("full"),
        key("early"),
        Line2D([], [], color="#555555", ls="--", marker="o", mfc="white", label="Raw"),
        Line2D([], [], color="#555555", ls="-", marker="o", label="24-h aligned"),
    ],
    columns=4,
    y_mm=11,
)
legend(
    fig,
    [
        Line2D([], [], color=color, marker=marker, label=label)
        for color, marker, label in [
            ("#3269B2", "o", "Available at lead"),
            ("#7142C4", "s", "Fixed 144 h"),
            ("#D97706", "^", "Fixed 168 h"),
        ]
    ],
    columns=3,
    y_mm=2,
)
save(fig, axes, "S13_performance", 3, 3)
fig, axes = layout(2, 3, top=8, bottom=18)
for i, ax in enumerate(axes):
    center, metric = (CENTERS[i // 3], ["ATE", "CTE", "distance"][i % 3])
    for source_method, method in [("full", "full"), ("selected", "early")]:
        curve(
            ax,
            choose(
                geometry_summary,
                cohort="common168",
                center=center,
                metric=source_method + "_" + metric,
            ),
            method,
        )
    zero(ax)
    ax.set(
        xlabel="Lead (h)",
        ylabel=("Mean signed " + metric if metric != "distance" else "Mean distance")
        + " (km)",
        xticks=[48, 96, 144, 168],
    )
    heading(ax, i, center + " · " + metric)
legend(fig, [key("full"), key("early")])
save(fig, axes, "S08_extended_errors", 3, 4)
for filename, cohort, number in [
    ("S02_count168", "common168", 5),
    ("S03_count_available", "available", 6),
    ("S04_count144", "common144", 7),
]:
    fig, axes = layout(1, 4, top=8, bottom=20)
    for i, ax in enumerate(axes):
        center, mode = (CENTERS[i // 2], ["raw", "translated"][i % 2])
        data = choose(counts, center=center, cohort=cohort, mode=mode)
        for lead, color in zip(sorted(data.lead.unique()), LEAD_COLORS):
            part = choose(data, lead=lead).sort_values("fraction")
            minimum = choose(
                minima, center=center, cohort=cohort, mode=mode, lead=lead
            ).iloc[0]
            assert np.isclose(part.ratio.min(), minimum.ratio)
            ax.plot(part.fraction, part.ratio, color=color, lw=1)
            ax.fill_between(
                part.fraction, part.lo, part.hi, color=color, alpha=0.045, lw=0
            )
            ax.plot(
                minimum.fraction,
                minimum.ratio,
                "v",
                color=color,
                ms=3.5,
                mec="white",
                mew=0.3,
                zorder=5,
            )
        ax.axhline(1, color="#777777", lw=0.55, ls=":")
        ax.set(
            xlim=(0, 102),
            xticks=[0, 50, 100],
            xlabel="Requested retention (%)",
            ylabel="Normalized position error",
        )
        heading(ax, i, center + " · " + ("raw" if mode == "raw" else "aligned"))
    for pair in [axes[:2], axes[2:]]:
        lo = min((a.get_ylim()[0] for a in pair))
        hi = max((a.get_ylim()[1] for a in pair))
        for ax in pair:
            ax.set_ylim(lo, hi)
    legend(
        fig,
        [
            Line2D([], [], color=color, marker="v", ms=3, label=f"{lead} h")
            for lead, color in zip(sorted(data.lead.unique()), LEAD_COLORS)
        ],
        columns=6,
        y_mm=2.5,
    )
    save(fig, axes, filename, 4, number)
fig, axes = layout(1, 4, top=16, bottom=20)
for i, ax in enumerate(axes):
    center = CENTERS[i % 2]
    if i < 2:
        date = ["2015-05-12T12:00:00", "2015-10-18T00:00:00"][i]
        data = choose(history, center=center, base=date, scoring="raw")
        assert not data.empty
        for member, part in data.groupby("member"):
            part = part.sort_values("lead")
            color = RANK_COLORS[min(4, int(part.early_fraction.iloc[0] * 5))]
            ax.plot(
                part.lead,
                100 * part.rank_fraction,
                color=color,
                lw=0.6,
                ls="-" if part.selected.iloc[0] else "--",
            )
        ax.set(
            xlim=(24, 120),
            ylim=(100, 0),
            xticks=[24, 48, 72, 96, 120],
            yticks=[0, 50, 100],
            xlabel="Lead (h)",
            ylabel="Raw member rank (%)",
        )
        ax.text(0.5, -0.32, date[:10], transform=ax.transAxes, fontsize=6, ha="center")
        title = "member ranks"
    else:
        for metric, method, marker in [
            ("best_stays", "full", "o"),
            ("worst_stays", "early", "s"),
        ]:
            data = choose(rank817, center=center, metric=metric).sort_values(
                "total_lead"
            )
            ax.errorbar(
                data.total_lead,
                100 * data.estimate,
                yerr=[
                    100 * (data.estimate - data.low),
                    100 * (data.high - data.estimate),
                ],
                color=COL[method],
                marker=marker,
                ms=3,
                lw=1,
                capsize=1.5,
            )
        ax.axhline(20, color="#888888", ls=":", lw=0.55)
        ax.set(
            xlabel="Lead (h)",
            ylabel="Same-fifth retention (%)",
            xticks=[48, 72, 96, 120],
        )
        title = "817-cohort groups"
    heading(ax, i, center + " · " + title)
legend(
    fig,
    [
        Line2D([], [], color=color, lw=1, label=label)
        for color, label in zip(RANK_COLORS, ["1 (best)", "2", "3", "4", "5 (worst)"])
    ],
    columns=5,
    y_mm=fig.get_figheight() * 25.4 - 7,
)
legend(
    fig,
    [
        Line2D([], [], color=COL["full"], marker="o", lw=1, ms=3, label="Best → best"),
        Line2D(
            [], [], color=COL["early"], marker="s", lw=1, ms=3, label="Worst → worst"
        ),
    ],
    columns=2,
    y_mm=1,
)
save(fig, axes, "S05_member_ranks", 4, 8)
fig, axes = layout(2, 3, top=8, bottom=21)
for i, ax in enumerate(axes):
    center, metric = (CENTERS[i // 3], ["ATE", "CTE", "distance"][i % 3])
    for method in ["full", "early", "path", "terminal"]:
        curve(ax, choose(fixed, center=center, metric=metric, method=method), method)
    zero(ax)
    ax.set(
        xlabel="Lead (h)",
        ylabel=(
            "Mean signed " + metric
            if metric != "distance"
            else "Mean projected distance"
        )
        + " (km)",
        xticks=[24, 48, 72, 96, 120],
    )
    heading(ax, i, center + " · " + metric)
legend(
    fig, [key(m) for m in ["full", "early", "path", "terminal"]], columns=2, y_mm=1.5
)
save(fig, axes, "S06_fixed_histories", 3, 9)
fig, axes = layout(2, 3, top=8, bottom=21)
for i, ax in enumerate(axes):
    center, metric = (CENTERS[i // 3], ["ATE", "CTE", "distance"][i % 3])
    for j, method in enumerate(["early", "path", "terminal"]):
        data = choose(
            storm, center=center, lead=120, method=method + "_minus_full"
        ).sort_values("sid")
        assert len(data) == (64 if center == "ECMWF" else 100)
        horizontal_distribution(ax, data[metric].to_numpy(), data.sid, 2 - j, method)
        summary = choose(
            fixed, center=center, lead=120, method=method + "_minus_full", metric=metric
        ).iloc[0]
        interval(ax, summary, 2 - j + 0.25, method, open_marker=method != "early")
    ax.axvline(0, color="#777777", lw=0.6)
    ax.set(
        yticks=[0, 1, 2],
        yticklabels=["Endpoint", "Path", "Early"],
        ylim=(-0.5, 2.6),
        xlabel="Method − full " + metric + " (km)",
    )
    heading(ax, i, center + " · paired " + metric)
legend(fig, [key(m) for m in ["early", "path", "terminal"]], columns=3, y_mm=2)
save(fig, axes, "S07_paired_changes", 3, 10)
fig, axes = layout(1, 4, top=8, bottom=19)
for i, ax in enumerate(axes):
    center, metric = (CENTERS[i // 2], ["ATE", "CTE"][i % 2])
    for lead, method, marker in [(48, "full", "o"), (120, "early", "s")]:
        data = choose(fractions, center=center, lead=lead, metric=metric).sort_values(
            "fraction"
        )
        ax.errorbar(
            data.fraction,
            data["mean"],
            yerr=[data["mean"] - data.ci_low, data.ci_high - data["mean"]],
            color=COL[method],
            marker=marker,
            ms=2.7,
            lw=0.9,
            capsize=1.5,
        )
    zero(ax)
    ax.set(
        xlabel="Requested retention (%)",
        ylabel="Mean signed " + metric + " (km)",
        xticks=[10, 50, 100],
    )
    heading(ax, i, center + " · " + metric)
legend(
    fig,
    [
        Line2D([], [], color=COL["full"], marker="o", lw=1, ms=3, label="48 h"),
        Line2D([], [], color=COL["early"], marker="s", lw=1, ms=3, label="120 h"),
    ],
    columns=2,
    y_mm=2,
)
save(fig, axes, "S12_fraction_components", 4, 11)
for filename, cohort, number in [
    ("S09_coverage120", "motion806", 12),
    ("S10_coverage168", "common168", 13),
]:
    fig, axes = layout(2, 3, top=8, bottom=21)
    for i, ax in enumerate(axes):
        center, metric = (CENTERS[i // 3], ["r70", "distance", "in70"][i % 3])
        for method in ["full", "selected", "random", "selected_recentered"]:
            data = choose(
                coverage, cohort=cohort, center=center, comparison=method, metric=metric
            ).copy()
            if metric == "in70":
                data[["mean", "ci_low", "ci_high"]] *= 100
            curve(ax, data, method)
        titles = {
            "r70": "member radius",
            "distance": "center error",
            "in70": "coverage",
        }
        labels = {
            "r70": "Member radius (km)",
            "distance": "Center distance (km)",
            "in70": "Empirical coverage (%)",
        }
        ax.set(
            xlabel="Lead (h)",
            ylabel=labels[metric],
            xticks=[48, 72, 96, 120] if cohort == "motion806" else [48, 96, 144, 168],
        )
        heading(ax, i, center + " · " + titles[metric])
    legend(
        fig,
        [key(m) for m in ["full", "selected", "random", "selected_recentered"]],
        columns=2,
        y_mm=1.5,
    )
    save(fig, axes, filename, 3, number)
fig, axes = layout(1, 2, top=8, bottom=21)
tail_colors = ["#7A5AF8", "#DDE3E9", "#F4C430"]
for i, (center, ax) in enumerate(zip(CENTERS, axes)):
    for y, method in [(1, "selected"), (0, "selected_recentered")]:
        left = 0
        for metric, color in zip(["behind10", "in_ATE80", "ahead90"], tail_colors):
            value = (
                100
                * choose(
                    coverage,
                    cohort="motion806",
                    center=center,
                    lead=120,
                    comparison=method,
                    metric=metric,
                )["mean"].iloc[0]
            )
            ax.barh(
                y, value, left=left, height=0.35, color=color, edgecolor="white", lw=0.5
            )
            ax.text(
                left + value / 2,
                y + 0.3,
                f"{value:.1f}",
                ha="center",
                va="center",
                fontsize=6.5,
            )
            left += value
        assert np.isclose(left, 100)
    ax.set(
        xlim=(-4, 104),
        ylim=(-0.6, 1.6),
        xticks=[0, 25, 50, 75, 100],
        yticks=[0, 1],
        yticklabels=["Adjusted", "Selected"],
        xlabel="Storm-weighted cases (%)",
    )
    heading(ax, i, center + " · along-track range")
legend(
    fig,
    [
        Line2D([], [], color=color, lw=4, label=label)
        for color, label in zip(
            tail_colors, ["Truth behind range", "Inside range", "Truth ahead of range"]
        )
    ],
    columns=3,
    y_mm=2,
)
save(fig, axes, "S11_along_track_tails", 2, 14)
assert len(checks) == 14
report = dict(
    status="PASS",
    forecast_scores_recomputed=0,
    new_models=0,
    new_member_selection=0,
    new_bootstrap_draws=0,
    source_sha256=inputs,
    shared_style="white; four borders; only panel letters bold; square plot rectangles",
    figures=checks,
)
(QA / "support_figures_data_and_layout.json").write_text(
    json.dumps(report, indent=2), encoding="utf-8"
)
print("Rebuilt all fourteen SI figures; main figures and S15 were not written.")
