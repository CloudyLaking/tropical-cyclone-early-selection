"""Render Figure 1: performance and retained-count curves from saved data."""

from pathlib import Path
import sys, json
import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

B = Path(__file__).resolve().parent.parent
D = B / "data"
O = B / "results"
perf = pd.read_csv(D / "performance_summary.csv")
counts = pd.read_csv(D / "count_curves.csv")
minimum = pd.read_csv(D / "count_minima.csv")
plt.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "DejaVu Sans"],
        "font.size": 7.5,
        "axes.labelsize": 7.5,
        "xtick.labelsize": 7.5,
        "ytick.labelsize": 7.5,
        "axes.linewidth": 0.65,
        "axes.edgecolor": "#262626",
        "axes.spines.top": True,
        "axes.spines.right": True,
        "pdf.fonttype": 42,
        "svg.fonttype": "none",
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "legend.frameon": False,
    }
)
cols = {"full": "#0095A6", "half": "#F18A24"}
markers = {"full": "s", "half": "D"}
leadcols = ["#7142C4", "#BC4FA0", "#EC527A", "#F48B3A", "#CEA526", "#078C98"]


def sel(t, **kw):
    m = np.ones(len(t), bool)
    for k, v in kw.items():
        m &= t[k] == v
    return t[m]


f = plt.figure(figsize=(170 / 25.4, 100 / 25.4))
axs = []
for r, c in enumerate(["ECMWF", "NCEP"]):
    for j in range(4):
        a = f.add_axes(
            [(12 + j * 42) / 170, (10 + (1 - r) * 44) / 100, 29 / 170, 29 / 100]
        )
        axs.append(a)
        letter = chr(97 + r * 4 + j)
        a.text(
            0,
            1.065,
            letter,
            transform=a.transAxes,
            weight="bold",
            fontsize=8.6,
            va="bottom",
        )
        title = c + " · " + ["error", "gain", "raw", "aligned"][j]
        a.text(0.125, 1.065, title, transform=a.transAxes, fontsize=7.5, va="bottom")
        a.tick_params(length=2.5, pad=1.8)
        if j == 0:
            handles = []
            for raw in [True, False]:
                for m, key in [("full", "full"), ("half", "selected")]:
                    t = sel(
                        perf,
                        center=c,
                        cohort="common168",
                        metric=("raw_" if raw else "translated_") + key,
                    ).sort_values("lead")
                    a.plot(
                        t.lead,
                        t.value,
                        color=cols[m],
                        ls="--" if raw else "-",
                        marker=markers[m],
                        mfc="white" if raw else cols[m],
                        ms=2.5,
                        lw=1.0,
                        mew=0.7,
                    )
                    handles.append(
                        Line2D(
                            [],
                            [],
                            color=cols[m],
                            ls="--" if raw else "-",
                            marker=markers[m],
                            mfc="white" if raw else cols[m],
                            ms=2.7,
                            lw=1,
                            label=("Raw " if raw else "Aligned ")
                            + ("full" if m == "full" else "half"),
                        )
                    )
            a.set(
                xlabel="Lead (h)",
                ylabel="Position error (km)",
                xticks=[48, 96, 144, 168],
                ylim=(0, None),
            )
            if r == 0:
                f.legend(
                    handles=handles,
                    title="a, e",
                    loc="upper center",
                    bbox_to_anchor=(26.5 / 170, 1.0),
                    ncol=2,
                    fontsize=7.5,
                    title_fontsize=7.5,
                    handlelength=1.2,
                    handletextpad=0.35,
                    columnspacing=0.6,
                    labelspacing=0.12,
                    borderaxespad=0,
                )
        elif j == 1:
            handles = []
            for metric, m, label in [
                ("translation_full_km", "full", "Alignment"),
                ("selection_after_translation_km", "half", "Added selection"),
            ]:
                t = sel(perf, center=c, cohort="common168", metric=metric).sort_values(
                    "lead"
                )
                a.fill_between(t.lead, t.lo, t.hi, color=cols[m], alpha=0.13, lw=0)
                a.plot(
                    t.lead, t.value, color=cols[m], marker=markers[m], ms=2.7, lw=1.1
                )
                handles.append(
                    Line2D(
                        [],
                        [],
                        color=cols[m],
                        marker=markers[m],
                        ms=2.7,
                        lw=1,
                        label=label,
                    )
                )
            a.axhline(0, color="#888888", lw=0.55, zorder=0)
            a.set(
                xlabel="Lead (h)",
                ylabel="Error reduction (km)",
                xticks=[48, 96, 144, 168],
            )
            if r == 0:
                f.legend(
                    handles=handles,
                    title="b, f",
                    loc="upper center",
                    bbox_to_anchor=(68.5 / 170, 1.0),
                    ncol=1,
                    fontsize=7.5,
                    title_fontsize=7.5,
                    handlelength=1.2,
                    handletextpad=0.35,
                    labelspacing=0.12,
                    borderaxespad=0,
                )
        else:
            mode = ["raw", "translated"][j - 2]
            for lead, col in zip([48, 72, 96, 120, 144, 168], leadcols):
                t = sel(
                    counts, center=c, cohort="common168", mode=mode, lead=lead
                ).sort_values("fraction")
                q = sel(
                    minimum, center=c, cohort="common168", mode=mode, lead=lead
                ).iloc[0]
                a.plot(t.fraction, t.ratio, color=col, lw=1)
                a.fill_between(t.fraction, t.lo, t.hi, color=col, alpha=0.035, lw=0)
                a.plot(
                    q.fraction,
                    q.ratio,
                    "v",
                    ms=3.1,
                    color=col,
                    mec="white",
                    mew=0.3,
                    zorder=5,
                )
            a.axhline(1, color="#888888", ls=":", lw=0.55)
            a.set(
                xlim=(0, 102),
                xticks=[0, 50, 100],
                xlabel="Retained members (%)",
                ylabel="Normalized error",
            )
handles = [
    Line2D([], [], color=leadcols[k], lw=1.1, label=str([48, 72, 96, 120, 144, 168][k]))
    for k in [0, 3, 1, 4, 2, 5]
]
leg = f.legend(
    handles=handles,
    title="c, d, g, h: lead (h)",
    loc="upper center",
    bbox_to_anchor=(131.5 / 170, 1.0),
    ncol=3,
    fontsize=7.5,
    title_fontsize=7.5,
    handlelength=1.2,
    columnspacing=0.8,
    labelspacing=0.12,
    borderaxespad=0,
)
f.add_artist(
    Line2D(
        [96 / 170, 167 / 170],
        [0.895, 0.895],
        transform=f.transFigure,
        color="#B6B6B6",
        lw=0.5,
    )
)
for x in [96, 167]:
    f.add_artist(
        Line2D(
            [x / 170, x / 170],
            [0.895, 0.89],
            transform=f.transFigure,
            color="#B6B6B6",
            lw=0.5,
        )
    )
f.canvas.draw()
f.savefig(O / "Figure1.pdf", dpi=600)
f.savefig(O / "Figure1.svg", dpi=600)
f.savefig(O / "Figure1.png", dpi=600)
(O / "Figure1.dimensions.json").write_text(
    json.dumps(
        {
            "width_mm": 170,
            "height_mm": 100,
            "axes_mm": [[29, 29]] * 8,
            "manuscript_width_mm": 160,
            "placed_height_mm": 100 * 160 / 170,
            "minimum_source_font_pt": 7.5,
            "minimum_placed_font_pt": 7.5 * 160 / 170,
            "legend_scope": {
                "a,e": "four explicit raw/aligned by full/half combinations",
                "b,f": "alignment and added selection; selection is incremental after alignment",
                "c,d,g,h": "forecast lead shared only across count columns",
            },
        },
        indent=2,
    ),
    encoding="utf8",
)
plt.close(f)
print("Figure 1: 170 x 100 mm source; 160 x 94.1 mm in manuscript; local legends.")
