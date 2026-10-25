"""Recompose fourteen supporting figures from frozen V23.2 evidence.

Four-panel figures use one horizontal row; six-panel figures use two rows;
all plot areas are square. No model fitting, member selection, resampling or
forecast scoring is performed. SI numbers and filenames retain their mapping.
"""
from pathlib import Path
import hashlib
import json
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Circle, Rectangle

sys.path.insert(0, str(Path(__file__).resolve().parent / 'helpers'))
from audit_panel_alignment import require_matplotlib_panel_alignment

import os
from plot_fonts import configure_font
ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT/'data/si_readability'
OUT = Path(os.environ['REPRO_OUTPUT'])/'si_readability'
QA = OUT
OUT.mkdir(parents=True,exist_ok=False)
WIDTH_MM = 180
COL = {'full': '#178F9F', 'early': '#F18A24', 'selected': '#F18A24',
       'path': '#9457BA', 'terminal': '#D0A018', 'random': '#777C86',
       'selected_recentered': '#7355B4'}
MARK = {'full': 's', 'early': 'D', 'selected': 'D', 'path': '^',
        'terminal': 'v', 'random': 'o', 'selected_recentered': 'o'}
LINE = {'full': '-', 'early': '-', 'selected': '-', 'path': '--',
        'terminal': '-.', 'random': '--', 'selected_recentered': '--'}
LABEL = {'full': 'Full ensemble', 'early': 'Early-selected half',
         'selected': 'Early-selected half', 'path': 'Hindsight fixed path half',
         'terminal': 'Hindsight endpoint half', 'random': 'Same-size random half',
         'selected_recentered': 'Retrospective center adjustment'}
LEAD_COLORS = ['#7142C4', '#BC4FA0', '#EC527A', '#F48B3A', '#CEA526', '#078C98']
RANK_COLORS = ['#EDAE28', '#EA7B72', '#D4539A', '#955BC4', '#55319B']
plt.rcParams.update({'font.family': configure_font(),
                     'font.size': 7.5, 'axes.labelsize': 7.5, 'axes.titlesize': 7.5,
                     'xtick.labelsize': 7.5, 'ytick.labelsize': 7.5,
                     'axes.spines.top': True, 'axes.spines.right': True,
                     'axes.linewidth': .65, 'axes.edgecolor': '#333333',
                     'font.weight': 'normal', 'axes.titleweight': 'normal',
                     'figure.facecolor': 'white', 'axes.facecolor': 'white',
                     'legend.frameon': False, 'legend.fontsize': 7.5,
                     'pdf.fonttype': 42, 'svg.fonttype': 'none'})
CENTERS = ['ECMWF', 'NCEP']
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
        side, gutter, left, row_gap = 32, 12, 13, 18
    elif cols == 3:
        side, gutter, left, row_gap = 43, 13, 20, 16
    else:
        side, gutter, left, row_gap = 61, 24, 20, 20
    height = top + bottom + rows * side + (rows - 1) * row_gap
    fig = plt.figure(figsize=(WIDTH_MM / 25.4, height / 25.4))
    axes = []
    for index in range(rows * cols):
        row, col = divmod(index, cols)
        rect = [(left + col * (side + gutter)) / WIDTH_MM,
                (bottom + (rows - row - 1) * (side + row_gap)) / height,
                side / WIDTH_MM, side / height]
        ax = fig.add_axes(rect)
        ax.tick_params(length=2.2, width=.55, pad=2)
        axes.append(ax)
    return fig, axes


def heading(ax, index, title):
    ax.text(0, 1.065, chr(97 + index), transform=ax.transAxes,
            ha='left', va='bottom', fontsize=8, fontweight='bold')
    ax.text(.11, 1.065, title, transform=ax.transAxes,
            ha='left', va='bottom', fontsize=7.5)


def key(method, label=None, line=None):
    return Line2D([], [], color=COL[method], marker=MARK[method],
                  ls=line or LINE[method], lw=1.1, ms=3,
                  label=label or LABEL[method])


def legend(fig, handles, columns=2, y_mm=2.5):
    height = fig.get_figheight() * 25.4
    fig.legend(handles=handles, loc='lower center', ncol=columns,
               bbox_to_anchor=(.52, y_mm / height), fontsize=7.5,
               columnspacing=1.25, handlelength=1.6, handletextpad=.5)


def save(fig, axes, filename, cols, si_number):
    ids = list('abcdefghijklmnop')[:len(axes)]
    fig.canvas.draw()
    require_matplotlib_panel_alignment(
        fig, axes=axes, panel_ids=ids,
        row_groups=[ids[i:i + cols] for i in range(0, len(ids), cols)],
        column_groups=[ids[i::cols] for i in range(cols) if len(ids[i::cols]) > 1],
        json_out=QA / (filename + '.alignment.json'),
        tolerance_pt=1.5, gutter_tolerance_pt=1.5, strict=True)
    fig.savefig(OUT / f'{filename}.pdf', facecolor='white')
    fig.savefig(OUT / f'{filename}.svg', facecolor='white')
    fig.savefig(OUT / f'{filename}.png', dpi=600, facecolor='white')
    dimensions = dict(si_number=si_number, filename=filename,
                      width_mm=fig.get_figwidth() * 25.4,
                      height_mm=fig.get_figheight() * 25.4,
                      plot_areas_mm=[dict(width=a.get_position().width * WIDTH_MM,
                                          height=a.get_position().height * fig.get_figheight() * 25.4)
                                     for a in axes])
    checks[filename] = dimensions
    (QA / (filename + '.dimensions.json')).write_text(json.dumps(dimensions, indent=2))
    plt.close(fig)


def interval(ax, row, y, method, scale=1, marker=None, open_marker=False):
    point, lower, upper = (float(row[k]) * scale for k in ['mean', 'ci_low', 'ci_high'])
    ax.errorbar(point, y, xerr=[[point - lower], [upper - point]],
                fmt=marker or MARK[method], color=COL[method],
                mfc='white' if open_marker else COL[method], ms=3,
                lw=.8, capsize=1.8, zorder=5)


def curve(ax, rows, method, x='lead'):
    rows = rows.sort_values(x)
    ax.fill_between(rows[x], rows.ci_low, rows.ci_high,
                    color=COL[method], alpha=.12, linewidth=0)
    ax.plot(rows[x], rows['mean'], color=COL[method],
            marker=MARK[method], ls=LINE[method], ms=3, lw=1.1)


def zero(ax):
    ax.axhline(0, color='#888888', lw=.55, zorder=0)


def horizontal_distribution(ax, values, sids, y, method):
    jitter = np.array([(int(hashlib.sha256(str(sid).encode()).hexdigest()[:8], 16)
                        % 1000 / 999 - .5) * .25 for sid in sids])
    ax.scatter(values, y + jitter, color=COL[method], alpha=.35,
               s=5, linewidth=0, rasterized=True, zorder=2)
    quartiles = np.quantile(values, [.25, .5, .75])
    ax.plot(quartiles[[0, 2]], [y, y], color=COL[method], alpha=.65,
            lw=2.8, solid_capstyle='butt')
    ax.plot(quartiles[1], y, '|', color='#222222', ms=7, mew=.7)



counts=read('count_curves.csv')
minima=read('count_minima.csv')
history=read('member_histories.csv.gz')
rank817=read('rank_817_summary.csv')
geometry_summary=read('geometry_summary.csv')
# S4: the fixed-168 cohort remains separate from the 806-case motion cohort.
fig, axes = layout(2, 3, top=8, bottom=18)
for i, ax in enumerate(axes):
    center, metric = CENTERS[i // 3], ['ATE', 'CTE', 'distance'][i % 3]
    for source_method, method in [('full', 'full'), ('selected', 'early')]:
        curve(ax, choose(geometry_summary, cohort='common168', center=center,
                         metric=source_method + '_' + metric), method)
    zero(ax)
    ax.set(xlabel='Lead (h)', ylabel=('Mean signed ' + metric if metric != 'distance' else 'Mean distance') + ' (km)',
           xticks=[48, 96, 144, 168])
    heading(ax, i, center + ' · ' + metric)
legend(fig, [key('full'), key('early')])
save(fig, axes, 'S08_extended_errors', 3, 4)

# S5–S7: every saved curve and retrospective minimum is retained.
for filename, cohort, number in [('S02_count168', 'common168', 5),
                                  ('S03_count_available', 'available', 6),
                                  ('S04_count144', 'common144', 7)]:
    fig, axes = layout(1, 4, top=8, bottom=20)
    for i, ax in enumerate(axes):
        center, mode = CENTERS[i // 2], ['raw', 'translated'][i % 2]
        data = choose(counts, center=center, cohort=cohort, mode=mode)
        for lead, color in zip(sorted(data.lead.unique()), LEAD_COLORS):
            part = choose(data, lead=lead).sort_values('fraction')
            minimum = choose(minima, center=center, cohort=cohort, mode=mode, lead=lead).iloc[0]
            assert np.isclose(part.ratio.min(), minimum.ratio)
            ax.plot(part.fraction, part.ratio, color=color, lw=1)
            ax.fill_between(part.fraction, part.lo, part.hi, color=color, alpha=.045, lw=0)
            ax.plot(minimum.fraction, minimum.ratio, 'v', color=color, ms=3.5,
                    mec='white', mew=.3, zorder=5)
        ax.axhline(1, color='#777777', lw=.55, ls=':')
        ax.set(xlim=(0, 102), xticks=[0, 50, 100], xlabel='Requested retention (%)', ylabel='Normalized position error')
        heading(ax, i, center + ' · ' + ('raw' if mode == 'raw' else 'aligned'))
    for pair in [axes[:2], axes[2:]]:
        lo = min(a.get_ylim()[0] for a in pair)
        hi = max(a.get_ylim()[1] for a in pair)
        for ax in pair:
            ax.set_ylim(lo, hi)
    legend(fig, [Line2D([], [], color=color, marker='v', ms=3, label=f'{lead} h')
                 for lead, color in zip(sorted(data.lead.unique()), LEAD_COLORS)], columns=6, y_mm=2.5)
    save(fig, axes, filename, 4, number)

