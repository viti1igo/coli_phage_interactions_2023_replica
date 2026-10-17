"""Shared plotting helpers for the Phase 10 panels: one style, PNG+SVG+CSV per figure, one stats table (paper vs ours)."""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
OUT = f"{ROOT}/results/figures"
os.makedirs(OUT, exist_ok=True)
SURF, INK, MUTED, GRID, AXIS = "#fcfcfb", "#0b0b0b", "#52514e", "#e6e5e0", "#c9c8c2"
BLUE, ORANGE = "#2a78d6", "#eb6834"  # categorical slots 1-2 of the reference palette
plt.rcParams.update({"font.size": 9, "axes.titlesize": 10, "axes.titleweight": "regular", "svg.fonttype": "none"})


def style(ax, grid="y"):
    ax.set_facecolor(SURF)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(AXIS)
    ax.tick_params(colors=MUTED, labelsize=8)
    ax.xaxis.label.set_color(MUTED); ax.yaxis.label.set_color(MUTED)
    ax.title.set_color(INK)
    if grid:
        ax.grid(axis=grid, color=GRID, lw=0.6); ax.set_axisbelow(True)


def save(fig, name, data):
    """data: DataFrame behind the figure (same numbers as plotted)."""
    fig.patch.set_facecolor(SURF)
    for ext in ("png", "svg"):
        fig.savefig(f"{OUT}/{name}.{ext}", dpi=160, facecolor=SURF, bbox_inches="tight")
    data.to_csv(f"{OUT}/{name}.csv", index=False)
    plt.close(fig)
    print("wrote", name)


STATS = []


def stat(figure, metric, paper, ours, source):
    """Record one number for the comparison tables (paper = value printed in the article, source = where)."""
    STATS.append(dict(figure=figure, metric=metric, paper=paper, ours=ours, paper_source=source))


def flush_stats(tag):
    pd.DataFrame(STATS).to_csv(f"{OUT}/stats_{tag}.csv", index=False)
