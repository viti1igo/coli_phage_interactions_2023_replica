"""Redrawn schematics (not data): Fig 2a-b (MLC scoring), Fig 6a (cocktail experiment), ED8 (four-step recommender), from the figure legends.
python scripts/fig_schematics.py -> results/figures/schematic_*.png/svg (+ a CSV holding the text of each box)"""
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
from figlib import *


def canvas(w, h, title):
    fig, ax = plt.subplots(figsize=(w, h)); ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.axis("off"); ax.set_title(title, loc="left", color=INK, fontsize=11)
    return fig, ax


def box(ax, x, y, w, h, text, fc="#bcd3f2", ec=BLUE, fs=9):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.4,rounding_size=1.5", fc=fc, ec=ec, lw=1.2))
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=fs, color=INK, wrap=True)


def arrow(ax, x0, y0, x1, y1):
    ax.annotate("", xy=(x1, y1), xytext=(x0, y0), arrowprops=dict(arrowstyle="->", color=MUTED, lw=1.3))


rows = []
# Fig 2a-b: plaque assay at three concentrations -> MLC
fig, ax = canvas(11, 4.2, "Fig 2a-b  Plaque assay and MLC score (redrawn from the legend)")
conc = [("5 x 10^8 pfu/ml\n(phage:bacteria ~ 10)", "MLC 1"), ("5 x 10^7 pfu/ml\n(~ 1)", "MLC 2"), ("5 x 10^6 pfu/ml\n(~ 0.1)", "MLC 3 (plaques) or\nMLC 4 (full clearing)")]
box(ax, 2, 70, 20, 22, "Phage x bacterial lawn\n3 replicates", fc="#f6f5f1", ec=AXIS)
for i, (c, m) in enumerate(conc):
    x = 30 + i * 23
    box(ax, x, 70, 20, 22, c, fc="#f6f5f1", ec=AXIS); box(ax, x, 22, 20, 28, "Lysis of the lawn at the\nLOWEST concentration\nwhere it is seen:\n" + m)
    arrow(ax, x + 10, 70, x + 10, 51); rows += [dict(figure="2a-b", box=c.replace("\n", " "), text=m)]
arrow(ax, 22, 81, 29, 81)
box(ax, 2, 22, 20, 28, "MLC 0: no lysis at\nany concentration", fc="#f6f5f1", ec=AXIS)
ax.text(2, 8, "Score = lowest phage concentration that lyses the lawn; a lower concentration means a higher MLC.", fontsize=9, color=MUTED)
save(fig, "schematic_fig2ab", pd.DataFrame(rows))

# Fig 6a: experiment
fig, ax = canvas(11, 3.6, "Fig 6a  Cocktail experiment (redrawn from the legend)")
steps = ["Interaction matrix\n403 strains x 96 phages", "Recommender trained on\nthe original matrix", "100 unseen E. coli\n(genomic traits only)", "3-phage cocktail\nper strain", "Tested at 6 concentrations\n(5x10^3 to 5x10^8 pfu/ml)\nvs baseline triplet"]
for i, t in enumerate(steps):
    box(ax, 1 + i * 20, 35, 17, 30, t, fs=8)
    if i: arrow(ax, 1 + i * 20 - 2.2, 50, 1 + i * 20 - 0.4, 50)
save(fig, "schematic_fig6a", pd.DataFrame({"step": range(1, 6), "text": [s.replace("\n", " ") for s in steps]}))

# ED8: four-step recommender
fig, ax = canvas(11, 4.4, "ED8  Four-step recommender (redrawn from the legend)")
st = [("1  Similar as the\nisolation strain", ""), ("2  Similar as strains\nin the Picard collection", ""), ("3  Trait matching", ""), ("4  Rescue", "")]  # box titles from the figure; the legend gives no further detail
for i, (t, d) in enumerate(st):
    box(ax, 2 + i * 24, 40, 21, 28, t, fs=9); 
    if i: arrow(ax, 2 + i * 24 - 2.6, 54, 2 + i * 24 - 0.6, 54)
ax.text(2, 10, "Rank phages by predicted probability, add top ones if p > 0.5; stop when k = 3; never two phages of the same genus and isolation strain.", fontsize=8.5, color=INK)
save(fig, "schematic_ed8", pd.DataFrame({"step": [s[0].replace("\n", " ") for s in st], "criterion": [s[1].replace("\n", " ") for s in st]}))
