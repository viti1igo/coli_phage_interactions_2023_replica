"""Fig 1a from OUR regenerated core-genome tree (PanACoTA + FastTree): rectangular tree with rings for ST, O-type, H-type and phylogroup
taken from our regenerated feature table.  python scripts/fig1a_tree.py -> results/figures/fig1a_tree.*"""
import pandas as pd, numpy as np
import matplotlib.pyplot as plt
from Bio import Phylo
from figlib import *

TREE = f"{ROOT}/replication/T14_panacota/tree/ESCO.nucl.grp.aln.fasttree_tree.nwk"
lst = pd.read_csv(f"{ROOT}/replication/T14_panacota/annot/LSTINFO-list.lst", sep="\t").assign(strain=lambda d: d.orig_name.str.replace(r"\.fasta$", "", regex=True)).set_index("gembase_name").strain
feat = pd.read_csv(f"{ROOT}/replication/T25_features/picard_collection_regenerated.csv", sep=";", index_col=0)
t = Phylo.read(TREE, "newick")
for l in t.get_terminals(): l.name = lst.get(l.name, l.name)
t.root_at_midpoint()
leaves = t.get_terminals(); y = {l: i for i, l in enumerate(leaves)}
depth = t.depths()
fig, ax = plt.subplots(figsize=(9, 11))
def draw(c):
    if c.is_terminal(): return y[c]
    ys = [draw(k) for k in c.clades]
    for k, yy in zip(c.clades, ys): ax.plot([depth[c], depth[k]], [yy, yy], color=MUTED, lw=0.4)
    ax.plot([depth[c]] * 2, [min(ys), max(ys)], color=MUTED, lw=0.4)
    return float(np.mean(ys))
draw(t.root)
xmax = max(depth.values())
ring_cols = [("Clermont_Phylo", "phylogroup"), ("ST_Warwick", "ST"), ("O-type", "O-type"), ("H-type", "H-type")]
cm = plt.get_cmap("tab20"); rows = []
for j, (col, lab) in enumerate(ring_cols):
    s = feat[col].astype(str).reindex([l.name for l in leaves]); top = s.value_counts().index[:12]
    for i, v in enumerate(s):
        ax.add_patch(plt.Rectangle((xmax * (1.03 + 0.045 * j), i - 0.5), xmax * 0.04, 1, color=cm(list(top).index(v)) if v in top else "#e6e5e0", lw=0))
    ax.text(xmax * (1.05 + 0.045 * j), -4, lab, rotation=90, fontsize=8, color=INK, ha="center", va="bottom")
    rows += [dict(ring=lab, leaf=n, value=v) for n, v in s.items()]
ax.set_ylim(len(leaves), -40); ax.set_xlim(0, xmax * 1.25); ax.axis("off")
ax.set_title("Fig 1a  Our regenerated core-genome tree (403 strains) with typing rings\n(rings: top 12 values coloured, others grey)", loc="left", fontsize=10, color=INK)
save(fig, "fig1a_tree", pd.DataFrame(rows))
