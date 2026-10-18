"""Fig 3b and ED3: genome maps of phages of one genus, protein-level links between neighbouring genomes and, for Fig 3b, the position-wise entropy of a
whole-genome multiple alignment (MAFFT).  The authors' annotation (Prokka) is not shipped: genes are called again with pyrodigal, and the RBP genes of
data/genomics/phages/RBP/RBP_list.csv are found by the number in their identifier (the n-th gene of the genome, as Prokka numbers them).
python scripts/fig_phage_maps.py  (needs mafft on PATH, or MAFFT=/path/to/mafft)  -> results/figures/fig3b_*.{png,svg,csv}, ed3_*.{png,svg,csv}"""
import os, subprocess, tempfile, warnings
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrow
import pyrodigal
from Bio import Align
from Bio.Seq import Seq
import phage_genomes as G
from figlib import *
from matrix import load_matrix
import supp

warnings.filterwarnings("ignore")  # pyrodigal warns for every genome shorter than 100 kb
RBP = pd.read_csv(f"{ROOT}/data/genomics/phages/RBP/RBP_list.csv", sep=";")
RBP = RBP[RBP.RBP.notna()]
COL = pd.read_csv(f"{ROOT}/data/genomics/phages/guelin_collection.csv", sep=";", index_col=0)
aligner = Align.PairwiseAligner(mode="global", substitution_matrix=Align.substitution_matrices.load("BLOSUM62"), open_gap_score=-10, extend_gap_score=-0.5)
_cache = {}


def genes(ph):
    """Genes as (start, end, strand, protein); single mode, or metagenomic mode when an RBP identifier exceeds the gene count (Prokka calls a few more genes)."""
    if ph in _cache: return _cache[ph]
    seq = G.read_genome(ph).encode()
    f = pyrodigal.GeneFinder(meta=False); f.train(seq, translation_table=11); gs = f.find_genes(seq)
    need = max([int(x.rsplit("_", 1)[1]) for x in RBP[RBP.phage == ph].RBP] or [0])
    if need > len(gs):
        gs = pyrodigal.GeneFinder(meta=True).find_genes(seq)
    out = [(g.begin, g.end, g.strand, g.translate().rstrip("*")) for g in gs]
    rbp = {int(x.rsplit("_", 1)[1]) for x in RBP[RBP.phage == ph].RBP}
    _cache[ph] = (len(seq), out, rbp)
    return _cache[ph]


def identity(a, b):
    aln = aligner.align(a, b)[0]
    cols = aln.length
    ident = sum(1 for (s1, e1), (s2, e2) in zip(*aln.aligned) for i, j in zip(range(s1, e1), range(s2, e2)) if a[i] == b[j])
    return ident / cols


def links(p1, p2, thr=0.3):
    """Reciprocal best protein matches between two genomes with identity above thr (clinker-like)."""
    g1, g2 = genes(p1)[1], genes(p2)[1]
    S = np.zeros((len(g1), len(g2)))
    for i, a in enumerate(g1):
        for j, b in enumerate(g2):
            if 0.6 < len(a[3]) / max(len(b[3]), 1) < 1 / 0.6:
                S[i, j] = identity(a[3], b[3])
    out = []
    for i in range(len(g1)):
        j = S[i].argmax()
        if S[i, j] >= thr and S[:, j].argmax() == i: out.append((i, j, S[i, j]))
    return out


def draw_map(ax, phages, y0=0, thr=0.3, label_hosts=True, color_map=None):
    rows = []
    L = max(genes(p)[0] for p in phages)
    for k, p in enumerate(phages):
        n, gs, rbp = genes(p); y = y0 - k * 1.6
        ax.plot([0, n], [y, y], color=AXIS, lw=0.8, zorder=1)
        for i, (b, e, s, prot) in enumerate(gs, 1):
            c = "#8fbcf2" if i in rbp else "#c9c8c2"
            ax.add_patch(FancyArrow(b if s == 1 else e, y, (e - b) * s, 0, width=0.45, head_width=0.45, head_length=min(300, (e - b) * 0.4), length_includes_head=True, color=c, lw=0, zorder=2))
            rows.append(dict(phage=p, gene=i, start=b, end=e, strand=s, aa=len(prot), rbp=i in rbp))
        ax.text(-L * 0.01, y, f"{p}" + (f" ({COL.loc[p, 'Phage_host']})" if label_hosts else ""), ha="right", va="center", fontsize=8, color=INK)
    lk = []
    for k in range(len(phages) - 1):
        p1, p2 = phages[k], phages[k + 1]
        g1, g2 = genes(p1)[1], genes(p2)[1]
        for i, j, idn in links(p1, p2, thr):
            x1, x2 = (g1[i][0] + g1[i][1]) / 2, (g2[j][0] + g2[j][1]) / 2; ya, yb = y0 - k * 1.6 - 0.3, y0 - (k + 1) * 1.6 + 0.3
            ax.plot([x1, x2], [ya, yb], color=BLUE, lw=0.6 + 1.2 * idn, alpha=0.25 + 0.5 * idn, zorder=0)
            lk.append(dict(phage_a=p1, gene_a=i + 1, phage_b=p2, gene_b=j + 1, identity=round(idn, 3)))
    ax.set_xlim(-L * 0.02, L * 1.01); ax.set_ylim(y0 - 1.6 * len(phages) + 0.8, y0 + 0.8); ax.axis("off")
    return pd.DataFrame(rows), pd.DataFrame(lk)


def mafft(phages):
    exe = os.environ.get("MAFFT") or "mafft"
    with tempfile.TemporaryDirectory() as d:
        open(f"{d}/in.fa", "w").write("".join(f">{p}\n{G.read_genome(p)}\n" for p in phages))
        out = subprocess.run([exe, "--auto", "--thread", "4", "--quiet", f"{d}/in.fa"], capture_output=True, text=True, check=True).stdout
    seqs, name = {}, None
    for l in out.splitlines():
        if l.startswith(">"): name = l[1:].strip(); seqs[name] = ""
        else: seqs[name] += l.strip().upper()
    return seqs


def entropy(seqs):
    A = np.array([list(s) for s in seqs.values()])
    H = np.zeros(A.shape[1])
    for nt in "ACGT":
        p = (A == nt).sum(0) / len(A)
        H -= np.where(p > 0, p * np.log2(np.where(p > 0, p, 1)), 0)
    return H


if __name__ == "__main__":
    # ---- Fig 3b: six Vectrevirus phages (the phages of the Source Data sheet Figure_3b: LI10_P4, LI10_P5, 409_P5, 409_P3, LF31_P3, 412_P2)
    six = ["LI10_P4", "LI10_P5", "409_P5", "409_P3", "LF31_P3", "412_P2"]
    M = load_matrix()
    order = [b for b in open(f"{ROOT}/data/genomics/bacteria/panacota/tree/370+host_leaf_order_root=albertii.txt").read().strip().split(",") if b in M.index]
    fig = plt.figure(figsize=(14, 7.5)); gs_ = fig.add_gridspec(2, 2, width_ratios=[2.2, 1], height_ratios=[1.6, 1], wspace=0.05, hspace=0.25)
    ax1 = fig.add_subplot(gs_[0, 0]); genes_df, links_df = draw_map(ax1, six); ax1.set_title("b  Six Vectrevirus phages: genes (RBP blue) and protein links (identity >= 30%)", loc="left", fontsize=9)
    seqs = mafft(six); H = entropy(seqs)
    roll = pd.Series(H).rolling(200, center=True, min_periods=1).mean()
    ax2 = fig.add_subplot(gs_[1, 0]); ax2.fill_between(range(len(H)), roll, color=BLUE, lw=0); ax2.set_xlim(0, len(H)); ax2.set_xlabel("Alignment position (nt)"); ax2.set_ylabel("Entropy (bits, 200-nt mean)"); style(ax2)
    ax3 = fig.add_subplot(gs_[:, 1]); ax3.imshow(M.loc[order, six].fillna(0).values.T, aspect="auto", cmap="Blues", vmin=0, vmax=4, interpolation="nearest"); ax3.set_yticks(range(6)); ax3.set_yticklabels(six, fontsize=8); ax3.set_xticks([]); ax3.set_xlabel("Strains (tree order)"); ax3.set_title("Host range (MLC 0-4)", loc="left", fontsize=9); style(ax3, None)
    save(fig, "fig3b_vectrevirus", genes_df)
    links_df.to_csv(f"{OUT}/fig3b_links.csv", index=False)
    pd.DataFrame({"position": range(len(H)), "entropy": H}).to_csv(f"{OUT}/fig3b_entropy.csv", index=False)
    # entropy inside the RBP genes against the other genes (gene coordinates -> alignment columns, per phage)
    inside, outside = [], []
    for p in six:
        cols = np.array([i for i, c in enumerate(seqs[p]) if c != "-"]); _, gl, rb = genes(p)
        for gi, (b, e, _, _) in enumerate(gl, 1):
            (inside if gi in rb else outside).append(H[cols[b:e]].mean())
    stat("Fig 3b", "mean entropy in RBP genes vs other genes (bits)", "RBP genes are hotspots (qualitative)", f"{np.mean(inside):.2f} vs {np.mean(outside):.2f}", "Fig 3b legend; ED Fig. 3")
    stat("Fig 3b", "mean entropy of the alignment (bits)", "not printed", round(float(H.mean()), 3), "Fig 3b")
    # ---- ED3: five genera
    genera = ["Vectrevirus", "Dhillonvirus", "Gamaleyavirus", "Felixounavirus", "Tequatrovirus"]
    have = set(RBP.phage)
    fig, axs = plt.subplots(len(genera), 1, figsize=(13, 3.2 * len(genera)))
    allg, alll = [], []
    for ax, gname in zip(axs, genera):
        cand = [p for p in COL.index[COL.Genus == gname] if p in have]
        pick, hosts = [], set()
        for p in sorted(cand):                         # up to four phages from different isolation hosts
            if COL.loc[p, "Phage_host"] not in hosts and len(pick) < 4: pick.append(p); hosts.add(COL.loc[p, "Phage_host"])
        gd, ld = draw_map(ax, pick); ax.set_title(f"{gname} ({len(cand)} phages with an RBP annotation; {len(pick)} shown)", loc="left", fontsize=9)
        allg.append(gd.assign(genus=gname)); alll.append(ld.assign(genus=gname))
    fig.tight_layout(); save(fig, "ed3_genus_maps", pd.concat(allg))
    pd.concat(alll).to_csv(f"{OUT}/ed3_links.csv", index=False)
    flush_stats("phage_maps")
