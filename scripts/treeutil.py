"""Patristic (path-length) distances between the leaves of a Newick tree."""
import numpy as np
from Bio import Phylo


def patristic(tree_path):
    """tree_path: a Newick file path or an already loaded Bio.Phylo tree."""
    tree = Phylo.read(tree_path, "newick") if isinstance(tree_path, str) else tree_path
    leaves = tree.get_terminals()
    paths = []  # per leaf: [(clade, cumulative root distance)]
    for l in leaves:
        d, p = 0.0, []
        for c in tree.get_path(l):
            d += c.branch_length or 0.0
            p.append((c, d))
        paths.append(p)
    n = len(leaves)
    D = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            k = 0  # deepest shared clade
            while k < min(len(paths[i]), len(paths[j])) and paths[i][k][0] is paths[j][k][0]:
                k += 1
            lca = paths[i][k - 1][1] if k else 0.0
            D[i, j] = D[j, i] = paths[i][-1][1] + paths[j][-1][1] - 2 * lca
    return [l.name for l in leaves], D
