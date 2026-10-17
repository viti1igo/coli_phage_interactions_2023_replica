"""Readers for the paper's Supplementary Tables (MOESM4) and Source Data workbooks (MOESM5-13) in docs/paper/supplement/ (downloaded from the
publisher; not redistributed in the GitHub repo).  sheet(file_no, sheet_name, text_first=True) -> DataFrame (title row dropped, first column as text unless text_first=False)."""
import os
import pandas as pd

SUPP = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "docs", "paper", "supplement")


def sheet(n, name, text_first=True):
    d = pd.read_excel(f"{SUPP}/MOESM{n}.xlsx", sheet_name=name, header=1)
    if text_first:
        d[d.columns[0]] = d[d.columns[0]].astype(str)
    return d


def available():
    return os.path.exists(f"{SUPP}/MOESM4.xlsx")


def phage_name_map():
    """New phage names of the paper (Supplementary Tables 2, 3, 7) -> public names of the authors' repository, found by identical interaction vectors."""
    import numpy as np
    root = os.path.join(SUPP, "..", "..", "..")
    st3 = sheet(4, "Supplementary Table 3").set_index("bacteria").apply(pd.to_numeric, errors="coerce")
    pub = pd.read_csv(f"{root}/data/interactions/interaction_matrix.csv", sep=";", index_col=0)
    pub.index = pub.index.astype(str)
    key = lambda v: np.nan_to_num(v, nan=-9)
    out = {}
    for c in st3.columns:
        hit = [m for m in pub.columns if np.array_equal(key(st3.loc[pub.index, c].values), key(pub[m].values))]
        assert len(hit) == 1, (c, hit)
        out[c] = hit[0]
    assert len(set(out.values())) == len(out)
    return out


def matrix403():
    """Interaction matrix of all 403 strains (Supplementary Table 3) with the public phage names; the public file lacks strain LF110."""
    m = phage_name_map()
    st3 = sheet(4, "Supplementary Table 3").set_index("bacteria").apply(pd.to_numeric, errors="coerce")
    return st3.rename(columns=m)
