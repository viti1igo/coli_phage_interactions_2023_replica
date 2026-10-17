"""Interaction matrix loader: the 403-strain matrix of the paper (Supplementary Table 3, public phage names) when the supplement has been
downloaded (scripts/download_supplement.sh), otherwise the authors' public 402-strain matrix (it lacks strain LF110)."""
import os
import pandas as pd
import supp

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")


def load_matrix():
    if supp.available():
        return supp.matrix403()
    m = pd.read_csv(f"{ROOT}/data/interactions/interaction_matrix.csv", sep=";", index_col=0)
    m.index = m.index.astype(str)
    return m
