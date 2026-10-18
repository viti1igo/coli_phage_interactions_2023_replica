"""Edits applied to the authors' notebooks in `run_notebooks.py patched`. Each one replaces a Windows/private path by its public equivalent or
fixes an incompatibility with current pandas/matplotlib; nothing changes an analysis. Listed in reports/10_repository_analyses.md."""
import re

BS = "\\" * 2  # two backslashes as they appear in the notebook source
PATCHES = {
    "visualize_predictions.ipynb": [
        (re.escape(BS), "/"),                                                       # Windows path separators written as "\\"
        (r'preds_dir = f"C:/Users/hvays/[^"]*"', 'preds_dir = f"{base_dir}/results/predictions"'),  # extracted predictions.rar
        (r'"D:/These/[^"]*ordered_tree_leaves_370\+host_origname\.txt"', '"../../data/genomics/bacteria/panacota/tree/370+host_leaf_order_root=albertii.txt"'),  # public leaf order
        (r'"D:/These/[^"]*/([^"/]+)"', r'"output/\1"'),                               # author's figure folders
        (r"(?m)^import os$", 'import os\nos.makedirs("output", exist_ok=True)'),
    ],
    "analyze_cocktail_interactions_tailored_vs_generic_vs_baseline.ipynb": [
        (r'plt\.savefig\("([^"]*)\.tsv"\)', r'plt.savefig("\1.png")'),              # the notebook saves a figure under a .tsv name
        (r'"C:/Users/hvays/[^"]*/([^"/]+)"', r'"output/\1"'),                       # author's output folders
        (r"(?m)^import os$", 'import os\nos.makedirs("output", exist_ok=True)'),
        (r'\.fillna\(""\)\.apply\(lambda x: "/"\.join\(x\), axis=1\)', '.fillna("").astype(str).apply(lambda x: "/".join(x), axis=1)'),     # numeric ST column
        (r'\.groupby\("bacteria"\)\["chosen_at_step"\]\.mean\(\)', '.groupby("bacteria")["chosen_at_step"].apply(lambda x: x.astype(float).mean())'),
        (r"nearest_neighbours\[0\]", "nearest_neighbours.iloc[0]"),                 # positional indexing on a labelled Series
        (r'bact_features\.loc\[bact, \["ST_Warwick", "O-type", "H-type"\]\]', 'bact_features.reindex(bact)[["ST_Warwick", "O-type", "H-type"]]'),  # test strain BCH857 is called BCH857-72 in the feature table
    ],
    "build_bacteria_arsenale_pairwise_distance_table.ipynb": [
        (r'\.\./panacota/370\+host_leaf_order_root=albertii\.txt', "../panacota/tree/370+host_leaf_order_root=albertii.txt"),   # file now lives in panacota/tree
    ],
}
