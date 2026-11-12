"""Phase 21 (T93): run folder for biophage-reproduce on the authors' shipped inputs (same inputs as T2.1).
Copies the biophage code unchanged to replication/T93_biophage/<env>/ and builds the pipeline_results/ files
run_training.py reads.   python scripts/biophage_inputs.py bp_asis bp_pinned"""
import shutil, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "external/biophage-reproduce-main"
B = ROOT / "data/genomics/bacteria"
for env in sys.argv[1:]:
    run = ROOT / "replication/T93_biophage" / env
    if run.exists():
        sys.exit(f"{run} exists; remove it first")
    shutil.copytree(SRC, run, ignore=shutil.ignore_patterns("tests", "__pycache__"))
    pr = run / "pipeline_results"
    (pr / "panacota/tree").mkdir(parents=True)
    (pr / "panacota/annotate").mkdir(parents=True)
    shutil.copy(B / "picard_collection.csv", pr / "picard_collection.csv")
    shutil.copy(B / "umap_phylogeny/coli_umap_8_dims.tsv", pr / "umap_coordinates.tsv")
    # tree leaves already carry strain names, so a header-only LSTINFO maps every name to itself
    shutil.copy(B / "isolation_strains/panacota/tree/370+host_tree_orignames.nwk",
                pr / "panacota/tree/ESCO.nucl.grp.aln.fasttree_tree.nwk")
    (pr / "panacota/annotate/LSTINFO-file_list.lst").write_text("gembase_name\torig_name\n")
    print("built", run)
