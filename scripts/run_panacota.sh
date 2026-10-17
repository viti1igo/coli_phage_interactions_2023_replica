#!/usr/bin/env bash
# T14: core-genome tree with PanACoTA 1.3.1 (annotate, pangenome, corepers, align, tree --soft fasttree) on all assemblies.
# QC thresholds are relaxed (--l90 1000 --nbcont 9999): the defaults (L90<=100, <999 contigs) drop 10 genomes that are in the
# published tree (incl. the E. albertii outgroup); the paper does not state its thresholds.
# Annotation uses --prodigal (Prokka under Rosetta needs ~10 min per genome, ~8-13 h for 403; see docs/COMPARISON_REPORT.md).
# Resumable: a step is skipped when its output exists.   Usage: scripts/run_panacota.sh [threads]   Output: replication/T14_panacota/
set -euo pipefail
cd "$(dirname "$0")/.."
export MAMBA_ROOT_PREFIX=${MAMBA_ROOT_PREFIX:-$HOME/micromamba-phage}
T=${1:-6}
G=data/genomics/bacteria/assemblies
D=replication/T14_panacota
N=ESCO
pc() { micromamba run -n pan PanACoTA "$@"; }
mkdir -p $D
find $G -name '._*' -delete

[ -s $D/list.txt ] || (cd $G && ls *.fasta | sed "s/\$/ :: $N.0622/") > $D/list.txt
echo "genomes: $(wc -l < $D/list.txt)"

[ -d $D/annot/Proteins ] || pc annotate -d $G -r $D/annot -l $D/list.txt -n $N --threads $T --prodigal --l90 1000 --nbcont 9999 > $D/annotate.log 2>&1
L=$D/annot/LSTINFO-list.lst
echo "annotated genomes (LSTINFO rows): $(($(wc -l < $L)-1))"

[ -d $D/pan ] || pc pangenome -l $L -n $N -d $D/annot/Proteins -o $D/pan --threads $T > $D/pangenome.log 2>&1
PAN=$(ls $D/pan/PanGenome-*.lst | head -1)
echo "pangenome: $PAN ($(wc -l < $PAN) families)"

[ -d $D/core ] || pc corepers -p $PAN -o $D/core > $D/corepers.log 2>&1
CORE=$(ls $D/core/PersGenome_*.lst | head -1)
echo "core: $CORE ($(wc -l < $CORE) families)"

[ -d $D/align ] || pc align -c $CORE -l $L -n $N -d $D/annot -o $D/align --threads $T > $D/align.log 2>&1
ALN=$(ls $D/align/Phylo-$N/$N.*.aln | head -1)
echo "alignment: $ALN"

[ -d $D/tree ] || pc tree -a $ALN -o $D/tree -s fasttree --threads $T > $D/tree.log 2>&1
echo "tree: $(ls $D/tree/*.nwk $D/tree/*tree* 2>/dev/null | head -3)"
