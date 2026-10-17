#!/usr/bin/env bash
# Remove regenerable intermediates after Phase 6. Keeps every result table and the inputs needed to re-run the later steps
# (assemblies, proteomes and genes, the pangenome family list, final alignment, tree, per-genome key outputs).
# Usage: DRY=1 scripts/cleanup_intermediates.sh   (only list sizes)      scripts/cleanup_intermediates.sh   (delete)
set -uo pipefail
cd "$(dirname "$0")/.."
R=replication; A=$R/T14_panacota
rm_() { for p in "$@"; do [ -e "$p" ] || continue; printf "%8s  %s\n" "$(du -sh "$p" 2>/dev/null | cut -f1)" "$p"; [ -n "${DRY:-}" ] || rm -rf "$p"; done; }

echo "== PanACoTA intermediates (annotate/pangenome/align work files)"
rm_ $A/annot/tmp_files $A/annot/Replicons $A/annot/gff3 $A/annot/LSTINFO $A/annot/Proteins_clean $A/align/Align-ESCO $A/align/List-ESCO $A/pan/tmp_ESCO.All.prt_0.8-mode1-th4 $A/pan/*.bin $A/pan/*quanti.txt $A/pan/*quali.txt
echo "== capsule: per-genome HMMER outputs (the *.tsv results are kept)"
n=0; for d in $R/T20_capsule/runs/*/hmmer_results; do [ -d "$d" ] && n=$((n+1)); done; echo "  $n directories of HMMER output"
[ -n "${DRY:-}" ] || rm -rf $R/T20_capsule/runs/*/hmmer_results
echo "== ClermonTyping: per-genome mash/blast files (all_phylogroups.txt is kept)"
echo "  $(find $R/T17_clermont/all -type f ! -name all_phylogroups.txt | wc -l | tr -d ' ') files"
[ -n "${DRY:-}" ] || find $R/T17_clermont/all -type f ! -name all_phylogroups.txt -delete
echo "== OMP BLAST index and Kaptive temp"
rm_ $R/T23_omp/db.* $R/T21_kaptive/kaptive_temp_*
echo "== macOS AppleDouble files (._*): metadata twins that the exFAT drive creates for every file"
echo "  $(find . -path ./.git -prune -o -name '._*' -type f -print | wc -l | tr -d ' ') files"
[ -n "${DRY:-}" ] || find . -path ./.git -prune -o -name '._*' -type f -delete
