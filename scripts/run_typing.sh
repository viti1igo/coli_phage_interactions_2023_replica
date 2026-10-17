#!/usr/bin/env bash
# T16 sequence type (mlst, scheme ecoli_achtman_4: the paper used SRST2 on reads, but only assemblies are public),
# T17 phylogroup (ClermonTyping), T21 Klebsiella-like capsule (Kaptive 2.0.6) on all assemblies. Resumable.
# Usage: [LIMIT=n] [OUT=dir] [THREADS=n] scripts/run_typing.sh      Output: replication/T16_mlst, T17_clermont, T21_kaptive
set -euo pipefail
cd "$(dirname "$0")/.."
export MAMBA_ROOT_PREFIX=${MAMBA_ROOT_PREFIX:-$HOME/micromamba-phage}
G=data/genomics/bacteria/assemblies; O=${OUT:-replication}; T=${THREADS:-2}
find $G -name '._*' -delete
FILES=$(ls $G/*.fasta | head -${LIMIT:-1000})
echo "genomes: $(echo $FILES | wc -w)"

mkdir -p $O/T16_mlst $O/T17_clermont $O/T21_kaptive
[ -s $O/T16_mlst/mlst.tsv ] || micromamba run -n mlst mlst --scheme ecoli_achtman_4 --threads $T $FILES > $O/T16_mlst/mlst.tsv 2> $O/T16_mlst/mlst.log
echo "T16 mlst rows: $(wc -l < $O/T16_mlst/mlst.tsv)"

if [ ! -s $O/T17_clermont/all/all_phylogroups.txt ]; then
  for f in $FILES; do echo "$PWD/$f"; done > $O/T17_clermont/files.txt
  (cd $O/T17_clermont && micromamba run -n clermont clermonTyping.sh --fastafile files.txt --name all --minimal > clermont.log 2>&1)
fi
echo "T17 clermont: $(ls $O/T17_clermont | tr '\n' ' ')"

# Kaptive: ~2 min per genome (Rosetta), so run KJ chunks in parallel and merge the tables.
# The figshare ROAR131 assembly has sequence glued to two header lines; makeblastdb (inside Kaptive) aborts on it and Kaptive then skips the
# genomes after it in its list. Rewrite such records into sanitized copies (original files untouched) and retry any genome still missing.
S=$O/T21_kaptive/sanitized; mkdir -p $S; KFILES=""
for f in $FILES; do
  if grep -q -E '^>.*[ACGTN]{20,}$' $f; then perl -pe 's/^(>\S*_cov_[0-9.]+)([ACGTN]{20,})$/$1\n$2/' $f > $S/$(basename $f); KFILES="$KFILES $S/$(basename $f)"; else KFILES="$KFILES $f"; fi
done
K=$O/T21_kaptive; KJ=${KJOBS:-3}; NF=$(echo $FILES | wc -w)
if [ ! -s $K/kaptive_results_table.txt ] || [ $(wc -l < $K/kaptive_results_table.txt) -le $NF ]; then
  rm -rf $K/chunk_* $K/lists $K/kaptive_results*; mkdir -p $K/lists
  echo $KFILES | tr ' ' '\n' | awk -v d=$K/lists -v n=$KJ '{print > (d "/chunk_" (NR % n))}'
  for l in $K/lists/chunk_*; do c=$(basename $l)
    micromamba run -n kaptive python tools/Kaptive/kaptive.py -a $(cat $l) -k tools/Kaptive/reference_database/Klebsiella_k_locus_primary_reference.gbk \
      -o $K/$c -t 1 --no_seq_out > $K/$c.log 2>&1 &
  done
  wait
  awk 'FNR==1 && NR!=1 {next} {print}' $K/chunk_*_table.txt > $K/kaptive_results_table.txt
  have=$(tail -n +2 $K/kaptive_results_table.txt | cut -f1 | sort); miss=""
  for f in $KFILES; do echo "$have" | grep -qx "$(basename $f .fasta)" || miss="$miss $f"; done
  if [ -n "$miss" ]; then
    micromamba run -n kaptive python tools/Kaptive/kaptive.py -a $miss -k tools/Kaptive/reference_database/Klebsiella_k_locus_primary_reference.gbk \
      -o $K/chunk_extra -t 3 --no_seq_out > $K/chunk_extra.log 2>&1
    awk 'FNR==1 && NR!=1 {next} {print}' $K/chunk_*_table.txt > $K/kaptive_results_table.txt
  fi
fi
echo "T21 kaptive: $(ls $O/T21_kaptive | tr '\n' ' ')"
