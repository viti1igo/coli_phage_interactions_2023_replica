#!/usr/bin/env bash
# T22 defense systems (DefenseFinder 1.2.0 + models 1.2.3, from github.com/mdmparis) and
# T20 ABC capsule typing (MacSyFinder with the authors' custom models in data/genomics/bacteria/capsules/ABC_capsules_types)
# on the PanACoTA proteomes (gembase format). DefenseFinder works in a temp dir on the internal disk: on the exFAT drive macOS adds
# ._* files that break its post-processing; only the small result TSVs are copied back. Resumable: a finished genome is skipped.
# Usage: [LIMIT=n] [OUT=dir] [JOBS=n] scripts/run_proteome_tools.sh [defense|capsule|all]
set -uo pipefail
cd "$(dirname "$0")/.."
export MAMBA_ROOT_PREFIX=${MAMBA_ROOT_PREFIX:-$HOME/micromamba-phage}
P=replication/T14_panacota/annot/Proteins_clean; O=${OUT:-replication}; J=${JOBS:-2}
WHAT=${1:-all}
# PanACoTA headers look like >ESCO.0622.00001.0001i_00002 <description>; the b/i suffix (border/inside) makes MacSyFinder's
# gembase parser see one contig as two replicons. Rewrite as >ESCO.0622.00001.0001_00002 (gene numbers are already ordered).
if [ ! -d $P ]; then mkdir -p $P; for f in $(ls replication/T14_panacota/annot/Proteins/*.prt | grep -v '\.All\.prt$'); do
  awk '/^>/ {id=$1; n=id; sub(/.*_/,"",n); sub(/_[0-9]+$/,"",id); sub(/[bi]$/,"",id); print id "_" n; next} {print}' $f > $P/$(basename $f); done; fi
PROTS=$(ls $P/*.prt | grep -v '\.All\.prt$' | head -${LIMIT:-1000})   # ESCO.All.prt = all proteomes concatenated by the pangenome step
echo "proteomes: $(echo $PROTS | wc -w)"

if [ $WHAT = defense ] || [ $WHAT = all ]; then
  M=$O/T22_defense/models; mkdir -p $O/T22_defense/runs $M
  [ -d $M/defense-finder-models ] || micromamba run -n defense macsydata install --org mdmparis --target $M defense-finder-models==1.2.3 > $O/T22_defense/models.log 2>&1
  # DefenseFinder 1.2.x also needs the CasFinder models (which `defense-finder update` normally installs)
  [ -d $M/CasFinder ] || micromamba run -n defense macsydata install --org macsy-models --target $M CasFinder==3.1.0 >> $O/T22_defense/models.log 2>&1
  find $M -name '._*' -delete   # AppleDouble junk from the exFAT drive breaks MacSyFinder's model parser
  echo "models: $(ls $M | tr '\n' ' ') | defense-finder-models $(grep -i '^vers' $M/defense-finder-models/metadata.yml) | CasFinder $(grep -i '^vers' $M/CasFinder/metadata.yml)"
  export O M
  for f in $PROTS; do n=$(basename $f .prt); [ -s $O/T22_defense/runs/$n/${n}_defense_finder_systems.tsv ] || echo $f; done |
    xargs -n 1 -P $J sh -c 'n=$(basename $0 .prt); w=$(mktemp -d); micromamba run -n defense defense-finder run $0 -o $w --db-type gembase --models-dir $M -w 1 > $O/T22_defense/runs_$n.log 2>&1; mkdir -p $O/T22_defense/runs/$n; cp $w/*.tsv $O/T22_defense/runs/$n/ 2>/dev/null; rm -rf $w'
  echo "T22 genomes done: $(ls $O/T22_defense/runs/*/*_defense_finder_systems.tsv 2>/dev/null | wc -l)"
fi

if [ $WHAT = capsule ] || [ $WHAT = all ]; then
  CM=$O/T20_capsule/models; mkdir -p $O/T20_capsule/runs $CM   # a clean models dir: MacSyFinder scans every sub-directory
  [ -d $CM/ABC_capsules_types ] || cp -R data/genomics/bacteria/capsules/ABC_capsules_types $CM/
  find $CM -name '._*' -delete
  export O CM
  for f in $PROTS; do n=$(basename $f .prt); [ -d $O/T20_capsule/runs/$n ] || echo $f; done |
    xargs -n 1 -P $J sh -c 'n=$(basename $0 .prt); micromamba run -n macsy macsyfinder --sequence-db $0 --db-type gembase --models-dir $CM --models ABC_capsules_types all -o $O/T20_capsule/runs/$n --worker 1 > $O/T20_capsule/runs_$n.log 2>&1'
  echo "T20 genomes done: $(ls -d $O/T20_capsule/runs/*/ 2>/dev/null | wc -l)"
fi
