#!/usr/bin/env bash
# One micromamba env per tool, at the paper's versions where bioconda has them (see config/paper_config.yaml).
# Tries native arm64 first, falls back to osx-64 (needs Rosetta). Only conda-forge and bioconda are used.
# Envs live on the internal disk (the external exFAT drive wastes ~10x space on conda's small files); outputs stay in the project.
# Usage: scripts/setup_envs.sh [env ...]   builds (or rebuilds) the named envs, default: all.  Always re-checks all envs
# and writes results/env_check.csv.  Written for bash 3.2 (macOS default): no associative arrays.
set -uo pipefail
cd "$(dirname "$0")/.."
export MAMBA_ROOT_PREFIX=${MAMBA_ROOT_PREFIX:-$HOME/micromamba-phage}
mkdir -p logs results
CH="--override-channels -c conda-forge -c bioconda"
ALL=(pan ectyper kaptive defense macsy clermont blast mlst)

# panacota 1.3.1 needs python-utils<3 (python 3.6); defense-finder 1.2.0 needs macsyfinder<=2.1.4 (module macsypy);
# kaptive 2.0.6 needs biopython<1.78 and is not on PyPI: run the authors' kaptive.py from tag v2.0.6 (cloned to tools/Kaptive; python 3.8, osx-64).
spec() { case $1 in
  pan) echo "panacota=1.3.1 python-utils<3.0";; ectyper) echo "ectyper=1.0.0";; kaptive) echo "python=3.8 blast";;
  defense) echo "defense-finder=1.2.0 macsyfinder=2.1.4";; macsy) echo "macsyfinder";; clermont) echo "clermontyping";;
  blast) echo "blast";; mlst) echo "mlst";; esac; }
check() { case $1 in
  pan) echo "PanACoTA --version";; ectyper) echo "ectyper --version";; kaptive) echo "python tools/Kaptive/kaptive.py --version";;
  defense) echo "defense-finder --help";; macsy) echo "macsyfinder --version";; clermont) echo "clermonTyping.sh -v";;
  blast) echo "blastp -version";; mlst) echo "mlst --version";; esac; }

ENVS=("$@"); [ ${#ENVS[@]} -eq 0 ] && ENVS=("${ALL[@]}")
for e in "${ENVS[@]}"; do
  micromamba env remove -y -n "$e" >/dev/null 2>&1
  s=$(spec "$e"); plat=arm64
  [ "$e" = kaptive ] && plat=osx-64
  if [ $plat = arm64 ]; then
    micromamba create -y -n "$e" $CH $s >"logs/env_$e.log" 2>&1 || plat=osx-64
  fi
  [ $plat = osx-64 ] && micromamba create -y -n "$e" --platform osx-64 $CH $s >>"logs/env_$e.log" 2>&1
  if [ "$e" = kaptive ]; then
    micromamba run -n kaptive pip install -q "biopython<1.78" >>"logs/env_$e.log" 2>&1
    [ -d tools/Kaptive ] || git clone -q --depth 1 --branch v2.0.6 https://github.com/katholt/Kaptive.git tools/Kaptive >>"logs/env_$e.log" 2>&1
  fi
  echo "$e built ($plat)"
  echo "$plat" > "logs/env_$e.platform"
done

echo "env,platform,status,version_output" > results/env_check.csv
for e in "${ALL[@]}"; do
  plat=$(cat "logs/env_$e.platform" 2>/dev/null || echo unknown)
  out=$(micromamba run -n "$e" $(check "$e") 2>&1); rc=$?
  st=FAIL; [ $rc -eq 0 ] && st=OK
  echo "$e,$plat,$st,$(echo "$out" | tail -1 | tr ',' ';' | cut -c1-110)" | tee -a results/env_check.csv
done
