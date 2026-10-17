#!/usr/bin/env bash
# Idempotent. Run from anywhere:  scripts/download_data.sh [--genomes]
#   (default)  repo data + unpack published feature importances / predictions
#   --genomes  also fetch the 403 Picard assemblies from figshare (~2 GB, md5-checked)
set -euo pipefail
cd "$(dirname "$0")/.."

[ -d data/genomics ] || { git clone --depth 1 https://github.com/mdmparis/coli_phage_interactions_2023 _tmp && cp -R _tmp/. . && rm -rf _tmp; }

for a in feature_importances predictions; do
  [ -d dev/predictions/results/$a ] || bsdtar -xf dev/predictions/results/$a.rar -C dev/predictions/results
done
rm -f dev/predictions/results/feature_importances/._* dev/predictions/results/predictions/._*

if [ "${1:-}" = "--genomes" ]; then
  mkdir -p data/genomics/bacteria/assemblies
  curl -s https://api.figshare.com/v2/articles/25941691 |
    python3 -c "import json,sys;[print(f['name'],f['download_url'],f['computed_md5']) for f in json.load(sys.stdin)['files']]" |
    while read -r name url md5; do
      out=data/genomics/bacteria/assemblies/$name
      [ -f "$out" ] && [ "$(md5 -q "$out")" = "$md5" ] && continue   # macOS md5; use md5sum on Linux
      curl -sL "$url" -o "$out"
      [ "$(md5 -q "$out")" = "$md5" ] || { echo "md5 mismatch: $name" >&2; exit 1; }
    done
fi
echo done
