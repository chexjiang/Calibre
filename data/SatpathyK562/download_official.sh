#!/usr/bin/env bash
set -euo pipefail

# Official public IGVF files for Satpathy et al. K562 multiome data.

root=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
out="$root/data/SatpathyK562/raw"
mkdir -p "$out"
base='https://api.data.igvf.org'

download() {
  local name="$1" url="$2"
  if [[ -s "$out/$name" ]]; then
    echo "exists: $name"
  else
    echo "download: $name"
    curl -L --fail --retry 3 --silent --show-error "$url" -o "$out/$name"
  fi
}

download IGVFFI2641LMFU.h5 "$base/matrix-files/IGVFFI2641LMFU/@@download/IGVFFI2641LMFU.h5"
download IGVFFI5156ENLL.tsv.gz "$base/tabular-files/IGVFFI5156ENLL/@@download/IGVFFI5156ENLL.tsv.gz"
download IGVFFI8318PQLH.bed.gz "$base/tabular-files/IGVFFI8318PQLH/@@download/IGVFFI8318PQLH.bed.gz"
download IGVFFI3930OIYP.tbi "$base/index-files/IGVFFI3930OIYP/@@download/IGVFFI3930OIYP.tbi"
download IGVFFI0627MDWY.tsv.gz "$base/tabular-files/IGVFFI0627MDWY/@@download/IGVFFI0627MDWY.tsv.gz"
if [[ ! -e "$out/IGVFFI0627MDWY.tsv.gz.tbi" ]]; then
  cp "$out/IGVFFI3930OIYP.tbi" "$out/IGVFFI0627MDWY.tsv.gz.tbi"
fi

echo "Downloaded IGVFDS2018URDY files into $out"
