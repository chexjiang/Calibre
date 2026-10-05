#!/usr/bin/env bash
set -euo pipefail

ROOT="${1:-$(pwd)}"
REF="$ROOT/results/validation/GWAS/reference"
FG="$REF/FinnGen"
mkdir -p "$FG/raw/summary" "$REF/RA" "$REF/T1D" "$REF/IBD" "$REF/liftover"

fetch() {
  local url="$1"
  local out="$2"
  if [[ -s "$out" ]]; then
    echo "present: $out"
  else
    curl -L --fail --retry 3 -o "$out" "$url"
  fi
}

fetch "https://storage.googleapis.com/finngen-public-data-r12/summary_stats/finngen_R12_manifest.tsv" \
  "$FG/finngen_R12_manifest.tsv"

# The full public standard SuSiE 95% credible-set collection is about 80 MB.
# Case-count filtering (>1,000) is performed by prepare_extended_references.R.
gsutil -m cp -n 'gs://finngen-public-data-r12/finemap/summary/*.SUSIE.snp.filter.tsv' \
  "$FG/raw/summary/"

fetch "https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41588-022-01213-w/MediaObjects/41588_2022_1213_MOESM4_ESM.xlsx" \
  "$REF/RA/Ishigaki_2022_supplementary_tables.xlsx"
fetch "https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41586-021-03552-w/MediaObjects/41586_2021_3552_MOESM4_ESM.xlsx" \
  "$REF/T1D/Chiou_2021_T1D_99pct_credible_sets.xlsx"
fetch "https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fnature22969/MediaObjects/41586_2017_BFnature22969_MOESM2_ESM.xlsx" \
  "$REF/IBD/Huang_2017_supplementary_table1.xlsx"
fetch "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz" \
  "$REF/liftover/hg19ToHg38.over.chain.gz"

echo "Reference downloads complete under $REF"
