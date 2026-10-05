#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
raw_dir="${repo_root}/reference/eQTL/shared/raw"
catalogue_dir="${raw_dir}/eqtl_catalogue_onek1k"
pgboost_dir="${raw_dir}/pgboost"
mkdir -p "${catalogue_dir}" "${pgboost_dir}"

for dataset_id in QTD000606 QTD000607 QTD000608 QTD000609 QTD000610 QTD000612 QTD000615 QTD000620; do
  url="https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/susie/QTS000038/${dataset_id}/${dataset_id}.credible_sets.tsv.gz"
  curl -L --fail --retry 3 -o "${catalogue_dir}/${dataset_id}.credible_sets.tsv.gz" "${url}"
done

curl -L --fail --retry 3 \
  -o "${raw_dir}/eqtl_catalogue_metadata.tsv" \
  "https://raw.githubusercontent.com/eQTL-Catalogue/eQTL-Catalogue-resources/master/tabix/tabix_ftp_paths.tsv"

tmp_clone="$(mktemp -d)"
trap 'rm -rf "${tmp_clone}"' EXIT
git clone --depth 1 https://github.com/elizabethdorans/pgBoost.git "${tmp_clone}/pgBoost"
cp "${tmp_clone}/pgBoost/resources/SNPs.bed.gz" "${pgboost_dir}/1000G_MAC5_SNPs.hg38.bed.gz"
cp "${tmp_clone}/pgBoost/resources/TSS.txt.gz" "${pgboost_dir}/TSS.hg38.tsv.gz"
cp "${tmp_clone}/pgBoost/resources/training_data.tsv.gz" "${pgboost_dir}/GTEx_training_PIPgt0.2_lt0.01.tsv.gz"

echo "Downloaded eQTL Catalogue OneK1K and pgBoost reference resources."
