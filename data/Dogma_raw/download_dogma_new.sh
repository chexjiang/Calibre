#!/usr/bin/env bash
set -euo pipefail

# Download the four DOGMA-seq control/stimulation libraries from the original
# Mimitou et al. GEO submission (GSE156478).  The modality files are stored
# under different GSM accessions but belong to the same library.

ROOT="$(cd "$(dirname "$0")" && pwd)"
RAW="$ROOT/Dogma_new_raw"
BASE="https://ftp.ncbi.nlm.nih.gov/geo/samples"
mkdir -p "$RAW"

download() {
  local name="$1"
  local accession="${name%%_*}"
  local series_dir="${accession:0:7}nnn"
  local url="$BASE/$series_dir/$accession/suppl/$name"
  local output="$RAW/$name"
  if [[ -s "$output" ]]; then
    echo "exists: $name"
    return
  fi
  echo "downloading: $name"
  curl -L --fail --retry 3 -C - -o "$output" "$url"
}

# GEX+ATAC filtered 10x matrices.  These contain the raw counts that will be
# jointly QC-filtered below.
for name in \
  GSM5065525_LLL_CTRL_GExp_ATAC_filtered_barcodes.tsv.gz \
  GSM5065525_LLL_CTRL_GExp_ATAC_filtered_features.tsv.gz \
  GSM5065525_LLL_CTRL_GExp_ATAC_filtered_matrix.mtx.gz \
  GSM5065528_LLL_STIM_GExp_ATAC_filtered_barcodes.tsv.gz \
  GSM5065528_LLL_STIM_GExp_ATAC_filtered_features.tsv.gz \
  GSM5065528_LLL_STIM_GExp_ATAC_filtered_matrix.mtx.gz \
  GSM5065531_DIG_CTRL_GExp_ATAC_filtered_barcodes.tsv.gz \
  GSM5065531_DIG_CTRL_GExp_ATAC_filtered_features.tsv.gz \
  GSM5065531_DIG_CTRL_GExp_ATAC_filtered_matrix.mtx.gz \
  GSM5065534_DIG_STIM_GExp_ATAC_filtered_barcodes.tsv.gz \
  GSM5065534_DIG_STIM_GExp_ATAC_filtered_features.tsv.gz \
  GSM5065534_DIG_STIM_GExp_ATAC_filtered_matrix.mtx.gz; do
  download "$name"
done

# Cell Ranger per-barcode metrics include the fragment/nucleosome/TSS metrics
# used by the paper's DOGMA QC. Fragment files are retained for provenance
# and optional independent Signac QC. Set DOGMA_DOWNLOAD_FRAGMENTS=0 if only
# the matrices and QC metrics are needed.
for name in \
  GSM5065524_LLL_CTRL_per_barcode_metrics.csv.gz \
  GSM5065527_LLL_STIM_per_barcode_metrics.csv.gz \
  GSM5065530_DIG_CTRL_per_barcode_metrics.csv.gz \
  GSM5065533_DIG_STIM_per_barcode_metrics.csv.gz; do
  download "$name"
done

if [[ "${DOGMA_DOWNLOAD_FRAGMENTS:-1}" == "1" ]]; then
  for name in \
    GSM5065524_LLL_CTRL_fragments.tsv.gz \
    GSM5065527_LLL_STIM_fragments.tsv.gz \
    GSM5065530_DIG_CTRL_fragments.tsv.gz \
    GSM5065533_DIG_STIM_fragments.tsv.gz; do
    download "$name"
  done
else
  echo "Skipping fragment files (DOGMA_DOWNLOAD_FRAGMENTS=0)"
fi

cat > "$RAW/README.txt" <<'EOF'
Source: Mimitou et al., Nat Biotechnol. 2021, GEO GSE156478.

Condition mapping:
  LLL_CTRL  = GSM5065524 metrics/fragments + GSM5065525 GEX/ATAC matrix
  LLL_STIM  = GSM5065527 metrics/fragments + GSM5065528 GEX/ATAC matrix
  DIG_CTRL  = GSM5065530 metrics/fragments + GSM5065531 GEX/ATAC matrix
  DIG_STIM  = GSM5065533 metrics/fragments + GSM5065534 GEX/ATAC matrix

Each matrix contains the filtered barcode set for the corresponding library;
the metrics/fragments GSM accession is the same library's ATAC record.
EOF

echo "DOGMA-seq raw files are in $RAW"
