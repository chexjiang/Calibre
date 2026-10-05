#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
OUT="$ROOT/reference/CRISPR_new"
TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT

mkdir -p "$OUT"

# The Nature Genetics paper used the ensemble K562 benchmark before the
# upstream repository reorganized it in 2026. Pin the last paper-era commit.
CRISPR_COMMIT="554443c7d840643473279e8d62514eacbc807b41"
PGBOOST_COMMIT="58a6ceb774d76d2523a2e5223681c065c2bfaa04"

git clone https://github.com/EngreitzLab/CRISPR_comparison.git "$TMP/CRISPR_comparison"
git -C "$TMP/CRISPR_comparison" show \
  "$CRISPR_COMMIT:resources/crispr_data/EPCrisprBenchmark_ensemble_data_GRCh38.tsv.gz" \
  > "$OUT/EPCrisprBenchmark_ensemble_data_GRCh38.tsv.gz"

git clone --depth 1 https://github.com/elizabethdorans/pgBoost.git "$TMP/pgBoost"
git -C "$TMP/pgBoost" fetch --depth 1 origin "$PGBOOST_COMMIT"
git -C "$TMP/pgBoost" checkout --detach "$PGBOOST_COMMIT"
cp "$TMP/pgBoost/resources/SNPs.bed.gz" "$OUT/SNPs.bed.gz"

pushd "$OUT" >/dev/null
sha256sum \
  EPCrisprBenchmark_ensemble_data_GRCh38.tsv.gz \
  SNPs.bed.gz > SHA256SUMS
popd >/dev/null

echo "Downloaded paper-era CRISPR benchmark and pgBoost SNP universe to $OUT"
