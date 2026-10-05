#!/usr/bin/env bash
set -euo pipefail

chunk="$1"
chr="$2"
chunk_dir="$SLDSC_LDSCORES/chunks/$chunk"
out="$chunk_dir/joint.$chr"

if [[ ! -s "${out}.l2.ldscore.gz" ]]; then
  "$SLDSC_PY" "$SLDSC_CODE/ldsc.py" \
    --l2 \
    --bfile "${SLDSC_BFILE_PREFIX}${chr}" \
    --ld-wind-cm 1 \
    --annot "$chunk_dir/joint.${chr}.annot.gz" \
    --thin-annot \
    --out "$out" \
    --print-snps "$SLDSC_HM3" > "$SLDSC_LOGS/l2_${chunk}_chr${chr}.log" 2>&1
fi
echo "LD scores complete: $chunk chr$chr"
