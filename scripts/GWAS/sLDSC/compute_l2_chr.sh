#!/usr/bin/env bash
set -euo pipefail

chr="$1"
if [[ ! -s "$SLDSC_LDSCORES/joint.${chr}.l2.ldscore.gz" ]]; then
  "$SLDSC_PY" "$SLDSC_CODE/ldsc.py" \
    --l2 \
    --bfile "${SLDSC_BFILE_PREFIX}${chr}" \
    --ld-wind-cm 1 \
    --annot "$SLDSC_LDSCORES/joint.${chr}.annot.gz" \
    --thin-annot \
    --out "$SLDSC_LDSCORES/joint.${chr}" \
    --print-snps "$SLDSC_HM3" > "$SLDSC_LOGS/l2_chr${chr}.log" 2>&1
fi
echo "LD scores complete: chr${chr}"
