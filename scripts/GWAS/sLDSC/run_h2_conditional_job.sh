#!/usr/bin/env bash
set -euo pipefail

annotation_id="$1"
background_id="$2"
trait="$3"

outdir="$SLDSC_RESULTS/$trait"
mkdir -p "$outdir"
out="$outdir/$annotation_id"
target="$SLDSC_LDSCORES/split/$annotation_id/$annotation_id."
background="$SLDSC_LDSCORES/split/$background_id/$background_id."

if [[ ! -s "${out}.results" ]]; then
  "$SLDSC_PY" "$SLDSC_CODE/ldsc.py" \
    --h2 "$SLDSC_SUMMUNGED/${trait}.sumstats.gz" \
    --ref-ld-chr "$SLDSC_BASE_PREFIX,$background,$target" \
    --w-ld-chr "$SLDSC_WEIGHT_PREFIX" \
    --overlap-annot \
    --frqfile-chr "$SLDSC_FRQ_PREFIX" \
    --print-coefficients \
    --out "$out" > "$SLDSC_LOGS/h2_${trait}_${annotation_id}.log" 2>&1
fi

echo "conditional S-LDSC complete: $trait $annotation_id"
