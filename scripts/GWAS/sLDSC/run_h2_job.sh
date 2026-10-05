#!/usr/bin/env bash
set -euo pipefail

annotation_id="$1"
panel="$2"  # retained in the job table for auditing
trait="$3"

outdir="$SLDSC_RESULTS/$trait"
mkdir -p "$outdir"
out="$outdir/$annotation_id"
target="$SLDSC_LDSCORES/split/$annotation_id/$annotation_id."

if [[ ! -s "${out}.results" ]]; then
  "$SLDSC_PY" "$SLDSC_CODE/ldsc.py" \
    --h2 "$SLDSC_SUMMUNGED/${trait}.sumstats.gz" \
    --ref-ld-chr "$SLDSC_BASE_PREFIX,$target" \
    --w-ld-chr "$SLDSC_WEIGHT_PREFIX" \
    --overlap-annot \
    --frqfile-chr "$SLDSC_FRQ_PREFIX" \
    --print-coefficients \
    --out "$out" > "$SLDSC_LOGS/h2_${trait}_${annotation_id}.log" 2>&1
fi

echo "sLDSC complete: $trait $annotation_id"
