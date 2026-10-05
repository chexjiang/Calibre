#!/usr/bin/env bash
set -euo pipefail

ROOT="${1:-$(pwd)}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE="$ROOT/results/validation/GWAS/sLDSC"
REF="$BASE/reference/processed"
LDSC="$BASE/software/ldsc"
ENV_DIR="$BASE/software/ldsc-env"
SUMMUNGED="$BASE/sumstats/munged"

if [[ ! -x "$ENV_DIR/bin/python" ]]; then
  bash "$SCRIPT_DIR/setup_ldsc_env.sh" "$ROOT"
fi
PY="$ENV_DIR/bin/python"
PY3="${PYTHON3:-python3}"

BIM1=$(find "$REF" -name '1000G.EUR.QC.1.bim' -print -quit)
BFILE_PREFIX="${BIM1%1.bim}"
HM3="$REF/hm3_no_MHC.list.txt"
BASE1=$(find "$REF" -name 'baselineLD.1.l2.ldscore.gz' -print -quit)
WEIGHT1=$(find "$REF" -name 'weights.hm3_noMHC.1.l2.ldscore.gz' -print -quit)
FRQ1=$(find "$REF" -name '1000G.EUR.QC.1.frq' -print -quit)
BASE_PREFIX="${BASE1%1.l2.ldscore.gz}"
WEIGHT_PREFIX="${WEIGHT1%1.l2.ldscore.gz}"
FRQ_PREFIX="${FRQ1%1.frq}"

if [[ "${SLDSC_SKIP_PREP:-0}" != "1" ]]; then
  Rscript "$SCRIPT_DIR/prepare_scmultimap_style_annotations.R" "$ROOT" 0.05
fi

for VIEW in celltype pooled; do
  VIEW_ROOT="$BASE/source_results/$VIEW"
  ANN="$VIEW_ROOT/annotations"
  BG="$VIEW_ROOT/backgrounds"
  LDSCORES="$VIEW_ROOT/ldscores"
  RESULTS="$VIEW_ROOT/results/raw"
  LOGS="$VIEW_ROOT/logs"
  mkdir -p "$LDSCORES" "$RESULTS" "$LOGS"

  if [[ ! -s "$LDSCORES/chunk_manifest.tsv" ]]; then
    "$PY3" "$SCRIPT_DIR/build_conditional_annotations.py" \
      --target-manifest "$VIEW_ROOT/manifest.tsv" \
      --background-manifest "$VIEW_ROOT/background_manifest.tsv" \
      --annotation-dir "$ANN" \
      --background-dir "$BG" \
      --bfile-prefix "$BFILE_PREFIX" \
      --out-dir "$LDSCORES" \
      --chunk-size "${SLDSC_ANNOTATION_CHUNK_SIZE:-20}"
  fi

  export SLDSC_PY="$PY"
  export SLDSC_CODE="$LDSC"
  export SLDSC_BFILE_PREFIX="$BFILE_PREFIX"
  export SLDSC_LDSCORES="$LDSCORES"
  export SLDSC_HM3="$HM3"
  export SLDSC_LOGS="$LOGS"
  find "$LDSCORES/chunks" -mindepth 1 -maxdepth 1 -type d -name 'chunk_*' -printf '%f\n' | sort | \
    while read -r chunk; do seq 1 22 | awk -v c="$chunk" '{print c "\t" $1}'; done | \
    xargs -r -n 2 -P "${SLDSC_N_JOBS:-4}" bash "$SCRIPT_DIR/compute_l2_chunk.sh"

  if [[ ! -s "$LDSCORES/split.complete" ]]; then
    "$PY3" "$SCRIPT_DIR/split_chunk_ldscores.py" \
      --chunk-root "$LDSCORES/chunks" \
      --out-dir "$LDSCORES"
    date -Iseconds > "$LDSCORES/split.complete"
  fi

  if [[ ! -s "$LDSCORES/annot.complete" ]]; then
    find "$LDSCORES/chunks" -mindepth 1 -maxdepth 1 -type d -name 'chunk_*' -printf '%f\n' | \
      while read -r chunk; do seq 1 22 | awk -v c="$chunk" '{print c "\t" $1}'; done | \
      xargs -r -n 2 -P "${SLDSC_ANNOT_N_JOBS:-12}" bash -c \
      'python3 "'"$SCRIPT_DIR"'/make_split_annotation_files.py" \
        --chunk-root "'"$LDSCORES"'/chunks" \
        --out-dir "'"$LDSCORES"'" \
        --bfile-prefix "'"$BFILE_PREFIX"'" \
        --chunk-name "$1" --chrom "$2"' _
    date -Iseconds > "$LDSCORES/annot.complete"
  fi

  Rscript "$SCRIPT_DIR/make_conditional_jobs.R" "$ROOT" "$VIEW"
  export SLDSC_RESULTS="$RESULTS"
  export SLDSC_SUMMUNGED="$SUMMUNGED"
  export SLDSC_BASE_PREFIX="$BASE_PREFIX"
  export SLDSC_WEIGHT_PREFIX="$WEIGHT_PREFIX"
  export SLDSC_FRQ_PREFIX="$FRQ_PREFIX"
  tail -n +2 "$LDSCORES/sldsc_jobs.tsv" | cut -f1,2,4 | \
    xargs -r -n 3 -P "${SLDSC_H2_N_JOBS:-8}" bash "$SCRIPT_DIR/run_h2_conditional_job.sh"

  Rscript "$SCRIPT_DIR/collect_conditional_results.R" "$ROOT" "$VIEW"
done

Rscript "$SCRIPT_DIR/plot_scmultimap_style.R" "$ROOT"
