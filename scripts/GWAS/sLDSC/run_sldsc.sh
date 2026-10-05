#!/usr/bin/env bash
set -euo pipefail

ROOT="${1:-$(pwd)}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE="$ROOT/results/validation/GWAS/sLDSC"
REF="$BASE/reference/processed"
LDSCORES="${SLDSC_LDSCORES_DIR:-$BASE/ldscores}"
SUMRAW="$BASE/sumstats/raw"
SUMMUNGED="${SLDSC_SUMMUNGED_DIR:-$BASE/sumstats/munged}"
RESULTS="${SLDSC_RESULTS_DIR:-$BASE/results}"
LOGS="${SLDSC_LOGS_DIR:-$BASE/logs}"
LDSC="$BASE/software/ldsc"
ENV_DIR="$BASE/software/ldsc-env"

mkdir -p "$LDSCORES" "$SUMMUNGED" "$RESULTS" "$LOGS"

if [[ ! -d "$LDSC/.git" || ! -d "$REF" ]]; then
  bash "$SCRIPT_DIR/download_references.sh" "$ROOT"
fi
if [[ ! -x "$ENV_DIR/bin/python" ]] || ! "$ENV_DIR/bin/python" -c 'import numpy, pandas, scipy, bitarray' >/dev/null 2>&1; then
  bash "$SCRIPT_DIR/setup_ldsc_env.sh" "$ROOT"
fi
PY="$ENV_DIR/bin/python"
PY3="${PYTHON3:-python3}"

MANIFEST="${SLDSC_MANIFEST:-$ROOT/results/validation/GWAS/SCENT_style/config/annotation_manifest.tsv}"
ANNOTATION_DIR="${SLDSC_ANNOTATION_DIR:-$BASE/annotations}"
if [[ ! -s "$MANIFEST" ]]; then
  Rscript "$ROOT/scripts/GWAS/prepare_annotations.R" "$ROOT"
fi

BIM1=$(find "$REF" -name '1000G.EUR.QC.1.bim' -print -quit)
if [[ -z "$BIM1" ]]; then
  echo "Cannot locate extracted 1000G PLINK files" >&2
  exit 1
fi
BFILE_PREFIX="${BIM1%1.bim}"
HM3="$REF/hm3_no_MHC.list.txt"
HM3_MERGE="$REF/hm3_no_MHC.1000G_EUR.snplist"

# The Zenodo hm3 file used by --print-snps contains rsIDs only, whereas
# munge_sumstats.py --merge-alleles requires SNP/A1/A2.  Recover the alleles
# from the matching 1000G EUR BIM files once and reuse the resulting list.
if [[ ! -s "$HM3_MERGE" ]]; then
  {
    printf 'SNP\tA1\tA2\n'
    awk 'NR==FNR { keep[$1]=1; next } ($2 in keep) { print $2 "\t" $5 "\t" $6 }' \
      "$HM3" "${BFILE_PREFIX}"{1..22}.bim
  } > "$HM3_MERGE"
fi

if [[ ! -s "$LDSCORES/annotations.complete" ]]; then
  "$PY3" "$SCRIPT_DIR/build_joint_annotations.py" \
    --manifest "$MANIFEST" \
    --annotation-dir "$ANNOTATION_DIR" \
    --bfile-prefix "$BFILE_PREFIX" \
    --out-dir "$LDSCORES" \
    --methods "${SLDSC_METHODS:-epilink_nb,epilink_poisson,scmultimap,signac_corrected,signac_uncorrected}" \
    --min-intervals "${SLDSC_MIN_INTERVALS:-200}"
fi

export SLDSC_PY="$PY"
export SLDSC_CODE="$LDSC"
export SLDSC_BFILE_PREFIX="$BFILE_PREFIX"
export SLDSC_LDSCORES="$LDSCORES"
export SLDSC_HM3="$HM3"
export SLDSC_LOGS="$LOGS"
seq 1 22 | xargs -n 1 -P "${SLDSC_N_JOBS:-4}" bash "$SCRIPT_DIR/compute_l2_chr.sh"

if [[ ! -s "$LDSCORES/split.complete" ]]; then
  "$PY3" "$SCRIPT_DIR/split_joint_ldscores.py" \
    --joint-prefix "$LDSCORES/joint." \
    --manifest "$LDSCORES/sldsc_annotation_manifest.tsv" \
    --out-dir "$LDSCORES"
  date -Iseconds > "$LDSCORES/split.complete"
fi

tail -n +2 "$ROOT/results/validation/GWAS/SCENT_style/config/traits.tsv" | cut -f1,5,6 | sort -u | while IFS=$'\t' read -r trait cases controls; do
  if [[ ! -s "$SUMMUNGED/${trait}.sumstats.gz" ]]; then
    n=$((cases + controls))
    "$PY" "$LDSC/munge_sumstats.py" \
      --sumstats "$SUMRAW/finngen_R12_${trait}.gz" \
      --out "$SUMMUNGED/${trait}" \
      --merge-alleles "$HM3_MERGE" \
      --snp rsids --a1 alt --a2 ref --p pval \
      --signed-sumstats beta,0 --N "$n" > "$LOGS/munge_${trait}.log" 2>&1
  fi
done

BASE1=$(find "$REF" -name 'baselineLD.1.l2.ldscore.gz' -print -quit)
WEIGHT1=$(find "$REF" -name 'weights.hm3_noMHC.1.l2.ldscore.gz' -print -quit)
FRQ1=$(find "$REF" -name '1000G.EUR.QC.1.frq' -print -quit)
BASE_PREFIX="${BASE1%1.l2.ldscore.gz}"
WEIGHT_PREFIX="${WEIGHT1%1.l2.ldscore.gz}"
FRQ_PREFIX="${FRQ1%1.frq}"

SLDSC_MANIFEST_FILE="$LDSCORES/sldsc_annotation_manifest.tsv" \
SLDSC_JOB_FILE="$LDSCORES/sldsc_jobs.tsv" Rscript "$SCRIPT_DIR/make_jobs.R" "$ROOT"
export SLDSC_RESULTS="$RESULTS"
export SLDSC_SUMMUNGED="$SUMMUNGED"
export SLDSC_BASE_PREFIX="$BASE_PREFIX"
export SLDSC_WEIGHT_PREFIX="$WEIGHT_PREFIX"
export SLDSC_FRQ_PREFIX="$FRQ_PREFIX"
tail -n +2 "$LDSCORES/sldsc_jobs.tsv" | \
  xargs -r -n 3 -P "${SLDSC_H2_N_JOBS:-8}" bash "$SCRIPT_DIR/run_h2_job.sh"

export SLDSC_LDSCORE_DIR="$LDSCORES"
export SLDSC_RESULTS_DIR="$RESULTS"
Rscript "$SCRIPT_DIR/collect_results.R" "$ROOT"
Rscript "$SCRIPT_DIR/plot_sldsc.R" "$ROOT"
