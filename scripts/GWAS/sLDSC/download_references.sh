#!/usr/bin/env bash
set -euo pipefail

ROOT="${1:-$(pwd)}"
GWAS_ROOT="$ROOT/results/validation/GWAS"
RAW="$GWAS_ROOT/sLDSC/reference/raw"
REF="$GWAS_ROOT/sLDSC/reference/processed"
SUM="$GWAS_ROOT/sLDSC/sumstats/raw"
SOFT="$GWAS_ROOT/sLDSC/software"
ZENODO="https://zenodo.org/records/10515792/files"
FG="https://storage.googleapis.com/finngen-public-data-r12"

mkdir -p "$RAW" "$REF" "$SUM" "$SOFT"

fetch() {
  local url="$1"
  local out="$2"
  if [[ -s "$out" ]]; then
    echo "present: $out"
  else
    curl -L --fail --retry 3 -o "$out" "$url"
  fi
}

fetch "$ZENODO/1000G_Phase3_plinkfiles.tgz?download=1" "$RAW/1000G_Phase3_plinkfiles.tgz"
fetch "$ZENODO/1000G_Phase3_frq.tgz?download=1" "$RAW/1000G_Phase3_frq.tgz"
fetch "$ZENODO/1000G_Phase3_baselineLD_v2.2_ldscores.tgz?download=1" "$RAW/1000G_Phase3_baselineLD_v2.2_ldscores.tgz"
fetch "$ZENODO/1000G_Phase3_weights_hm3_no_MHC.tgz?download=1" "$RAW/1000G_Phase3_weights_hm3_no_MHC.tgz"
fetch "$ZENODO/hm3_no_MHC.list.txt?download=1" "$RAW/hm3_no_MHC.list.txt"

tar -xzf "$RAW/1000G_Phase3_plinkfiles.tgz" -C "$REF"
tar -xzf "$RAW/1000G_Phase3_frq.tgz" -C "$REF"
tar -xzf "$RAW/1000G_Phase3_baselineLD_v2.2_ldscores.tgz" -C "$REF"
tar -xzf "$RAW/1000G_Phase3_weights_hm3_no_MHC.tgz" -C "$REF"
cp "$RAW/hm3_no_MHC.list.txt" "$REF/hm3_no_MHC.list.txt"

for trait in M13_RHEUMA T1D K11_IBD_STRICT G6_MS G6_ALZHEIMER D3_ANAEMIA; do
  fetch "$FG/summary_stats/release/finngen_R12_${trait}.gz" "$SUM/finngen_R12_${trait}.gz"
done

if [[ ! -d "$SOFT/ldsc/.git" ]]; then
  git clone https://github.com/bulik/ldsc.git "$SOFT/ldsc"
fi

echo "References downloaded. Run scripts/GWAS/sLDSC/setup_ldsc_env.sh once, then scripts/GWAS/sLDSC/run_sldsc.sh."
